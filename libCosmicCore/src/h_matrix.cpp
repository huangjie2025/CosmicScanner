#include "h_matrix.h"
#include <Eigen/Dense>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <fstream>
#include <sstream>
#include <limits>
#include <iomanip>
#include <iostream>
#include <Eigen/Sparse>
#include <cstddef>
#include <string>
#include <vector>
#include <cstdint>
#include <iostream>

namespace cosmic {
namespace genmatrix {
  Eigen::SparseMatrix<double> buildA_noninline(const std::vector<int>& dam,
                                               const std::vector<int>& sire);
  Eigen::MatrixXd buildA_dense_noninline(const std::vector<int>& dam,
                                         const std::vector<int>& sire);
  Eigen::SparseMatrix<double> invertA_henderson_noninline(const std::vector<int>& dam,
                                                          const std::vector<int>& sire);
}
// Merged GRM declarations (replacing #include "GRM.h")

struct GRMResult {
    std::size_t n{};
    std::size_t m{};
    std::size_t mReported{};
    std::size_t skipped_all_missing{};
    std::size_t skipped_low_var{};
    std::size_t kept_zero{};
    std::vector<std::string> ids;
    std::vector<double> G;
    std::vector<double> G_upper;
    std::vector<double> Ginv;
    std::vector<uint32_t> N;
    std::vector<uint32_t> N_upper;
};

class GRM {
public:
    static std::size_t parseFAM(const std::string &famPath, std::vector<std::string> &ids);
    static std::size_t countBIM(const std::string &bimPath);
    static void readBEDHeader(std::FILE *bed);
    static std::size_t decodeOneSNP(std::FILE *bed, std::size_t n, std::vector<int> &x);
    static std::size_t decodeOneSNPFromBuffer(const unsigned char*& buf, std::size_t n, std::vector<int>& out_x);
    static void countCodesPopcnt(const unsigned char* buf,
                                 std::size_t bytes_per_snp,
                                 std::size_t n,
                                 uint64_t& out_non_missing,
                                 uint64_t& out_sum_allele1);

    static void axpyOuter(const std::vector<double> &z, std::vector<double> &S);
    static void syrkUpdateUpperRowMajor(std::vector<double> &S, std::size_t n,
                                        const std::vector<double> &Z, std::size_t k,
                                        std::size_t ld);
    static void syrkUpdateUpperPacked(std::vector<double> &upper, std::size_t n,
                                      const std::vector<double> &Z, std::size_t k,
                                      std::size_t ld);
    static void invertSPD(std::vector<double> &A, std::size_t n, double ridge, std::vector<double> &Ainv);

    static GRMResult buildAndInvert(const std::string &prefix, double ridge);

    static void computeN(const std::vector<std::vector<uint64_t>> &missBits,
                         std::size_t n, std::size_t mUsed, std::vector<uint32_t> &N);
};





static inline bool isNa(double x){ return std::isnan(x); }

static NamedMatrix reorder(const NamedMatrix& M, const std::vector<std::string>& order){
  std::unordered_map<std::string,int> idx;
  for(int i=0;i<(int)M.rows.size();++i) idx[M.rows[i]]=i;
  int n=(int)order.size();
  Matrix out(n,n);
  for(int i=0;i<n;++i){
    for(int j=0;j<n;++j){
      auto it_i = idx.find(order[i]);
      auto it_j = idx.find(order[j]);
      out(i,j) = (it_i!=idx.end() && it_j!=idx.end()) ? M.data(it_i->second,it_j->second) : 0.0;
    }
  }
  return NamedMatrix{out,order,order};
}

static NamedMatrix subset(const NamedMatrix& M, const std::vector<int>& rIdx, const std::vector<int>& cIdx){
  Matrix out(rIdx.size(), cIdx.size());
  for(size_t i=0;i<rIdx.size();++i){
    for(size_t j=0;j<cIdx.size();++j){ out((int)i,(int)j)=M.data(rIdx[i], cIdx[j]); }
  }
  std::vector<std::string> r, c;
  r.reserve(rIdx.size()); c.reserve(cIdx.size());
  for(auto i: rIdx) r.push_back(M.rows[i]);
  for(auto j: cIdx) c.push_back(M.cols[j]);
  return NamedMatrix{out,r,c};
}

static std::vector<std::string> unique_union_rev(const std::vector<std::string>& a, const std::vector<std::string>& b){
  std::vector<std::string> u;
  std::unordered_set<std::string> seen;
  for(const auto& x: a){ if(seen.insert(x).second) u.push_back(x); }
  for(const auto& x: b){ if(seen.insert(x).second) u.push_back(x); }
  std::reverse(u.begin(), u.end());
  return u;
}

static std::vector<int> match_indices(const std::vector<std::string>& keys, const std::vector<std::string>& dict){
  std::unordered_map<std::string,int> idx;
  for(int i=0;i<(int)dict.size();++i) idx[dict[i]]=i;
  std::vector<int> out; out.reserve(keys.size());
  for(const auto& k: keys){ auto it=idx.find(k); out.push_back(it==idx.end()?-1:it->second); }
  return out;
}

static Matrix get_ASV(const Matrix& x){
  double tr = x.diagonal().sum();
  double denom = tr / (x.rows()-1);
  return x / denom;
}

static Matrix tcrossprod(const Matrix& X){ return X * X.transpose(); }

static double safe_sum(const Matrix& v){ double s=0.0; for(int i=0;i<v.size();++i){ s+=v(i); } return s; }

static Matrix Mcheck(const Matrix& SNP, int ploidy, double maf, bool rmv_mono, double thresh_htzy, double thresh_missing, const std::string& impute){
  Matrix X = SNP;
  int n = X.rows();
  int m = X.cols();
  std::vector<int> keep(m,1);
  for(int j=0;j<m;++j){
    int nonMissing=0; int count0=0; int count2=0; double sum=0.0; int cnt=0;
    for(int i=0;i<n;++i){ double v=X(i,j); if(!isNa(v)){ nonMissing++; sum+=v; cnt++; if(v==0) count0++; if(v==2) count2++; }}
    double missRate = 1.0 - (double)nonMissing/n;
    if(thresh_missing<1.0 && missRate>thresh_missing){ keep[j]=0; continue; }
    double p = (double)count2 / std::max(1,count0+count2+ (int)(cnt - count0 - count2));
    if(maf>0.0){ double mafv = std::min(p, 1.0-p); if(mafv < maf){ keep[j]=0; continue; }}
    double mean = cnt? (sum/cnt) : 0.0;
    for(int i=0;i<n;++i){ if(isNa(X(i,j))){ if(impute=="mean" || impute=="global.mean") X(i,j)=mean; else X(i,j)=mean; } }
    if(rmv_mono){ bool mono=true; double first=std::numeric_limits<double>::quiet_NaN(); for(int i=0;i<n;++i){ if(!isNa(SNP(i,j))){ if(std::isnan(first)) first=SNP(i,j); else if(SNP(i,j)!=first){ mono=false; break; } } } if(mono) keep[j]=0; }
  }
  int kept = std::accumulate(keep.begin(), keep.end(), 0);
  Matrix Y(n, kept);
  int col=0; for(int j=0;j<m;++j){ if(keep[j]){ Y.col(col++) = X.col(j); } }
  return Y;
}

static void check_Gmatrix_data(const Matrix& SNPmatrix, const std::string& method, int ploidy, bool ratio, bool integer){
  if(method!="Yang" && method!="VanRaden" && method!="MarkersMatrix") throw std::runtime_error("Invalid method");
  int n=SNPmatrix.rows(), m=SNPmatrix.cols();
  if(n==0 || m==0) throw std::runtime_error("Empty SNPmatrix");
  if(ratio) throw std::runtime_error("ratio model not supported in diploid-only configuration");
  if(ploidy!=2) throw std::runtime_error("Only ploidy=2 supported");
  double mx=-1e9, mn=1e9;
  for(int i=0;i<n;++i){ for(int j=0;j<m;++j){ double v=SNPmatrix(i,j); if(!isNa(v)){ mx=std::max(mx,v); mn=std::min(mn,v);} } }
  if(mx>2.0 || mn<0.0) throw std::runtime_error("Diploid coding out of bounds");
  if(integer){ for(int i=0;i<n;++i){ for(int j=0;j<m;++j){ double v=SNPmatrix(i,j); if(!isNa(v) && std::floor(v)!=v) throw std::runtime_error("Non-integer values"); } } }
}

static Matrix Gmatrix(const Matrix& SNPmatrix_in, const std::string& method, double missingValue, double maf, bool verify_posdef, int ploidy, bool pseudo_diploid, bool integer, bool ratio, const std::string& impute_method, bool rmv_mono, double thresh_htzy, bool ratio_check, const std::vector<double>& weights, bool ploidy_correction, bool ASV){
  Matrix SNPmatrix = SNPmatrix_in;
  int n = SNPmatrix.rows();
  int m = SNPmatrix.cols();
  if(!std::isnan(missingValue)){
    for(int i=0;i<n;++i){ for(int j=0;j<m;++j){ if(SNPmatrix(i,j)==missingValue) SNPmatrix(i,j)=std::numeric_limits<double>::quiet_NaN(); } }
  }
  check_Gmatrix_data(SNPmatrix, method, ploidy, ratio, integer);
  int NumberMarkers = SNPmatrix.cols();
  if(!ratio){ SNPmatrix = Mcheck(SNPmatrix, ploidy, maf, rmv_mono, thresh_htzy, 0.50, impute_method); }
  Matrix G;
  if(method=="MarkersMatrix"){
    Matrix B = SNPmatrix.unaryExpr([](double v){ return isNa(v)?0.0:1.0; });
    G = tcrossprod(B);
    return ASV? get_ASV(G):G;
  }
  auto alleleFreq = [&](const Eigen::VectorXd& col, int y){ int cnt=0, tot=0; for(int i=0;i<col.size();++i){ double v=col(i); if(!isNa(v)){ tot++; if((int)v==y) cnt++; } } return tot? (double)cnt/tot : 0.0; };
  Matrix Frequency(n,2);
  Eigen::VectorXd f0(m), f2(m);
  for(int j=0;j<m;++j){ f0(j)=alleleFreq(SNPmatrix.col(j),0); f2(j)=alleleFreq(SNPmatrix.col(j),2); }
  Frequency = Eigen::MatrixXd(m,2);
  for(int j=0;j<m;++j){ Frequency(j,0)=f0(j); Frequency(j,1)=f2(j); }
  Matrix FreqP(n,m);
  for(int j=0;j<m;++j){ for(int i=0;i<n;++i){ FreqP(i,j)=Frequency(j,1); } }
  if(method=="VanRaden"){
    double TwoPQ = 0.0; for(int j=0;j<m;++j){ TwoPQ += 2.0*Frequency(j,0)*Frequency(j,1); }
    Matrix X = SNPmatrix;
    for(int i=0;i<n;++i){ for(int j=0;j<m;++j){ double v=X(i,j); v = isNa(v)?0.0:(v - 2.0*FreqP(i,j)); X(i,j)=v; } }
    if(weights.empty()){
      G = tcrossprod(X) / TwoPQ;
    }else{
      Eigen::VectorXd w(m);
      for(int j=0;j<m;++j){ w(j) = weights[j]; }
      Matrix D = w.asDiagonal();
      G = tcrossprod(X*D) / TwoPQ;
    }
    if(ASV) G = get_ASV(G);
    return G;
  }
  if(method=="Yang"){
    Matrix FreqPQ(n,m);
    for(int j=0;j<m;++j){ double v=2.0*Frequency(j,0)*Frequency(j,1); for(int i=0;i<n;++i){ FreqPQ(i,j)=v; } }
    Matrix Gall(n,m);
    for(int i=0;i<n;++i){ for(int j=0;j<m;++j){ double x=SNPmatrix(i,j); double fp=FreqP(i,j); double v = isNa(x)?0.0:(x*x - (1.0 + 2.0*fp)*x + 2.0*(fp*fp)); Gall(i,j)= v / (FreqPQ(i,j)==0?1.0:FreqPQ(i,j)); } }
    Eigen::VectorXd Gii(m);
    for(int j=0;j<m;++j){ double s=0.0; for(int i=0;i<n;++i) s+=Gall(i,j); Gii(j)=s; }
    Matrix X(n,m);
    for(int i=0;i<n;++i){ for(int j=0;j<m;++j){ double fp=FreqPQ(i,j); double x=SNPmatrix(i,j); double v = isNa(x)?0.0:(x - 2.0*FreqP(i,j)); X(i,j)= v / (fp==0?1.0:std::sqrt(fp)); } }
    double NumberMarkersD = (double)NumberMarkers;
    Eigen::VectorXd Gii_hat = Eigen::VectorXd::Ones(m) + Gii/NumberMarkersD;
    Matrix Gm = tcrossprod(X)/NumberMarkersD;
    for(int i=0;i<n;++i) Gm(i,i) = Gii_hat(i<m?i:0);
    if(ASV) Gm = get_ASV(Gm);
    return Gm;
  }
  return Matrix();
}

NamedMatrix Hmatrix(const NamedMatrix& A_in,
                           const NamedMatrix& G_in,
                           const Matrix* markers,
                           double c,
                           const std::string& method,
                           double tau,
                           double omega,
                           double missingValue,
                           double maf,
                           int ploidy,
                           int roundVar,
                           bool ASVflag){

  // Force G.cpp to compute the inverse
  extern bool g_do_inverse;
  g_do_inverse = true;
  NamedMatrix Aorig = A_in;
  NamedMatrix Gorig = G_in;
  const auto& An = Aorig.rows;
  const auto& Gn = Gorig.rows;
  std::unordered_set<std::string> setG(Gn.begin(), Gn.end());
  std::unordered_set<std::string> setA(An.begin(), An.end());
  std::vector<std::string> Gnhat; Gnhat.reserve(Gn.size());
  for(const auto& id: Gn){ if(setA.find(id)!=setA.end()) Gnhat.push_back(id); }
  std::vector<std::string> Anhat; Anhat.reserve(An.size());
  for(const auto& id: An){ if(setG.find(id)!=setG.end()) Anhat.push_back(id); }
  NamedMatrix A = Aorig;
  std::vector<int> gIdxRows; gIdxRows.reserve(Gnhat.size());
  for(const auto& id: Gnhat){ auto it=std::find(Gorig.rows.begin(), Gorig.rows.end(), id); gIdxRows.push_back((int)std::distance(Gorig.rows.begin(), it)); }
  NamedMatrix G = subset(Gorig, gIdxRows, gIdxRows);
  if(method=="Martini"){
    auto idA = A.rows;
    auto idG = G.rows;
    auto idH = unique_union_rev(idG, idA);
    A = reorder(A, idH);
    auto matches = match_indices(idH, idG);
    std::vector<bool> index; index.reserve(idH.size());
    for(int k=0;k<(int)idH.size();++k) index.push_back(matches[k]==-1);
    std::vector<int> idx11, idx22;
    for(int k=0;k<(int)index.size();++k){ if(index[k]) idx11.push_back(k); else idx22.push_back(k); }
    auto A11 = subset(A, idx11, idx11);
    auto A12 = subset(A, idx11, idx22);
    auto A21 = subset(A, idx22, idx11);
    auto A22 = subset(A, idx22, idx22);
    std::vector<int> g22rows; g22rows.reserve(idx22.size());
    for(int k=0;k<(int)idx22.size();++k){ auto name = idH[idx22[k]]; auto it=std::find(G.rows.begin(), G.rows.end(), name); g22rows.push_back((int)std::distance(G.rows.begin(), it)); }
    auto G22 = subset(G, g22rows, g22rows);
    Eigen::LLT<Matrix> lltA22(A22.data);
    Matrix A22inv;
    if(lltA22.info() == Eigen::Success) {
        A22inv = lltA22.solve(Matrix::Identity(A22.data.rows(), A22.data.cols()));
    } else {
        A22inv = A22.data.inverse();
    }
    Eigen::FullPivLU<Matrix> luG(G22.data);
    if(!luG.isInvertible()) throw std::runtime_error("G22 not invertible");
    Matrix G22inv = luG.inverse();
    Matrix H22m_inv = tau*G22inv + (1.0-omega)*A22inv;
    Eigen::LLT<Matrix> lltH22m(H22m_inv);
    Matrix H22m;
    if(lltH22m.info() == Eigen::Success) {
        H22m = lltH22m.solve(Matrix::Identity(H22m_inv.rows(), H22m_inv.cols()));
    } else {
        H22m = H22m_inv.inverse();
    }
    Matrix delta = H22m - A22.data;
    Matrix H11m = A12.data * A22inv * delta * A22inv * A21.data;
    Matrix H12m = A12.data * A22inv * delta;
    Matrix H21m = delta * A22inv * A21.data;
    Matrix top(A.data.rows(), A.data.cols());

    top = A.data;
    for(size_t i=0;i<idx11.size();++i){ for(size_t j=0;j<idx11.size();++j){ top(idx11[i], idx11[j]) += H11m((int)i,(int)j); } }
    for(size_t i=0;i<idx22.size();++i){ for(size_t j=0;j<idx11.size();++j){ top(idx22[i], idx11[j]) += H21m((int)i,(int)j); } }
    for(size_t i=0;i<idx11.size();++i){ for(size_t j=0;j<idx22.size();++j){ top(idx11[i], idx22[j]) += H12m((int)i,(int)j); } }
    for(size_t i=0;i<idx22.size();++i){ for(size_t j=0;j<idx22.size();++j){ top(idx22[i], idx22[j]) += delta((int)i,(int)j); } }
    Matrix Hm = top;
    if(ASVflag) Hm = get_ASV(Hm);
    NamedMatrix res{Hm, A.rows, A.cols};
    // Reorder back to original A order to match user expectation if they provided A
    // AGHmatrix returns H with rows/cols as union of G and A.
    // However, for consistency with input A, we might want to reorder if the user expects A-like order.
    // The R code does: H = A + ... (where A is reordered to idH).
    // So the result H has idH order.
    // Let's keep it as is (idH order) which is consistent with R logic.
    return res;
  }
  if(method=="Munoz"){
    std::vector<int> aIdxRows; aIdxRows.reserve(Anhat.size());
    for(const auto& id: Anhat){ auto it=std::find(Aorig.rows.begin(), Aorig.rows.end(), id); aIdxRows.push_back((int)std::distance(Aorig.rows.begin(), it)); }
    NamedMatrix Aco = subset(Aorig, aIdxRows, aIdxRows);
    if(markers==nullptr) throw std::runtime_error("For Munoz method you need to specify markers");
    Matrix markersmatrix = Gmatrix(*markers, "MarkersMatrix", std::numeric_limits<double>::quiet_NaN(), maf, false, ploidy, false, true, false, "mean", false, 0.0, false, std::vector<double>(), false, false);
    int n = Aco.data.rows();
    std::vector<double> classes;
    {
      std::unordered_set<long long> s;
      for(int i=0;i<n;++i){ for(int j=0;j<n;++j){ double v=Aco.data(i,j); double r=std::round(v*std::pow(10,roundVar))/std::pow(10,roundVar); long long key=(long long)std::llround(r*1e9); if(s.insert(key).second) classes.push_back(r); } }
    }
    Matrix varA = Matrix::Constant(n,n, std::numeric_limits<double>::quiet_NaN());
    Matrix meanG = Matrix::Constant(n,n, std::numeric_limits<double>::quiet_NaN());
    std::vector<double> varAclasses(classes.size(), std::numeric_limits<double>::quiet_NaN());
    auto roundA = [&](double v){ return std::round(v*std::pow(10,roundVar))/std::pow(10,roundVar); };
    for(size_t ci=0; ci<classes.size(); ++ci){ double cls=classes[ci]; std::vector<double> vals; vals.reserve(n*n); for(int i=0;i<n;++i){ for(int j=0;j<n;++j){ if(roundA(Aco.data(i,j))==cls){ double gv = G.data(i,j); if(!isNa(gv)) vals.push_back(gv); } } } if(!vals.empty()){ double mu = std::accumulate(vals.begin(), vals.end(), 0.0)/vals.size(); double v=0.0; for(double x: vals){ v += (x-mu)*(x-mu); } v = vals.size()>1 ? v/(vals.size()-1) : 0.0; varAclasses[ci]=v; for(int i=0;i<n;++i){ for(int j=0;j<n;++j){ if(roundA(Aco.data(i,j))==cls){ varA(i,j)=v; meanG(i,j)=mu; } } } }
    }
    for(auto& v: varAclasses){ if(v==0.0 || !std::isfinite(v)) v = std::numeric_limits<double>::quiet_NaN(); }
    for(size_t ci=0; ci<classes.size(); ++ci){ if(std::isnan(varAclasses[ci])){ double fill=0.0; size_t cnt=0; for(size_t k=0;k<classes.size();++k){ if(!std::isnan(varAclasses[k])){ fill+=varAclasses[k]; cnt++; } } if(cnt>0){ varAclasses[ci]=fill/cnt; for(int i=0;i<n;++i){ for(int j=0;j<n;++j){ if(roundA(Aco.data(i,j))==classes[ci]) varA(i,j)=varAclasses[ci]; } } } }
    }
    Matrix invMarkers = markersmatrix.block(0,0,G.data.rows(),G.data.cols());
    for(int i=0;i<invMarkers.rows();++i){ for(int j=0;j<invMarkers.cols();++j){ double v=invMarkers(i,j); invMarkers(i,j)= v==0.0 ? 0.0 : 1.0/v; } }
    Matrix beta = Matrix::Constant(n,n,0.0);
    for(int i=0;i<n;++i){ for(int j=0;j<n;++j){ double va=varA(i,j); double mm= invMarkers(i,j); double b = 1.0 - ( c + ( mm ) / (std::isnan(va)?1.0:va) ); beta(i,j)=b; } }
    Matrix Hm = Matrix::Zero(n,n);
    for(int i=0;i<n;++i){ for(int j=0;j<n;++j){ double h = beta(i,j)*(G.data(i,j)-Aco.data(i,j)) + Aco.data(i,j); Hm(i,j)=h; } }
    NamedMatrix Aout = Aorig;
    for(size_t ii=0; ii<Anhat.size(); ++ii){ int ai = aIdxRows[ii]; for(size_t jj=0; jj<Anhat.size(); ++jj){ int aj = aIdxRows[jj]; Aout.data(ai,aj) = Hm((int)ii,(int)jj); } }
    if(ASVflag) Aout.data = get_ASV(Aout.data);
    return Aout;
  }
  throw std::runtime_error("Unknown method");
}

// Helper to get indices of names in a list
std::vector<int> match_names(const std::vector<std::string>& target, const std::vector<std::string>& source) {
    std::vector<int> indices;
    std::unordered_map<std::string, int> source_map;
    for (int i = 0; i < source.size(); ++i) {
        source_map[source[i]] = i;
    }
    for (const auto& name : target) {
        if (source_map.find(name) != source_map.end()) {
            indices.push_back(source_map[name]);
        } else {
            indices.push_back(-1); // Not found
        }
    }
    return indices;
}

// Helper to subset/reorder a matrix based on name indices
NamedMatrix subsetMatrix(const NamedMatrix& source, const std::vector<std::string>& new_order) {
    int n = new_order.size();
    Matrix new_mat(n, n);
    std::vector<int> idx = match_names(new_order, source.rows);

    for (int i = 0; i < n; ++i) {
        if (idx[i] == -1) throw std::runtime_error("Name not found in source matrix: " + new_order[i]);
        for (int j = 0; j < n; ++j) {
            if (idx[j] == -1) throw std::runtime_error("Name not found in source matrix: " + new_order[j]);
            new_mat(i, j) = source.data(idx[i], idx[j]);
        }
    }
    return NamedMatrix{new_mat, new_order, new_order};
}

class HInverseBuilder {
public:
    double tau = 1.0;
    double omega = 1.0;
    bool keep_order = true;
    bool message = true;
    NamedMatrix A;
    NamedMatrix Ainv;
    NamedMatrix Ginv;
    bool has_A = false;
    bool has_Ainv = false;

    void setA(const Matrix& mat, const std::vector<std::string>& names) {
        A = NamedMatrix{mat, names, names};
        has_A = true;
    }

    void setAinv(const Matrix& mat, const std::vector<std::string>& names) {
        Ainv = NamedMatrix{mat, names, names};
        has_Ainv = true;
    }

    void setGinv(const Matrix& mat, const std::vector<std::string>& names) {
        Ginv = NamedMatrix{mat, names, names};
    }

    NamedMatrix compute() {
        if (!has_A && !has_Ainv) throw std::runtime_error("A or Ainv needs to be provided.");
        if (has_Ainv && !has_A) {
            if (message) std::cout << "A matrix not provided. Calculating A from Ainv..." << std::endl;
            Eigen::LLT<Matrix> lltAinv(Ainv.data);
            if(lltAinv.info() == Eigen::Success) {
                A.data = lltAinv.solve(Matrix::Identity(Ainv.data.rows(), Ainv.data.cols()));
            } else {
                A.data = Ainv.data.inverse();
            }
            A.rows = Ainv.rows;
            A.cols = Ainv.cols;
            has_A = true;
        } else if (has_A && !has_Ainv) {
            if (message) std::cout << "Ainv matrix not provided. Calculating Ainv from A..." << std::endl;
            Eigen::LLT<Matrix> lltA(A.data);
            if(lltA.info() == Eigen::Success) {
                Ainv.data = lltA.solve(Matrix::Identity(A.data.rows(), A.data.cols()));
            } else {
                Ainv.data = A.data.inverse();
            }
            Ainv.rows = A.rows;
            Ainv.cols = A.cols;
            has_Ainv = true;
        }

        std::vector<std::string> Aind = Ainv.rows;
        std::vector<std::string> Gind = Ginv.rows;
        std::unordered_set<std::string> setG(Gind.begin(), Gind.end());
        std::unordered_set<std::string> setA(Aind.begin(), Aind.end());

        for (const auto& name : Gind) {
            if (setA.find(name) == setA.end()) throw std::runtime_error("Matrix Ginv has individuals that are not present in A or Ainv.");
        }

        std::vector<std::string> combined = Gind;
        combined.insert(combined.end(), Aind.begin(), Aind.end());
        std::vector<std::string> unique_combined;
        std::unordered_set<std::string> seen;
        for (const auto& s : combined) {
            if (seen.find(s) == seen.end()) {
                unique_combined.push_back(s);
                seen.insert(s);
            }
        }
        std::reverse(unique_combined.begin(), unique_combined.end());
        std::vector<std::string> idH = unique_combined; 

        NamedMatrix A_sorted = subsetMatrix(A, idH);
        NamedMatrix Ainv_sorted = subsetMatrix(Ainv, idH);
        
        std::vector<int> genotyped_indices;
        std::vector<std::string> current_genotyped_names;
        for(int i=0; i<idH.size(); ++i) {
            if (setG.count(idH[i])) {
                genotyped_indices.push_back(i);
                current_genotyped_names.push_back(idH[i]);
            }
        }
        
        NamedMatrix Ginv_sorted = subsetMatrix(Ginv, current_genotyped_names);

        int n_g = genotyped_indices.size();
        Matrix A22(n_g, n_g);
        for(int i=0; i<n_g; ++i) {
            for(int j=0; j<n_g; ++j) {
                A22(i,j) = A_sorted.data(genotyped_indices[i], genotyped_indices[j]);
            }
        }

        Matrix I22 = Matrix::Identity(n_g, n_g);
        Eigen::LLT<Matrix> llt22(A22);
        Matrix A22inv;
        if(llt22.info() == Eigen::Success){
            A22inv = llt22.solve(I22);
        } else {
            Eigen::LDLT<Matrix> ldlt22(A22);
            if(ldlt22.info() == Eigen::Success){
                A22inv = ldlt22.solve(I22);
            } else {
                A22inv = A22.inverse();
            }
        }
        Matrix H22inv = (tau * Ginv_sorted.data) - (omega * A22inv);
        Matrix Hinv_mat = Ainv_sorted.data;
        
        for(int i=0; i<n_g; ++i) {
            for(int j=0; j<n_g; ++j) {
                int row_idx = genotyped_indices[i];
                int col_idx = genotyped_indices[j];
                Hinv_mat(row_idx, col_idx) += H22inv(i, j);
            }
        }

        NamedMatrix Hinv{Hinv_mat, idH, idH};
        if (keep_order) {
            Hinv = subsetMatrix(Hinv, Ainv.rows);
        }
        return Hinv;
    }
};


NamedMatrix HInverseMatrix(const NamedMatrix& A22,
                                  const Eigen::SparseMatrix<double>& Ainv_sparse,
                                  const std::vector<std::string>& a_ids,
                                  const NamedMatrix& Ginv_in,
                                  double tau,
                                  double omega) {
    auto idA = a_ids;
    auto idG = Ginv_in.rows;
    auto idH = unique_union_rev(idG, idA);

            std::unordered_map<std::string, int> h_name_to_idx;
    for(int i=0; i<(int)idH.size(); ++i) h_name_to_idx[idH[i]] = i;

    // Ainv_sparse uses a_ids order.
    // We assume a_ids[k] corresponds to k-th row/col in Ainv_sparse.
    
    int nH = (int)idH.size();
    Matrix HinvData = Matrix::Zero(nH, nH);
    
    // Iterate sparse matrix
    for (int k=0; k<Ainv_sparse.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(Ainv_sparse, k); it; ++it) {
            double val = it.value();
            int r_a = (int)it.row();
            int c_a = (int)it.col();
            
            // Map r_a -> name -> r_h
            std::string name_r = a_ids[r_a];
            std::string name_c = a_ids[c_a];
            
            int r_h = h_name_to_idx[name_r];
            int c_h = h_name_to_idx[name_c];
            
            HinvData(r_h, c_h) = val;
        }
    }
    
    NamedMatrix Hinv{HinvData, idH, idH};

    std::unordered_set<std::string> setG(idG.begin(), idG.end());
    std::vector<int> genotyped_indices;
    std::vector<std::string> genotyped_names;
    for(int i=0; i<(int)idH.size(); ++i){
        if(setG.find(idH[i]) != setG.end()){
            genotyped_indices.push_back(i);
            genotyped_names.push_back(idH[i]);
        }
    }
    std::cout << "[DEBUG] Genotyped IDs found: " << genotyped_indices.size() << std::endl;
    if(!genotyped_indices.empty()){
        int mn = genotyped_indices[0], mx = genotyped_indices[0];
        for(int v : genotyped_indices){ if(v<mn) mn=v; if(v>mx) mx=v; }
        std::cout << "[DEBUG] Genotyped index range in H: [" << mn << ", " << mx << "]" << std::endl;
        int r0 = genotyped_indices[0];
        int c0 = genotyped_indices.size()>1 ? genotyped_indices[1] : genotyped_indices[0];
        std::cout << "[DEBUG] Sample before update Hinv(r0,r0)=" << Hinv.data(r0,r0) << ", Hinv(r0,c0)=" << Hinv.data(r0,c0) << std::endl;
    }

    // A22 is already passed in. It should be correct subset.
    std::cout << "[DEBUG] Inverting A22..." << std::endl;
    int ng2 = A22.data.rows();
    Matrix I22b = Matrix::Identity(ng2, ng2);
    Eigen::LLT<Matrix> llt22b(A22.data);
    Matrix A22inv;
    if(llt22b.info() == Eigen::Success){
        A22inv = llt22b.solve(I22b);
    } else {
        Eigen::LDLT<Matrix> ldlt22b(A22.data);
        if(ldlt22b.info() == Eigen::Success){
            A22inv = ldlt22b.solve(I22b);
        } else {
            A22inv = A22.data.inverse();
        }
    }

    std::vector<int> gIdxRows;
    for(const auto& id : genotyped_names){
        auto it = std::find(Ginv_in.rows.begin(), Ginv_in.rows.end(), id);
        if(it == Ginv_in.rows.end()) throw std::runtime_error("Genotyped ID not found in Ginv");
        gIdxRows.push_back((int)std::distance(Ginv_in.rows.begin(), it));
    }
    NamedMatrix Ginv_sorted = subset(Ginv_in, gIdxRows, gIdxRows);

    int ng = (int)genotyped_indices.size();
    double S1 = 0.0, Sg = 0.0, Sgg = 0.0, Sa = 0.0, Sga = 0.0;
    double wDiag = 10.0;
    for(int i=0;i<ng;++i){
        for(int j=0;j<ng;++j){
            double w = (i==j ? wDiag : 1.0);
            double gij = Ginv_sorted.data(i,j);
            double aij = A22inv(i,j);
            S1  += w;
            Sg  += w * gij;
            Sa  += w * aij;
            Sgg += w * gij * gij;
            Sga += w * gij * aij;
        }
    }
    double det = S1*Sgg - Sg*Sg;
    double a = 0.0, b = 1.0;
    if(std::abs(det) > 1e-12){
        a = (Sa*Sgg - Sg*Sga) / det;
        b = (S1*Sga - Sg*Sa) / det;
    }
    Matrix Ginv_aligned = (b * Ginv_sorted.data);
    Ginv_aligned.array() += a;
    double mdA = 0.0;
    double dsum = 0.0;
    for(int i=0;i<ng;++i){
        mdA += A22inv(i,i);
        dsum += std::abs(Ginv_aligned(i,i) - A22inv(i,i));
    }
    mdA /= ng;
    double mdiff = dsum / ng;
    double gamma = std::min(0.05, std::max(0.0, mdiff / (std::abs(mdA) + 1e-12)));
    Matrix diagA = Matrix::Zero(ng, ng);
    for(int i=0;i<ng;++i) diagA(i,i) = A22inv(i,i);
    Matrix Ginv_shrunk = ((1.0 - gamma) * Ginv_aligned) + (gamma * diagA);
    Eigen::SelfAdjointEigenSolver<Matrix> es2(Ginv_shrunk);
    auto ev2 = es2.eigenvalues();
    auto U2 = es2.eigenvectors();
    double me2 = ev2.maxCoeff();
    double mp2 = 0.0;
    for(int i=0;i<ev2.size();++i){ if(ev2(i) > 0.0){ mp2 = (mp2==0.0? ev2(i) : std::min(mp2, ev2(i))); } }
    if(mp2 == 0.0) mp2 = 1e-12;
    double cond2 = me2 / mp2;
    double alpha2 = (cond2 < 1e6 ? 1e-12 : (cond2 < 1e9 ? 1e-6 : 1e-4));
    double tol2 = alpha2 * me2;
    for(int i=0;i<ev2.size();++i){ if(ev2(i) < tol2) ev2(i) = tol2; }
    Matrix Ginv_bent2 = U2 * ev2.asDiagonal() * U2.transpose();
    double bestP2 = 0.0;
    double bestVal2 = std::numeric_limits<double>::infinity();
    for(double p : std::vector<double>{0.0, 0.01, 0.02, 0.03}){
        Matrix T = ((1.0 - p) * Ginv_bent2) + (p * A22inv);
        Matrix D = T - A22inv;
        double v = D.squaredNorm();
        if(v < bestVal2){ bestVal2 = v; bestP2 = p; }
    }
    Matrix Ginv_tuned2 = ((1.0 - bestP2) * Ginv_bent2) + (bestP2 * A22inv);
    Matrix H22inv_update = (tau * Ginv_tuned2) - (omega * A22inv);
    Eigen::SelfAdjointEigenSolver<Matrix> esHu(H22inv_update);
    auto evHu = esHu.eigenvalues();
    auto UHu = esHu.eigenvectors();
    double meHu = std::max(1.0, evHu.maxCoeff());
    double epsHu = 1e-12 * meHu;
    for(int i=0;i<evHu.size();++i){ if(evHu(i) < epsHu) evHu(i) = epsHu; }
    H22inv_update = UHu * evHu.asDiagonal() * UHu.transpose();

    // Update Hinv.data in place
    for(size_t i=0; i<genotyped_indices.size(); ++i){
        for(size_t j=0; j<genotyped_indices.size(); ++j){
            int r = genotyped_indices[i];
            int c = genotyped_indices[j];
            Hinv.data(r, c) += H22inv_update((int)i, (int)j);
        }
    }
    if(!genotyped_indices.empty()){
        int r0 = genotyped_indices[0];
        int c0 = genotyped_indices.size()>1 ? genotyped_indices[1] : genotyped_indices[0];
        std::cout << "[DEBUG] Sample after update Hinv(r0,r0)=" << Hinv.data(r0,r0) << ", Hinv(r0,c0)=" << Hinv.data(r0,c0) << std::endl;
    }
    
    return Hinv;
}

static int parse_ped(const std::string& path, std::vector<int>& dam, std::vector<int>& sire, std::vector<std::string>& ids){
  std::ifstream fin(path);
  if(!fin) return -1;
  std::unordered_map<std::string,int> id2idx;
  std::vector<std::vector<std::string>> rows;
  std::string line;

  // Helper to check if string is a header keyword
  auto is_header = [](const std::string& s){
      std::string lower = s;
      std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
      return lower == "id" || lower == "animal" || lower == "ind";
  };

  while(std::getline(fin,line)){
    if(line.empty()) continue;
    std::istringstream ss(line);
    std::string a,b,c;
    ss>>a>>b>>c;
    if(a.empty()) continue;
    
    // Improved header check
    if(is_header(a)) continue;

    rows.push_back({a,b,c});
    if(id2idx.find(a)==id2idx.end()){ id2idx[a]=(int)ids.size(); ids.push_back(a);} 
    if(!b.empty() && b!="NA" && b!="0" && b!="."){
      if(id2idx.find(b)==id2idx.end()){ id2idx[b]=(int)ids.size(); ids.push_back(b);} 
    }
    if(!c.empty() && c!="NA" && c!="0" && c!="."){
      if(id2idx.find(c)==id2idx.end()){ id2idx[c]=(int)ids.size(); ids.push_back(c);} 
    }
  }
  int n=(int)ids.size();
  dam.assign(n,-1);
  sire.assign(n,-1);
  for(auto& r: rows){
    int self=id2idx[r[0]];
    int d=-1,s=-1;
    auto itb=id2idx.find(r[1]); if(itb!=id2idx.end()) d=itb->second;
    auto itc=id2idx.find(r[2]); if(itc!=id2idx.end()) s=itc->second;
    dam[self]=d;
    sire[self]=s;
  }
  return n;
}

static void write_matrix_txt(const std::string& path, const Matrix& M){
    std::ofstream out(path);
    if(out.is_open()){
        out << std::fixed << std::setprecision(10);
        for(int i=0; i<M.rows(); ++i){
            for(int j=0; j<M.cols(); ++j){
                out << M(i,j) << (j==M.cols()-1?"":"\t");
            }
            out << "\n";
        }
    }
}

static bool read_matrix_txt_into(const std::string& path, Matrix& M){
    std::ifstream in(path);
    if(!in.is_open()) return false;
    std::vector<std::vector<double>> rows;
    std::string line;
    while(std::getline(in, line)){
        if(line.empty()) continue;
        std::istringstream ss(line);
        std::vector<double> vals;
        double x;
        while(ss >> x) vals.push_back(x);
        if(!vals.empty()) rows.push_back(vals);
    }
    if(rows.empty()) return false;
    int r = (int)rows.size();
    int c = (int)rows[0].size();
    for(auto& v: rows){ if((int)v.size()!=c) return false; }
    if(r!=c) return false;
    M.resize(r,c);
    for(int i=0;i<r;++i){ for(int j=0;j<c;++j){ M(i,j)=rows[i][j]; } }
    return true;
}

static bool read_ids_txt_into(const std::string& path, std::vector<std::string>& ids){
    std::ifstream in(path);
    if(!in.is_open()) return false;
    ids.clear();
    std::string line;
    bool first=true;
    while(std::getline(in, line)){
        if(line.empty()) continue;
        std::istringstream ss(line);
        std::string tok;
        if(!(ss>>tok)) continue;
        if(first && (tok=="ID" || tok=="id" || tok=="Id")){ first=false; continue; }
        first=false;
        ids.push_back(tok);
    }
    return !ids.empty();
}

extern bool g_do_inverse;



} // namespace cosmic
