#include "kinship_operator.h"
#include <filesystem>
#include <fstream>
#include <iostream>
#include <random>
#include <numeric>
#include <stdexcept>

void require(bool value, const char* message) {
    if (!value) throw std::runtime_error(message);
}
void compare(const Eigen::MatrixXd& actual, const Eigen::MatrixXd& expected) {
    const double rel = (actual-expected).cwiseAbs().maxCoeff() /
        std::max(1e-30, expected.cwiseAbs().maxCoeff());
    const double scaled = ((actual-expected).array().abs() /
                          (1.0+expected.array().abs())).maxCoeff();
    std::cout << "max_abs_difference=" << (actual-expected).cwiseAbs().maxCoeff()
              << " max_relative_error=" << rel << " max_scaled_error=" << scaled << '\n';
    require(rel < 1e-10 && scaled < 1e-10, "batch equivalence failed");
}
int main(int argc, char** argv) {
    try {
        require(argc == 2, "expected scratch directory");
        std::filesystem::create_directories(argv[1]);
        const std::string prefix = std::string(argv[1])+"/fixture";
        const int n=37, m=1031; // partial BED byte and multiple partial SNP blocks
        std::mt19937_64 rng(8128);
        {
            std::ofstream bed(prefix+".bed", std::ios::binary), bim(prefix+".bim"), fam(prefix+".fam");
            const unsigned char magic[] = {0x6c,0x1b,0x01};
            bed.write(reinterpret_cast<const char*>(magic),3);
            for (int i=0; i<n; ++i) fam << "f i" << i << " 0 0 0 1\n";
            for (int j=0; j<m; ++j) {
                bim << 1+j%3 << " s" << j << " 0 " << j+1 << " A C\n";
                for (int b=0; b<(n+3)/4; ++b) {
                    unsigned char byte=0;
                    for (int k=0; k<4; ++k) {
                        unsigned char code = j%31 == 0 ? 3 : j%37 == 0 ? 1 : rng()%4;
                        byte |= code << (2*k);
                    }
                    bed.put(static_cast<char>(byte));
                }
            }
        }
        cosmic::PlinkReader reader;
        reader.load(prefix);
        std::vector<int> subset{35,1,20,3,4,31,7,8,10,12,14,16,18,22,24,26,28,30,32,34};
        std::vector<int> identity(n), permuted(n), repeated(n);
        std::iota(identity.begin(), identity.end(), 0);
        for (int i=0; i<n; ++i) {
            permuted[i] = n-1-i;
            repeated[i] = i/2; // explicit repeated rows, not an identity selection
        }
        std::vector<int> variants;
        for (int j=0; j<m; ++j) if (j%3 != 1) variants.push_back(j); // LOCO-like mask
        for (int selection=0; selection<6; ++selection) {
            const std::vector<int> keep = selection==1 ? subset : selection==3 ? identity :
                selection==4 ? permuted : selection==5 ? repeated : std::vector<int>{};
            const std::vector<int> snps = selection==0 ? std::vector<int>{} : variants;
            cosmic::PlinkKinshipOperator op(&reader, keep, snps);
            Eigen::MatrixXd z(op.rows(), op.numVariants());
            for (int j=0; j<z.cols(); ++j) {
                const int source = snps.empty() ? j : snps[j];
                Eigen::VectorXf raw;
                reader.getSnpGenotypes(source, raw);
                for (int i=0; i<z.rows(); ++i) {
                    const float dosage=raw[keep.empty()?i:keep[i]];
                    z(i,j)=std::isnan(dosage)?0.0:
                        (dosage-reader.getMeans()[source])*reader.getInvStds()[source];
                }
            }
            const Eigen::MatrixXd kernel = z*z.transpose()/op.numVariants();
            require(std::abs(op.trace()-kernel.trace()) < 1e-10, "mapped trace mismatch");
            for (int width : {1,7,41}) {
                std::normal_distribution<double> normal;
                Eigen::MatrixXd omega(op.rows(),width), old(op.rows(),width);
                for (int j=0; j<width; ++j) {
                    for (int i=0; i<op.rows(); ++i) omega(i,j)=normal(rng);
                    Eigen::VectorXd v=omega.col(j), result;
                    op.apply(v,result,1);
                    old.col(j)=result;
                    for (int threads : {2,4,8}) {
                        Eigen::VectorXd parallel, transpose, alias=v;
                        op.apply(v,parallel,threads);
                        compare(parallel,result);
                        compare(parallel,kernel*v);
                        op.applyTranspose(v,transpose,threads);
                        compare(transpose,parallel);
                        op.apply(alias,alias,threads);
                        compare(alias,parallel);
                        Eigen::VectorXd again;
                        op.apply(v,again,threads);
                        require((again.array()==parallel.array()).all(),
                                "fixed-team vector repeat is not deterministic");
                    }
                }
                Eigen::MatrixXd oracle = z*(z.transpose()*omega)/op.numVariants();
                compare(old,oracle);
                for (int threads : {1,4}) {
                    Eigen::MatrixXd result;
                    cosmic::PlinkKinshipOperator::BlockProfile p;
                    op.apply_block(omega,result,threads,&p);
                    compare(result,old);
                    require(p.blocks == (op.numVariants()+511)/512, "wrong block count");
                    require(p.packed_bytes == size_t(op.numVariants())*reader.getBytesPerSnp(),
                            "more than one packed pass");
                    std::uint64_t active=0;
                    for (int j=0; j<op.numVariants(); ++j)
                        if (reader.getInvStds()[snps.empty()?j:snps[j]] != 0) active+=op.rows();
                    require(p.decoded_genotypes == active, "SNP decoded more than once");
                    Eigen::MatrixXd alias=omega;
                    op.apply(alias,alias,threads);
                    compare(alias,old);
                }
            }
            Eigen::MatrixXd empty(op.rows(),0), out;
            op.apply(empty,out,1);
            require(out.rows()==op.rows() && out.cols()==0, "empty RHS shape mismatch");
            bool rejected=false;
            try { Eigen::MatrixXd wrong(op.rows()+1,1); op.apply(wrong,out,1); }
            catch (const std::runtime_error&) { rejected=true; }
            require(rejected,"bad dimensions not rejected");
        }
        std::cout << "BATCH_KOMEGA_PASS\n";
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n'; return 1;
    }
}
