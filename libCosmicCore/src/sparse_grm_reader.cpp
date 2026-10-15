#include "sparse_grm_reader.h"
#include <stdexcept>

namespace cosmic {

SparseGrmReader::Meta SparseGrmReader::readMeta(const std::string& prefix) {
    Meta meta;
    
    // Read ID file
    std::string id_file = prefix + ".grm.id";
    std::ifstream ifs_id(id_file);
    if (!ifs_id) throw std::runtime_error("Cannot open " + id_file);
    
    std::string fid, iid;
    while (ifs_id >> fid >> iid) {
        meta.fids.push_back(fid);
        meta.iids.push_back(iid);
    }
    meta.n_samples = (int)meta.fids.size();
    
    // Read Header from .sp file
    std::string sp_file = prefix + ".grm.sp";
    std::ifstream ifs_sp(sp_file, std::ios::binary);
    if (!ifs_sp) throw std::runtime_error("Cannot open " + sp_file);
    
    char magic[6] = {0};
    ifs_sp.read(magic, 5);
    if (std::string(magic) != "SPGRM") throw std::runtime_error("Invalid Sparse GRM format");
    
    int32_t version;
    ifs_sp.read((char*)&version, sizeof(int32_t));
    if (version != 1) throw std::runtime_error("Unsupported Sparse GRM version");
    
    int32_t n_check;
    ifs_sp.read((char*)&n_check, sizeof(int32_t));
    
    if (n_check != meta.n_samples) {
        // Warning or Error?
        // Usually should match.
        // std::cerr << "Warning: N in binary (" << n_check << ") != N in .id (" << meta.n_samples << ")\n";
    }
    
    ifs_sp.read((char*)&meta.threshold, sizeof(float));
    
    return meta;
}

template<typename Scalar>
void read_impl(const std::string& prefix, Eigen::SparseMatrix<Scalar>& G) {
    auto meta = SparseGrmReader::readMeta(prefix);
    int n = meta.n_samples;
    
    std::string sp_file = prefix + ".grm.sp";
    std::ifstream ifs(sp_file, std::ios::binary);
    if (!ifs) throw std::runtime_error("Cannot open " + sp_file);
    
    // Skip header
    ifs.seekg(5 + sizeof(int32_t) + sizeof(int32_t) + sizeof(float));
    
    std::vector<Eigen::Triplet<Scalar>> triplets;
    // Estimate non-zeros? Hard to guess.
    
    // Read triplets
    // Buffer for reading
    const size_t buffer_size = 1024 * 1024;
    std::vector<char> buffer(buffer_size * (sizeof(int32_t)*2 + sizeof(float)));
    
    // We read one triplet at a time or block?
    // Format is row(i32), col(i32), val(float)
    // Let's read in chunks.
    
    struct TripletRaw {
        int32_t r, c;
        float v;
    };
    
    const size_t chunk_size = 100000;
    std::vector<TripletRaw> raw_buf(chunk_size);
    
    while (ifs) {
        ifs.read((char*)raw_buf.data(), chunk_size * sizeof(TripletRaw));
        std::streamsize bytes_read = ifs.gcount();
        size_t count = bytes_read / sizeof(TripletRaw);
        
        for (size_t k = 0; k < count; ++k) {
            int r = raw_buf[k].r;
            int c = raw_buf[k].c;
            Scalar v = (Scalar)raw_buf[k].v;
            
            if (r >= n || c >= n) continue; // Safety
            
            triplets.emplace_back(r, c, v);
            if (r != c) {
                triplets.emplace_back(c, r, v);
            }
        }
    }
    
    G.resize(n, n);
    G.setFromTriplets(triplets.begin(), triplets.end());
}

void SparseGrmReader::read(const std::string& prefix, Eigen::SparseMatrix<float>& G) {
    read_impl(prefix, G);
}

void SparseGrmReader::read(const std::string& prefix, Eigen::SparseMatrix<double>& G) {
    read_impl(prefix, G);
}

} // namespace cosmic
