#include "repeatability_helpers.h"

#include <cctype>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

namespace cosmic {

std::string sanitizeChromLabel(std::string chrom) {
    for (char& c : chrom) {
        if (!(std::isalnum(static_cast<unsigned char>(c)) || c == '_' || c == '-')) c = '_';
    }
    return chrom;
}

std::string getCacheFilePath(const std::string& prefix, const std::string& label) {
    return prefix + "." + label + ".f32grm.bin";
}

std::string getCacheMetaPath(const std::string& prefix, const std::string& label) {
    return prefix + "." + label + ".meta.txt";
}

namespace {

struct MappedFloat32Matrix {
    int fd = -1;
    size_t mapped_bytes = 0;
    void* mapped_ptr = MAP_FAILED;
    const float* data = nullptr;
    uint32_t rows = 0, cols = 0;

    ~MappedFloat32Matrix() { close(); }
    void close() {
        if (mapped_ptr != MAP_FAILED) { munmap(mapped_ptr, mapped_bytes); mapped_ptr = MAP_FAILED; }
        if (fd >= 0) { ::close(fd); fd = -1; }
        data = nullptr; mapped_bytes = 0; rows = 0; cols = 0;
    }
    bool open(const std::string& file) {
        close();
        fd = ::open(file.c_str(), O_RDONLY);
        if (fd < 0) return false;
        struct stat st;
        if (fstat(fd, &st) != 0 || st.st_size < (off_t)(sizeof(uint32_t) * 2)) { close(); return false; }
        mapped_bytes = (size_t)st.st_size;
        mapped_ptr = mmap(nullptr, mapped_bytes, PROT_READ, MAP_PRIVATE, fd, 0);
        if (mapped_ptr == MAP_FAILED) { close(); return false; }
        const uint32_t* header = reinterpret_cast<const uint32_t*>(mapped_ptr);
        rows = header[0]; cols = header[1];
        if (rows == 0 || cols == 0) { close(); return false; }
        size_t expected = sizeof(uint32_t) * 2 + sizeof(float) * (size_t)rows * (size_t)cols;
        if (mapped_bytes < expected) { close(); return false; }
        data = reinterpret_cast<const float*>(header + 2);
        return true;
    }
};

} // namespace

bool readMatrixFloat32Mmap(const std::string& file, Eigen::MatrixXf& out) {
    MappedFloat32Matrix mapped;
    if (!mapped.open(file)) return false;
    Eigen::Map<const Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>> view(mapped.data, (int)mapped.rows, (int)mapped.cols);
    out = view;
    return true;
}

void writeMatrixFloat32(const std::string& file, const Eigen::MatrixXf& mat) {
    std::filesystem::path p(file);
    if (!p.parent_path().empty()) std::filesystem::create_directories(p.parent_path());
    std::ofstream ofs(file, std::ios::binary);
    if (!ofs) throw std::runtime_error("Cannot write GRM cache file: " + file);
    uint32_t rows = (uint32_t)mat.rows(), cols = (uint32_t)mat.cols();
    ofs.write(reinterpret_cast<const char*>(&rows), sizeof(rows));
    ofs.write(reinterpret_cast<const char*>(&cols), sizeof(cols));
    ofs.write(reinterpret_cast<const char*>(mat.data()), sizeof(float) * rows * cols);
}

bool readCacheMeta(const std::string& file, double* scaling_factor, std::set<std::string>* observed_chroms) {
    std::ifstream ifs(file);
    if (!ifs) return false;
    std::string key;
    while (ifs >> key) {
        if (key == "scale" && scaling_factor) { ifs >> *scaling_factor; }
        else if (key == "chroms" && observed_chroms) {
            observed_chroms->clear();
            std::string chrom_line; std::getline(ifs, chrom_line);
            std::istringstream iss(chrom_line); std::string chrom;
            while (iss >> chrom) observed_chroms->insert(chrom);
        } else { std::string dummy; std::getline(ifs, dummy); }
    }
    return true;
}

void writeCacheMeta(const std::string& file, double scaling_factor, const std::set<std::string>& observed_chroms) {
    std::filesystem::path p(file);
    if (!p.parent_path().empty()) std::filesystem::create_directories(p.parent_path());
    std::ofstream ofs(file);
    if (!ofs) throw std::runtime_error("Cannot write GRM cache meta file: " + file);
    ofs << "scale " << scaling_factor << "\nchroms";
    for (const auto& chrom : observed_chroms) ofs << " " << chrom;
    ofs << "\n";
}

} // namespace cosmic
