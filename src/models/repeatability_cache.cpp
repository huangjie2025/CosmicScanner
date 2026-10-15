#include "repeatability_helpers.h"

#include <cctype>
#include <cstdint>
#include <cmath>
#include <iomanip>
#include <limits>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include "mmap_matrix.h"

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
    MMapFile mapping;
    const float* data = nullptr;
    uint32_t rows = 0, cols = 0;

    void close() {
        mapping.close();
        data = nullptr; rows = 0; cols = 0;
    }
    bool open(const std::string& file) {
        close();
        try { mapping.open(file); }
        catch (const std::runtime_error&) { close(); return false; }
        const size_t header_bytes = sizeof(uint32_t) * 2;
        if (!mapping.get_data() || mapping.get_size() < header_bytes) {
            close(); return false;
        }
        const auto* header = static_cast<const uint32_t*>(mapping.get_data());
        rows = header[0]; cols = header[1];
        if (!rows || !cols || rows > static_cast<uint32_t>(std::numeric_limits<int>::max()) ||
            cols > static_cast<uint32_t>(std::numeric_limits<int>::max()) ||
            static_cast<size_t>(rows) >
                (std::numeric_limits<size_t>::max() - header_bytes) / sizeof(float) / cols) {
            close(); return false;
        }
        const size_t expected = header_bytes + sizeof(float) * static_cast<size_t>(rows) * cols;
        if (mapping.get_size() < expected) { close(); return false; }
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

bool readCacheMeta(const std::string& file, double* scaling_factor, std::set<std::string>* observed_chroms,
                   const std::string& expected_identity) {
    std::ifstream ifs(file);
    if (!ifs) return false;
    bool have_scale = false, have_chroms = false;
    int version = 0;
    std::string identity;
    double scale = 0.0;
    std::set<std::string> chroms;
    std::string key;
    while (ifs >> key) {
        if (key == "version") { if (!(ifs >> version)) return false; }
        else if (key == "identity") { if (!(ifs >> identity)) return false; }
        else if (key == "scale") {
            if (!(ifs >> scale) || !std::isfinite(scale) || scale < 0.0) return false;
            have_scale = true;
        } else if (key == "chroms") {
            std::string chrom_line; std::getline(ifs, chrom_line);
            std::istringstream iss(chrom_line); std::string chrom;
            while (iss >> chrom) chroms.insert(chrom);
            have_chroms = true;
        } else { std::string dummy; std::getline(ifs, dummy); }
    }
    // Legacy caches lack an input identity and must be recomputed once.
    if (ifs.bad() || version != 2 || expected_identity.empty() ||
        identity != expected_identity || !have_scale || !have_chroms) return false;
    if (scaling_factor) *scaling_factor = scale;
    if (observed_chroms) *observed_chroms = std::move(chroms);
    return true;
}

void writeCacheMeta(const std::string& file, double scaling_factor, const std::set<std::string>& observed_chroms,
                    const std::string& identity) {
    std::filesystem::path p(file);
    if (!p.parent_path().empty()) std::filesystem::create_directories(p.parent_path());
    std::ofstream ofs(file);
    if (!ofs) throw std::runtime_error("Cannot write GRM cache meta file: " + file);
    ofs << std::setprecision(std::numeric_limits<double>::max_digits10)
        << "version 2\nidentity " << identity << "\nscale " << scaling_factor << "\nchroms";
    for (const auto& chrom : observed_chroms) ofs << " " << chrom;
    ofs << "\n";
    ofs.flush();
    if (!ofs) throw std::runtime_error("Failed writing GRM cache meta file: " + file);
}

} // namespace cosmic
