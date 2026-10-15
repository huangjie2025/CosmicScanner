#ifndef COSMIC_VCF_READER_H
#define COSMIC_VCF_READER_H

#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <Eigen/Dense>
#include "matrix_adapter.h"

namespace cosmic {

struct VcfVariant {
    std::string chrom;
    std::string id;
    uint32_t pos = 0;
    std::string ref;
    std::string alt;
    std::string qual;
    std::string filter;
    std::string info;
    std::string format;
    std::vector<std::string> genotypes;  // genotype strings for each sample
};

class VcfReader : public GenotypeMatrix {
public:
    VcfReader();
    ~VcfReader();

    void open(const std::string& vcf_file, const std::string& tbi_file = "");
    void close();
    void reset() const;

    size_t getNumSamples() const { return num_samples_; }
    size_t getNumVariants() const { return num_variants_; }
    const std::vector<std::string>& getSampleIds() const { return sample_ids_; }
    int rows() const override { return static_cast<int>(num_samples_); }
    int cols() const override { return static_cast<int>(num_variants_); }

    bool readVariant(std::string& chrom, std::string& id, uint32_t& pos,
                     std::string& ref, std::string& alt,
                     Eigen::VectorXf& dosage) const;

    void computeStats(int n_threads = 1);
    void multiply_Z_v(const Eigen::VectorXd& v, Eigen::VectorXd& y, int n_threads = 1) override;
    void multiply_Zt_v(const Eigen::VectorXd& v, Eigen::VectorXd& y, int n_threads = 1) override;

    const std::vector<double>& getMeans() const { return snp_means_; }
    const std::vector<double>& getInvStds() const { return snp_inv_stds_; }

    void loadVariantInfo() const;
    void seekVariant(uint32_t idx) const;
    std::pair<uint32_t, uint32_t> getChromRange(const std::string& chrom) const;
    const std::vector<std::string>& getChromLabels() const { return chrom_labels_; }
    const std::vector<VcfVariant>& getVariants() const { return variants_; }

    void releaseVariantRange(uint32_t start_idx, uint32_t count) const;

private:
    std::string vcf_file_;
    std::string tbi_file_;
    mutable std::ifstream ifs_;
    std::ifstream tbi_ifs_;
    
    uint32_t num_samples_ = 0;
    mutable uint32_t num_variants_ = 0;
    std::vector<std::string> sample_ids_;
    mutable std::vector<VcfVariant> variants_;
    mutable std::vector<std::string> chrom_labels_;
    mutable std::vector<std::pair<uint32_t, uint32_t>> chrom_ranges_;
    
    mutable size_t current_offset_ = 0;
    mutable size_t variant_data_offset_ = 0;
    size_t header_length_ = 0;
    
    std::vector<double> snp_means_;
    std::vector<double> snp_inv_stds_;
    std::vector<double> snp_mafs_;
    std::vector<int> snp_missing_counts_;
    mutable bool stats_computed_ = false;
    
    std::string parseHeader() const;
    void parseSampleLine(const std::string& line);
    bool parseVariantLine(const std::string& line, VcfVariant& variant, Eigen::VectorXf& dosage) const;
    void buildChromIndex() const;
    void computeDosageFromGT(const std::string& gt, const std::string& ref, const std::string& alt, float& dosage) const;
    bool isMissing(const std::string& gt) const;

public:
    // Export genotype data to various formats
    void exportToText(const std::string& geno_file, const std::string& var_info_file, const std::string& sample_info_file) const;
    void exportToNpy(const std::string& npy_file, const std::string& var_info_file, const std::string& sample_info_file) const;
    void exportToPedMap(const std::string& ped_file, const std::string& map_file) const;
    void exportToVcf(const std::string& vcf_file) const;
};

} // namespace cosmic

#endif // COSMIC_VCF_READER_H