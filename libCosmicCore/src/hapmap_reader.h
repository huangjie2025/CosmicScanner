#pragma once
#include <string>
#include <vector>
#include <Eigen/Dense>
#include "matrix_adapter.h"

namespace cosmic {

struct HapMapVariant {
    std::string id;
    std::string chrom;
    uint32_t pos = 0;
    std::string allele1;
    std::string allele2;
    std::vector<uint8_t> genotypes; // 0=hom allele1, 1=het, 2=hom allele2, 3=missing
};

class HapMapReader : public SequentialGenotypeMatrix {
public:
    HapMapReader() = default;
    ~HapMapReader() = default;

    // Standard HapMap/HMP: 11 metadata columns followed by one column per sample.
    void open(const std::string& hmp_file);
    void close();

    uint32_t getNumSamples() const override { return num_samples_; }
    uint32_t getNumVariants() const override { return num_variants_; }
    const std::vector<std::string>& getSampleIds() const override { return sample_ids_; }
    int rows() const override { return static_cast<int>(num_samples_); }
    int cols() const override { return static_cast<int>(num_variants_); }

    bool readVariant(std::string& chrom, std::string& id, uint32_t& pos,
                     std::string& ref, std::string& alt,
                     Eigen::VectorXf& dosage) const override;

    void computeStats(int n_threads = 1);
    void multiply_Z_v(const Eigen::VectorXd& v, Eigen::VectorXd& out, int n_threads = 1) override;
    void multiply_Zt_v(const Eigen::VectorXd& v, Eigen::VectorXd& out, int n_threads = 1) override;

    const std::vector<double>& getMeans() const { return snp_means_; }
    const std::vector<double>& getInvStds() const { return snp_inv_stds_; }

    void seekVariant(uint32_t idx) const;
    void reset() const override { current_variant_idx_ = 0; }
    std::pair<uint32_t, uint32_t> getChromRange(const std::string& chrom) const;
    const std::vector<std::string>& getChromLabels() const { return chrom_labels_; }
    const std::vector<HapMapVariant>& getVariants() const { return variants_; }

public:
    void exportToText(const std::string& geno_file, const std::string& var_info_file, const std::string& sample_info_file) const;
    void exportToNpy(const std::string& npy_file, const std::string& var_info_file, const std::string& sample_info_file) const;
    void exportToPedMap(const std::string& ped_file, const std::string& map_file) const;
    void exportToVcf(const std::string& vcf_file) const;

private:
    std::string hmp_file_;
    
    uint32_t num_samples_ = 0;
    uint32_t num_variants_ = 0;
    std::vector<std::string> sample_ids_;
    std::vector<HapMapVariant> variants_;
    std::vector<std::string> chrom_labels_;
    std::vector<std::pair<uint32_t, uint32_t>> chrom_ranges_;
    
    std::vector<double> snp_means_;
    std::vector<double> snp_inv_stds_;
    bool stats_computed_ = false;
    
    mutable size_t current_variant_idx_ = 0;
};

} // namespace cosmic
