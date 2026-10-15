#pragma once
#include <string>
#include <vector>
#include <Eigen/Dense>
#include "matrix_adapter.h"

namespace cosmic {

struct PedSample {
    std::string fid;
    std::string iid;
    std::string father;
    std::string mother;
    int sex = 0;
    float phenotype = -9.0f;
};

struct MapVariant {
    std::string chrom;
    std::string id;
    double cm_pos = 0.0;
    uint32_t bp_pos = 0;
    std::string allele1;
    std::string allele2;
};

class PedMapReader : public SequentialGenotypeMatrix {
public:
    PedMapReader() = default;
    ~PedMapReader() = default;

    // Open PED/MAP files
    // ped_file: PLINK .ped format (FID IID PAT MAT SEX PHENO genotype...)
    // map_file: PLINK .map format (CHR SNP CM BP A1 A2)
    void open(const std::string& ped_file, const std::string& map_file);
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
    const std::vector<PedSample>& getSamples() const { return samples_; }
    const std::vector<MapVariant>& getVariants() const { return variants_; }

    // Export to PED/MAP format
    void exportToPedMap(const std::string& ped_file, const std::string& map_file) const;

private:
    std::string ped_file_;
    std::string map_file_;
    
    uint32_t num_samples_ = 0;
    uint32_t num_variants_ = 0;
    std::vector<PedSample> samples_;
    std::vector<MapVariant> variants_;
    std::vector<std::string> sample_ids_;
    std::vector<std::string> chrom_labels_;
    std::vector<std::pair<uint32_t, uint32_t>> chrom_ranges_;
    
    // Genotype data: [variant][sample] = {allele1, allele2} encoded as 0,1,2,3
    // 0 = hom allele1 (A1/A1), 1 = het (A1/A2), 2 = hom allele2 (A2/A2), 3 = missing
    std::vector<std::vector<uint8_t>> genotypes_;
    
    std::vector<double> snp_means_;
    std::vector<double> snp_inv_stds_;
    bool stats_computed_ = false;
    
    mutable size_t current_variant_idx_ = 0;
    
    void parsePedFile();
    void parseMapFile();
    void buildChromIndex();
    uint8_t encodeGenotype(const std::string& a1, const std::string& a2, const std::string& ref, const std::string& alt) const;
    void decodeGenotype(uint8_t code, const std::string& ref, const std::string& alt, std::string& a1, std::string& a2) const;
};

} // namespace cosmic