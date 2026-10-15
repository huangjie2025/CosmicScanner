#include "association/cojo/cojo_engine.h"
#include "association/sumstats/sumstats.h"
#include "plink_reader.h"
#include "bgen_reader.h"
#include "pgen_reader.h"
#include "vcf/vcf_reader.h"
#include "genetics/ld/ld_calculator.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <unordered_set>

namespace cosmic::association::cojo {

namespace {

double normal_cdf(double x) {
    return 0.5 * std::erfc(-x * M_SQRT1_2);
}

// Parse summary line from canonical .cosmic.sumstats format
SummaryRecord parseSummaryLine(const std::string& line) {
    std::istringstream iss(line);
    SummaryRecord s;
    std::string dummy;
    iss >> s.snp_id >> s.chrom >> s.position >> s.a1 >> s.a2
        >> s.beta >> s.se >> s.p_value >> s.n >> s.freq;
    if (s.se > 0.0) s.z_score = s.beta / s.se;
    if (s.se > 0.0) s.weight = 1.0 / (s.se * s.se);
    return s;
}

// Parse from multiple possible formats
SummaryRecord parseSummaryLineFlexible(const std::string& line) {
    std::istringstream iss(line);
    std::vector<std::string> fields;
    std::string field;
    while (iss >> field) fields.push_back(field);
    
    if (fields.size() < 9) return SummaryRecord{};
    
    SummaryRecord s;
    // Detect format: canonical sumstats has CHR (numeric) first, old format has SNP_ID first
    bool is_canonical = false;
    if (!fields[0].empty() && std::isdigit(fields[0][0])) {
        // First field looks like a chromosome number - canonical format
        // CHR, SNP, POS, A1, A2, BETA, SE, P, N, AF, ...
        is_canonical = true;
    }
    
    if (is_canonical) {
        // Canonical: CHR, SNP, POS, A1, A2, BETA, SE, P, N, AF, ...
        s.chrom = fields[0];
        s.snp_id = fields[1];
        s.position = std::stoll(fields[2]);
        s.a1 = fields[3];
        s.a2 = fields[4];
        s.beta = std::stod(fields[5]);
        s.se = std::stod(fields[6]);
        s.p_value = std::stod(fields[7]);
        s.n = std::stoi(fields[8]);
        if (fields.size() > 9) s.freq = std::stod(fields[9]);
    } else {
        // Old format: SNP, CHR, POS, A1, A2, BETA, SE, P, N, FREQ
        s.snp_id = fields[0];
        s.chrom = fields[1];
        s.position = std::stoll(fields[2]);
        s.a1 = fields[3];
        s.a2 = fields[4];
        s.beta = std::stod(fields[5]);
        s.se = std::stod(fields[6]);
        s.p_value = std::stod(fields[7]);
        s.n = std::stoi(fields[8]);
        if (fields.size() > 9) s.freq = std::stod(fields[9]);
    }
    if (s.se > 0.0) s.z_score = s.beta / s.se;
    if (s.se > 0.0) s.weight = 1.0 / (s.se * s.se);
    return s;
}

std::vector<std::string> readIdList(const std::string& filepath) {
    std::ifstream file(filepath);
    std::vector<std::string> ids;
    std::string line;
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream iss(line);
        std::string id;
        iss >> id;
        if (!id.empty()) ids.push_back(id);
    }
    return ids;
}

void normalizeDosage(Eigen::MatrixXd& mat) {
    for (int j = 0; j < mat.cols(); ++j) {
        Eigen::VectorXd col = mat.col(j);
        double mean = col.mean();
        double std = col.array().square().mean();
        std = std::sqrt(std - mean * mean);
        if (std > 1e-8) {
            mat.col(j) = (col.array() - mean) / std;
        }
    }
}

} // namespace

CoJoEngine::CoJoEngine(const CoJoOptions& opts) : options_(opts) {}

void CoJoEngine::loadSummaryStats(const std::string& filepath) {
    summary_stats_.clear();
    summary_idx_.clear();
    if (!options_.canonical_summary_only) {
        std::ifstream file(filepath);
        if (!file) throw std::runtime_error("Cannot open summary stats file: " + filepath);
        std::string line;
        while (std::getline(file, line)) {
            if (line.empty() || line[0] == '#') continue;
            std::istringstream header(line);
            std::string first_field;
            header >> first_field;
            if (first_field == "CHR" || first_field == "chr" || first_field == "Chromosome" ||
                first_field == "SNP" || first_field == "rsid") continue;
            SummaryRecord s = parseSummaryLineFlexible(line);
            if (s.se > 0.0 && std::isfinite(s.beta) && std::isfinite(s.se) && s.n >= 1 &&
                s.freq >= options_.maf_filter && s.freq <= 1.0 - options_.maf_filter)
                summary_stats_.push_back(std::move(s));
        }
        for (size_t i = 0; i < summary_stats_.size(); ++i)
            summary_idx_[summary_stats_[i].snp_id] = static_cast<int>(i);
        std::cout << "Loaded " << summary_stats_.size() << " summary statistics records\n";
        return;
    }

    using cosmic::association::sumstats::SumstatsReader;
    using cosmic::association::sumstats::SumstatsRecord;
    SumstatsReader reader(filepath);
    SumstatsRecord record;
    while (reader.readRecord(record)) {
        if (record.se <= 0.0 || !std::isfinite(record.beta) || !std::isfinite(record.se) || record.n < 1)
            continue;
        if (record.freq < options_.maf_filter || record.freq > 1.0 - options_.maf_filter) continue;
        if (!options_.region_chrom.empty() && record.chrom != options_.region_chrom) continue;
        if (options_.region_start >= 0 && record.pos < options_.region_start) continue;
        if (options_.region_end >= 0 && record.pos > options_.region_end) continue;
        SummaryRecord s;
        s.snp_id = record.snp_id;
        s.chrom = record.chrom;
        s.position = record.pos;
        s.a1 = record.a1;
        s.a2 = record.a2;
        s.beta = record.beta;
        s.se = record.se;
        s.p_value = record.p_value;
        s.freq = record.freq;
        s.n = static_cast<int>(record.n);
        s.z_score = record.beta / record.se;
        s.weight = 1.0 / (record.se * record.se);
        summary_stats_.push_back(std::move(s));
    }

    for (size_t i = 0; i < summary_stats_.size(); ++i) {
        summary_idx_[summary_stats_[i].snp_id] = static_cast<int>(i);
    }

    std::cout << "Loaded " << summary_stats_.size() << " summary statistics records\n";
}
void CoJoEngine::loadLDReference(const std::string& filepath) {
    ld_ref_.snp_ids.clear();
    ld_ref_.chroms.clear();
    ld_ref_.positions.clear();
    ld_ref_.a1s.clear();
    ld_ref_.a2s.clear();
    ld_ref_.snp_to_idx.clear();
    ld_ref_.ld_matrix.resize(0, 0);

    std::ifstream file(filepath);
    if (!file) throw std::runtime_error("Cannot open LD reference file: " + filepath);

    std::string line;
    std::vector<std::string> header;
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream input(line);
        std::string field;
        while (input >> field) header.push_back(field);
        break;
    }
    if (header.size() < 6 || header[0] != "SNP" || header[1] != "CHR" ||
        header[2] != "BP" || header[3] != "A1" || header[4] != "A2")
        throw std::runtime_error("External LD header must be: SNP CHR BP A1 A2 <ordered SNP IDs>");

    const int m = static_cast<int>(header.size()) - 5;
    ld_ref_.snp_ids.assign(header.begin() + 5, header.end());
    std::unordered_set<std::string> seen_columns;
    for (int i = 0; i < m; ++i) {
        if (ld_ref_.snp_ids[i].empty() || !seen_columns.insert(ld_ref_.snp_ids[i]).second)
            throw std::runtime_error("External LD matrix contains duplicate SNP columns");
        ld_ref_.snp_to_idx[ld_ref_.snp_ids[i]] = i;
    }
    ld_ref_.ld_matrix = Eigen::MatrixXd::Zero(m, m);

    int row = 0;
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream input(line);
        std::vector<std::string> fields;
        std::string field;
        while (input >> field) fields.push_back(field);
        if (fields.size() != static_cast<size_t>(m + 5))
            throw std::runtime_error("External LD matrix row has the wrong number of fields");
        if (row >= m) throw std::runtime_error("External LD matrix has too many rows");
        if (fields[0] != ld_ref_.snp_ids[row])
            throw std::runtime_error("External LD row order does not match SNP column order");
        ld_ref_.chroms.push_back(fields[1]);
        try {
            ld_ref_.positions.push_back(std::stoll(fields[2]));
        } catch (...) {
            throw std::runtime_error("External LD matrix contains an invalid BP value");
        }
        ld_ref_.a1s.push_back(fields[3]);
        ld_ref_.a2s.push_back(fields[4]);
        for (int column = 0; column < m; ++column) {
            try {
                ld_ref_.ld_matrix(row, column) = std::stod(fields[column + 5]);
            } catch (...) {
                throw std::runtime_error("External LD matrix contains an invalid correlation value");
            }
        }
        ++row;
    }
    if (row != m) throw std::runtime_error("External LD matrix is not square");
    validateLDReference();
    std::cout << "LD reference loaded: " << m << " SNPs with a validated signed-correlation matrix\n";
}

void CoJoEngine::loadLDReferenceFromGenotype(const std::string& summary_stats_file,
                                             PlinkReader* plink_reader,
                                             BgenReader* bgen_reader,
                                             PgenReader* pgen_reader,
                                             VcfReader* vcf_reader,
                                             bool is_bgen, bool is_pgen, bool is_vcf) {
    // First load summary stats to get the SNP list
    loadSummaryStats(summary_stats_file);
    
    // Build set of SNP IDs from summary stats
    std::unordered_set<std::string> target_snps;
    for (const auto& s : summary_stats_) {
        target_snps.insert(s.snp_id);
    }
    
    if (target_snps.empty()) {
        throw std::runtime_error("No SNPs in summary statistics to build LD reference");
    }
    
    std::cout << "Building LD reference from genotype data for " << target_snps.size() << " target SNPs...\n";
    
    // Collect genotype data for target SNPs
    std::vector<std::string> snp_ids;
    std::vector<std::string> chroms;
    std::vector<long long> positions;
    std::vector<std::string> a1s;
    std::vector<std::string> a2s;
    
    int n_samples = 0;
    if (is_bgen && bgen_reader) {
        n_samples = bgen_reader->getNumSamples();
    } else if (is_pgen && pgen_reader) {
        n_samples = pgen_reader->getNumSamples();
    } else if (is_vcf && vcf_reader) {
        n_samples = vcf_reader->getNumSamples();
    } else if (plink_reader) {
        n_samples = plink_reader->getNumSamples();
    } else {
        throw std::runtime_error("No genotype reader available for LD computation");
    }
    
    if (n_samples < 2) {
        throw std::runtime_error("Insufficient samples for LD computation");
    }
    
    std::vector<Eigen::VectorXd> dosages_list;
    int valid_snp_count = 0;
    
    std::string chrom, rsid, ref, alt;
    uint32_t pos = 0;
    Eigen::VectorXf dosage_f;
    
    auto process_variant = [&](const std::string& chrom, const std::string& rsid,
                               uint32_t pos, const std::string& ref, const std::string& alt,
                               const Eigen::VectorXf& dosage_f) -> bool {
        if (!target_snps.count(rsid)) return false;
        
        // Convert to double and check for valid dosages
        Eigen::VectorXd dosage(n_samples);
        int obs = 0;
        double sum = 0.0;
        for (int i = 0; i < n_samples; ++i) {
            float val = dosage_f[i];
            if (std::isfinite(val) && val >= 0.0 && val <= 2.0) {
                dosage[i] = val;
                sum += val;
                obs++;
            } else {
                dosage[i] = std::numeric_limits<double>::quiet_NaN();
            }
        }
        
        if (obs < 3) return false;
        double freq = sum / (2.0 * obs);
        double maf = std::min(freq, 1.0 - freq);
        if (maf < options_.maf_filter) return false;
        
        // Check for monomorphic
        bool all_same = true;
        double first_val = std::numeric_limits<double>::quiet_NaN();
        for (int i = 0; i < n_samples; ++i) {
            if (std::isfinite(dosage[i])) {
                if (std::isnan(first_val)) first_val = dosage[i];
                else if (dosage[i] != first_val) { all_same = false; break; }
            }
        }
        if (all_same) return false;
        
        snp_ids.push_back(rsid);
        chroms.push_back(chrom);
        positions.push_back(pos);
        a1s.push_back(alt);
        a2s.push_back(ref);
        dosages_list.push_back(std::move(dosage));
        valid_snp_count++;
        return true;
    };
    
    if (is_bgen && bgen_reader) {
        bgen_reader->reset();
        while (bgen_reader->readVariant(chrom, rsid, pos, ref, alt, dosage_f)) {
            process_variant(chrom, rsid, pos, ref, alt, dosage_f);
        }
    } else if (is_pgen && pgen_reader) {
        pgen_reader->reset();
        while (pgen_reader->readVariant(chrom, rsid, pos, ref, alt, dosage_f)) {
            process_variant(chrom, rsid, pos, ref, alt, dosage_f);
        }
    } else if (is_vcf && vcf_reader) {
        vcf_reader->reset();
        while (vcf_reader->readVariant(chrom, rsid, pos, ref, alt, dosage_f)) {
            process_variant(chrom, rsid, pos, ref, alt, dosage_f);
        }
    } else if (plink_reader) {
        // PLINK reader uses a different interface - need to iterate through all SNPs
        size_t m = plink_reader->getNumSnps();
        const auto& bim = plink_reader->getBimInfo();
        for (size_t j = 0; j < m; ++j) {
            if (!target_snps.count(bim[j].id)) continue;
            
            Eigen::VectorXd dosage(n_samples);
            plink_reader->readSnp(j, dosage);
            
            const std::string& rsid = bim[j].id;
            const std::string& chrom = bim[j].chrom;
            uint32_t pos = bim[j].bp_pos;
            const std::string& ref = bim[j].ref;
            const std::string& alt = bim[j].alt;
            
            // Convert to float vector for process_variant
            Eigen::VectorXf dosage_f = dosage.cast<float>();
            if (!process_variant(chrom, rsid, pos, ref, alt, dosage_f)) continue;
        }
    }
    
    if (valid_snp_count == 0) {
        throw std::runtime_error("No valid SNPs found in genotype data matching summary statistics");
    }
    
    std::cout << "Found " << valid_snp_count << " valid SNPs for LD computation\n";
    
    // Build dosage matrix (samples x variants)
    Eigen::MatrixXd dosage_matrix(n_samples, valid_snp_count);
    for (int j = 0; j < valid_snp_count; ++j) {
        dosage_matrix.col(j) = dosages_list[j];
    }
    
    // Compute LD matrix using existing LDCalculator
    std::cout << "Computing LD matrix (" << valid_snp_count << " x " << valid_snp_count << ")...\n";
    Eigen::MatrixXd ld_matrix = cosmic::genetics::LDCalculator::computePairwiseLD(dosage_matrix);
    
    std::vector<int> kept_indices(valid_snp_count);
    std::iota(kept_indices.begin(), kept_indices.end(), 0);
    
    // Build final LD reference with kept SNPs
    ld_ref_.snp_ids.clear();
    ld_ref_.chroms.clear();
    ld_ref_.positions.clear();
    ld_ref_.a1s.clear();
    ld_ref_.a2s.clear();
    ld_ref_.snp_to_idx.clear();
    
    int kept_count = kept_indices.size();
    Eigen::MatrixXd final_ld = Eigen::MatrixXd::Zero(kept_count, kept_count);
    
    for (int i = 0; i < kept_count; ++i) {
        int orig_i = kept_indices[i];
        ld_ref_.snp_ids.push_back(snp_ids[orig_i]);
        ld_ref_.chroms.push_back(chroms[orig_i]);
        ld_ref_.positions.push_back(positions[orig_i]);
        ld_ref_.a1s.push_back(a1s[orig_i]);
        ld_ref_.a2s.push_back(a2s[orig_i]);
        ld_ref_.snp_to_idx[snp_ids[orig_i]] = i;
        for (int j = 0; j < kept_count; ++j) {
            int orig_j = kept_indices[j];
            final_ld(i, j) = ld_matrix(orig_i, orig_j);
        }
    }
    
    ld_ref_.ld_matrix = std::move(final_ld);
    
    std::cout << "LD reference built: " << kept_count << " SNPs\n";
    std::cout << "LD matrix condition number check: ";
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(ld_ref_.ld_matrix);
    double min_eig = eig.eigenvalues().minCoeff();
    double max_eig = eig.eigenvalues().maxCoeff();
    if (max_eig > 0) {
        std::cout << "cond = " << (max_eig / std::max(min_eig, 1e-12)) << "\n";
    } else {
        std::cout << "degenerate\n";
    }
    validateLDReference();
}

void CoJoEngine::loadConditionSnps(const std::string& filepath) {
    condition_snps_ = readIdList(filepath);
    std::cout << "Loaded " << condition_snps_.size() << " condition SNPs\n";
}

void CoJoEngine::alignAlleles() {
    using namespace cosmic::association::sumstats;
    
    // Build reference allele map from LD reference
    std::unordered_map<std::string, std::pair<std::string, std::string>> ref_alleles;
    for (size_t i = 0; i < ld_ref_.snp_ids.size(); ++i) {
        ref_alleles[ld_ref_.snp_ids[i]] = {ld_ref_.a1s[i], ld_ref_.a2s[i]};
    }
    
    // Convert summary stats to sumstats records for harmonization
    std::vector<SumstatsRecord> sumstats_records;
    sumstats_records.reserve(summary_stats_.size());
    for (const auto& s : summary_stats_) {
        SumstatsRecord r;
        r.snp_id = s.snp_id;
        r.chrom = s.chrom;
        r.pos = s.position;
        r.a1 = s.a1;
        r.a2 = s.a2;
        r.beta = s.beta;
        r.se = s.se;
        r.p_value = s.p_value;
        r.n = s.n;
        r.freq = s.freq;
        sumstats_records.push_back(r);
    }
    
    // Harmonize alleles
    auto harmonized = harmonizeSumstats(sumstats_records, ref_alleles);
    
    std::vector<SummaryRecord> aligned;
    aligned.reserve(summary_stats_.size());
    
    size_t harmonized_count = 0;
    size_t mismatched_count = 0;
    size_t ambiguous_count = 0;
    std::unordered_set<std::string> aligned_ids;
    
    for (const auto& h : harmonized) {
        if (h.harmonized) {
            if (!aligned_ids.insert(h.record.snp_id).second)
                throw std::runtime_error("Duplicate SNP in COJO summary statistics: " + h.record.snp_id);
            SummaryRecord s;
            s.snp_id = h.record.snp_id;
            s.chrom = h.record.chrom;
            s.position = h.record.pos;
            s.a1 = h.record.a1;
            s.a2 = h.record.a2;
            s.beta = h.record.beta;
            s.se = h.record.se;
            s.p_value = h.record.p_value;
            s.n = h.record.n;
            s.freq = h.record.freq;
            s.z_score = s.beta / s.se;
            s.weight = 1.0 / (s.se * s.se);
            aligned.push_back(s);
            harmonized_count++;
        } else {
            if (h.status == AlleleMatchStatus::AMBIGUOUS) {
                ambiguous_count++;
            } else {
                mismatched_count++;
            }
            std::cerr << "Warning: " << h.reason << " for " << h.record.snp_id << "\n";
        }
    }

    // Conditional routines index directly into the LD matrix, so preserve its
    // genotype-reference order after dropping unmatched summary records.
    std::stable_sort(aligned.begin(), aligned.end(), [&](const SummaryRecord& left, const SummaryRecord& right) {
        return ld_ref_.snp_to_idx.at(left.snp_id) < ld_ref_.snp_to_idx.at(right.snp_id);
    });

    summary_stats_ = std::move(aligned);
    summary_idx_.clear();
    for (size_t i = 0; i < summary_stats_.size(); ++i) {
        summary_idx_[summary_stats_[i].snp_id] = static_cast<int>(i);
    }
    
    std::cout << "After allele alignment: " << summary_stats_.size() << " SNPs "
              << "(harmonized: " << harmonized_count 
              << ", mismatched: " << mismatched_count 
              << ", ambiguous: " << ambiguous_count << ")\n";
}

void CoJoEngine::validateLDReference() const {
    const Eigen::Index m = static_cast<Eigen::Index>(ld_ref_.snp_ids.size());
    if (m == 0 || ld_ref_.ld_matrix.rows() != m || ld_ref_.ld_matrix.cols() != m)
        throw std::runtime_error("LD reference must contain a non-empty square matrix matching its SNP list");
    if (ld_ref_.chroms.size() != static_cast<size_t>(m) ||
        ld_ref_.positions.size() != static_cast<size_t>(m) ||
        ld_ref_.a1s.size() != static_cast<size_t>(m) ||
        ld_ref_.a2s.size() != static_cast<size_t>(m))
        throw std::runtime_error("LD reference metadata dimensions do not match its SNP list");
    if (!ld_ref_.ld_matrix.allFinite())
        throw std::runtime_error("LD reference contains non-finite correlation values");
    if ((ld_ref_.ld_matrix.array().abs() > 1.0 + 1e-8).any())
        throw std::runtime_error("LD reference correlation lies outside [-1,1]");
    if ((ld_ref_.ld_matrix - ld_ref_.ld_matrix.transpose()).cwiseAbs().maxCoeff() > 1e-8)
        throw std::runtime_error("LD reference matrix is not symmetric");
    if ((ld_ref_.ld_matrix.diagonal().array() - 1.0).abs().maxCoeff() > 1e-6)
        throw std::runtime_error("LD reference diagonal must be approximately one");
    std::unordered_set<std::string> seen;
    for (int i = 0; i < m; ++i) {
        if (!seen.insert(ld_ref_.snp_ids[i]).second)
            throw std::runtime_error("LD reference contains duplicate SNP IDs");
        const auto found = ld_ref_.snp_to_idx.find(ld_ref_.snp_ids[i]);
        if (found == ld_ref_.snp_to_idx.end() || found->second != i)
            throw std::runtime_error("LD reference SNP index map is inconsistent");
    }
}

int CoJoEngine::ldIndexForSummary(int summary_index) const {
    if (summary_index < 0 || summary_index >= static_cast<int>(summary_stats_.size()))
        throw std::out_of_range("COJO summary index is out of range");
    const auto found = ld_ref_.snp_to_idx.find(summary_stats_[summary_index].snp_id);
    if (found == ld_ref_.snp_to_idx.end() || found->second < 0 ||
        found->second >= ld_ref_.ld_matrix.rows())
        throw std::runtime_error("COJO summary SNP is not mapped to the LD matrix");
    return found->second;
}

std::vector<int> CoJoEngine::conditionIndices() const {
    std::vector<int> indices;
    std::unordered_set<int> seen;
    for (const auto& id : condition_snps_) {
        const auto found = summary_idx_.find(id);
        if (found == summary_idx_.end())
            throw std::runtime_error("Condition SNP was absent after COJO harmonization: " + id);
        if (seen.insert(found->second).second) indices.push_back(found->second);
    }
    return indices;
}

Eigen::MatrixXd CoJoEngine::extractLDMatrix(const std::vector<int>& indices) const {
    int k = static_cast<int>(indices.size());
    Eigen::MatrixXd R = Eigen::MatrixXd::Identity(k, k);

    if (options_.use_reference_ld && !ld_ref_.ld_matrix.isZero(0)) {
        for (int i = 0; i < k; ++i) {
            for (int j = i + 1; j < k; ++j) {
                const int idx1 = ldIndexForSummary(indices[i]);
                const int idx2 = ldIndexForSummary(indices[j]);
                const double r = ld_ref_.ld_matrix(idx1, idx2);
                R(i, j) = R(j, i) = r;
            }
        }
    } else {
        for (int i = 0; i < k; ++i) {
            for (int j = i + 1; j < k; ++j) {
                R(i, j) = R(j, i) = 0.0;
            }
        }
    }
    return R;
}

Eigen::VectorXd CoJoEngine::computeConditionalEffects(const std::vector<int>& indices) const {
    int k = static_cast<int>(indices.size());
    if (k == 0) return Eigen::VectorXd();

    Eigen::VectorXd z(k);
    Eigen::VectorXd w(k);
    for (int i = 0; i < k; ++i) {
        const auto& s = summary_stats_[indices[i]];
        z[i] = s.z_score;
        w[i] = std::sqrt(s.weight);
    }

    Eigen::MatrixXd R = extractLDMatrix(indices);
    Eigen::LLT<Eigen::MatrixXd> llt(R);
    if (llt.info() != Eigen::Success) {
        Eigen::MatrixXd R_reg = R;
        R_reg.diagonal().array() += 1e-6;
        llt.compute(R_reg);
        if (llt.info() != Eigen::Success) {
            return z;
        }
    }

    Eigen::VectorXd z_cond = llt.solve(z);
    return z_cond;
}

Eigen::MatrixXd CoJoEngine::computeConditionalCovariance(const std::vector<int>& indices) const {
    int k = static_cast<int>(indices.size());
    if (k == 0) return Eigen::MatrixXd(0, 0);
    Eigen::MatrixXd R = extractLDMatrix(indices);
    Eigen::LLT<Eigen::MatrixXd> llt(R);
    if (llt.info() != Eigen::Success) {
        Eigen::MatrixXd R_reg = R;
        R_reg.diagonal().array() += 1e-6;
        llt.compute(R_reg);
        if (llt.info() != Eigen::Success) {
            return Eigen::MatrixXd::Identity(k, k);
        }
    }
    Eigen::MatrixXd R_inv = llt.solve(Eigen::MatrixXd::Identity(k, k));
    return R_inv;
}

std::vector<CoJoResult> CoJoEngine::runConditional() const {
    std::vector<CoJoResult> results;
    results.reserve(summary_stats_.size());

    const std::vector<int> cond_indices = conditionIndices();

    if (cond_indices.empty()) {
        for (const auto& s : summary_stats_) {
            CoJoResult r;
            r.snp_id = s.snp_id;
            r.chrom = s.chrom;
            r.position = s.position;
            r.a1 = s.a1;
            r.a2 = s.a2;
            r.beta = s.beta;
            r.se = s.se;
            r.p_value = s.p_value;
            r.freq = s.freq;
            r.n = s.n;
            r.conditional_beta = s.beta;
            r.conditional_se = s.se;
            r.conditional_p = s.p_value;
            r.valid = true;
            results.push_back(r);
        }
        return results;
    }

    const Eigen::VectorXd z_cond = computeConditionalEffects(cond_indices);
    const Eigen::MatrixXd cov_cond = computeConditionalCovariance(cond_indices);

    for (size_t i = 0; i < summary_stats_.size(); ++i) {
        const auto& s = summary_stats_[i];
        CoJoResult r;
        r.snp_id = s.snp_id;
        r.chrom = s.chrom;
        r.position = s.position;
        r.a1 = s.a1;
        r.a2 = s.a2;
        r.freq = s.freq;
        r.n = s.n;

        bool is_condition = std::find(cond_indices.begin(), cond_indices.end(), i) != cond_indices.end();
        if (is_condition) {
            int cond_pos = std::distance(cond_indices.begin(),
                std::find(cond_indices.begin(), cond_indices.end(), i));
            r.conditional_beta = z_cond[cond_pos] / std::sqrt(s.weight);
            r.conditional_se = std::sqrt(cov_cond(cond_pos, cond_pos) / s.weight);
            r.selected = true;
            r.selected_order = cond_pos + 1;
        } else {
            std::vector<int> joint_indices = cond_indices;
            joint_indices.push_back(static_cast<int>(i));
            const Eigen::VectorXd joint_effects = computeConditionalEffects(joint_indices);
            const Eigen::MatrixXd joint_covariance = computeConditionalCovariance(joint_indices);
            const int candidate = static_cast<int>(joint_indices.size()) - 1;
            r.conditional_beta = joint_effects[candidate] / std::sqrt(s.weight);
            r.conditional_se = std::sqrt(joint_covariance(candidate, candidate) / s.weight);
        }

        r.conditional_p = 2.0 * normal_cdf(-std::abs(r.conditional_beta / r.conditional_se));
        r.beta = s.beta;
        r.se = s.se;
        r.p_value = s.p_value;
        r.valid = true;
        results.push_back(r);
    }

    return results;
}

std::vector<CoJoResult> CoJoEngine::runStepwise() const {
    std::vector<CoJoResult> results;
    results.reserve(summary_stats_.size());

    std::vector<int> selected;
    std::vector<int> candidates;

    for (size_t i = 0; i < summary_stats_.size(); ++i) {
        if (summary_stats_[i].p_value <= options_.threshold) {
            candidates.push_back(static_cast<int>(i));
        }
    }

    std::sort(candidates.begin(), candidates.end(),
              [&](int a, int b) {
                  return summary_stats_[a].p_value < summary_stats_[b].p_value;
              });

    for (int iter = 0; iter < options_.max_iter && !candidates.empty(); ++iter) {
        int best_idx = -1;
        double best_p = 1.0;

        for (int idx : candidates) {
            if (std::find(selected.begin(), selected.end(), idx) != selected.end()) continue;

            bool collinear = false;
            for (int sel : selected) {
                if (checkCollinearity(idx, sel)) {
                    collinear = true;
                    break;
                }
            }
            if (collinear) continue;

            std::vector<int> test_selected = selected;
            test_selected.push_back(idx);
            Eigen::VectorXd z_cond = computeConditionalEffects(test_selected);
            Eigen::MatrixXd cov_cond = computeConditionalCovariance(test_selected);

            int pos = static_cast<int>(test_selected.size()) - 1;
            double beta_cond = z_cond[pos] / std::sqrt(summary_stats_[idx].weight);
            double se_cond = std::sqrt(cov_cond(pos, pos) / summary_stats_[idx].weight);
            double p_cond = 2.0 * normal_cdf(-std::abs(beta_cond / se_cond));

            if (p_cond < best_p) {
                best_p = p_cond;
                best_idx = idx;
            }
        }

        if (best_idx >= 0 && best_p <= options_.threshold) {
            selected.push_back(best_idx);
            candidates.erase(std::remove(candidates.begin(), candidates.end(), best_idx),
                           candidates.end());
        } else {
            break;
        }
    }

    const Eigen::VectorXd selected_effects = computeConditionalEffects(selected);
    const Eigen::MatrixXd selected_covariance = computeConditionalCovariance(selected);
    for (size_t summary_index = 0; summary_index < summary_stats_.size(); ++summary_index) {
        const auto& s = summary_stats_[summary_index];
        CoJoResult r;
        r.snp_id = s.snp_id;
        r.chrom = s.chrom;
        r.position = s.position;
        r.a1 = s.a1;
        r.a2 = s.a2;
        r.freq = s.freq;
        r.n = s.n;
        r.beta = s.beta;
        r.se = s.se;
        r.p_value = s.p_value;
        const auto selected_position = std::find(selected.begin(), selected.end(), static_cast<int>(summary_index));
        if (selected_position != selected.end()) {
            const int position = std::distance(selected.begin(), selected_position);
            r.conditional_beta = selected_effects[position] / std::sqrt(s.weight);
            r.conditional_se = std::sqrt(selected_covariance(position, position) / s.weight);
            r.selected = true;
            r.selected_order = position + 1;
        } else if (!selected.empty()) {
            std::vector<int> joint_indices = selected;
            joint_indices.push_back(static_cast<int>(summary_index));
            const auto effects = computeConditionalEffects(joint_indices);
            const auto covariance = computeConditionalCovariance(joint_indices);
            const int position = static_cast<int>(joint_indices.size()) - 1;
            r.conditional_beta = effects[position] / std::sqrt(s.weight);
            r.conditional_se = std::sqrt(covariance(position, position) / s.weight);
        } else {
            r.conditional_beta = s.beta;
            r.conditional_se = s.se;
        }
        r.conditional_p = 2.0 * normal_cdf(-std::abs(r.conditional_beta / r.conditional_se));
        r.valid = true;
        results.push_back(r);
    }

    return results;
}

std::vector<CoJoResult> CoJoEngine::runJoint() const {
    const std::vector<int> selected = condition_snps_.empty() ? selectIndependentSignals() : conditionIndices();
    const Eigen::VectorXd effects = computeConditionalEffects(selected);
    const Eigen::MatrixXd covariance = computeConditionalCovariance(selected);
    std::vector<CoJoResult> results;
    results.reserve(summary_stats_.size());
    for (size_t i = 0; i < summary_stats_.size(); ++i) {
        const auto& s = summary_stats_[i];
        CoJoResult r;
        r.snp_id = s.snp_id;
        r.chrom = s.chrom;
        r.position = s.position;
        r.a1 = s.a1;
        r.a2 = s.a2;
        r.freq = s.freq;
        r.n = s.n;
        r.beta = s.beta;
        r.se = s.se;
        r.p_value = s.p_value;
        const auto found = std::find(selected.begin(), selected.end(), static_cast<int>(i));
        if (found != selected.end()) {
            const int position = std::distance(selected.begin(), found);
            r.conditional_beta = effects[position] / std::sqrt(s.weight);
            r.conditional_se = std::sqrt(covariance(position, position) / s.weight);
            r.conditional_p = 2.0 * normal_cdf(-std::abs(r.conditional_beta / r.conditional_se));
            r.selected = true;
            r.selected_order = position + 1;
        } else {
            r.conditional_beta = s.beta;
            r.conditional_se = s.se;
            r.conditional_p = s.p_value;
        }
        r.valid = true;
        results.push_back(r);
    }
    return results;
}

std::vector<int> CoJoEngine::selectIndependentSignals() const {
    std::vector<int> candidates;
    for (size_t i = 0; i < summary_stats_.size(); ++i) {
        if (summary_stats_[i].p_value <= options_.threshold) {
            candidates.push_back(static_cast<int>(i));
        }
    }

    std::sort(candidates.begin(), candidates.end(),
              [&](int a, int b) {
                  return summary_stats_[a].p_value < summary_stats_[b].p_value;
              });

    std::vector<int> selected;
    for (int idx : candidates) {
        bool collinear = false;
        for (int sel : selected) {
            if (checkCollinearity(idx, sel)) {
                collinear = true;
                break;
            }
        }
        if (!collinear) selected.push_back(idx);
    }

    return selected;
}

bool CoJoEngine::checkCollinearity(int idx1, int idx2) const {
    if (!options_.use_reference_ld || ld_ref_.ld_matrix.rows() == 0) return false;
    const auto& left = summary_stats_.at(idx1);
    const auto& right = summary_stats_.at(idx2);
    if (left.chrom != right.chrom || std::llabs(left.position - right.position) > options_.window_size) return false;
    const double r = ld_ref_.ld_matrix(ldIndexForSummary(idx1), ldIndexForSummary(idx2));
    return r * r >= options_.collinearity;
}

std::vector<CoJoResult> CoJoEngine::run() {
    validateLDReference();
    alignAlleles();

    switch (options_.mode) {
        case CoJoMode::CONDITIONAL:
            return runConditional();
        case CoJoMode::STEPWISE:
            return runStepwise();
        case CoJoMode::JOINT:
            return runJoint();
    }
    return runConditional();
}

void writeCoJoResults(
    const std::vector<CoJoResult>& results,
    const std::string& out_prefix,
    bool write_jma) {

    std::ofstream out(out_prefix + ".cojo.assoc");
    out << "SNP\tCHR\tBP\tA1\tA2\tFREQ\tN\tBETA\tSE\tP\tCOND_BETA\tCOND_SE\tCOND_P\tSELECTED\tSELECTED_ORDER\n";
    out << std::setprecision(10);
    for (const auto& r : results) {
        out << r.snp_id << '\t' << r.chrom << '\t' << r.position << '\t'
            << r.a1 << '\t' << r.a2 << '\t' << r.freq << '\t'
            << r.n << '\t' << r.beta << '\t' << r.se << '\t' << r.p_value << '\t'
            << r.conditional_beta << '\t' << r.conditional_se << '\t' << r.conditional_p << '\t'
            << (r.selected ? 1 : 0) << '\t' << r.selected_order << '\n';
    }

    if (write_jma) {
        std::ofstream jout(out_prefix + ".jma");
        jout << "SNP\tCHR\tBP\tA1\tA2\tFREQ\tN\tBETA\tSE\tP\tSELECTED\n";
        jout << std::setprecision(10);
        for (const auto& r : results) {
            if (r.selected) {
                jout << r.snp_id << '\t' << r.chrom << '\t' << r.position << '\t'
                     << r.a1 << '\t' << r.a2 << '\t' << r.freq << '\t'
                     << r.n << '\t' << r.conditional_beta << '\t'
                     << r.conditional_se << '\t' << r.conditional_p << '\t' << 1 << '\n';
            }
        }
    }
}

} // namespace cosmic::association::cojo
