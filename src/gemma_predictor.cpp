#include "gemma_predictor.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

#include <Eigen/Dense>

#include "plink_reader.h"
#include "bgen_reader.h"
#include "pgen_reader.h"
#include "vcf/vcf_reader.h"

namespace cosmic {

namespace {

std::vector<std::string> split_ws(const std::string& line) {
    std::vector<std::string> out;
    std::stringstream ss(line);
    std::string tok;
    while (ss >> tok) out.push_back(tok);
    return out;
}

int find_col(const std::vector<std::string>& header,
             const std::vector<std::string>& candidates,
             const std::string& requested,
             const std::string& label) {
    if (!requested.empty()) {
        for (int i = 0; i < (int)header.size(); ++i) {
            if (header[i] == requested) return i;
        }
        throw std::runtime_error("Prediction effect file is missing requested " + label + " column: " + requested);
    }
    for (const auto& cand : candidates) {
        for (int i = 0; i < (int)header.size(); ++i) {
            if (header[i] == cand) return i;
        }
    }
    throw std::runtime_error("Prediction effect file is missing " + label + " column");
}

std::unordered_map<std::string, double> load_effects(const GemmaPredictOptions& options) {
    std::ifstream in(options.effect_file);
    if (!in) throw std::runtime_error("Cannot open prediction effect file: " + options.effect_file);

    std::string line;
    while (std::getline(in, line)) {
        if (!line.empty() && line[0] != '#') break;
    }
    if (line.empty()) throw std::runtime_error("Prediction effect file is empty: " + options.effect_file);

    const auto header = split_ws(line);
    const int snp_col = find_col(header, {"SNP", "rs", "id", "ID", "variant_id"}, options.snp_col, "SNP ID");
    const int eff_col = find_col(header, {"BETA", "beta", "alpha", "effect", "EFFECT", "EB_BETA"}, options.effect_col, "effect");

    std::unordered_map<std::string, double> effects;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_ws(line);
        if ((int)fields.size() <= std::max(snp_col, eff_col)) continue;
        try {
            double beta = std::stod(fields[eff_col]);
            if (std::isfinite(beta)) effects[fields[snp_col]] = beta;
        } catch (...) {
        }
    }
    if (effects.empty()) throw std::runtime_error("No usable SNP effects loaded from: " + options.effect_file);
    return effects;
}

void add_dosage_score(Eigen::VectorXd& score,
                      const Eigen::VectorXd& dosage,
                      double effect,
                      int& variants_used) {
    double sum = 0.0;
    int n_valid = 0;
    for (int i = 0; i < dosage.size(); ++i) {
        if (std::isfinite(dosage[i])) {
            sum += dosage[i];
            ++n_valid;
        }
    }
    if (n_valid == 0) return;
    const double mean = sum / n_valid;
    for (int i = 0; i < dosage.size(); ++i) {
        score[i] += (std::isfinite(dosage[i]) ? dosage[i] : mean) * effect;
    }
    ++variants_used;
}

} // namespace

GemmaPredictor::GemmaPredictor(GemmaPredictOptions options) : options_(std::move(options)) {}

void GemmaPredictor::run(const std::vector<std::string>& fids,
                         const std::vector<std::string>& iids,
                         PlinkReader& plink_reader,
                         BgenReader& bgen_reader,
                         PgenReader& pgen_reader,
                         VcfReader& vcf_reader,
                         bool is_bgen,
                         bool is_pgen,
                         bool is_vcf) {
    auto effects = load_effects(options_);
    const int n = (int)fids.size();
    Eigen::VectorXd score = Eigen::VectorXd::Zero(n);
    int variants_used = 0;

    if (is_pgen) {
        pgen_reader.reset();
        std::string chrom, rsid, ref, alt;
        uint32_t pos = 0;
        Eigen::VectorXf dosage_f;
        while (pgen_reader.readVariant(chrom, rsid, pos, ref, alt, dosage_f)) {
            auto it = effects.find(rsid);
            if (it == effects.end()) continue;
            add_dosage_score(score, dosage_f.cast<double>(), it->second, variants_used);
        }
    } else if (is_bgen) {
        bgen_reader.reset();
        std::string chrom, rsid, a1, a2;
        uint32_t pos = 0;
        Eigen::VectorXf dosage_f;
        while (bgen_reader.readVariant(chrom, rsid, pos, a1, a2, dosage_f)) {
            auto it = effects.find(rsid);
            if (it == effects.end()) continue;
            add_dosage_score(score, dosage_f.cast<double>(), it->second, variants_used);
        }
    } else if (is_vcf) {
        vcf_reader.reset();
        std::string chrom, rsid, ref, alt;
        uint32_t pos = 0;
        Eigen::VectorXf dosage_f;
        while (vcf_reader.readVariant(chrom, rsid, pos, ref, alt, dosage_f)) {
            auto it = effects.find(rsid);
            if (it == effects.end()) continue;
            add_dosage_score(score, dosage_f.cast<double>(), it->second, variants_used);
        }
    } else {
        const auto& bim = plink_reader.getBimInfo();
        const size_t block_size = 512;
        for (size_t start = 0; start < bim.size(); start += block_size) {
            size_t end = std::min(bim.size(), start + block_size);
            Eigen::MatrixXd raw;
            plink_reader.getSnpBlock(start, end, raw);
            for (size_t k = 0; k < end - start; ++k) {
                const auto& snp = bim[start + k];
                auto it = effects.find(snp.id);
                if (it == effects.end()) continue;
                add_dosage_score(score, raw.col((int)k), it->second, variants_used);
            }
        }
    }

    const std::string out = options_.output_file.empty() ? "cosmic_predict.txt" : options_.output_file;
    std::ofstream ofs(out);
    if (!ofs) throw std::runtime_error("Cannot write prediction output: " + out);
    ofs << "FID\tIID\tPRED\tN_VARIANTS_USED\n";
    ofs.setf(std::ios::scientific);
    ofs.precision(12);
    for (int i = 0; i < n; ++i) {
        ofs << fids[i] << '\t' << iids[i] << '\t' << score[i] << '\t' << variants_used << '\n';
    }
    std::cout << "GEMMA-style prediction complete. Loaded " << effects.size()
              << " effects, matched " << variants_used << " variants. Results written to "
              << out << "\n";
}

} // namespace cosmic
