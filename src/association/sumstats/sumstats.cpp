#include "association/sumstats/sumstats.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <ctime>

namespace cosmic::association::sumstats {

// --- SumstatsWriter ---

SumstatsWriter::SumstatsWriter(const std::string& filepath, const SumstatsMetadata& meta)
    : meta_(meta) {
    file_.open(filepath);
    if (!file_) {
        throw std::runtime_error("Cannot open sumstats file for writing: " + filepath);
    }
    
    if (meta_.date.empty()) {
        auto now = std::chrono::system_clock::now();
        auto in_time_t = std::chrono::system_clock::to_time_t(now);
        std::stringstream ss;
        ss << std::put_time(std::gmtime(&in_time_t), "%Y-%m-%d %H:%M:%S UTC");
        meta_.date = ss.str();
    }
}

SumstatsWriter::~SumstatsWriter() {
    close();
}

void SumstatsWriter::writeHeader() {
    if (header_written_) return;
    
    // Write metadata as comments
    file_ << "# CosmicScanner Summary Statistics v" << meta_.version << "\n";
    file_ << "# Created: " << meta_.date << "\n";
    if (!meta_.phenotype.empty()) file_ << "# Phenotype: " << meta_.phenotype << "\n";
    if (!meta_.model.empty()) file_ << "# Model: " << meta_.model << "\n";
    if (!meta_.reference_genome.empty()) file_ << "# Reference: " << meta_.reference_genome << "\n";
    if (meta_.n_samples > 0) file_ << "# N_Samples: " << meta_.n_samples << "\n";
    if (!meta_.command_line.empty()) file_ << "# Command: " << meta_.command_line << "\n";
    file_ << "#\n";
    
    // Column header
    file_ << "CHR\tSNP\tPOS\tA1\tA2\tBETA\tSE\tP\tN\tAF\tLOG10P\tMODEL";
    file_ << "\tZ\tWEIGHT\tINFO\tWARNING\n";
    
    header_written_ = true;
}

void SumstatsWriter::writeRecord(const SumstatsRecord& record) {
    writeHeader();
    
    file_ << record.chrom << '\t'
          << record.snp_id << '\t'
          << record.pos << '\t'
          << record.a1 << '\t'
          << record.a2 << '\t'
          << std::setprecision(10) << record.beta << '\t'
          << std::setprecision(10) << record.se << '\t'
          << std::setprecision(10) << record.p_value << '\t'
          << record.n << '\t'
          << std::setprecision(6) << record.freq << '\t'
          << std::setprecision(6) << record.log10_p << '\t'
          << record.model << '\t'
          << std::setprecision(6) << record.z_score << '\t'
          << std::setprecision(6) << record.weight << '\t'
          << std::setprecision(4) << record.info << '\t'
          << (record.warning.empty() ? "." : record.warning) << '\n';
}

void SumstatsWriter::close() {
    if (file_.is_open()) {
        file_.close();
    }
}

// --- SumstatsReader ---

SumstatsReader::SumstatsReader(const std::string& filepath) {
    file_.open(filepath);
    if (!file_) {
        throw std::runtime_error("Cannot open sumstats file for reading: " + filepath);
    }
}

SumstatsReader::~SumstatsReader() {
    if (file_.is_open()) file_.close();
}

bool SumstatsReader::readHeader() {
    if (header_read_) return true;
    
    std::string line;
    while (std::getline(file_, line)) {
        if (line.empty()) continue;
        if (line[0] == '#') {
            // Parse metadata
            if (line.rfind("# Phenotype: ", 0) == 0) meta_.phenotype = line.substr(13);
            else if (line.rfind("# Model: ", 0) == 0) meta_.model = line.substr(9);
            else if (line.rfind("# Reference: ", 0) == 0) meta_.reference_genome = line.substr(13);
            else if (line.rfind("# N_Samples: ", 0) == 0) meta_.n_samples = std::stoi(line.substr(13));
            else if (line.rfind("# Command: ", 0) == 0) meta_.command_line = line.substr(11);
            else if (line.rfind("# Created: ", 0) == 0) meta_.date = line.substr(11);
            else if (line.rfind("# CosmicScanner Summary Statistics v", 0) == 0) {
                meta_.version = line.substr(36);
            }
            continue;
        }
        
        // This is the column header
        std::istringstream iss(line);
        std::string col;
        while (iss >> col) {
            col_idx_[col] = static_cast<int>(column_names_.size());
            column_names_.push_back(col);
        }
        header_read_ = true;
        return true;
    }
    return false;
}

bool SumstatsReader::readRecord(SumstatsRecord& record) {
    if (!header_read_) readHeader();
    
    std::string line;
    while (std::getline(file_, line)) {
        if (line.empty() || line[0] == '#') continue;
        
        std::istringstream iss(line);
        std::vector<std::string> fields;
        std::string field;
        while (iss >> field) fields.push_back(field);
        
        if (fields.size() < 11) continue; // Minimum required columns
        
        auto get = [&](const std::string& name, const std::string& def = "") -> std::string {
            auto it = col_idx_.find(name);
            if (it != col_idx_.end() && it->second < (int)fields.size()) {
                return fields[it->second];
            }
            return def;
        };
        
        auto get_d = [&](const std::string& name, double def = 0.0) -> double {
            auto s = get(name);
            if (s.empty() || s == "." || s == "NA") return def;
            try { return std::stod(s); } catch (...) { return def; }
        };
        
        auto get_ll = [&](const std::string& name, long long def = 0) -> long long {
            auto s = get(name);
            if (s.empty() || s == "." || s == "NA") return def;
            try { return std::stoll(s); } catch (...) { return def; }
        };
        
        record.chrom = get("CHR");
        record.snp_id = get("SNP");
        record.pos = get_ll("POS", get_ll("BP"));
        record.a1 = get("A1");
        record.a2 = get("A2");
        record.beta = get_d("BETA_A", get_d("BETA_D", get_d("BETA")));
        record.se = get_d("SE_A", get_d("SE_D", get_d("SE")));
        record.p_value = get_d("P_A", get_d("P_D", get_d("P")));
        record.n = get_ll("N");
        record.freq = get_d("AF", get_d("FREQ"));
        record.log10_p = get_d("LOG10P_A", get_d("LOG10P_D", get_d("LOG10P")));
        record.model = get("MODEL");
        record.z_score = get_d("Z");
        record.weight = get_d("WEIGHT");
        record.info = get_d("INFO", 1.0);
        record.warning = get("WARNING", "");
        
        // Compute missing fields
        if (record.se > 0) {
            record.z_score = record.beta / record.se;
            record.weight = 1.0 / (record.se * record.se);
        }
        if (record.p_value > 0 && record.log10_p == 0) {
            record.log10_p = -std::log10(record.p_value);
        }
        
        return true;
    }
    return false;
}

bool SumstatsReader::hasColumn(const std::string& col) const {
    return col_idx_.find(col) != col_idx_.end();
}

// --- Allele harmonization ---

std::string complementAllele(const std::string& allele) {
    std::string result = allele;
    for (char& c : result) {
        switch (std::toupper(c)) {
            case 'A': c = 'T'; break;
            case 'T': c = 'A'; break;
            case 'C': c = 'G'; break;
            case 'G': c = 'C'; break;
            default: c = c; break;
        }
    }
    return result;
}

bool isAmbiguousSNP(const std::string& a1, const std::string& a2) {
    std::string u1 = a1, u2 = a2;
    std::transform(u1.begin(), u1.end(), u1.begin(), ::toupper);
    std::transform(u2.begin(), u2.end(), u2.begin(), ::toupper);
    
    // A/T or T/A or C/G or G/C
    return (u1 == "A" && u2 == "T") || (u1 == "T" && u2 == "A") ||
           (u1 == "C" && u2 == "G") || (u1 == "G" && u2 == "C");
}

AlleleMatchStatus checkAlleleMatch(const std::string& a1, const std::string& a2,
                                   const std::string& ref_a1, const std::string& ref_a2) {
    auto uppercase = [](std::string allele) {
        std::transform(allele.begin(), allele.end(), allele.begin(),
                       [](unsigned char value) { return static_cast<char>(std::toupper(value)); });
        return allele;
    };
    const std::string u_a1 = uppercase(a1);
    const std::string u_a2 = uppercase(a2);
    const std::string u_ref_a1 = uppercase(ref_a1);
    const std::string u_ref_a2 = uppercase(ref_a2);

    // Without a trusted frequency reference, strand-ambiguous SNPs cannot be
    // oriented safely even when their allele strings appear to match.
    if (isAmbiguousSNP(u_a1, u_a2) || isAmbiguousSNP(u_ref_a1, u_ref_a2))
        return AlleleMatchStatus::AMBIGUOUS;

    if (u_a1 == u_ref_a1 && u_a2 == u_ref_a2) return AlleleMatchStatus::MATCH;
    if (u_a1 == u_ref_a2 && u_a2 == u_ref_a1) return AlleleMatchStatus::SWAPPED;
    
    // Check strand complement
    std::string c_a1 = complementAllele(u_a1);
    std::string c_a2 = complementAllele(u_a2);
    
    if (c_a1 == u_ref_a1 && c_a2 == u_ref_a2) return AlleleMatchStatus::STRAND_FLIP;
    if (c_a1 == u_ref_a2 && c_a2 == u_ref_a1) return AlleleMatchStatus::STRAND_FLIP_SWAPPED;
    
    return AlleleMatchStatus::MISMATCH;
}

AlleleHarmonizationResult harmonizeAlleles(const SumstatsRecord& record,
                                           const std::string& ref_a1,
                                           const std::string& ref_a2) {
    AlleleHarmonizationResult result;
    result.record = record;
    result.original_a1 = record.a1;
    result.original_a2 = record.a2;
    result.original_beta = record.beta;
    
    if (ref_a1.empty() || ref_a2.empty()) {
        result.status = AlleleMatchStatus::MISSING_REF;
        result.reason = "Reference alleles not provided";
        return result;
    }
    
    result.status = checkAlleleMatch(record.a1, record.a2, ref_a1, ref_a2);
    
    switch (result.status) {
        case AlleleMatchStatus::MATCH:
            result.record.a1 = ref_a1;
            result.record.a2 = ref_a2;
            result.harmonized = true;
            result.reason = "Alleles match reference";
            break;
            
        case AlleleMatchStatus::SWAPPED:
            result.record.a1 = ref_a1;
            result.record.a2 = ref_a2;
            result.record.beta = -record.beta;
            result.record.freq = 1.0 - record.freq;
            result.harmonized = true;
            result.reason = "Alleles swapped, beta negated";
            break;
            
        case AlleleMatchStatus::STRAND_FLIP:
            result.record.a1 = ref_a1;
            result.record.a2 = ref_a2;
            result.harmonized = true;
            result.reason = "Strand flipped to match reference";
            break;

        case AlleleMatchStatus::STRAND_FLIP_SWAPPED:
            result.record.a1 = ref_a1;
            result.record.a2 = ref_a2;
            result.record.beta = -record.beta;
            result.record.freq = 1.0 - record.freq;
            result.harmonized = true;
            result.reason = "Strand flipped and alleles swapped, beta negated";
            break;
            
        case AlleleMatchStatus::AMBIGUOUS:
            result.harmonized = false;
            result.reason = "Ambiguous A/T or C/G SNP - cannot resolve strand";
            break;
            
        case AlleleMatchStatus::MISMATCH:
            result.harmonized = false;
            result.reason = "Allele mismatch: " + record.a1 + "/" + record.a2 + " vs " + ref_a1 + "/" + ref_a2;
            break;
            
        case AlleleMatchStatus::MISSING_REF:
            result.harmonized = false;
            break;
    }

    if (result.harmonized && result.record.se > 0.0) {
        result.record.z_score = result.record.beta / result.record.se;
        result.record.weight = 1.0 / (result.record.se * result.record.se);
    }
    
    return result;
}

std::vector<AlleleHarmonizationResult> harmonizeSumstats(
    const std::vector<SumstatsRecord>& records,
    const std::unordered_map<std::string, std::pair<std::string, std::string>>& ref_alleles) {
    
    std::vector<AlleleHarmonizationResult> results;
    results.reserve(records.size());
    
    for (const auto& record : records) {
        auto it = ref_alleles.find(record.snp_id);
        if (it == ref_alleles.end()) {
            AlleleHarmonizationResult r;
            r.record = record;
            r.status = AlleleMatchStatus::MISSING_REF;
            r.reason = "SNP not in reference panel";
            results.push_back(r);
            continue;
        }
        results.push_back(harmonizeAlleles(record, it->second.first, it->second.second));
    }
    
    return results;
}

// --- Convenience function ---

void writeCosmicSumstats(const std::string& filepath,
                         const std::vector<SumstatsRecord>& records,
                         const SumstatsMetadata& meta) {
    SumstatsWriter writer(filepath, meta);
    for (const auto& record : records) {
        writer.writeRecord(record);
    }
    writer.close();
}

} // namespace cosmic::association::sumstats
