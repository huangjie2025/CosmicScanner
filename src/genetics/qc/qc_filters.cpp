#include "qc_filters.h"
#include "plink_reader.h"
#include "qc_utils.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace cosmic {
namespace genetics {

QCResult QCFilter::preprocess(std::string input_prefix, std::string output_prefix,
                              QCOptions options) {
    namespace fs = std::filesystem;
    const auto validate = [](double value, double maximum, const char* name) {
        if (!std::isfinite(value) || value < 0 || value > maximum)
            throw std::invalid_argument(std::string("QC --") + name +
                                        " must be finite and in [0, " +
                                        std::to_string(maximum) + "]");
    };
    validate(options.geno, 1, "geno");
    validate(options.mind, 1, "mind");
    validate(options.maf, .5, "maf");
    validate(options.hwe, 1, "hwe");
    if (input_prefix.empty() || output_prefix.empty())
        throw std::invalid_argument("QC input and output prefixes must not be empty");
    if (fs::weakly_canonical(input_prefix) == fs::weakly_canonical(output_prefix))
        throw std::invalid_argument("QC input and output prefixes must differ");

    const fs::path destination(output_prefix);
    const fs::path parent = destination.has_parent_path() ? destination.parent_path() : fs::path(".");
    const std::string basename = destination.filename().string();
    if (basename.empty()) throw std::invalid_argument("QC output prefix needs a filename");
    for (const auto& entry : fs::directory_iterator(parent)) {
        const auto name = entry.path().filename().string();
        if (name == basename || name.compare(0, basename.size() + 1, basename + ".") == 0)
            throw std::runtime_error("QC destination prefix already exists: " + output_prefix);
    }

    // Preserve metadata verbatim, while validating rows before the permissive reader.
    const auto metadata = [](const std::string& path, bool fam) {
        std::ifstream in(path);
        if (!in) throw std::runtime_error("Cannot open QC metadata: " + path);
        std::vector<std::string> lines;
        std::string line;
        while (std::getline(in, line)) {
            if (line.empty()) continue;
            std::istringstream row(line);
            std::string a, b, c, d, extra;
            double real;
            long long integer;
            bool valid;
            if (fam) {
                int sex;
                valid = static_cast<bool>(row >> a >> b >> c >> d >> sex >> real);
            } else {
                valid = static_cast<bool>(row >> a >> b >> real >> integer >> c >> d);
            }
            if (!valid || !std::isfinite(real) || (row >> extra))
                throw std::runtime_error("Invalid six-column QC metadata row: " + path);
            lines.push_back(line);
        }
        if (in.bad()) throw std::runtime_error("Cannot read QC metadata: " + path);
        if (lines.empty() || lines.size() > static_cast<size_t>(std::numeric_limits<int>::max()))
            throw std::runtime_error("QC metadata must contain 1..INT_MAX rows: " + path);
        return lines;
    };
    const auto fam_lines = metadata(input_prefix + ".fam", true);
    const auto bim_lines = metadata(input_prefix + ".bim", false);
    const size_t stride = (fam_lines.size() + 3) / 4;
    if (bim_lines.size() > (std::numeric_limits<size_t>::max() - 3) / stride)
        throw std::runtime_error("QC BED dimensions overflow");
    const size_t expected_size = 3 + stride * bim_lines.size();
    {
        std::ifstream bed(input_prefix + ".bed", std::ios::binary | std::ios::ate);
        if (!bed || bed.tellg() != static_cast<std::streamoff>(expected_size))
            throw std::runtime_error("QC BED file length does not match BIM/FAM");
        bed.seekg(0);
        std::array<unsigned char, 3> header{};
        if (!bed.read(reinterpret_cast<char*>(header.data()), header.size()) ||
            header != std::array<unsigned char, 3>{0x6c, 0x1b, 0x01})
            throw std::runtime_error("QC requires a valid SNP-major BED header");
    }

    // Sample lists require exactly FID IID, never IID-only or token-wise matching.
    using ID = std::pair<std::string, std::string>;
    const auto load_list = [](const std::string& path, bool paired) {
        std::set<ID> ids;
        if (path.empty()) return ids;
        std::ifstream in(path);
        if (!in) throw std::runtime_error("Cannot open QC selection list: " + path);
        std::string line;
        while (std::getline(in, line)) {
            std::istringstream row(line);
            std::string first, second, extra;
            if (!(row >> first)) continue;
            if ((paired && !(row >> second)) || (row >> extra))
                throw std::runtime_error("QC list requires " +
                    std::string(paired ? "exact FID IID pairs: " : "one variant ID per row: ") + path);
            ids.emplace(first, second);
        }
        if (in.bad() || ids.empty())
            throw std::runtime_error("Empty or unreadable QC selection list: " + path);
        return ids;
    };
    const auto extract = load_list(options.extract_file, false);
    const auto exclude = load_list(options.exclude_file, false);
    const auto keep = load_list(options.keep_file, true);
    const auto remove = load_list(options.remove_file, true);

    PlinkReader reader;
    reader.load(input_prefix); // Header and exact file length were checked before mmap.
    if (reader.getNumSamples() != fam_lines.size() || reader.getNumSnps() != bim_lines.size())
        throw std::runtime_error("QC input metadata changed during loading");
    std::vector<int> candidates, selected_samples;
    size_t fail_extract = 0, fail_exclude = 0, fail_keep = 0, fail_remove = 0;
    for (size_t j = 0; j < bim_lines.size(); ++j) {
        const ID id{reader.getBimInfo()[j].id, ""};
        if (exclude.count(id)) ++fail_exclude;
        else if (!options.extract_file.empty() && !extract.count(id)) ++fail_extract;
        else candidates.push_back(static_cast<int>(j));
    }
    if (candidates.empty()) throw std::runtime_error("QC selected zero candidate variants");
    for (size_t i = 0; i < fam_lines.size(); ++i) {
        const auto& sample = reader.getFamInfo()[i];
        const ID id{sample.fid, sample.iid};
        if (remove.count(id)) ++fail_remove;
        else if (!options.keep_file.empty() && !keep.count(id)) ++fail_keep;
        else selected_samples.push_back(static_cast<int>(i));
    }
    if (selected_samples.empty()) throw std::runtime_error("QC selected zero samples");

    Eigen::VectorXf genotypes;
    std::vector<size_t> missing(selected_samples.size(), 0);
    for (int j : candidates) {
        reader.getSnpGenotypes(j, genotypes);
        for (size_t i = 0; i < selected_samples.size(); ++i)
            if (std::isnan(genotypes[selected_samples[i]])) ++missing[i];
    }
    QCResult result{output_prefix, {}, {}};
    for (size_t i = 0; i < selected_samples.size(); ++i)
        if (static_cast<double>(missing[i]) / candidates.size() <= options.mind)
            result.samples.push_back(selected_samples[i]);
    if (result.samples.empty()) throw std::runtime_error("QC --mind removed all samples");

    size_t fail_geno = 0, fail_maf = 0, fail_hwe = 0, fail_no_calls = 0;
    for (int j : candidates) {
        reader.getSnpGenotypes(j, genotypes);
        long long hom1 = 0, het = 0, hom2 = 0, absent = 0;
        for (int i : result.samples) {
            const float g = genotypes[i];
            if (std::isnan(g)) ++absent;
            else if (g == 2) ++hom1;
            else if (g == 1) ++het;
            else ++hom2;
        }
        if (static_cast<double>(absent) / result.samples.size() > options.geno) ++fail_geno;
        else if (hom1 + het + hom2 == 0) ++fail_no_calls;
        else if (cosmic::calculateMAF(hom1, het, hom2) < options.maf) ++fail_maf;
        else if (options.hwe > 0 && cosmic::calculateHWE(hom1, het, hom2) < options.hwe) ++fail_hwe;
        else result.variants.push_back(j);
    }
    if (result.variants.empty()) throw std::runtime_error("QC removed all variants");

    const std::array<std::string, 4> paths{output_prefix + ".bed", output_prefix + ".bim",
        output_prefix + ".fam", output_prefix + ".qc_report.txt"};
    using Stream = std::unique_ptr<FILE, decltype(&std::fclose)>;
    std::array<Stream, 4> outputs{Stream(nullptr, &std::fclose), Stream(nullptr, &std::fclose),
        Stream(nullptr, &std::fclose), Stream(nullptr, &std::fclose)};
    size_t created = 0;
    try {
        // 'x' is exclusive creation: a racing writer or symlink cannot be overwritten.
        for (size_t k = 0; k < paths.size(); ++k) {
            outputs[k].reset(std::fopen(paths[k].c_str(), "wbx"));
            if (!outputs[k]) throw std::runtime_error("Cannot exclusively create QC output: " + paths[k]);
            ++created;
        }
        const auto write = [](FILE* out, const void* data, size_t size) {
            if (std::fwrite(data, 1, size, out) != size)
                throw std::runtime_error("Failed writing QC output");
        };
        const std::array<unsigned char, 3> header{0x6c, 0x1b, 0x01};
        write(outputs[0].get(), header.data(), header.size());
        std::vector<unsigned char> packed((result.samples.size() + 3) / 4);
        for (int j : result.variants) {
            std::fill(packed.begin(), packed.end(), 0);
            const uint8_t* raw = reader.getData() + 3 + static_cast<size_t>(j) * stride;
            for (size_t k = 0; k < result.samples.size(); ++k) {
                const size_t i = result.samples[k];
                const unsigned char call = (raw[i / 4] >> (2 * (i % 4))) & 3;
                packed[k / 4] |= call << (2 * (k % 4));
            }
            write(outputs[0].get(), packed.data(), packed.size());
            const auto line = bim_lines[j] + "\n";
            write(outputs[1].get(), line.data(), line.size());
        }
        for (int i : result.samples) {
            const auto line = fam_lines[i] + "\n";
            write(outputs[2].get(), line.data(), line.size());
        }
        std::ostringstream report;
        report << std::setprecision(17)
               << "QC Report (BED hardcalls; sequential exclusion counts)\n"
               << "input_prefix\t" << input_prefix << "\noutput_prefix\t" << output_prefix
               << "\ngeno\t" << options.geno << "\nmind\t" << options.mind
               << "\nmaf\t" << options.maf << "\nhwe\t" << options.hwe
               << "\nsamples_total\t" << fam_lines.size()
               << "\nsamples_fail_keep\t" << fail_keep << "\nsamples_fail_remove\t" << fail_remove
               << "\nsamples_fail_mind\t" << selected_samples.size() - result.samples.size()
               << "\nsamples_pass\t" << result.samples.size()
               << "\nvariants_total\t" << bim_lines.size()
               << "\nvariants_fail_extract\t" << fail_extract << "\nvariants_fail_exclude\t" << fail_exclude
               << "\nvariants_candidates\t" << candidates.size()
               << "\nvariants_fail_geno\t" << fail_geno << "\nvariants_fail_no_calls\t" << fail_no_calls
               << "\nvariants_fail_maf\t" << fail_maf << "\nvariants_fail_hwe\t" << fail_hwe
               << "\nvariants_pass\t" << result.variants.size() << '\n';
        const auto text = report.str();
        write(outputs[3].get(), text.data(), text.size());
        for (auto& out : outputs)
            if (std::fclose(out.release()) != 0) throw std::runtime_error("Failed closing QC output");
    } catch (...) {
        for (auto& out : outputs) out.reset();
        for (size_t k = 0; k < created; ++k) {
            std::error_code ignored;
            fs::remove(paths[k], ignored);
        }
        throw;
    }
    return result;
}

} // namespace genetics
} // namespace cosmic
