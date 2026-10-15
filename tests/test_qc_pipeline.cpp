#include "genetics/qc/qc_filters.h"
#include "plink_reader.h"
#include "qc_utils.h"

#include <cmath>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;
using cosmic::genetics::QCFilter;
using cosmic::genetics::QCOptions;
using cosmic::genetics::QCResult;

static void check(bool ok, const std::string& message) {
    if (!ok) throw std::runtime_error(message);
}

static void write(const std::string& path, const std::string& content) {
    std::ofstream out(path, std::ios::binary);
    out.exceptions(std::ios::badbit | std::ios::failbit);
    out.write(content.data(), content.size());
    out.close();
}

static std::string read(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    check(static_cast<bool>(in), "Cannot read " + path);
    return {std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>()};
}

template<class F> static void rejects(F action, const std::string& message) {
    bool rejected = false;
    try { action(); } catch (const std::exception&) { rejected = true; }
    check(rejected, "Expected rejection: " + message);
}

int main(int argc, char** argv) {
    try {
        check(argc == 2, "Usage: test_qc_pipeline <temporary-directory>");
        fs::create_directories(argv[1]);
        // Use a dedicated new directory, never overwrite a previous fixture run.
        const auto run_id = std::chrono::steady_clock::now().time_since_epoch().count();
        const fs::path dir = fs::path(argv[1]) / ("qc_pipeline_" + std::to_string(run_id));
        check(fs::create_directory(dir), "Fixture directory already exists: " + dir.string());
        const auto path = [&](const std::string& name) { return (dir / name).string(); };
        const std::string input = path("input");
        const std::vector<std::string> fam{
            "F1 same 0 0 1 -9", "F2 same 0 0 2 1.2500", "F3 i3 dad mom 1 2",
            "F4 i4 0 0 2 -9", "F5 i5 0 0 0 0", "F6 i6 0 0 1 3"};
        const std::vector<std::string> bim{
            "1 v0 0.000 100 A G", "1 v1 0.1250 200 C T", "2 v2 0 300 AT A",
            "2 v3 1.0 400 G C", "X v4 2.50 500 T C", "3 v5 0 600 A C"};
        // Raw BED codes: 00=2, 01=missing, 10=1, 11=0. Six samples cross byte boundaries.
        const std::vector<std::vector<unsigned char>> calls{
            {0, 1, 2, 3, 0, 3}, {1, 1, 2, 3, 0, 3}, {0, 0, 0, 0, 0, 0},
            {0, 0, 0, 3, 3, 3}, {0, 2, 2, 3, 0, 3}, {1, 1, 1, 1, 1, 1}};
        std::string fam_text, bim_text, bed_text("\x6c\x1b\x01", 3);
        for (const auto& row : fam) fam_text += row + '\n';
        for (const auto& row : bim) bim_text += row + '\n';
        for (const auto& variant : calls) {
            unsigned char bytes[2]{0, 0};
            for (size_t i = 0; i < variant.size(); ++i) bytes[i / 4] |= variant[i] << (2 * (i % 4));
            // Nonzero input padding must never leak into real output sample calls.
            bytes[1] |= 0xf0;
            bed_text.append(reinterpret_cast<char*>(bytes), 2);
        }
        write(input + ".fam", fam_text);
        write(input + ".bim", bim_text);
        write(input + ".bed", bed_text);
        const auto list = [&](const std::string& name, const std::string& text) {
            const auto file = path(name);
            write(file, text);
            return file;
        };
        const auto verify = [&](const QCResult& result, const std::vector<int>& samples,
                                const std::vector<int>& variants) {
            check(result.samples == samples, "Wrong retained samples: " + result.prefix);
            check(result.variants == variants, "Wrong retained variants: " + result.prefix);
            std::string expected_fam, expected_bim;
            for (int i : samples) expected_fam += fam[i] + '\n';
            for (int j : variants) expected_bim += bim[j] + '\n';
            check(read(result.prefix + ".fam") == expected_fam, "FAM metadata/order changed");
            check(read(result.prefix + ".bim") == expected_bim, "BIM alleles/metadata/order changed");
            const auto bytes = read(result.prefix + ".bed");
            const size_t stride = (samples.size() + 3) / 4;
            check(bytes.size() == 3 + stride * variants.size(), "Output BED length");
            cosmic::PlinkReader reader;
            reader.load(result.prefix);
            Eigen::VectorXf genotypes;
            const float dosage[]{2, NAN, 1, 0};
            for (size_t j = 0; j < variants.size(); ++j) {
                reader.getSnpGenotypes(j, genotypes);
                for (size_t i = 0; i < samples.size(); ++i) {
                    const unsigned char raw = (static_cast<unsigned char>(bytes[3 + j * stride + i / 4]) >>
                                               (2 * (i % 4))) & 3;
                    const auto expected = calls[variants[j]][samples[i]];
                    check(raw == expected, "Raw BED call changed");
                    check(expected == 1 ? std::isnan(genotypes[i]) : genotypes[i] == dosage[expected],
                          "Reader dosage mapping changed");
                }
                if (samples.size() % 4)
                    check((static_cast<unsigned char>(bytes[3 + (j + 1) * stride - 1]) >>
                           (2 * (samples.size() % 4))) == 0, "Output padding not zero");
            }
            check(read(result.prefix + ".qc_report.txt").find("samples_pass\t" +
                  std::to_string(samples.size())) != std::string::npos, "Missing QC report counts");
        };
        const std::vector<int> all{0, 1, 2, 3, 4, 5};
        QCOptions options;
        check(options.geno == 1 && options.mind == 1 && options.maf == .01 && options.hwe == 0,
              "API defaults");
        verify(QCFilter::preprocess(input, path("default.qc"), options), all, {0, 1, 3, 4});
        options.maf = 0;
        verify(QCFilter::preprocess(input, path("raw.qc"), options), all, {0, 1, 2, 3, 4});

        options.keep_file = list("keep", "F6 i6\nF2 same\nF1 same\nF3 i3\nF5 i5\n");
        options.remove_file = list("remove", "F2 same\n");
        options.extract_file = list("extract", "v4\nv3\nv1\nv0\n");
        options.exclude_file = list("exclude", "v3\n");
        verify(QCFilter::preprocess(input, path("pairs.qc"), options), {0, 2, 4, 5}, {0, 1, 4});
        options = QCOptions{};
        options.maf = 0;
        options.keep_file = list("collision", "F2 same\n");
        options.extract_file = list("only_v3", "v3\n");
        verify(QCFilter::preprocess(input, path("collision.qc"), options), {1}, {3});

        options = QCOptions{};
        options.maf = 0;
        options.extract_file = list("candidates", "v1\nv0\n");
        options.mind = .5;
        options.geno = .2;
        verify(QCFilter::preprocess(input, path("mind_boundary.qc"), options), {0, 2, 3, 4, 5}, {0, 1});
        options.geno = std::nextafter(.2, 0.0);
        verify(QCFilter::preprocess(input, path("geno_below.qc"), options), {0, 2, 3, 4, 5}, {0});
        options.geno = 0;
        options.mind = std::nextafter(.5, 0.0);
        verify(QCFilter::preprocess(input, path("mind_below.qc"), options), {2, 3, 4, 5}, {0, 1});
        options.mind = 0;
        verify(QCFilter::preprocess(input, path("mind_zero.qc"), options), {2, 3, 4, 5}, {0, 1});
        options.exclude_file = list("exclude_v1", "v1\n");
        verify(QCFilter::preprocess(input, path("exclude_before_mind.qc"), options), {0, 2, 3, 4, 5}, {0});
        options.exclude_file.clear();
        options.mind = 1;
        options.geno = 1.0 / 6;
        verify(QCFilter::preprocess(input, path("geno_boundary.qc"), options), all, {0});
        options.geno = std::nextafter(1.0 / 6, 0.0);
        rejects([&] { QCFilter::preprocess(input, path("geno_empty"), options); }, "geno zero variants");

        options = QCOptions{};
        options.keep_file = list("two_samples", "F1 same\nF3 i3\n");
        options.extract_file = list("only_v0", "v0\n");
        options.maf = .25;
        verify(QCFilter::preprocess(input, path("maf_boundary.qc"), options), {0, 2}, {0});
        options.maf = std::nextafter(.25, 1.0);
        rejects([&] { QCFilter::preprocess(input, path("maf_empty"), options); }, "MAF on surviving samples");
        options = QCOptions{};
        options.extract_file = path("only_v3");
        options.maf = .5;
        const double exact_hwe = cosmic::calculateHWE(3, 0, 3);
        check(std::abs(exact_hwe - 5.0 / 231) < 1e-14, "Exact HWE fixture value");
        options.hwe = exact_hwe;
        verify(QCFilter::preprocess(input, path("hwe_boundary.qc"), options), all, {3});
        options.hwe = std::nextafter(exact_hwe, 1.0);
        rejects([&] { QCFilter::preprocess(input, path("hwe_empty"), options); }, "Exact HWE inclusive boundary");
        options.hwe = 1;
        options.maf = 0;
        options.keep_file = list("hom_samples", "F1 same\nF2 same\nF3 i3\n");
        verify(QCFilter::preprocess(input, path("hwe_subset.qc"), options), {0, 1, 2}, {3});

        for (auto field : {&QCOptions::geno, &QCOptions::mind, &QCOptions::maf, &QCOptions::hwe}) {
            for (double invalid : {-0.01, std::numeric_limits<double>::infinity(),
                                   -std::numeric_limits<double>::infinity(),
                                   std::numeric_limits<double>::quiet_NaN(),
                                   field == &QCOptions::maf ? .50001 : 1.00001}) {
                options = QCOptions{};
                options.*field = invalid;
                rejects([&] { QCFilter::preprocess(input, path("invalid"), options); }, "threshold range");
            }
        }
        for (auto field : {&QCOptions::keep_file, &QCOptions::remove_file,
                           &QCOptions::extract_file, &QCOptions::exclude_file}) {
            options = QCOptions{};
            options.*field = list("empty_list", " \n\t\n");
            rejects([&] { QCFilter::preprocess(input, path("empty_list_output"), options); }, "empty list");
            options.*field = path("nonexistent_list");
            rejects([&] { QCFilter::preprocess(input, path("missing_list_output"), options); }, "missing list");
        }
        for (auto field : {&QCOptions::keep_file, &QCOptions::remove_file}) {
            for (const auto& text : {"same\n", "F1 same extra\n"}) {
                options = QCOptions{};
                options.*field = list("bad_pair", text);
                rejects([&] { QCFilter::preprocess(input, path("bad_pair_output"), options); }, "pair columns");
            }
        }
        options = QCOptions{};
        options.keep_file = list("no_match", "F9 same\n");
        rejects([&] { QCFilter::preprocess(input, path("no_samples"), options); }, "keep zero matches");
        options.keep_file = path("collision");
        options.remove_file = path("collision");
        rejects([&] { QCFilter::preprocess(input, path("remove_all"), options); }, "remove wins over keep");
        options = QCOptions{};
        options.extract_file = list("no_variant", "not_present\n");
        rejects([&] { QCFilter::preprocess(input, path("no_variants"), options); }, "extract zero matches");
        options.extract_file = path("only_v0");
        options.exclude_file = path("only_v0");
        rejects([&] { QCFilter::preprocess(input, path("exclude_all"), options); }, "exclude wins");
        options.exclude_file.clear();
        options.extract_file = list("only_missing", "v5\n");
        options.mind = 0;
        rejects([&] { QCFilter::preprocess(input, path("mind_empty"), options); }, "mind zero samples");
        options.mind = 1;
        options.maf = 0;
        rejects([&] { QCFilter::preprocess(input, path("all_missing"), options); }, "no called genotypes");

        options = QCOptions{};
        rejects([&] { QCFilter::preprocess(input, input, options); }, "same prefix");
        rejects([&] { QCFilter::preprocess(input, (dir / "." / "input").string(), options); }, "same canonical prefix");
        for (const auto& suffix : {"", ".bed", ".bim", ".fam", ".qc_report.txt", ".log"}) {
            const auto prefix = path(std::string("existing") + std::to_string(std::string(suffix).size()) + suffix);
            write(prefix + suffix, "sentinel");
            rejects([&] { QCFilter::preprocess(input, prefix, options); }, "existing output prefix");
            check(read(prefix + suffix) == "sentinel", "Existing output overwritten");
        }
        fs::create_symlink(path("absent_target"), path("symlink.bed"));
        rejects([&] { QCFilter::preprocess(input, path("symlink"), options); }, "dangling output symlink");
        rejects([&] { QCFilter::preprocess(input, path("default.qc"), options); }, "second QC run");
        check(read(input + ".bed") == bed_text && read(input + ".bim") == bim_text &&
              read(input + ".fam") == fam_text, "Input overwritten");

        const auto malformed = path("malformed");
        write(malformed + ".fam", fam_text);
        write(malformed + ".bim", bim_text);
        for (const auto& bad : {std::string{}, bed_text.substr(0, 2), bed_text.substr(0, bed_text.size() - 1),
                                bed_text + 'x', std::string("bad") + bed_text.substr(3),
                                std::string("\x6c\x1b\x00", 3) + bed_text.substr(3)}) {
            write(malformed + ".bed", bad);
            rejects([&] { QCFilter::preprocess(malformed, path("bad_bed_output"), options); }, "BED validation before mmap");
            check(!fs::exists(path("bad_bed_output.bed")), "Malformed BED created output");
        }
        write(malformed + ".bed", bed_text);
        write(malformed + ".fam", "");
        rejects([&] { QCFilter::preprocess(malformed, path("empty_fam"), options); }, "zero input samples");
        write(malformed + ".fam", fam_text);
        write(malformed + ".bim", "");
        rejects([&] { QCFilter::preprocess(malformed, path("empty_bim"), options); }, "zero input variants");
        std::cout << "QC pipeline tests passed\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "QC pipeline test failed: " << e.what() << '\n';
        return 1;
    }
}
