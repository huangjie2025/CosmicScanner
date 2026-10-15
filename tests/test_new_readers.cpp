#include "hapmap_reader.h"
#include "text_geno_reader.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>

namespace {

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

void write_file(const std::filesystem::path& path, const std::string& contents) {
    std::ofstream output(path);
    if (!output) throw std::runtime_error("Cannot create test fixture: " + path.string());
    output << contents;
}

template <typename Reader>
void check_matvec_adjoint(Reader& reader) {
    Eigen::VectorXd marker_weights(reader.cols());
    Eigen::VectorXd sample_weights(reader.rows());
    for (int i = 0; i < marker_weights.size(); ++i) marker_weights[i] = 0.25 + i;
    for (int i = 0; i < sample_weights.size(); ++i) sample_weights[i] = 0.5 - 0.2 * i;
    Eigen::VectorXd zv;
    Eigen::VectorXd ztv;
    reader.multiply_Z_v(marker_weights, zv, 1);
    reader.multiply_Zt_v(sample_weights, ztv, 1);
    require(zv.size() == reader.rows(), "Z*v returned the wrong number of samples");
    require(ztv.size() == reader.cols(), "Z'*v returned the wrong number of variants");
    require(std::abs(sample_weights.dot(zv) - marker_weights.dot(ztv)) < 1e-10,
            "Z and Z' matrix products are not adjoints");
}

} // namespace

int main(int argc, char** argv) {
    try {
        const std::filesystem::path output_dir = argc > 1 ? argv[1] : "/tmp/cosmic-reader-test";
        std::filesystem::create_directories(output_dir);

        const auto hmp = output_dir / "fixture.hmp.txt";
        write_file(hmp,
            "rs# alleles chrom pos strand assembly# center protLSID assayLSID panelLSID QCcode S1 S2 S3\n"
            "rs1 A/C 1 100 + NA NA NA NA NA NA AA AC CC\n"
            "rs2 G/T 1 200 + NA NA NA NA NA NA GT GG TT\n"
            "rs3 C/G 2 300 + NA NA NA NA NA NA NN CG GG\n");

        cosmic::HapMapReader hapmap;
        hapmap.open(hmp.string());
        require(hapmap.getNumSamples() == 3, "HapMap sample count mismatch");
        require(hapmap.getNumVariants() == 3, "HapMap variant count mismatch");
        require(hapmap.getChromLabels() == std::vector<std::string>({"1", "2"}), "HapMap chromosome labels mismatch");
        std::string chrom, id, ref, alt;
        uint32_t pos = 0;
        Eigen::VectorXf dosage;
        require(hapmap.readVariant(chrom, id, pos, ref, alt, dosage), "Cannot read first HapMap variant");
        require(id == "rs1" && chrom == "1" && pos == 100 && ref == "A" && alt == "C", "HapMap metadata mismatch");
        require(dosage[0] == 0.0f && dosage[1] == 1.0f && dosage[2] == 2.0f, "HapMap dosage conversion mismatch");
        hapmap.seekVariant(2);
        require(hapmap.readVariant(chrom, id, pos, ref, alt, dosage) && std::isnan(dosage[0]), "HapMap missing call mismatch");
        check_matvec_adjoint(hapmap);

        const auto variant_matrix = output_dir / "variants_by_samples.geno";
        const auto sample_matrix = output_dir / "samples_by_variants.geno";
        const auto variant_info = output_dir / "fixture.var.info";
        write_file(variant_matrix,
            "VariantID S1 S2 S3\n"
            "rs1 0 1 2\n"
            "rs2 1 0 2\n"
            "rs3 NA 1 2\n");
        write_file(sample_matrix,
            "SampleID rs1 rs2 rs3\n"
            "S1 0 1 NA\n"
            "S2 1 0 1\n"
            "S3 2 2 2\n");
        write_file(variant_info,
            "ID CHROM POS REF ALT\n"
            "rs1 1 100 A C\n"
            "rs2 1 200 G T\n"
            "rs3 2 300 C G\n");

        cosmic::TextGenoReader by_variant;
        by_variant.open(variant_matrix.string(), cosmic::TextGenoOrientation::VariantsBySamples,
                        variant_info.string());
        require(by_variant.getNumSamples() == 3, "Variants-by-samples sample count includes header");
        require(by_variant.getNumVariants() == 3, "Variants-by-samples variant count includes header");
        require(by_variant.getSampleIds() == std::vector<std::string>({"S1", "S2", "S3"}), "Text sample IDs mismatch");
        check_matvec_adjoint(by_variant);

        cosmic::TextGenoReader by_sample;
        by_sample.open(sample_matrix.string(), cosmic::TextGenoOrientation::SamplesByVariants,
                       variant_info.string());
        require(by_sample.getNumSamples() == 3, "Samples-by-variants sample count includes header");
        require(by_sample.getNumVariants() == 3, "Samples-by-variants variant count includes header");
        check_matvec_adjoint(by_sample);

        const auto exported = output_dir / "exported.geno";
        by_variant.exportToText(exported.string(), (output_dir / "exported.var.info").string());
        cosmic::TextGenoReader round_trip;
        round_trip.open(exported.string(), cosmic::TextGenoOrientation::VariantsBySamples,
                        (output_dir / "exported.var.info").string());
        require(round_trip.getNumSamples() == 3 && round_trip.getNumVariants() == 3,
                "Text genotype round trip changed dimensions");

        const auto ragged = output_dir / "ragged.geno";
        write_file(ragged, "VariantID S1 S2\nrs1 0\n");
        bool rejected = false;
        try {
            cosmic::TextGenoReader invalid;
            invalid.open(ragged.string());
        } catch (const std::runtime_error&) {
            rejected = true;
        }
        require(rejected, "Ragged text genotype row was not rejected");

        std::cout << "PASS HapMap/text genotype readers\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "FAIL HapMap/text genotype readers: " << error.what() << '\n';
        return 1;
    }
}
