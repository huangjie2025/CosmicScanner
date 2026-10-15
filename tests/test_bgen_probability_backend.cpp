// Phase E2.4: BGEN genotype probability backend unit tests.
//
// Covers: probability decode, range/sum validation, dosage parity against the
// legacy production dosage path, allele-frequency contract, posterior-mean
// Vitezica hardcall-limit equivalence, missing-data semantics, and
// unsupported-record fail-fast guards.
//
// BGEN dominance / additive-dominance remain NO-GO; nothing here is wired to
// a production dominance or GRM path.

#include <cmath>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <zlib.h>

#include "bgen_reader.h"
#include "bgen_dominance_encoding.h"

namespace {

int g_checks = 0;
int g_failures = 0;

void check(bool cond, const std::string& what) {
    ++g_checks;
    if (cond) {
        std::cout << "  ok: " << what << "\n";
    } else {
        ++g_failures;
        std::cout << "  FAIL: " << what << "\n";
    }
    std::cout.flush();
}

void checkClose(double a, double b, double tol, const std::string& what) {
    const bool ok = std::fabs(a - b) <= tol;
    ++g_checks;
    if (ok) {
        std::cout << "  ok: " << what << " (" << a << " vs " << b << ")\n";
    } else {
        ++g_failures;
        std::cout << "  FAIL: " << what << " (" << a << " vs " << b
                  << ", |diff|=" << std::fabs(a - b) << " > " << tol << ")\n";
    }
    std::cout.flush();
}

template <typename T>
void put(std::ofstream& out, T value) {
    out.write(reinterpret_cast<const char*>(&value), sizeof(T));
}

void putStr2(std::ofstream& out, const std::string& s) {
    put<uint16_t>(out, static_cast<uint16_t>(s.size()));
    out.write(s.data(), static_cast<std::streamsize>(s.size()));
}

void pushBytes16(std::vector<char>& v, uint16_t x) {
    v.push_back(static_cast<char>(x & 0xFF));
    v.push_back(static_cast<char>((x >> 8) & 0xFF));
}

void pushU32(std::vector<char>& v, uint32_t x) {
    for (int i = 0; i < 4; ++i) v.push_back(static_cast<char>((x >> (8 * i)) & 0xFF));
}

size_t sampleBlockLen(const std::vector<std::string>& sids) {
    // 4 bytes for the block length field, 4 bytes for the sample count, then
    // a uint16-length-prefixed identifier per sample.
    size_t s = 8;
    for (const auto& x : sids) s += 2 + x.size();
    return s;
}

// A single biallelic diploid unphased Layout-2 record.
struct ProbTriplet {
    double a0a0;  // stored probability 0
    double a0a1;  // stored probability 1
    bool missing;
};

// Uncompressed-but-zlib-framed genotype payload with `bits` probability depth.
std::vector<char> zlibWrap(const std::vector<char>& raw) {
    uLongf bound = compressBound(static_cast<uLong>(raw.size()));
    std::vector<char> comp(bound);
    if (compress2(reinterpret_cast<Bytef*>(comp.data()), &bound,
                  reinterpret_cast<const Bytef*>(raw.data()),
                  static_cast<uLong>(raw.size()), 6) != Z_OK) {
        throw std::runtime_error("zlib compress failed");
    }
    comp.resize(bound);
    return comp;
}

// Minimal BGEN v1.2 Layout-2 zlib writer. `bits` selects 8 or 16 bit
// probabilities. Each entry is one physical record for one variant.
void writeBgen(const std::string& path,
               const std::vector<std::string>& sample_ids,
               const std::vector<std::string>& variant_ids,
               const std::vector<std::vector<ProbTriplet>>& records,
               uint8_t bits = 8,
               bool phased = false,
               uint8_t ploidy_override = 0) {
    std::ofstream out(path, std::ios::binary);
    if (!out) throw std::runtime_error("cannot write " + path);

    const uint32_t M = static_cast<uint32_t>(variant_ids.size());
    const uint32_t N = static_cast<uint32_t>(sample_ids.size());

    // Header block: 4-byte offset, then header of header_length bytes:
    // [4 header_length][4 M][4 N][4 magic][4 flags]
    const uint32_t header_length = 20;
    const uint32_t flags = 1u | (2u << 2) | (1u << 31);  // zlib, Layout 2, sample identifiers
    // The leading BGEN offset field is relative to byte 4 of the file, and the
    // reader seeks to (offset + 4), so store header_length + sample block size.
    const uint32_t variant_data_offset =
        header_length + static_cast<uint32_t>(sampleBlockLen(sample_ids));
    const uint32_t sample_block_len = variant_data_offset - header_length;

    put<uint32_t>(out, variant_data_offset);
    put<uint32_t>(out, header_length);
    put<uint32_t>(out, M);
    put<uint32_t>(out, N);
    out.write("bgen", 4);
    put<uint32_t>(out, flags);
    put<uint32_t>(out, sample_block_len);
    put<uint32_t>(out, N);
    for (const auto& sid : sample_ids) putStr2(out, sid);

    for (size_t j = 0; j < M; ++j) {
        const std::string& id = variant_ids[j];
        putStr2(out, id);   // id
        putStr2(out, id);   // rsid
        putStr2(out, "1");  // chrom
        put<uint32_t>(out, static_cast<uint32_t>((j + 1) * 100));
        put<uint16_t>(out, 2);  // n_alleles
        put<uint32_t>(out, 1);
        out.write("A", 1);
        put<uint32_t>(out, 1);
        out.write("C", 1);

        std::vector<ProbTriplet> rec = records[j];
        std::vector<char> body;
        {
            const uint32_t NN = static_cast<uint32_t>(rec.size());
            pushU32(body, NN);
            pushBytes16(body, 2);
            const uint8_t pl = ploidy_override ? ploidy_override : 2;
            body.push_back(static_cast<char>(pl));
            body.push_back(static_cast<char>(pl));
            for (const auto& t : rec) {
                body.push_back(static_cast<char>(t.missing ? (0x80 | (pl & 0x3F)) : (pl & 0x3F)));
            }
            body.push_back(static_cast<char>(phased ? 1 : 0));
            body.push_back(static_cast<char>(bits));
            const uint32_t denom = (1u << bits) - 1u;
            for (const auto& t : rec) {
                const uint32_t a = static_cast<uint32_t>(std::lround(t.a0a0 * denom));
                const uint32_t b = static_cast<uint32_t>(std::lround(t.a0a1 * denom));
                if (bits == 8) {
                    body.push_back(static_cast<char>(a & 0xFF));
                    body.push_back(static_cast<char>(b & 0xFF));
                } else {
                    pushBytes16(body, static_cast<uint16_t>(a));
                    pushBytes16(body, static_cast<uint16_t>(b));
                }
            }
        }
        // C = 4 (D) + compressed payload length; D = UNCOMPRESSED length.
        const std::vector<char> comp = zlibWrap(body);
        put<uint32_t>(out, static_cast<uint32_t>(comp.size() + 4));
        put<uint32_t>(out, static_cast<uint32_t>(body.size()));
        out.write(comp.data(), static_cast<std::streamsize>(comp.size()));
    }
}

} // namespace

int main(int argc, char** argv) {
    const std::string dir = (argc > 1) ? argv[1] : "bgen_prob_backend";
    std::filesystem::create_directories(dir);

    std::vector<std::string> sids;
    for (int i = 0; i < 6; ++i) sids.push_back("S" + std::to_string(i + 1));

    // --- J.1 probability decode + D dosage parity + H missing ---
    // Deterministic mix: confident genotypes, a maximally ambiguous genotype,
    // probabilities that round to dosage 1 from two different (p0,p1) shapes,
    // and two hard-missing samples.
    std::vector<ProbTriplet> rec1 = {
        {1.0, 0.0, false},   // dosage 2, certain
        {0.0, 1.0, false},   // dosage 1, certain het
        {0.0, 0.0, false},   // dosage 0, certain
        {0.5, 0.0, false},   // dosage 1, ambiguous
        {0.25, 0.5, false},  // dosage 1, different ambiguous
        {0.0, 1.0, true},    // missing
    };
    std::vector<ProbTriplet> rec2 = {
        {0.0, 1.0, false},   // all het
        {0.0, 1.0, false},
        {0.0, 1.0, false},
        {1.0, 0.0, false},   // hom A
        {1.0, 0.0, true},    // missing
        {0.0, 0.0, false},   // hom C
    };

    const std::string f1 = dir + "/basic.bgen";
    writeBgen(f1, sids, {"rs1", "rs2"}, {rec1, rec2}, 8);

    {
        cosmic::BgenReader r;
        r.open(f1);
        r.initialize();
        check(r.getNumSamples() == 6, "sample count read from header");
        check(r.getNumVariants() == 2, "variant count read from header");
        check(r.getSampleIds().size() == 6, "sample identifiers decoded");

        std::string chrom, rsid, a1, a2;
        uint32_t pos;
        cosmic::BgenGenotypeProbs probs;

        // --- Legacy dosage path first (production behavior) ---
        Eigen::VectorXf legacy;
        r.reset();
        check(r.readVariant(chrom, rsid, pos, a1, a2, legacy), "legacy readVariant rs1");
        Eigen::VectorXf legacy_copy = legacy;
        std::string legacy_chrom = chrom, legacy_rsid = rsid;
        uint32_t legacy_pos = pos;
        std::string legacy_a1 = a1, legacy_a2 = a2;

        // --- Probability path ---
        r.reset();
        check(r.readVariantProbs(chrom, rsid, pos, a1, a2, probs), "probability readVariantProbs rs1");
        check(probs.size() == 6, "probability block size matches sample count");
        check(probs.n_alleles == 2, "n_alleles reported as 2");
        check(probs.ploidy == 2, "ploidy reported as 2");
        check(probs.prob_bits == 8, "probability bit depth reported as 8");
        check(!probs.phased, "record reported unphased");

        // Variant/sample ordering and allele labels must match legacy exactly.
        check(legacy_rsid == rsid && legacy_chrom == chrom && legacy_pos == pos,
              "probability API preserves variant ordering/identity");
        check(legacy_a1 == a1 && legacy_a2 == a2, "probability API preserves allele labels");

        // --- J.1 decode: confident genotypes ---
        checkClose(probs.p2[0], 1.0, 1e-6, "p2 = P(A0A0) decoded");
        checkClose(probs.p1[1], 1.0, 1e-6, "p1 = P(het) decoded");
        checkClose(probs.p0[2], 1.0, 1e-6, "p0 = P(A1A1) decoded");
        // Ambiguous genotypes: the third class is the omitted remainder.
        // 8-bit probabilities are quantized to steps of 1/255, so decoding
        // tolerances must be quantization-aware, not machine precision.
        const double Q8 = 1.0 / 255.0 + 1e-7;
        checkClose(probs.p2[3], 0.5, Q8, "ambiguous genotype p2 decoded");
        checkClose(probs.p0[3], 0.5, Q8, "ambiguous genotype p0 = remainder");
        checkClose(probs.p1[4], 0.5, Q8, "second ambiguous shape p1 decoded");
        checkClose(probs.p0[4], 0.25, Q8, "second ambiguous shape p0 decoded");
        // Quantization-aware sum validation still holds.
        for (int i = 0; i < 6; ++i) {
            if (probs.missing[i] != 0.0f) continue;
            check(std::fabs((probs.p0[i] + probs.p1[i] + probs.p2[i]) - 1.0) <= Q8,
                  "decoded probabilities sum to 1 within 8-bit quantization, sample " +
                  std::to_string(i));
        }

        // --- D. dosage parity: probability API vs legacy production path ---
        for (int i = 0; i < 6; ++i) {
            if (std::isnan(legacy_copy[i])) {
                check(std::isnan(probs.dosage[i]),
                      "dosage NaN parity for missing sample index " + std::to_string(i));
            } else {
                check(probs.dosage[i] == legacy_copy[i],
                      "bit-exact dosage parity sample " + std::to_string(i) +
                      " (prob=" + std::to_string(probs.dosage[i]) +
                      " legacy=" + std::to_string(legacy_copy[i]) + ")");
            }
        }
        // E[G] = p1 + 2 p2 identity.
        for (int i = 0; i < 6; ++i) {
            if (probs.missing[i] != 0.0f) continue;
            checkClose(probs.dosage[i], probs.p1[i] + 2.0 * probs.p2[i], 1e-6,
                       "E[G] = P(AB) + 2 P(BB) identity sample " + std::to_string(i));
        }
        // Dosage near 0/1/2 covered above; assert the extremes explicitly.
        checkClose(probs.dosage[0], 2.0, 1e-6, "dosage near 2");
        checkClose(probs.dosage[1], 1.0, 1e-6, "dosage near 1");
        checkClose(probs.dosage[2], 0.0, 1e-6, "dosage near 0");

        // --- H. missing semantics ---
        check(probs.missing[5] == 1.0f, "missing flag set for hard-missing sample");
        check(probs.missing[0] == 0.0f, "missing flag clear for observed sample");
        check(probs.anyMissing(), "anyMissing() reports the variant has missing data");
        check(probs.observedCount() == 5, "observed count excludes hard-missing samples");
        check(probs.p0[5] == 0.0f && probs.p1[5] == 0.0f && probs.p2[5] == 0.0f,
              "missing sample probability vectors are zeroed");

        // --- G. allele frequency contract ---
        // rec1 observed dosages are the quantized values of {2,1,0,1,1}.
        // p = mean(E[G])/2 over non-missing samples must equal the dosage-derived
        // frequency exactly, and must be stable under quantization.
        checkClose(probs.alleleFrequency(), probs.meanDosage() / 2.0, 1e-12,
                   "p = mean(E[G])/2 over non-missing samples");
        check(std::fabs(probs.alleleFrequency() - 0.5) <= Q8,
              "p close to the exact 0.5 design value under 8-bit quantization");
        check(probs.observedCount() == 5, "observed count used by allele frequency");

        // --- F. posterior-mean Vitezica: one-hot hardcall limit ---
        {
            const double p = probs.alleleFrequency();
            const Eigen::VectorXd mu =
                cosmic::bgen::posteriorMeanVitezicaVector(probs.p0, probs.p1, probs.p2,
                                                          probs.missing, p);
            // Samples 0,1,2 are one-hot; missing -> 0.
            checkClose(mu[0], cosmic::bgen::vitezicaCodeHardcall(2, p), 1e-6,
                       "posterior-mean Vitezica reduces to hardcall code G=2");
            checkClose(mu[1], cosmic::bgen::vitezicaCodeHardcall(1, p), 1e-6,
                       "posterior-mean Vitezica reduces to hardcall code G=1");
            checkClose(mu[2], cosmic::bgen::vitezicaCodeHardcall(0, p), 1e-6,
                       "posterior-mean Vitezica reduces to hardcall code G=0");
            check(mu[5] == 0.0, "posterior-mean Vitezica sets missing to 0");
            // Ambiguous sample 3: p0=p2=0.5, p1=0 -> mean of codes G=0 and G=2.
            const double expect = 0.5 * cosmic::bgen::vitezicaCodeHardcall(0, p) +
                                  0.5 * cosmic::bgen::vitezicaCodeHardcall(2, p);
            checkClose(mu[3], expect, Q8,
                       "posterior-mean Vitezica averages codes under ambiguity");
            // No posterior-variance diagonal correction: decomposition identity.
            const double q = 1.0 - p;
            const double h = probs.p1[3], x = probs.p1[3] + 2.0 * probs.p2[3];
            checkClose(mu[3], (h - 2.0 * p * q) + (2.0 * p - 1.0) * (x - 2.0 * p), 1e-6,
                       "posterior-mean Vitezica satisfies frozen decomposition");
            // Denominator: 4 p^2 q^2.
            checkClose(cosmic::bgen::vitezicaDenominator(p), 4.0 * p * p * q * q, 1e-12,
                       "Vitezica denominator is 4 p^2 q^2");
            // Centered heterozygosity marker contract.
            const Eigen::VectorXd xd = cosmic::bgen::centeredHeterozygosity(probs.p1, probs.missing);
            checkClose(xd.sum(), 0.0, 1e-6, "Xd = P(AB) - mean(P(AB)) sums to zero");
            check(xd[5] == 0.0, "Xd assigns 0 to missing samples after mean imputation");
        }

        // --- G. second variant, all-het + HWD shape ---
        check(r.readVariantProbs(chrom, rsid, pos, a1, a2, probs), "probability readVariantProbs rs2");
        check(rsid == "rs2", "second variant identity in order");
        check(probs.observedCount() == 5, "rs2 observed count");
        // dosage: 1,1,1,2,NaN,0 -> 5 observed
        check(probs.observedCount() == 5, "rs2 observed count");
        checkClose(probs.alleleFrequency(), probs.meanDosage() / 2.0, 1e-12,
                   "rs2 allele frequency");
    }

    // --- J.10 16-bit probability depth ---
    {
        const std::string f16 = dir + "/b16.bgen";
        writeBgen(f16, sids, {"rs1", "rs2"}, {rec1, rec2}, 16);
        cosmic::BgenReader r;
        r.open(f16);
        r.initialize();
        std::string chrom, rsid, a1, a2;
        uint32_t pos;
        cosmic::BgenGenotypeProbs probs;
        check(r.readVariantProbs(chrom, rsid, pos, a1, a2, probs), "16-bit readVariantProbs");
        check(probs.prob_bits == 16, "16-bit depth reported");
        checkClose(probs.p1[1], 1.0, 1e-9, "16-bit one-hot het decoded exactly");
        checkClose(probs.p0[4], 0.25, 1e-4, "16-bit fractional probability decoded");
    }

    // --- I. unsupported record guards ---
    {
        // phased
        const std::string fp = dir + "/phased.bgen";
        writeBgen(fp, sids, {"rs1", "rs2"}, {rec1, rec2}, 8, /*phased=*/true);
        cosmic::BgenReader r;
        r.open(fp);
        r.initialize();
        std::string chrom, rsid, a1, a2;
        uint32_t pos;
        cosmic::BgenGenotypeProbs probs;
        bool threw = false;
        std::string msg;
        try {
            r.readVariantProbs(chrom, rsid, pos, a1, a2, probs);
        } catch (const std::exception& e) {
            threw = true;
            msg = e.what();
        }
        check(threw && msg.find("Phased") != std::string::npos,
              "phased record is REJECTED by the probability API");

        // non-diploid (haploid Pmin/Pmax = 1)
        const std::string fh = dir + "/haploid.bgen";
        writeBgen(fh, sids, {"rs1", "rs2"}, {rec1, rec2}, 8, false, /*ploidy_override=*/1);
        cosmic::BgenReader rh;
        rh.open(fh);
        rh.initialize();
        threw = false;
        msg.clear();
        try {
            rh.readVariantProbs(chrom, rsid, pos, a1, a2, probs);
        } catch (const std::exception& e) {
            threw = true;
            msg = e.what();
        }
        check(threw && msg.find("diploid") != std::string::npos,
              "non-diploid record is REJECTED by the probability API");

        // polyploid Pmax = 3
        const std::string fp3 = dir + "/polyploid.bgen";
        writeBgen(fp3, sids, {"rs1", "rs2"}, {rec1, rec2}, 8, false, 3);
        cosmic::BgenReader rp;
        rp.open(fp3);
        rp.initialize();
        threw = false;
        try {
            rp.readVariantProbs(chrom, rsid, pos, a1, a2, probs);
        } catch (const std::exception&) {
            threw = true;
        }
        check(threw, "polyploid record is REJECTED by the probability API");
    }

    // --- I. multiallelic: legacy silently skips, probability API must reject ---
    {
        const std::string fm = dir + "/multi.bgen";
        // Write a biallelic record, then a multiallelic one, then biallelic.
        std::ofstream out(fm, std::ios::binary);
        const uint32_t M = 3, N = 6;
        const uint32_t header_length = 20;
        const uint32_t flags = 1u | (2u << 2) | (1u << 31);
        const uint32_t sample_region = static_cast<uint32_t>(sampleBlockLen(sids));
        const uint32_t variant_data_offset = header_length + sample_region;
        const uint32_t sample_block_len = sample_region - 4;
        put<uint32_t>(out, variant_data_offset);
        put<uint32_t>(out, header_length);
        put<uint32_t>(out, M);
        put<uint32_t>(out, N);
        out.write("bgen", 4);
        put<uint32_t>(out, flags);
        put<uint32_t>(out, static_cast<uint32_t>(sample_block_len));
        put<uint32_t>(out, N);
        for (const auto& s : sids) putStr2(out, s);

        auto emit = [&](const std::string& id, uint16_t n_alleles) {
            putStr2(out, id);
            putStr2(out, id);
            putStr2(out, "1");
            put<uint32_t>(out, 1234);
            put<uint16_t>(out, n_alleles);
            for (uint16_t a = 0; a < n_alleles; ++a) {
                put<uint32_t>(out, 1);
                out.write(a == 0 ? "A" : (a == 1 ? "C" : "G"), 1);
            }
            std::vector<char> body;
            pushU32(body, N);
            pushBytes16(body, static_cast<uint16_t>(n_alleles));
            body.push_back(2);
            body.push_back(2);
            for (uint32_t i = 0; i < N; ++i) body.push_back(2);
            body.push_back(0);
            body.push_back(8);
            const uint32_t nprob = (n_alleles == 2) ? 2 : 3;
            for (uint32_t i = 0; i < N * nprob; ++i) body.push_back(0);
            const auto comp = zlibWrap(body);
            put<uint32_t>(out, static_cast<uint32_t>(comp.size() + 4));
            put<uint32_t>(out, static_cast<uint32_t>(body.size()));
            out.write(comp.data(), static_cast<std::streamsize>(comp.size()));
        };
        emit("rs1", 2);
        emit("rs2multi", 3);
        emit("rs3", 2);
        out.close();

        cosmic::BgenReader r;
        r.open(fm);
        r.initialize();
        std::string chrom, rsid, a1, a2;
        uint32_t pos;
        cosmic::BgenGenotypeProbs probs;
        // First physical record is biallelic and must succeed.
        bool first_ok = r.readVariantProbs(chrom, rsid, pos, a1, a2, probs);
        check(first_ok && rsid == "rs1", "probability API reads the leading biallelic record");
        // The next physical record is multiallelic and must be rejected.
        bool threw = false;
        std::string msg;
        try {
            r.readVariantProbs(chrom, rsid, pos, a1, a2, probs);
        } catch (const std::exception& e) {
            threw = true;
            msg = e.what();
        }
        check(threw && msg.find("multiallelic") != std::string::npos,
              "multiallelic record is REJECTED (no silent skip) by the probability API");

        // Legacy path still skips it silently -> unchanged additive behavior.
        r.reset();
        Eigen::VectorXf dosage;
        check(r.readVariant(chrom, rsid, pos, a1, a2, dosage), "legacy readVariant biallelic rs1");
        check(r.readVariant(chrom, rsid, pos, a1, a2, dosage), "legacy readVariant after silent multiallelic skip");
        check(rsid == "rs3", "legacy path silently skipped the multiallelic record");
    }

    // --- C. probability sum validation ---
    {
        // Craft probabilities whose stored values do not sum to 1 within the
        // default quantization tolerance: 8-bit max sum is 1.0 when both stored
        // bytes are 255 -> p2=1, p1=1, p0=-1 => out of range AND sum error.
        const std::string fb = dir + "/badprob.bgen";
        std::vector<ProbTriplet> bad = {
            {1.0, 1.0, false},  // p0 = -1 : invalid
            {0.0, 1.0, false},
            {0.0, 0.0, false},
            {0.0, 1.0, false},
            {0.0, 0.0, false},
            {0.0, 1.0, false},
        };
        writeBgen(fb, sids, {"rs1", "rs2"}, {bad, bad}, 8);
        cosmic::BgenReader r;
        r.open(fb);
        r.initialize();
        std::string chrom, rsid, a1, a2;
        uint32_t pos;
        cosmic::BgenGenotypeProbs probs;
        bool threw = false;
        std::string msg;
        try {
            r.readVariantProbs(chrom, rsid, pos, a1, a2, probs);
        } catch (const std::exception& e) {
            threw = true;
            msg = e.what();
        }
        check(threw, "out-of-range probability is REJECTED");
        check(msg.find("[0,1]") != std::string::npos || msg.find("sum") != std::string::npos,
              "rejection message identifies range or sum violation");

        // Layout 2 omits the last genotype probability, so p0 is derived as the
        // remainder 1 - s0 - s1 and the sum-to-one invariant holds by
        // construction (up to floating point). The effective guard is therefore
        // the range check on the derived remainder, which rejects records whose
        // stored probabilities exceed 1 in total. A valid quantized record must
        // be accepted by the default tolerance.
        const std::string fq = dir + "/tol.bgen";
        std::vector<ProbTriplet> valid = {
            {0.25, 0.5, false},  // p0 = 0.25, p1 = 0.5, p2 = 0.25
            {0.0, 1.0, false},
            {0.0, 0.0, false},
            {0.0, 1.0, false},
            {0.0, 0.0, false},
            {0.0, 1.0, false},
        };
        writeBgen(fq, sids, {"rs1", "rs2"}, {valid, valid}, 8);
        cosmic::BgenReader rq;
        rq.open(fq);
        rq.initialize();
        bool threw_default = false;
        try {
            rq.reset();
            rq.readVariantProbs(chrom, rsid, pos, a1, a2, probs);
        } catch (const std::exception&) {
            threw_default = true;
        }
        check(!threw_default, "default tolerance accepts valid quantized probabilities");

        // Derived-remainder range guard: s0 + s1 > 1 must be rejected.
        const std::string fn = dir + "/negremainder.bgen";
        std::vector<ProbTriplet> over = {
            {0.75, 0.5, false},  // 0.75 + 0.5 > 1 -> p0 < 0
            {0.0, 1.0, false},
            {0.0, 0.0, false},
            {0.0, 1.0, false},
            {0.0, 0.0, false},
            {0.0, 1.0, false},
        };
        writeBgen(fn, sids, {"rs1", "rs2"}, {over, over}, 8);
        cosmic::BgenReader rn;
        rn.open(fn);
        rn.initialize();
        bool threw_range = false;
        msg.clear();
        try {
            rn.reset();
            rn.readVariantProbs(chrom, rsid, pos, a1, a2, probs);
        } catch (const std::exception& e) {
            threw_range = true;
            msg = e.what();
        }
        check(threw_range && msg.find("[0,1]") != std::string::npos,
              "record whose stored probabilities exceed 1 in total is REJECTED");
    }

    // BGEN dominance production wiring and its AD rejection are exercised by
    // the cross-format, QC, and production gate integration tests.

    std::cout << "\nBGEN probability backend unit tests: " << (g_checks - g_failures) << "/"
              << g_checks << " checks passed\n";
    if (g_failures > 0) {
        std::cout << "BGEN probability backend unit tests FAILED (" << g_failures << ")\n";
        return 1;
    }
    std::cout << "BGEN probability backend unit tests PASSED.\n";
    return 0;
}
