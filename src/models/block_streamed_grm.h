#pragma once
#include <Eigen/Dense>
#include <atomic>
#include <functional>
#include <memory>
#include <string>
#include <unordered_set>
#include <vector>
#include "plink_reader.h"

namespace cosmic {

// ---------------------------------------------------------------------------
// Gcn Block Iterator - streaming interface for Gcn blocks.
// Allows iterating over Gcn in blocks without materializing the full matrix.
// ---------------------------------------------------------------------------
struct GcnBlock {
    Eigen::MatrixXd block;      // n_core x block_size
    std::vector<int> noncore_indices;  // global noncore indices for this block
    int start_col = 0;          // starting column in the full Gcn
    int end_col = 0;            // ending column (exclusive)
};

// E3-D3 stream-audit counters: measure the repeated work of the two streaming
// passes (computeStats + computeBlocks) independently of the downstream APY /
// PEV consumers.  A "pass" is one sweep over the SNP array reading genotype
// bytes (full pass = all reader samples; genotyped pass = genotyped rows only).
struct StreamAudit {
    long long compute_blocks_calls = 0;
    double compute_stats_time_s = 0.0;
    double compute_blocks_time_s = 0.0;
    long long full_passes = 0;        // sweeps over all samples
    long long genotyped_passes = 0;   // sweeps over genotyped rows only
    long long snp_bytes_read = 0;     // genotype bytes touched
    int n_samples = 0;
    int n_snps = 0;
    int n_genotyped = 0;
    int n_core = 0;
    int n_noncore = 0;
};

// Abstract interface for Gcn block iterator
class IGcnBlockIterator {
public:
    virtual ~IGcnBlockIterator() = default;
    virtual bool next(GcnBlock& out_block) = 0;
    virtual bool hasNext() const = 0;
    virtual int totalBlocks() const = 0;
    virtual int blockSize() const = 0;
};

using GcnBlockIteratorPtr = std::unique_ptr<IGcnBlockIterator>;
using GcnBlockIteratorFactory = std::function<GcnBlockIteratorPtr()>;

enum class GcnBlockView { Raw, Tuned, Final };

// Forward declaration
class BlockStreamedGrmBuilder;

// ---------------------------------------------------------------------------
// GcnBlockIteratorImpl - concrete implementation of IGcnBlockIterator
// ---------------------------------------------------------------------------
class GcnBlockIteratorImpl : public IGcnBlockIterator {
    friend class BlockStreamedGrmBuilder;
public:
    GcnBlockIteratorImpl(const BlockStreamedGrmBuilder* builder, int block_size,
                          GcnBlockView view);

     bool next(GcnBlock& out_block) override;
    bool hasNext() const override { return current_chunk_ < n_chunks_; }
    int totalBlocks() const override { return n_chunks_; }
     int blockSize() const override { return block_size_; }

 private:
    const BlockStreamedGrmBuilder* builder_ = nullptr;
    int block_size_ = 4096;
    int n_chunks_ = 0;
    int current_chunk_ = 0;
    GcnBlockView view_ = GcnBlockView::Final;
};

// ---------------------------------------------------------------------------
// Block-streamed construction of the APY-relevant blocks of the VanRaden
// genomic relationship matrix G = M M' / S (M = centered dosage, missing -> 0;
// S = sum_j 2p_j(1-p_j)) WITHOUT materializing the full n_g x n_g G.
//
// Only the blocks consumed by the APY block inverse (ApyGinv, Leite et al.
// 2024 Eqs. 5-6) are derived from the genotype stream:
//   Gcc        (n_core x n_core)   = Zc Zc' / S
//   Gcn        (n_core x n_noncore)= Zc Zn' / S
//   diag(Gnn)  (n_noncore)         = rowsum(Zn .* Zn) / S
// plus the scalar raw trace / off-diagonal sums of the genotyped block (used
// for tuneG alignment with A22 and for the diagnostics that currently scan
// the dense G).
//
// All centering conventions reproduce the production dense VanRaden path
// (GrmBuilder, algorithm 1): the per-SNP allele frequency is the GLOBAL
// frequency over ALL reader samples (missing excluded), the per-genotype value
// is dosage - 2*frequency (missing -> 0), and the SNP 2p(1-p) divisor is
// stored as float exactly as the dense builder does.
//
// Two streaming passes:
//   computeStats()   global allele stats + filters + scalar raw sums (the
//                    per-SNP frequency / validity are cached for pass 2);
//   computeBlocks()  rank-1 outer-product accumulation of Gcc / Gcn / diag(Gnn)
//                    over the genotyped rows only, reusing the cached stats.
// tuneG and blend are applied to the raw blocks (never to a full G):
//   G_tuned = b*G_raw + a*J ;  G_final = (1-alpha)*G_tuned + alpha*T ,
// T = I or the full A22 block (blend-a22), with the J (= all-ones) / T block
// contributions computed block-wise (Jcn = 1, Tcn = A22_cn, Tnn diag, ...).
// ---------------------------------------------------------------------------
class BlockStreamedGrmBuilder {
    friend class GcnBlockIteratorImpl;
public:
    struct Options {
        float min_maf = 0.01f;
        float min_hwe = 0.0f;
        std::unordered_set<std::string> extract_set;
        std::unordered_set<std::string> exclude_set;
        int n_threads = 1;
        int snp_block_size = 2048;   // stats pass SNP blocking (GrmBuilder default)
        int gcn_chunk_size = 4096;   // noncore chunk for the rank-1 updates
        bool materialize_gcn = true; // reference path; false generates Gcn from genotypes on demand
    };

    // keep_indices: genotyped sample indices (sample order), i.e. the rows of
    // the genotyped block; global allele statistics are over the full reader.
    BlockStreamedGrmBuilder(const PlinkReader& reader,
                            const std::vector<int>& keep_indices,
                            const Options& opts);

    // Pass 1: compute S_2pq, per-SNP frequency/validity, and the scalar raw
    // trace / off-diagonal sums over the genotyped block.  Required before
    // tuneG (and hence before the auto-core Lanczos) and before computeBlocks.
    void computeStats();

    // Store the tuneG scalars (G_tuned = b*G_raw + a*J).  These are produced
    // by the caller from the raw diagonal/off-diagonal means (see stats()).
    void setTune(double a, double b);

    // Pass 2: accumulate Gcc/Gcn/diag(Gnn) from the cached stats over the
    // genotyped rows.  core_idx / noncore_idx are indices into the genotyped
    // order (ascending, disjoint, covering all genotyped animals).
    void computeBlocks(const std::vector<int>& core_idx,
                       const std::vector<int>& noncore_idx);

    // Apply tuneG + blend to the raw blocks and finalize the output blocks.
    //   blend_alpha: alpha of G = (1-alpha)G_tuned + alpha T.
    //   a22_full:    full n_g x n_g A22 of the genotyped animals (only needed
    //                when blend_alpha > 0 AND blending with A22; pass nullptr
    //                for the identity blend).  When a22_full is provided the
    //                tuned-unblended blocks are retained for the effect
    //                back-solve operator (HIBLUP/preGSf90 convention).
    void finalize(double blend_alpha, const Eigen::MatrixXd* a22_full);

    // Create an iterator to stream Gcn blocks without materializing the full matrix.
    // block_size: number of noncore columns per block (default from options).
    GcnBlockIteratorPtr makeGcnIterator(int block_size = -1,
                                        GcnBlockView view = GcnBlockView::Final) const;

     // --- outputs ---
    double S2pq() const { return S_2pq_; }
    double rawDiagSum() const { return raw_diag_sum_; }
    double rawOffSum() const { return raw_off_sum_; }
    double tuneA() const { return tune_a_; }
    double tuneB() const { return tune_b_; }
    double finalDiagSum() const { return final_diag_sum_; }   // trace of final G
    double finalOffSum() const { return final_off_sum_; }     // off-diag sum of final G
    int nGenotyped() const { return n_geno_; }
    int nCore() const { return n_core_; }
    int nNoncore() const { return n_noncore_; }
    bool fullGcnStored() const { return gcn_raw_.size() != 0; }
    // Phase E3-D2D diagnostic: number of on-demand full K x Nn Gcn block
    // generations.  With the matrix-free cross operator active the APY MME /
    // PEV applies never regenerate the block, so the count stays strictly at
    // its construction-time value during the per-SNP scan.
    long long gcnRegenCalls() const { return gcn_regen_calls_.load(); }
    const std::vector<char>& snpValidMask() const { return snp_valid_; }
    const std::vector<double>& snpFrequencies() const { return snp_freq_; }

    // E3-D3 stream-audit counters (COMPUTE_BLOCKS_CALLS / COMPUTE_BLOCKS_TIME /
    // GENOTYPE_PASSES).
    const StreamAudit& streamAudit() const { return audit_; }
     std::string streamAuditSummary() const;

     // Accessors for the partition indices
    const std::vector<int>& coreIndices() const { return core_idx_; }
    const std::vector<int>& noncoreIndices() const { return noncore_idx_; }

    // Final blended blocks (for the MME APY inverse).
    const Eigen::MatrixXd& gccFinal() const { return gcc_final_; }
    const Eigen::MatrixXd& gcnFinal() const { return gcn_final_; }
    const Eigen::VectorXd& gnnDiagFinal() const { return gnn_diag_final_; }
    // Tuned-unblended blocks (for the effect APY inverse under blend-a22;
    // equal to the final blocks when no A22 blend is applied).
    const Eigen::MatrixXd& gccTuned() const { return gcc_tuned_; }
    const Eigen::MatrixXd& gcnTuned() const { return gcn_tuned_; }
    const Eigen::VectorXd& gnnDiagTuned() const { return gnn_diag_tuned_; }

    // Raw (pre-tune, pre-blend) blocks, held for validation/debugging.
    const Eigen::MatrixXd& gccRaw() const { return gcc_raw_; }
    const Eigen::MatrixXd& gcnRaw() const { return gcn_raw_; }
    const Eigen::VectorXd& gnnDiagRaw() const { return gnn_diag_raw_; }

private:
    const PlinkReader& reader_;
    std::vector<int> keep_indices_;   // genotyped sample indices
    std::vector<char> sample_is_geno_; // size n_samples_; 1 if genotyped
    Options opts_;
    int n_samples_ = 0;
    int n_snps_ = 0;
    int n_geno_ = 0;

    // Pass-1 cached statistics.
    std::vector<double> snp_freq_;    // global frequency (all samples)
    std::vector<char> snp_valid_;
    std::vector<double> snp_2pq_;     // float-per-SNP divisor, stored as double
    double S_2pq_ = 0.0;
    double raw_diag_sum_ = 0.0;       // sum_i G_raw(i,i) over genotyped block
    double raw_off_sum_ = 0.0;        // sum_{i!=j} G_raw(i,j) over genotyped block

    double tune_a_ = 0.0;
    double tune_b_ = 1.0;

    // Pass-2 raw blocks (divided by S).
    Eigen::MatrixXd gcc_raw_;
    Eigen::MatrixXd gcn_raw_;
    Eigen::VectorXd gnn_diag_raw_;

    // Final / tuned block sets (after block-wise tune + blend).
    Eigen::MatrixXd gcc_tuned_, gcn_tuned_;
    Eigen::VectorXd gnn_diag_tuned_;
    Eigen::MatrixXd gcc_final_, gcn_final_;
    Eigen::VectorXd gnn_diag_final_;
    double final_diag_sum_ = 0.0;
    double final_off_sum_ = 0.0;
    double blend_alpha_ = 0.0;
    const Eigen::MatrixXd* blend_target_ = nullptr;
    bool blocks_computed_ = false;
    bool finalized_ = false;

    int n_core_ = 0;
    int n_noncore_ = 0;
    std::vector<int> core_idx_;      // partition passed to computeBlocks
    std::vector<int> noncore_idx_;
    mutable std::atomic<long long> gcn_regen_calls_{0};
    StreamAudit audit_;

    void generateGcnBlock(int start_col, int end_col, GcnBlockView view,
                          Eigen::MatrixXd& block) const;
};

} // namespace cosmic
