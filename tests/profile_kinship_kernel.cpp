#include "kinship_operator.h"
#include <chrono>
#include <iostream>
#include <iomanip>
#include <random>
#include <sys/resource.h>

using Clock = std::chrono::steady_clock;
double seconds(Clock::time_point start) {
    return std::chrono::duration<double>(Clock::now() - start).count();
}

// Independent probes, NOT an additive decomposition of the fused old loops.
int main(int argc, char** argv) {
    if (argc != 4) {
        std::cerr << "usage: profile_kinship_kernel BED_PREFIX WIDTH THREADS\n";
        return 1;
    }
    cosmic::PlinkReader reader;
    reader.load(argv[1]);
    const int width = std::stoi(argv[2]), threads = std::stoi(argv[3]);
    if (width < 1 || threads < 1) return 1;
    const int n = reader.getNumSamples(), m = reader.getNumSnps();
    const auto* bed = reader.getData() + 3;
    const size_t bytes = reader.getBytesPerSnp() * m;
    rusage before{}, after{};
    getrusage(RUSAGE_SELF, &before);
    auto start = Clock::now();
    uint64_t checksum = 0;
    for (size_t i = 0; i < bytes; i += 4096) checksum += bed[i];
    const double page_touch = seconds(start);
    getrusage(RUSAGE_SELF, &after);
    const long major = after.ru_majflt - before.ru_majflt;
    const long minor = after.ru_minflt - before.ru_minflt;
    start = Clock::now();
    cosmic::PlinkKinshipOperator op(&reader);
    const double stats_seconds = seconds(start);
    start = Clock::now();
    double decode_checksum = 0;
    constexpr double dosage[] = {2, 0, 1, 0};
    Eigen::MatrixXd block(n, std::min(256, m));
    for (int first = 0; first < m; first += block.cols()) {
        const int count = std::min<int>(block.cols(), m - first);
        for (int j = 0; j < count; ++j) {
            const auto* col = bed + size_t(first + j) * reader.getBytesPerSnp();
            const double mean = reader.getMeans()[first+j];
            const double inv = reader.getInvStds()[first+j];
            for (int i = 0; i < n; ++i) {
                const int code = (col[i/4] >> (2*(i%4))) & 3;
                block(i,j) = code == 1 ? 0 : (dosage[code]-mean)*inv;
            }
        }
        decode_checksum += block.leftCols(count).squaredNorm();
    }
    const double decode_seconds = seconds(start);
    std::mt19937_64 rng(42);
    std::normal_distribution<double> normal;
    Eigen::MatrixXd omega(n, width);
    for (int j=0; j<width; ++j)
        for (int i=0; i<n; ++i) omega(i,j)=normal(rng);
    Eigen::MatrixXd result;
    getrusage(RUSAGE_SELF, &before);
    start = Clock::now();
    op.apply(omega, result, threads);
    const double kernel_seconds = seconds(start);
    getrusage(RUSAGE_SELF, &after);
    cosmic::PlinkKinshipOperator::BlockProfile profile;
    Eigen::MatrixXd profiled;
    op.apply_block(omega, profiled, threads, &profile);
    Eigen::VectorXd v = omega.col(0), serial, parallel, repeat;
    op.apply(v, serial, 1);
    op.apply(v, parallel, threads);
    op.apply(v, repeat, threads);
    std::cout << std::setprecision(15)
              << "N=" << n << " M=" << m << " WIDTH=" << width << " THREADS=" << threads
              << "\nPAGE_TOUCH_SECONDS=" << page_touch << " MAJOR_FAULTS=" << major
              << " MINOR_FAULTS=" << minor << " CHECKSUM=" << checksum
              << "\nSTATS_SECONDS=" << stats_seconds
              << "\nDECODE_ONLY_SERIAL_SECONDS=" << decode_seconds
              << " DECODE_CHECKSUM=" << decode_checksum
              << "\nKERNEL_SECONDS=" << kernel_seconds
              << " KERNEL_MAJOR_FAULTS=" << after.ru_majflt-before.ru_majflt
              << " RESULT_NORM=" << result.norm()
              << "\nVECTOR_PARALLEL_REL_ERROR=" << (parallel-serial).norm()/serial.norm()
              << " VECTOR_REPEAT_REL_ERROR=" << (repeat-parallel).norm()/serial.norm()
              << "\nBATCH_BLOCKS=" << profile.blocks
              << " PACKED_BYTES=" << profile.packed_bytes
              << " DECODED_GENOTYPES=" << profile.decoded_genotypes
              << "\nBATCH_PACKED_READ_SECONDS=" << profile.packed_read_seconds
              << " BATCH_DECODE_SECONDS=" << profile.decode_seconds
              << " BATCH_ZT_OMEGA_SECONDS=" << profile.transpose_product_seconds
              << " BATCH_Z_PROJECTED_SECONDS=" << profile.forward_product_seconds
              << "\nBATCH_FIRST_COLUMN_REL_ERROR=" << (result.col(0)-serial).norm()/serial.norm()
              << " BATCH_REPEAT_REL_ERROR=" << (result-profiled).norm()/result.norm()
              << "\n";
}
