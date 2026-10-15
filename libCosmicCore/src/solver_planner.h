#pragma once

#include <algorithm>
#include <cerrno>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace cosmic {

enum class PlannerSolverKind {
    Unknown,
    Auto,
    DenseDirect,
    SparseDirect,
    PCG,
    MatrixFreePCG,
    STCG,
    APY
};

enum class PlannerMatrixKind {
    Unknown,
    Dense,
    Sparse,
    MatrixFreeOperator
};

enum class PlannerSparsity {
    Unknown,
    Dense,
    Sparse
};

enum class MemoryBudgetSource {
    Unknown,
    UserSpecified,
    ConfiguredDefault
};

enum class EstimateQuality {
    Unknown,
    Known,
    Approximate,
    LowerBound
};

struct MemoryBudget {
    bool specified = false;
    std::uint64_t bytes = 0;
    MemoryBudgetSource source = MemoryBudgetSource::Unknown;
};

struct MemoryBudgetParseResult {
    bool ok = false;
    MemoryBudget budget;
    std::string error;
};

struct SolverEstimate {
    PlannerSolverKind solver = PlannerSolverKind::Unknown;
    bool estimate_available = false;
    bool feasible = false;
    bool overflow = false;
    EstimateQuality quality = EstimateQuality::Unknown;
    std::size_t estimated_memory_bytes = 0;
    std::size_t estimated_working_memory_bytes = 0;
    std::string reason;
};

struct SolverRequest {
    std::string model_type = "unknown";
    std::size_t n_samples = 0;
    std::size_t n_markers = 0;
    std::size_t system_dimension = 0;
    PlannerMatrixKind matrix_type = PlannerMatrixKind::Unknown;
    PlannerSparsity matrix_sparsity = PlannerSparsity::Unknown;
    PlannerSolverKind requested_solver = PlannerSolverKind::Auto;
    std::string requested_solver_label = "auto";
    std::string current_solver_label = "unknown";
    bool memory_budget_available = false;
    std::uint64_t memory_budget_bytes = 0;
    MemoryBudget memory_budget;
    int threads = 1;
    double tolerance = 1e-8;

    bool sparse_matrix_available = false;
    std::size_t sparse_nonzeros = 0;
    bool matrix_free_operator_available = false;
    bool genotype_operator_available = false;
    bool stcg_available = false;
    bool apy_available = false;
    std::size_t apy_core_size = 0;
    bool diagnostic_only = true;
};

struct SolverPlan {
    PlannerSolverKind requested_solver = PlannerSolverKind::Auto;
    PlannerSolverKind recommended_solver = PlannerSolverKind::Unknown;
    PlannerSolverKind fallback_solver = PlannerSolverKind::Unknown;
    std::string requested_solver_label = "auto";
    std::string current_solver_label = "unknown";
    std::string recommendation_reason;
    std::size_t estimated_memory_bytes = 0;
    bool estimate_available = false;
    EstimateQuality estimate_quality = EstimateQuality::Unknown;
    bool memory_budget_available = false;
    std::uint64_t memory_budget_bytes = 0;
    MemoryBudgetSource memory_budget_source = MemoryBudgetSource::Unknown;
    std::size_t system_dimension = 0;
    std::size_t n_samples = 0;
    std::size_t n_markers = 0;
    bool diagnostic_only = true;
    std::vector<SolverEstimate> candidates;
};

inline const char* to_string(PlannerSolverKind solver) {
    switch (solver) {
        case PlannerSolverKind::Auto: return "AUTO";
        case PlannerSolverKind::DenseDirect: return "DENSE_DIRECT";
        case PlannerSolverKind::SparseDirect: return "SPARSE_DIRECT";
        case PlannerSolverKind::PCG: return "PCG";
        case PlannerSolverKind::MatrixFreePCG: return "MATRIX_FREE_PCG";
        case PlannerSolverKind::STCG: return "STCG";
        case PlannerSolverKind::APY: return "APY";
        case PlannerSolverKind::Unknown:
        default: return "unknown";
    }
}

inline const char* to_string(MemoryBudgetSource source) {
    switch (source) {
        case MemoryBudgetSource::UserSpecified: return "USER_SPECIFIED";
        case MemoryBudgetSource::ConfiguredDefault: return "CONFIGURED_DEFAULT";
        case MemoryBudgetSource::Unknown:
        default: return "UNKNOWN";
    }
}

inline const char* to_string(EstimateQuality quality) {
    switch (quality) {
        case EstimateQuality::Known: return "KNOWN";
        case EstimateQuality::Approximate: return "APPROXIMATE";
        case EstimateQuality::LowerBound: return "LOWER_BOUND";
        case EstimateQuality::Unknown:
        default: return "UNKNOWN";
    }
}

inline std::string bytes_to_string(std::uint64_t bytes) {
    const double kib = 1024.0;
    const double mib = kib * 1024.0;
    const double gib = mib * 1024.0;
    std::ostringstream os;
    os << std::fixed << std::setprecision(2);
    if (bytes >= static_cast<std::uint64_t>(gib)) {
        os << static_cast<double>(bytes) / gib << " GiB";
    } else if (bytes >= static_cast<std::uint64_t>(mib)) {
        os << static_cast<double>(bytes) / mib << " MiB";
    } else if (bytes >= static_cast<std::uint64_t>(kib)) {
        os << static_cast<double>(bytes) / kib << " KiB";
    } else {
        os << bytes << " bytes";
    }
    return os.str();
}

inline MemoryBudget unknown_memory_budget() {
    return MemoryBudget{};
}

inline MemoryBudget make_memory_budget(std::uint64_t bytes, MemoryBudgetSource source) {
    MemoryBudget budget;
    budget.specified = source != MemoryBudgetSource::Unknown;
    budget.bytes = bytes;
    budget.source = source;
    return budget;
}

inline MemoryBudgetParseResult parse_memory_gib(const std::string& text,
                                                MemoryBudgetSource source = MemoryBudgetSource::UserSpecified) {
    MemoryBudgetParseResult result;
    if (text.empty()) {
        result.error = "memory budget is empty";
        return result;
    }

    char* end = nullptr;
    errno = 0;
    const long double value = std::strtold(text.c_str(), &end);
    if (end == text.c_str() || *end != '\0' || errno == ERANGE || !std::isfinite(static_cast<double>(value))) {
        result.error = "invalid memory budget '" + text + "'; expected a positive GiB value";
        return result;
    }
    if (value <= 0.0L) {
        result.error = "memory budget must be positive GiB";
        return result;
    }

    constexpr long double gib_bytes = 1024.0L * 1024.0L * 1024.0L;
    const long double bytes_ld = value * gib_bytes;
    const long double max_bytes = static_cast<long double>(std::numeric_limits<std::uint64_t>::max());
    if (bytes_ld > max_bytes || bytes_ld < 1.0L) {
        result.error = "memory budget overflows uint64 byte storage";
        return result;
    }
    const auto bytes = static_cast<std::uint64_t>(std::ceil(bytes_ld));

    result.ok = true;
    result.budget = make_memory_budget(bytes, source);
    return result;
}

inline PlannerSolverKind planner_solver_from_string(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });
    if (value == "auto" || value.empty()) return PlannerSolverKind::Auto;
    if (value == "dense" || value == "direct" || value == "ldlt" || value == "cholesky") return PlannerSolverKind::DenseDirect;
    if (value == "sparse" || value == "sparseldlt" || value == "sparse_ldlt") return PlannerSolverKind::SparseDirect;
    if (value == "pcg") return PlannerSolverKind::PCG;
    if (value == "matrix-free" || value == "matrix_free" || value == "matrixfree" || value == "mfpcg") return PlannerSolverKind::MatrixFreePCG;
    if (value == "stcg") return PlannerSolverKind::STCG;
    if (value == "apy") return PlannerSolverKind::APY;
    return PlannerSolverKind::Unknown;
}

namespace solver_planner_detail {

inline bool checked_mul(std::size_t a, std::size_t b, std::size_t& out) {
    if (a != 0 && b > std::numeric_limits<std::size_t>::max() / a) {
        out = std::numeric_limits<std::size_t>::max();
        return false;
    }
    out = a * b;
    return true;
}

inline bool checked_add(std::size_t a, std::size_t b, std::size_t& out) {
    if (b > std::numeric_limits<std::size_t>::max() - a) {
        out = std::numeric_limits<std::size_t>::max();
        return false;
    }
    out = a + b;
    return true;
}

inline bool checked_scale(std::size_t value, std::size_t numerator, std::size_t denominator, std::size_t& out) {
    std::size_t scaled = 0;
    if (!checked_mul(value, numerator, scaled)) {
        out = std::numeric_limits<std::size_t>::max();
        return false;
    }
    out = scaled / denominator;
    return true;
}

inline std::size_t effective_dimension(const SolverRequest& request) {
    if (request.system_dimension > 0) return request.system_dimension;
    if (request.n_samples > 0) return request.n_samples;
    return request.n_markers;
}

inline MemoryBudget effective_memory_budget(const SolverRequest& request) {
    if (request.memory_budget.specified && request.memory_budget.source != MemoryBudgetSource::Unknown) {
        return request.memory_budget;
    }
    if (request.memory_budget_available) {
        return make_memory_budget(static_cast<std::uint64_t>(request.memory_budget_bytes),
                                  MemoryBudgetSource::ConfiguredDefault);
    }
    return unknown_memory_budget();
}

inline void apply_budget(SolverEstimate& estimate, const SolverRequest& request) {
    if (!estimate.estimate_available) {
        estimate.feasible = false;
        if (estimate.reason.empty()) estimate.reason = "memory estimate is unavailable";
        return;
    }
    if (estimate.overflow) {
        estimate.feasible = false;
        estimate.reason += " Estimate overflowed size_t.";
        return;
    }
    const MemoryBudget budget = effective_memory_budget(request);
    if (budget.specified) {
        estimate.feasible = static_cast<std::uint64_t>(estimate.estimated_memory_bytes) <= budget.bytes;
        estimate.reason += estimate.feasible ? " Estimated memory is within configured budget."
                                             : " Estimated memory exceeds configured budget.";
    } else {
        estimate.feasible = true;
        estimate.reason += " Memory budget is unknown; feasibility is not proven.";
    }
}

inline SolverEstimate dense_estimate(const SolverRequest& request) {
    SolverEstimate estimate;
    estimate.solver = PlannerSolverKind::DenseDirect;
    const std::size_t d = effective_dimension(request);
    if (d == 0) {
        estimate.reason = "system dimension is unknown";
        return estimate;
    }

    // Conservative dense working-set model:
    //   dense coefficient/covariance matrix: 8*d*d bytes
    //   factorization/copy/workspace reserve: 16*d*d bytes
    //   RHS and iterative/refinement vectors: 8 vectors of d doubles
    //   final safety margin: +25%
    std::size_t d2 = 0;
    bool ok = checked_mul(d, d, d2);
    std::size_t matrix = 0;
    ok = checked_mul(d2, sizeof(double), matrix) && ok;
    std::size_t factor = 0;
    ok = checked_mul(d2, 2 * sizeof(double), factor) && ok;
    std::size_t vectors = 0;
    ok = checked_mul(d, 8 * sizeof(double), vectors) && ok;

    std::size_t subtotal = 0;
    ok = checked_add(matrix, factor, subtotal) && ok;
    ok = checked_add(subtotal, vectors, subtotal) && ok;
    std::size_t total = 0;
    ok = checked_scale(subtotal, 5, 4, total) && ok;

    estimate.estimate_available = true;
    estimate.overflow = !ok;
    estimate.quality = EstimateQuality::Known;
    estimate.estimated_working_memory_bytes = subtotal;
    estimate.estimated_memory_bytes = total;
    estimate.reason = "dense model uses matrix + factorization workspace + RHS/vector reserve + 25% safety margin.";
    apply_budget(estimate, request);
    return estimate;
}

inline SolverEstimate sparse_estimate(const SolverRequest& request) {
    SolverEstimate estimate;
    estimate.solver = PlannerSolverKind::SparseDirect;
    const std::size_t d = effective_dimension(request);
    if (!request.sparse_matrix_available && request.matrix_sparsity != PlannerSparsity::Sparse) {
        estimate.reason = "sparse matrix was not reported by caller";
        return estimate;
    }
    if (d == 0 || request.sparse_nonzeros == 0) {
        estimate.reason = "sparse nnz or dimension is unknown; sparse direct memory cannot be estimated precisely";
        return estimate;
    }

    // Sparse direct estimate is intentionally conservative:
    //   compressed sparse storage: values + inner indices + outer pointers
    //   factorization fill/workspace reserve: 3x compressed storage
    //   RHS/vector reserve: 8 dense vectors
    //   final safety margin: +25%
    std::size_t per_nnz = sizeof(double) + sizeof(int);
    std::size_t storage = 0;
    bool ok = checked_mul(request.sparse_nonzeros, per_nnz, storage);
    std::size_t outer = 0;
    std::size_t outer_len = 0;
    ok = checked_add(d, 1, outer_len) && ok;
    ok = checked_mul(outer_len, sizeof(int), outer) && ok;
    ok = checked_add(storage, outer, storage) && ok;
    std::size_t factor = 0;
    ok = checked_mul(storage, 3, factor) && ok;
    std::size_t vectors = 0;
    ok = checked_mul(d, 8 * sizeof(double), vectors) && ok;
    std::size_t subtotal = 0;
    ok = checked_add(storage, factor, subtotal) && ok;
    ok = checked_add(subtotal, vectors, subtotal) && ok;
    std::size_t total = 0;
    ok = checked_scale(subtotal, 5, 4, total) && ok;

    estimate.estimate_available = true;
    estimate.overflow = !ok;
    estimate.quality = EstimateQuality::Approximate;
    estimate.estimated_working_memory_bytes = subtotal;
    estimate.estimated_memory_bytes = total;
    estimate.reason = "sparse direct estimate uses compressed storage + 3x fill/workspace reserve + vector reserve + 25% safety margin.";
    apply_budget(estimate, request);
    return estimate;
}

inline SolverEstimate pcg_estimate(const SolverRequest& request) {
    SolverEstimate estimate;
    estimate.solver = PlannerSolverKind::PCG;
    const std::size_t d = effective_dimension(request);
    if (d == 0) {
        estimate.reason = "system dimension is unknown";
        return estimate;
    }

    // Explicit PCG estimate:
    //   sparse matrix storage if nnz is known
    //   PCG vectors and diagonal/Jacobi preconditioner reserve: 10 vectors
    //   final safety margin: +25%
    bool ok = true;
    std::size_t storage = 0;
    if (request.sparse_nonzeros > 0) {
        ok = checked_mul(request.sparse_nonzeros, sizeof(double) + sizeof(int), storage);
        std::size_t outer = 0;
        std::size_t outer_len = 0;
        ok = checked_add(d, 1, outer_len) && ok;
        ok = checked_mul(outer_len, sizeof(int), outer) && ok;
        ok = checked_add(storage, outer, storage) && ok;
    }
    std::size_t vectors = 0;
    ok = checked_mul(d, 10 * sizeof(double), vectors) && ok;
    std::size_t subtotal = 0;
    ok = checked_add(storage, vectors, subtotal) && ok;
    std::size_t total = 0;
    ok = checked_scale(subtotal, 5, 4, total) && ok;

    estimate.estimate_available = true;
    estimate.overflow = !ok;
    estimate.quality = request.sparse_nonzeros > 0 ? EstimateQuality::Approximate : EstimateQuality::LowerBound;
    estimate.estimated_working_memory_bytes = subtotal;
    estimate.estimated_memory_bytes = total;
    estimate.reason = request.sparse_nonzeros > 0
        ? "PCG estimate uses sparse storage + PCG vector/preconditioner reserve + 25% safety margin."
        : "PCG estimate is vector/preconditioner-only because sparse nnz is unknown.";
    apply_budget(estimate, request);
    return estimate;
}

inline SolverEstimate matrix_free_estimate(const SolverRequest& request) {
    SolverEstimate estimate;
    estimate.solver = PlannerSolverKind::MatrixFreePCG;
    const std::size_t d = effective_dimension(request);
    if (!request.matrix_free_operator_available && !request.genotype_operator_available) {
        estimate.reason = "matrix-free or genotype operator was not reported by caller";
        return estimate;
    }
    if (d == 0) {
        estimate.reason = "system dimension is unknown";
        return estimate;
    }

    // Matrix-free estimate is a lower-bound working set:
    //   PCG vector/preconditioner reserve: 12 vectors
    //   genotype block buffers are model-owned today and not yet exposed here
    //   final safety margin: +25%
    std::size_t vectors = 0;
    bool ok = checked_mul(d, 12 * sizeof(double), vectors);
    std::size_t total = 0;
    ok = checked_scale(vectors, 5, 4, total) && ok;

    estimate.estimate_available = true;
    estimate.overflow = !ok;
    estimate.quality = EstimateQuality::LowerBound;
    estimate.estimated_working_memory_bytes = vectors;
    estimate.estimated_memory_bytes = total;
    estimate.reason = "matrix-free estimate covers PCG vectors/preconditioner only; model-owned genotype block buffers are not exposed yet.";
    apply_budget(estimate, request);
    return estimate;
}

inline SolverEstimate stcg_estimate(const SolverRequest& request) {
    SolverEstimate estimate;
    estimate.solver = PlannerSolverKind::STCG;
    if (!request.stcg_available) {
        estimate.reason = "STCG infrastructure was not reported by caller";
        return estimate;
    }
    const std::size_t d = effective_dimension(request);
    if (d == 0) {
        estimate.reason = "system dimension is unknown";
        return estimate;
    }

    // STCG estimate is currently a lower-bound working set:
    //   Krylov/STCG vector reserve: 12 vectors
    //   stochastic trace samples and genotype block buffers are model-owned today
    //   final safety margin: +25%
    std::size_t vectors = 0;
    bool ok = checked_mul(d, 12 * sizeof(double), vectors);
    std::size_t total = 0;
    ok = checked_scale(vectors, 5, 4, total) && ok;

    estimate.estimate_available = true;
    estimate.overflow = !ok;
    estimate.quality = EstimateQuality::LowerBound;
    estimate.estimated_working_memory_bytes = vectors;
    estimate.estimated_memory_bytes = total;
    estimate.reason = "STCG estimate uses matrix-free PCG working-set lower bound; stochastic trace samples are not exposed to planner yet.";
    apply_budget(estimate, request);
    return estimate;
}

inline SolverEstimate apy_estimate(const SolverRequest& request) {
    SolverEstimate estimate;
    estimate.solver = PlannerSolverKind::APY;
    const std::size_t d = effective_dimension(request);
    const std::size_t c = request.apy_core_size;
    if (!request.apy_available || c == 0) {
        estimate.reason = "APY core was not configured by caller";
        return estimate;
    }

    // APY estimate:
    //   dense core matrix + factor/workspace: 24*c*c bytes
    //   full-length diagonal/vector reserve: 8 vectors over system dimension
    //   final safety margin: +25%
    std::size_t c2 = 0;
    bool ok = checked_mul(c, c, c2);
    std::size_t core = 0;
    ok = checked_mul(c2, 3 * sizeof(double), core) && ok;
    std::size_t vectors = 0;
    ok = checked_mul(d, 8 * sizeof(double), vectors) && ok;
    std::size_t subtotal = 0;
    ok = checked_add(core, vectors, subtotal) && ok;
    std::size_t total = 0;
    ok = checked_scale(subtotal, 5, 4, total) && ok;

    estimate.estimate_available = true;
    estimate.overflow = !ok;
    estimate.quality = EstimateQuality::Approximate;
    estimate.estimated_working_memory_bytes = subtotal;
    estimate.estimated_memory_bytes = total;
    estimate.reason = "APY estimate uses dense core workspace + full-length vector reserve + 25% safety margin.";
    apply_budget(estimate, request);
    return estimate;
}

inline const SolverEstimate* find_candidate(const std::vector<SolverEstimate>& candidates, PlannerSolverKind solver) {
    for (const auto& candidate : candidates) {
        if (candidate.solver == solver) return &candidate;
    }
    return nullptr;
}

inline bool is_usable(const SolverEstimate& estimate) {
    return estimate.estimate_available && estimate.feasible && !estimate.overflow;
}

inline PlannerSolverKind first_usable(const std::vector<SolverEstimate>& candidates,
                                      const std::vector<PlannerSolverKind>& preference) {
    for (auto solver : preference) {
        const SolverEstimate* candidate = find_candidate(candidates, solver);
        if (candidate && is_usable(*candidate)) return solver;
    }
    return PlannerSolverKind::Unknown;
}

} // namespace solver_planner_detail

inline SolverPlan plan_solver(const SolverRequest& request) {
    SolverPlan plan;
    const MemoryBudget budget = solver_planner_detail::effective_memory_budget(request);
    plan.requested_solver = request.requested_solver;
    plan.requested_solver_label = request.requested_solver_label.empty()
        ? to_string(request.requested_solver)
        : request.requested_solver_label;
    plan.current_solver_label = request.current_solver_label.empty() ? "unknown" : request.current_solver_label;
    plan.memory_budget_available = budget.specified;
    plan.memory_budget_bytes = budget.bytes;
    plan.memory_budget_source = budget.source;
    plan.system_dimension = solver_planner_detail::effective_dimension(request);
    plan.n_samples = request.n_samples;
    plan.n_markers = request.n_markers;
    plan.diagnostic_only = request.diagnostic_only;

    plan.candidates.push_back(solver_planner_detail::dense_estimate(request));
    if (request.sparse_matrix_available || request.matrix_sparsity == PlannerSparsity::Sparse) {
        plan.candidates.push_back(solver_planner_detail::sparse_estimate(request));
        plan.candidates.push_back(solver_planner_detail::pcg_estimate(request));
    }
    if (request.apy_available) {
        plan.candidates.push_back(solver_planner_detail::apy_estimate(request));
    }
    if (request.matrix_free_operator_available || request.genotype_operator_available) {
        plan.candidates.push_back(solver_planner_detail::matrix_free_estimate(request));
    }
    if (request.stcg_available) {
        plan.candidates.push_back(solver_planner_detail::stcg_estimate(request));
    }

    if (request.requested_solver != PlannerSolverKind::Auto &&
        request.requested_solver != PlannerSolverKind::Unknown) {
        plan.recommended_solver = request.requested_solver;
        const SolverEstimate* requested = solver_planner_detail::find_candidate(plan.candidates, request.requested_solver);
        if (requested) {
            plan.estimate_available = requested->estimate_available;
            plan.estimate_quality = requested->quality;
            plan.estimated_memory_bytes = requested->estimated_memory_bytes;
            plan.recommendation_reason = "Explicit solver request is respected in diagnostic-only mode. " + requested->reason;
        } else {
            plan.recommendation_reason = "Explicit solver request is respected in diagnostic-only mode, but this planner lacks enough problem metadata to estimate that solver.";
        }
        plan.fallback_solver = solver_planner_detail::first_usable(
            plan.candidates,
            {PlannerSolverKind::DenseDirect, PlannerSolverKind::APY, PlannerSolverKind::SparseDirect,
             PlannerSolverKind::MatrixFreePCG, PlannerSolverKind::STCG, PlannerSolverKind::PCG});
        return plan;
    }

    const bool dense_sparse = request.matrix_sparsity == PlannerSparsity::Dense ||
                              request.matrix_type == PlannerMatrixKind::Dense ||
                              request.matrix_sparsity == PlannerSparsity::Unknown;
    std::vector<PlannerSolverKind> preference;
    if (dense_sparse) preference.push_back(PlannerSolverKind::DenseDirect);
    preference.push_back(PlannerSolverKind::APY);
    preference.push_back(PlannerSolverKind::SparseDirect);
    preference.push_back(PlannerSolverKind::MatrixFreePCG);
    preference.push_back(PlannerSolverKind::STCG);
    preference.push_back(PlannerSolverKind::PCG);
    if (!dense_sparse) preference.push_back(PlannerSolverKind::DenseDirect);

    plan.recommended_solver = solver_planner_detail::first_usable(plan.candidates, preference);
    if (plan.recommended_solver == PlannerSolverKind::Unknown && !plan.candidates.empty()) {
        plan.recommendation_reason = "No feasible solver candidate was found within the available memory estimates.";
        plan.fallback_solver = PlannerSolverKind::Unknown;
        return plan;
    }

    const SolverEstimate* recommended = solver_planner_detail::find_candidate(plan.candidates, plan.recommended_solver);
    if (recommended) {
        plan.estimate_available = recommended->estimate_available;
        plan.estimate_quality = recommended->quality;
        plan.estimated_memory_bytes = recommended->estimated_memory_bytes;
        plan.recommendation_reason = recommended->reason;
    }

    if (plan.recommended_solver != PlannerSolverKind::DenseDirect) {
        const SolverEstimate* dense = solver_planner_detail::find_candidate(plan.candidates, PlannerSolverKind::DenseDirect);
        if (dense && dense->estimate_available && !dense->feasible) {
            plan.recommendation_reason = "Dense direct was not selected because " + dense->reason + " " + plan.recommendation_reason;
        }
    }

    plan.fallback_solver = solver_planner_detail::first_usable(
        plan.candidates,
        {PlannerSolverKind::SparseDirect, PlannerSolverKind::MatrixFreePCG, PlannerSolverKind::STCG,
         PlannerSolverKind::PCG, PlannerSolverKind::DenseDirect});
    if (plan.fallback_solver == plan.recommended_solver) {
        plan.fallback_solver = PlannerSolverKind::Unknown;
    }
    return plan;
}

inline std::string format_solver_plan(const SolverPlan& plan) {
    auto value_or_unknown = [](std::size_t value) {
        return value == 0 ? std::string("unknown") : std::to_string(value);
    };

    std::ostringstream os;
    os << "Solver planner mode : " << (plan.diagnostic_only ? "diagnostic-only" : "active") << "\n";
    os << "Requested solver    : " << plan.requested_solver_label << "\n";
    os << "Current solver      : " << plan.current_solver_label << "\n";
    os << "Recommended solver  : " << to_string(plan.recommended_solver) << "\n";
    os << "System dimension    : " << value_or_unknown(plan.system_dimension) << "\n";
    os << "N                    : " << value_or_unknown(plan.n_samples) << "\n";
    os << "M                    : " << value_or_unknown(plan.n_markers) << "\n";
    os << "Estimated memory     : "
       << (plan.estimate_available ? bytes_to_string(plan.estimated_memory_bytes) : std::string("unknown")) << "\n";
    os << "Estimate quality     : " << to_string(plan.estimate_quality) << "\n";
    os << "Memory budget        : "
       << (plan.memory_budget_available ? bytes_to_string(plan.memory_budget_bytes) : std::string("unknown")) << "\n";
    os << "Memory budget source : " << to_string(plan.memory_budget_source) << "\n";
    os << "Fallback solver      : " << to_string(plan.fallback_solver) << "\n";
    os << "Reason               : " << (plan.recommendation_reason.empty() ? "unknown" : plan.recommendation_reason) << "\n";
    return os.str();
}

} // namespace cosmic
