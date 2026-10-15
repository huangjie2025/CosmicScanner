#include "association/multilocus/multilocus_engine.h"
#include "association/multilocus/multilocus_internal.h"
#include "genetics/ld/ld_calculator.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

namespace cosmic::association {

MultiLocusOptions farmCpuOptions(int max_iterations, int max_pseudo_qtn) {
    MultiLocusOptions options;
    options.method = MultiLocusMethod::FARMCPU;
    options.max_iterations = max_iterations;
    options.max_pseudo_qtn = max_pseudo_qtn;
    options.ld_r2 = 0.7;
    return options;
}

MultiLocusEngine::MultiLocusEngine(MultiLocusOptions options) : options_(std::move(options)) {
    if (options_.max_iterations < 1 || options_.max_pseudo_qtn < 0 ||
        !std::isfinite(options_.ld_r2) || options_.ld_r2 < 0.0 || options_.ld_r2 > 1.0)
        throw std::invalid_argument("Invalid multi-locus options");
}

AssociationModel MultiLocusEngine::getModelType() const {
    if (options_.method == MultiLocusMethod::BLINK) return AssociationModel::BLINK;
    if (options_.method == MultiLocusMethod::CONDITIONAL) return AssociationModel::CONDITIONAL;
    if (options_.method == MultiLocusMethod::CONDITIONAL_JOINT) return AssociationModel::CONDITIONAL_JOINT;
    return AssociationModel::FARMCPU;
}

AssociationCapabilities MultiLocusEngine::capabilities() const { return {true, true, false}; }

void MultiLocusEngine::prepare(const ScanContext& context) {
    if (context.phenotypes.cols() != 1 || context.phenotypes.rows() < 3 ||
        context.covariates.rows() != context.phenotypes.rows())
        throw std::invalid_argument("Multi-locus analysis requires one phenotype and aligned covariates");
    if (!context.phenotypes.allFinite() || !context.covariates.allFinite())
        throw std::invalid_argument("Multi-locus phenotype/covariates must be finite");
    context_ = &context;
}

std::vector<AssociationResult> MultiLocusEngine::scanBlock(const GenotypeBlock& block) {
    if (!context_) throw std::logic_error("Multi-locus engine was not prepared");
    if (block.dosages.rows() <= 0 || block.dosages.cols() != static_cast<int>(block.variant_ids.size()) ||
        !block.dosages.allFinite()) throw std::invalid_argument("Invalid multi-locus genotype block");
    const int n = context_->phenotypes.rows();
    if (static_cast<int>(context_->observation_to_sample.size()) != n)
        throw std::invalid_argument("Invalid observation-to-sample map");
    Eigen::MatrixXd genotypes(n, block.dosages.cols());
    for (int i = 0; i < n; ++i) {
        const Eigen::Index row = context_->observation_to_sample[i];
        if (row < 0 || row >= block.dosages.rows()) throw std::invalid_argument("Genotype sample index is out of range");
        genotypes.row(i) = block.dosages.row(row);
    }
    const Eigen::VectorXd y = context_->phenotypes.col(0);
    const Eigen::MatrixXd& base = context_->covariates;
    const int m = genotypes.cols();
    std::unordered_map<std::string, int> by_id;
    for (int j = 0; j < m; ++j) {
        if (block.variant_ids[j].empty() || !by_id.emplace(block.variant_ids[j], j).second)
            throw std::invalid_argument("Duplicate or empty multi-locus variant ID");
    }

    selected_.clear();
    std::vector<int> condition_indices;
    if (options_.method == MultiLocusMethod::CONDITIONAL ||
        options_.method == MultiLocusMethod::CONDITIONAL_JOINT) {
        std::unordered_set<int> unique;
        for (const auto& id : options_.condition_ids) {
            const auto it = by_id.find(id);
            if (it == by_id.end()) throw std::invalid_argument("Condition SNP not found: " + id);
            if (!unique.insert(it->second).second) throw std::invalid_argument("Duplicate condition SNP: " + id);
            condition_indices.push_back(it->second);
        }
        if (condition_indices.empty()) throw std::invalid_argument("Conditional GWAS needs at least one condition SNP");
        if (options_.method == MultiLocusMethod::CONDITIONAL) selected_ = condition_indices;
    }

    std::vector<detail::MarkerFit> fits(m);
    auto scan = [&](const std::vector<int>& conditioning) {
        for (int j = 0; j < m; ++j) fits[j] = detail::fitMarker(y, base, genotypes, conditioning, j);
    };
    scan(selected_);
    int completed_iterations = 1;

    if (options_.method == MultiLocusMethod::CONDITIONAL_JOINT) {
        Eigen::MatrixXd lead_dosages(n, condition_indices.size());
        for (int j = 0; j < static_cast<int>(condition_indices.size()); ++j)
            lead_dosages.col(j) = genotypes.col(condition_indices[j]);
        const Eigen::MatrixXd lead_ld = genetics::LDCalculator::computePairwiseLD(std::move(lead_dosages));
        std::unordered_set<int> selected_set;
        while (selected_.size() < condition_indices.size()) {
            int best_marker = -1;
            double best_p = 1.0;
            for (int candidate_position = 0; candidate_position < static_cast<int>(condition_indices.size());
                 ++candidate_position) {
                const int candidate = condition_indices[candidate_position];
                if (selected_set.count(candidate) || !fits[candidate].tested) continue;
                bool linked = false;
                for (int retained : selected_) {
                    const auto retained_it = std::find(condition_indices.begin(), condition_indices.end(), retained);
                    const int retained_position = static_cast<int>(retained_it - condition_indices.begin());
                    const double r = lead_ld(candidate_position, retained_position);
                    if (r * r >= 0.99) { linked = true; break; }
                }
                if (!linked && (best_marker < 0 || fits[candidate].p < best_p)) {
                    best_marker = candidate;
                    best_p = fits[candidate].p;
                }
            }
            if (best_marker < 0) break;
            std::vector<int> trial = selected_;
            trial.push_back(best_marker);
            const auto bic_selected = detail::selectByBic(y, base, genotypes, trial,
                                                           static_cast<int>(trial.size()));
            if (bic_selected.size() != trial.size()) break;
            selected_ = std::move(trial);
            selected_set.insert(best_marker);
            completed_iterations = static_cast<int>(selected_.size());
            scan(selected_);
        }
    } else if (options_.method != MultiLocusMethod::CONDITIONAL) {
        std::vector<double> ranking_p(m);
        if (options_.initial_p_values.empty()) {
            for (int j = 0; j < m; ++j) ranking_p[j] = fits[j].tested ? fits[j].p : 1.0;
        } else {
            if (static_cast<int>(options_.initial_p_values.size()) != m)
                throw std::invalid_argument("Initial LMM ranking does not match genotype block");
            ranking_p = options_.initial_p_values;
        }
        for (int iteration = 1; iteration <= options_.max_iterations; ++iteration) {
            std::vector<int> ranked(m);
            std::iota(ranked.begin(), ranked.end(), 0);
            std::stable_sort(ranked.begin(), ranked.end(), [&](int a, int b) {
                return ranking_p[a] < ranking_p[b] || (ranking_p[a] == ranking_p[b] && a < b);
            });
            if (options_.method == MultiLocusMethod::BLINK) {
                std::vector<int> eligible;
                for (int index : ranked) if (ranking_p[index] <= 0.01 || eligible.size() < 2) eligible.push_back(index);
                ranked.swap(eligible);
            }
            const int pool = std::max(options_.max_pseudo_qtn, options_.max_pseudo_qtn * 3);
            const auto independent = detail::ldPruneRanked(genotypes, ranked, pool, options_.ld_r2);
            const auto next = detail::selectByBic(y, base, genotypes, independent, options_.max_pseudo_qtn);
            completed_iterations = iteration;
            if (next == selected_) break;
            selected_ = next;
            scan(selected_);
            for (int j = 0; j < m; ++j) ranking_p[j] = fits[j].tested ? fits[j].p : 1.0;
        }
    }

    if (options_.method == MultiLocusMethod::CONDITIONAL_JOINT) {
        for (const int marker : selected_) {
            std::vector<int> other_selected;
            other_selected.reserve(selected_.size() - 1);
            for (const int candidate : selected_)
                if (candidate != marker) other_selected.push_back(candidate);
            fits[marker] = detail::fitMarker(y, base, genotypes, other_selected, marker);
        }
    }

    std::unordered_set<int> selected_set(selected_.begin(), selected_.end());
    results_.clear();
    std::vector<AssociationResult> generic;
    results_.reserve(m);
    generic.reserve(m);
    for (int j = 0; j < m; ++j) {
        const auto& fit = fits[j];
        MultiLocusResult detail_result;
        detail_result.variant_id = block.variant_ids[j];
        detail_result.n_obs = n;
        detail_result.iteration = completed_iterations;
        detail_result.selected_pseudo_qtn = selected_set.count(j) != 0;
        AssociationResult result;
        result.test_id = block.variant_ids[j];
        result.n_obs = n;
        if (fit.tested) {
            detail_result.beta = fit.beta;
            detail_result.se = fit.se;
            detail_result.p_value = fit.p;
            result.status = AssociationStatus::SUCCESS;
            result.effects = Eigen::VectorXd::Constant(1, fit.beta);
            result.effect_covariance = Eigen::MatrixXd::Constant(1, 1, fit.se * fit.se);
            result.statistic = fit.beta * fit.beta / (fit.se * fit.se);
            result.degrees_of_freedom = 1.0;
            result.p_value = fit.p;
        } else {
            result.status = AssociationStatus::FILTERED;
            result.message = "Marker is constant or collinear with fixed effects";
        }
        results_.push_back(detail_result);
        generic.push_back(std::move(result));
    }
    return generic;
}

} // namespace cosmic::association
