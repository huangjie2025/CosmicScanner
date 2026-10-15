#pragma once

#include "analysis_dataset.h"
#include "solver_planner.h"
#include "validation.h"

#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#ifndef _WIN32
#include <sys/resource.h>
#include <unistd.h>
#endif

namespace cosmic {

class JsonValue {
public:
    enum class Type { Null, Bool, Int, Double, String, Array, Object };

    JsonValue() = default;
    static JsonValue null() { return JsonValue(); }
    static JsonValue boolean(bool v) { JsonValue j; j.type_ = Type::Bool; j.bool_ = v; return j; }
    static JsonValue integer(std::int64_t v) { JsonValue j; j.type_ = Type::Int; j.int_ = v; return j; }
    static JsonValue uinteger(std::uint64_t v) {
        if (v > static_cast<std::uint64_t>(std::numeric_limits<std::int64_t>::max())) {
            JsonValue j;
            j.type_ = Type::String;
            j.string_ = std::to_string(v);
            return j;
        }
        return integer(static_cast<std::int64_t>(v));
    }
    static JsonValue number(double v) {
        if (!std::isfinite(v)) return JsonValue::null();
        JsonValue j; j.type_ = Type::Double; j.double_ = v; return j;
    }
    static JsonValue string(std::string v) { JsonValue j; j.type_ = Type::String; j.string_ = std::move(v); return j; }
    static JsonValue array() { JsonValue j; j.type_ = Type::Array; return j; }
    static JsonValue object() { JsonValue j; j.type_ = Type::Object; return j; }

    void push(JsonValue value) {
        ensure(Type::Array);
        array_.push_back(std::move(value));
    }

    void set(std::string key, JsonValue value) {
        ensure(Type::Object);
        object_.push_back({std::move(key), std::move(value)});
    }

    std::string dump(int indent = 2) const {
        std::ostringstream os;
        write(os, 0, indent);
        os << "\n";
        return os.str();
    }

private:
    void ensure(Type type) {
        if (type_ != type) throw std::runtime_error("JSON value has incompatible type.");
    }

    static std::string escape(const std::string& input) {
        std::ostringstream os;
        for (unsigned char c : input) {
            switch (c) {
                case '"': os << "\\\""; break;
                case '\\': os << "\\\\"; break;
                case '\b': os << "\\b"; break;
                case '\f': os << "\\f"; break;
                case '\n': os << "\\n"; break;
                case '\r': os << "\\r"; break;
                case '\t': os << "\\t"; break;
                default:
                    if (c < 0x20) {
                        os << "\\u" << std::hex << std::setw(4) << std::setfill('0')
                           << static_cast<int>(c) << std::dec << std::setfill(' ');
                    } else {
                        os << static_cast<char>(c);
                    }
                    break;
            }
        }
        return os.str();
    }

    static void write_indent(std::ostream& os, int depth, int indent) {
        for (int i = 0; i < depth * indent; ++i) os.put(' ');
    }

    void write(std::ostream& os, int depth, int indent) const {
        switch (type_) {
            case Type::Null:
                os << "null";
                break;
            case Type::Bool:
                os << (bool_ ? "true" : "false");
                break;
            case Type::Int:
                os << int_;
                break;
            case Type::Double:
                os << std::setprecision(17) << double_;
                break;
            case Type::String:
                os << '"' << escape(string_) << '"';
                break;
            case Type::Array:
                os << "[";
                if (!array_.empty()) os << "\n";
                for (std::size_t i = 0; i < array_.size(); ++i) {
                    write_indent(os, depth + 1, indent);
                    array_[i].write(os, depth + 1, indent);
                    if (i + 1 < array_.size()) os << ",";
                    os << "\n";
                }
                if (!array_.empty()) write_indent(os, depth, indent);
                os << "]";
                break;
            case Type::Object:
                os << "{";
                if (!object_.empty()) os << "\n";
                for (std::size_t i = 0; i < object_.size(); ++i) {
                    write_indent(os, depth + 1, indent);
                    os << '"' << escape(object_[i].first) << "\": ";
                    object_[i].second.write(os, depth + 1, indent);
                    if (i + 1 < object_.size()) os << ",";
                    os << "\n";
                }
                if (!object_.empty()) write_indent(os, depth, indent);
                os << "}";
                break;
        }
    }

    Type type_ = Type::Null;
    bool bool_ = false;
    std::int64_t int_ = 0;
    double double_ = 0.0;
    std::string string_;
    std::vector<JsonValue> array_;
    std::vector<std::pair<std::string, JsonValue>> object_;
};

inline std::string iso8601_utc(std::chrono::system_clock::time_point tp) {
    const std::time_t t = std::chrono::system_clock::to_time_t(tp);
    std::tm tm{};
#ifdef _WIN32
    gmtime_s(&tm, &t);
#else
    gmtime_r(&t, &tm);
#endif
    std::ostringstream os;
    os << std::put_time(&tm, "%Y-%m-%dT%H:%M:%SZ");
    return os.str();
}

inline JsonValue input_file_manifest(const std::string& path) {
    JsonValue j = JsonValue::object();
    j.set("path", JsonValue::string(path));
    if (path.empty()) {
        j.set("exists", JsonValue::boolean(false));
        return j;
    }
    std::error_code ec;
    const bool exists = std::filesystem::exists(path, ec);
    j.set("exists", JsonValue::boolean(exists && !ec));
    if (exists && !ec && std::filesystem::is_regular_file(path, ec) && !ec) {
        const auto size = std::filesystem::file_size(path, ec);
        if (!ec) j.set("size_bytes", JsonValue::uinteger(static_cast<std::uint64_t>(size)));
    }
    return j;
}

inline std::uint64_t peak_rss_bytes() {
#ifndef _WIN32
    struct rusage usage;
    if (getrusage(RUSAGE_SELF, &usage) == 0) {
#if defined(__APPLE__)
        return static_cast<std::uint64_t>(usage.ru_maxrss);
#else
        return static_cast<std::uint64_t>(usage.ru_maxrss) * 1024ULL;
#endif
    }
#endif
    return 0;
}

inline JsonValue alignment_manifest(const SampleAlignmentDiagnostics& d) {
    JsonValue j = JsonValue::object();
    j.set("genotype_samples", JsonValue::uinteger(d.genotype_samples));
    j.set("phenotype_rows", JsonValue::uinteger(d.phenotype_rows));
    j.set("covariate_rows", JsonValue::uinteger(d.covariate_rows));
    j.set("pedigree_rows", JsonValue::uinteger(d.pedigree_rows));
    j.set("analysis_samples", JsonValue::uinteger(d.analysis_samples));
    j.set("matched_phenotypes", JsonValue::uinteger(d.matched_phenotypes));
    j.set("matched_covariates", JsonValue::uinteger(d.matched_covariates));
    j.set("matched_pedigree", JsonValue::uinteger(d.matched_pedigree));
    j.set("duplicate_genotype_ids", JsonValue::uinteger(d.duplicate_genotype_ids));
    j.set("duplicate_phenotype_ids", JsonValue::uinteger(d.duplicate_phenotype_ids));
    j.set("duplicate_covariate_ids", JsonValue::uinteger(d.duplicate_covariate_ids));
    j.set("genotype_without_phenotype", JsonValue::uinteger(d.genotype_without_phenotype));
    j.set("genotype_without_covariate", JsonValue::uinteger(d.genotype_without_covariate));
    j.set("phenotype_without_genotype", JsonValue::uinteger(d.phenotype_without_genotype));
    j.set("covariate_without_genotype", JsonValue::uinteger(d.covariate_without_genotype));
    return j;
}

inline JsonValue analysis_dataset_manifest(const AnalysisDataset& dataset) {
    JsonValue j = JsonValue::object();
    j.set("n_samples_analysis", JsonValue::uinteger(dataset.samples().size()));
    j.set("sample_alignment", alignment_manifest(dataset.samples().diagnostics()));
    const auto genotype = dataset.genotype();
    if (genotype) {
        const auto& meta = genotype->metadata();
        JsonValue g = JsonValue::object();
        g.set("format", JsonValue::string(to_string(meta.format)));
        g.set("n_samples", JsonValue::uinteger(meta.n_samples));
        g.set("n_markers", JsonValue::uinteger(meta.n_markers));
        g.set("source_label", JsonValue::string(meta.source_label));
        JsonValue caps = JsonValue::object();
        caps.set("supports_streaming", JsonValue::boolean(meta.capabilities.supports_streaming));
        caps.set("supports_random_access", JsonValue::boolean(meta.capabilities.supports_random_access));
        caps.set("supports_prefetch", JsonValue::boolean(meta.capabilities.supports_prefetch));
        caps.set("supports_sample_subset", JsonValue::boolean(meta.capabilities.supports_sample_subset));
        caps.set("materializes_full_matrix", JsonValue::boolean(meta.capabilities.materializes_full_matrix));
        g.set("capabilities", std::move(caps));
        j.set("genotype", std::move(g));
    } else {
        j.set("genotype", JsonValue::null());
    }
    if (const auto* p = dataset.phenotype()) {
        JsonValue traits = JsonValue::array();
        for (const auto& t : p->trait_names) traits.push(JsonValue::string(t));
        j.set("n_phenotype_rows", JsonValue::uinteger(p->rows()));
        j.set("traits", std::move(traits));
    }
    if (const auto* c = dataset.covariates()) {
        JsonValue covars = JsonValue::array();
        for (const auto& name : c->numeric_names) covars.push(JsonValue::string(name));
        j.set("n_covariate_rows", JsonValue::uinteger(c->rows()));
        j.set("covariates", std::move(covars));
    }
    if (const auto* obs = dataset.observations()) {
        j.set("n_observations", JsonValue::uinteger(obs->observations.size()));
    }
    if (const auto* ped = dataset.pedigree()) {
        j.set("n_pedigree_rows", JsonValue::uinteger(ped->records.size()));
    }
    return j;
}

inline JsonValue solver_plan_manifest(const SolverPlan& plan) {
    JsonValue j = JsonValue::object();
    j.set("requested", JsonValue::string(plan.requested_solver_label));
    j.set("current_actual", JsonValue::string(plan.current_solver_label));
    j.set("recommended", JsonValue::string(to_string(plan.recommended_solver)));
    j.set("recommendation_reason", JsonValue::string(plan.recommendation_reason));
    j.set("planner_mode", JsonValue::string(plan.diagnostic_only ? "diagnostic-only" : "dispatch"));
    j.set("diagnostic_only", JsonValue::boolean(plan.diagnostic_only));
    j.set("system_dimension", plan.system_dimension ? JsonValue::uinteger(plan.system_dimension) : JsonValue::null());
    j.set("n_samples", plan.n_samples ? JsonValue::uinteger(plan.n_samples) : JsonValue::null());
    j.set("n_markers", plan.n_markers ? JsonValue::uinteger(plan.n_markers) : JsonValue::null());
    j.set("estimated_memory_bytes", plan.estimate_available ? JsonValue::uinteger(plan.estimated_memory_bytes) : JsonValue::null());
    j.set("estimate_quality", JsonValue::string(to_string(plan.estimate_quality)));
    j.set("memory_budget_bytes", plan.memory_budget_available ? JsonValue::uinteger(plan.memory_budget_bytes) : JsonValue::null());
    j.set("memory_budget_source", JsonValue::string(to_string(plan.memory_budget_source)));
    j.set("fallback_solver", JsonValue::string(to_string(plan.fallback_solver)));
    return j;
}

inline JsonValue validation_report_manifest(const ValidationReport& report, bool include_pass_issues = false) {
    const auto counts = report.counts();
    JsonValue j = JsonValue::object();
    j.set("status", JsonValue::string(to_string(report.overall_status())));
    j.set("passes", JsonValue::uinteger(counts.passed));
    j.set("warnings", JsonValue::uinteger(counts.warnings));
    j.set("failures", JsonValue::uinteger(counts.failures));
    JsonValue issues = JsonValue::array();
    for (const auto& issue : report.issues()) {
        if (issue.severity == ValidationSeverity::Pass && !include_pass_issues) continue;
        JsonValue item = JsonValue::object();
        item.set("code", JsonValue::string(issue.code));
        item.set("stage", JsonValue::string(to_string(issue.stage)));
        item.set("severity", JsonValue::string(to_string(issue.severity)));
        item.set("context", JsonValue::string(issue.context));
        item.set("message", JsonValue::string(issue.message));
        if (!issue.field.empty()) item.set("field", JsonValue::string(issue.field));
        if (!issue.observed_value.empty()) item.set("observed_value", JsonValue::string(issue.observed_value));
        if (!issue.expected_value.empty()) item.set("expected_value", JsonValue::string(issue.expected_value));
        issues.push(std::move(item));
    }
    j.set("issues", std::move(issues));
    return j;
}

inline JsonValue output_record_manifest(const std::string& type, const std::string& path) {
    JsonValue j = JsonValue::object();
    j.set("type", JsonValue::string(type));
    j.set("path", JsonValue::string(path));
    std::error_code ec;
    if (!path.empty() && std::filesystem::exists(path, ec) && !ec) {
        j.set("exists", JsonValue::boolean(true));
        if (std::filesystem::is_regular_file(path, ec) && !ec) {
            const auto size = std::filesystem::file_size(path, ec);
            if (!ec) j.set("size_bytes", JsonValue::uinteger(static_cast<std::uint64_t>(size)));
        }
    } else {
        j.set("exists", JsonValue::boolean(false));
    }
    return j;
}

inline void write_json_atomic(const std::string& final_path, const JsonValue& value) {
    const std::filesystem::path target(final_path);
    if (target.has_parent_path()) {
        std::filesystem::create_directories(target.parent_path());
    }
#ifndef _WIN32
    const auto pid = static_cast<long long>(getpid());
#else
    const auto pid = static_cast<long long>(GetCurrentProcessId());
#endif
    const std::filesystem::path tmp = target.string() + ".tmp." + std::to_string(pid);
    {
        std::ofstream out(tmp, std::ios::binary);
        if (!out) throw std::runtime_error("Cannot write temporary JSON manifest: " + tmp.string());
        const std::string body = value.dump(2);
        out.write(body.data(), static_cast<std::streamsize>(body.size()));
        out.flush();
        if (!out) throw std::runtime_error("Failed while writing temporary JSON manifest: " + tmp.string());
    }
    std::error_code ec;
    std::filesystem::rename(tmp, target, ec);
    if (ec) {
        std::filesystem::remove(tmp);
        throw std::runtime_error("Cannot atomically replace JSON manifest '" + target.string() + "': " + ec.message());
    }
}

} // namespace cosmic
