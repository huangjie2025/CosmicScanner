#include "input/input_preparer.h"
#include "input/input_common.h"
#include "input/input_source_reader.h"
#include "input/time_normalizer.h"
#include "input/wide_to_long.h"

#include <algorithm>
#include <cerrno>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#ifdef _WIN32
#include <direct.h>
#else
#include <sys/stat.h>
#endif
#include <unordered_map>
#include <unordered_set>

namespace cosmic {
namespace {

static int makePreparedDirectory(const std::string& path) {
#ifdef _WIN32
    return ::_mkdir(path.c_str());
#else
    return ::mkdir(path.c_str(), 0755);
#endif
}

static void ensureDir(const std::string& path) {
    if (path.empty()) throw std::runtime_error("MISSING_REQUIRED_COLUMN: --prepared-input-dir is required with --input-file");
    if (makePreparedDirectory(path) != 0 && errno != EEXIST) {
        throw std::runtime_error("OUTPUT_DIR_NOT_WRITABLE: Cannot create prepared input directory: " + path);
    }
    std::string traits = path + "/traits";
    if (makePreparedDirectory(traits) != 0 && errno != EEXIST) {
        throw std::runtime_error("OUTPUT_DIR_NOT_WRITABLE: Cannot create prepared trait directory: " + traits);
    }
    std::string rrm = path + "/rrm";
    if (makePreparedDirectory(rrm) != 0 && errno != EEXIST) {
        throw std::runtime_error("OUTPUT_DIR_NOT_WRITABLE: Cannot create prepared RRM directory: " + rrm);
    }
    std::string subject = path + "/subject";
    if (makePreparedDirectory(subject) != 0 && errno != EEXIST) {
        throw std::runtime_error("OUTPUT_DIR_NOT_WRITABLE: Cannot create prepared subject directory: " + subject);
    }
}

static void writeSchemaPreview(const std::string& path,
                               const Table& table,
                               const InputPrepareOptions& options,
                               int matched_rows,
                               int unmatched_rows,
                               int duplicate_keys,
                               const std::map<std::string, int>& nonmissing) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("OUTPUT_DIR_NOT_WRITABLE: Cannot write schema preview: " + path);
    out << "section\titem\tvalue\n";
    out << "input\traw_rows\t" << table.rows.size() << "\n";
    out << "input\tcolumns\t";
    for (size_t i = 0; i < table.header.size(); ++i) out << (i ? "," : "") << table.header[i];
    out << "\n";
    out << "mapping\tid_col\t" << options.id_col << "\n";
    out << "mapping\tiid_col\t" << options.iid_col << "\n";
    out << "mapping\ttime_col\t" << options.time_col << "\n";
    out << "mapping\twide_to_long\t" << (options.wide_to_long ? "true" : "false") << "\n";
    out << "match\tmatched_rows\t" << matched_rows << "\n";
    out << "match\tunmatched_rows\t" << unmatched_rows << "\n";
    out << "time\tduplicate_id_time_keys\t" << duplicate_keys << "\n";
    for (const auto& kv : nonmissing) out << "trait\tnonmissing_" << kv.first << "\t" << kv.second << "\n";
}

static std::string lowerString(const std::string& s) {
    std::string out;
    out.reserve(s.size());
    for (char c : s) out.push_back(std::tolower(static_cast<unsigned char>(c)));
    return out;
}

static bool containsAny(const std::string& lower, const std::vector<std::string>& keys) {
    for (const auto& k : keys) if (lower.find(k) != std::string::npos) return true;
    return false;
}

static std::string jsonEscape(const std::string& s);

static void writeColumnInference(const std::string& path,
                                 const Table& table,
                                 const std::unordered_map<std::string, std::pair<std::string, std::string>>& raw_to_sample) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("OUTPUT_DIR_NOT_WRITABLE: Cannot write column inference report: " + path);
    out << "column\tsuggested_role\tconfidence\treason\tnonmissing\tmissing_rate\tunique_values\tnumeric_fraction\tid_match_rate\n";
    const std::vector<std::string> id_keys = {"fid", "iid", "id", "animal", "pig", "sample", "new_id", "ear", "tag"};
    const std::vector<std::string> time_keys = {"time", "date", "day", "batch", "month", "year", "parity", "collect"};
    const std::vector<std::string> covar_keys = {"age", "farm", "breed", "season", "year", "month", "batch", "line", "sex", "temp", "interval"};
    for (const auto& col : table.header) {
        int total = 0, miss = 0, numeric = 0, id_match = 0, date_like = 0;
        std::set<std::string> unique;
        for (const auto& row : table.rows) {
            std::string v = getCell(table, row, col);
            ++total;
            if (isMissing(v)) { ++miss; continue; }
            unique.insert(v);
            double d = 0.0;
            if (parseDouble(v, d)) ++numeric;
            int y = 0, m = 0, dd = 0;
            if (parseDate(v, y, m, dd)) ++date_like;
            if (raw_to_sample.count(v)) ++id_match;
        }
        int nonmissing = total - miss;
        double missing_rate = total ? (double)miss / total : 0.0;
        double numeric_frac = nonmissing ? (double)numeric / nonmissing : 0.0;
        double id_rate = nonmissing ? (double)id_match / nonmissing : 0.0;
        double date_frac = nonmissing ? (double)date_like / nonmissing : 0.0;
        std::string low = lowerString(col);
        std::string role = "unresolved";
        double conf = 0.1;
        std::string reason = "no strong heuristic match";
        if (id_rate > 0.5 || containsAny(low, id_keys)) {
            role = "id_col";
            conf = std::max(0.6, id_rate);
            reason = "name/id-match heuristic";
        } else if (date_frac > 0.5 || containsAny(low, time_keys)) {
            role = "time_col";
            conf = std::max(0.6, date_frac);
            reason = "name/date heuristic";
        } else if (containsAny(low, covar_keys)) {
            bool categorical = unique.size() <= 30 && numeric_frac < 0.9;
            role = categorical ? "covariate_c" : "covariate_q";
            conf = 0.6;
            reason = "known covariate name heuristic";
        } else if (numeric_frac > 0.9 && unique.size() > 5) {
            role = "trait_or_qcovar";
            conf = 0.55;
            reason = "mostly numeric with multiple values";
        } else if (unique.size() <= 30 && nonmissing > 0) {
            role = "categorical_candidate";
            conf = 0.45;
            reason = "low-cardinality nonnumeric column";
        }
        out << col << "\t" << role << "\t" << std::setprecision(4) << conf << "\t" << reason
            << "\t" << nonmissing << "\t" << missing_rate << "\t" << unique.size()
            << "\t" << numeric_frac << "\t" << id_rate << "\n";
    }
}

static void writeInferredConfigDraft(const std::string& path,
                                     const Table& table,
                                     const std::unordered_map<std::string, std::pair<std::string, std::string>>& raw_to_sample,
                                     const std::string& input_file,
                                     const std::string& output_dir) {
    std::string best_id, best_time;
    double best_id_rate = -1.0, best_time_score = -1.0;
    std::vector<std::string> trait_candidates, qcovar_candidates, ccovar_candidates;
    const std::vector<std::string> covar_keys = {"age", "farm", "breed", "season", "year", "month", "batch", "line", "sex", "temp", "interval"};
    for (const auto& col : table.header) {
        int nonmissing = 0, numeric = 0, id_match = 0, date_like = 0;
        std::set<std::string> unique;
        for (const auto& row : table.rows) {
            std::string v = getCell(table, row, col);
            if (isMissing(v)) continue;
            ++nonmissing;
            unique.insert(v);
            double d = 0.0;
            if (parseDouble(v, d)) ++numeric;
            int yy = 0, mm = 0, dd = 0;
            if (parseDate(v, yy, mm, dd)) ++date_like;
            if (raw_to_sample.count(v)) ++id_match;
        }
        if (nonmissing == 0) continue;
        double id_rate = (double)id_match / nonmissing;
        if (id_rate > best_id_rate) { best_id_rate = id_rate; best_id = col; }
        double time_score = std::max((double)date_like / nonmissing, lowerString(col).find("date") != std::string::npos || lowerString(col).find("time") != std::string::npos || lowerString(col).find("day") != std::string::npos ? 0.7 : 0.0);
        if (time_score > best_time_score) { best_time_score = time_score; best_time = col; }
        double numeric_frac = (double)numeric / nonmissing;
        std::string low = lowerString(col);
        bool known_covar = containsAny(low, covar_keys);
        if (col == best_id || col == best_time) continue;
        if (known_covar) {
            if (unique.size() <= 30 && numeric_frac < 0.9) ccovar_candidates.push_back(col);
            else qcovar_candidates.push_back(col);
        } else if (numeric_frac > 0.9 && unique.size() > 5) {
            trait_candidates.push_back(col);
        }
    }
    std::ofstream out(path);
    if (!out) throw std::runtime_error("OUTPUT_DIR_NOT_WRITABLE: Cannot write inferred config draft: " + path);
    out << "{\n";
    out << "  \"input_file\": \"" << jsonEscape(input_file) << "\",\n";
    out << "  \"input_format\": \"auto\",\n";
    out << "  \"input_model\": \"repeatability\",\n";
    out << "  \"prepared_input_dir\": \"" << jsonEscape(output_dir) << "\",\n";
    out << "  \"id_col\": \"" << jsonEscape(best_id) << "\",\n";
    out << "  \"time_col\": \"" << jsonEscape(best_time) << "\",\n";
    auto write_array = [&](const std::string& name, const std::vector<std::string>& vals, bool comma) {
        out << "  \"" << name << "\": [";
        for (size_t i = 0; i < vals.size(); ++i) out << (i ? ", " : "") << "\"" << jsonEscape(vals[i]) << "\"";
        out << "]" << (comma ? "," : "") << "\n";
    };
    write_array("trait_cols", trait_candidates, true);
    write_array("record_qcovars", qcovar_candidates, true);
    write_array("record_ccovars", ccovar_candidates, true);
    out << "  \"time_parser\": \"auto\",\n";
    out << "  \"time_unique_policy\": \"within-time-index\",\n";
    out << "  \"time_index_multiplier\": \"auto\"\n";
    out << "}\n";
}

static void writeUnresolvedColumns(const std::string& path,
                                   const Table& table,
                                   const InputPrepareOptions& options,
                                   const std::vector<WideSpec>& wide_specs) {
    std::set<std::string> used;
    used.insert(options.id_col);
    used.insert(options.iid_col);
    used.insert(options.time_col);
    for (const auto& c : options.trait_cols) used.insert(c);
    for (const auto& c : options.record_qcovars) used.insert(c);
    for (const auto& c : options.record_ccovars) used.insert(c);
    for (const auto& c : options.subject_qcovars) used.insert(c);
    for (const auto& c : options.subject_ccovars) used.insert(c);
    for (const auto& w : wide_specs) used.insert(w.column);
    std::ofstream out(path);
    if (!out) throw std::runtime_error("OUTPUT_DIR_NOT_WRITABLE: Cannot write unresolved column report: " + path);
    out << "column\tstatus\n";
    for (const auto& col : table.header) if (!used.count(col)) out << col << "\tunmapped\n";
}

static std::string jsonEscape(const std::string& s) {
    std::string out;
    for (char c : s) {
        if (c == '"' || c == '\\') { out.push_back('\\'); out.push_back(c); }
        else if (c == '\n') out += "\\n";
        else out.push_back(c);
    }
    return out;
}

static void writeHeaderAndRows(const std::string& path,
                               const std::vector<std::string>& header,
                               const std::vector<PreparedRecord>& records,
                               const std::vector<std::string>& cols,
                               bool with_time) {
    if (cols.empty()) return;
    std::ofstream out(path);
    if (!out) throw std::runtime_error("OUTPUT_DIR_NOT_WRITABLE: Cannot write prepared file: " + path);
    for (size_t i = 0; i < header.size(); ++i) out << (i ? "\t" : "") << header[i];
    for (const auto& c : cols) out << "\t" << c;
    out << "\n";
    std::set<std::string> seen_subject;
    for (const auto& rec : records) {
        if (!with_time) {
            std::string k = key2(rec.fid, rec.iid);
            if (seen_subject.count(k)) continue;
            seen_subject.insert(k);
        }
        out << rec.fid << "\t" << rec.iid;
        if (with_time) out << "\t" << rec.time;
        for (const auto& c : cols) {
            auto it = rec.values.find(c);
            out << "\t" << (it == rec.values.end() || isMissing(it->second) ? "NA" : it->second);
        }
        out << "\n";
    }
}

static void writeSubjectModelCovars(const std::string& path,
                                    const std::vector<PreparedRecord>& records,
                                    const std::vector<std::string>& cols) {
    if (cols.empty()) return;
    std::ofstream out(path);
    if (!out) throw std::runtime_error("OUTPUT_DIR_NOT_WRITABLE: Cannot write subject model covariate file: " + path);
    std::set<std::string> seen;
    for (const auto& rec : records) {
        std::string k = key2(rec.fid, rec.iid);
        if (seen.count(k)) continue;
        seen.insert(k);
        out << rec.fid << "\t" << rec.iid;
        for (const auto& c : cols) {
            auto it = rec.values.find(c);
            out << "\t" << (it == rec.values.end() || isMissing(it->second) ? "NA" : it->second);
        }
        out << "\n";
    }
}

} // namespace

PreparedInputResult InputPreparer::loadPrepared(const std::string& prepared_dir,
                                                const std::string& input_trait,
                                                const std::string& input_model) {
    if (prepared_dir.empty()) throw std::runtime_error("MISSING_REQUIRED_COLUMN: --prepared-input-dir is empty");
    PreparedInputResult res;
    res.output_dir = prepared_dir;
    res.selected_trait = input_trait.empty() ? "trait" : input_trait;
    res.pheno_file = prepared_dir + "/traits/" + res.selected_trait + ".pheno.tsv";
    res.rrm_pheno_file = prepared_dir + "/rrm/" + res.selected_trait + ".pheno.tsv";
    res.subject_pheno_file = prepared_dir + "/subject/" + res.selected_trait + ".pheno.tsv";
    res.multivariate_pheno_file = prepared_dir + "/subject/multivariate.pheno.tsv";
    res.qcovar_file = prepared_dir + "/model_qcovar.tsv";
    res.covar_file = prepared_dir + "/model_covar.tsv";
    res.subject_qcovar_file = prepared_dir + "/model_subject_qcovar.tsv";
    res.subject_covar_file = prepared_dir + "/model_subject_covar.tsv";
    std::string model = input_model.empty() ? "repeatability" : input_model;
    std::string pheno_path = res.pheno_file;
    if (model == "rrm") pheno_path = res.rrm_pheno_file;
    else if (model == "lmm") pheno_path = res.subject_pheno_file;
    else if (model == "multivariate") pheno_path = res.multivariate_pheno_file;
    std::ifstream p(pheno_path);
    if (!p) throw std::runtime_error("PREPARED_PHENOTYPE_NOT_FOUND: " + pheno_path + " (use --input-trait for multi-trait prepared input)");
    std::ifstream q(res.qcovar_file);
    if (!q) res.qcovar_file.clear();
    std::ifstream c(res.covar_file);
    if (!c) res.covar_file.clear();
    return res;
}

PreparedInputResult InputPreparer::prepare(const InputPrepareOptions& options,
                                           const std::vector<std::string>& genotype_fids,
                                           const std::vector<std::string>& genotype_iids) {
    if (options.input_model != "repeatability" && options.input_model != "rrm" && options.input_model != "lmm" && options.input_model != "multivariate") {
        throw std::runtime_error("UNSUPPORTED_INPUT_MODEL: Input adaptation currently supports --input-model repeatability, rrm, lmm, or multivariate");
    }
    if (options.input_format == "xlsx" || options.input_format == "xls") {
        throw std::runtime_error("UNSUPPORTED_INPUT_FORMAT: Excel input is not supported in the MVP; export to csv/tsv first");
    }
    ensureDir(options.output_dir);
    Table table = readTable(options.input_file, options.input_format);
    std::unordered_map<std::string, std::pair<std::string, std::string>> raw_to_sample;
    std::unordered_set<std::string> fid_iid_keys;
    for (size_t i = 0; i < genotype_iids.size(); ++i) {
        raw_to_sample.emplace(genotype_iids[i], std::make_pair(genotype_fids[i], genotype_iids[i]));
        raw_to_sample.emplace(genotype_fids[i] + ":" + genotype_iids[i], std::make_pair(genotype_fids[i], genotype_iids[i]));
        fid_iid_keys.insert(key2(genotype_fids[i], genotype_iids[i]));
    }

    if (options.infer_columns) {
        writeColumnInference(options.output_dir + "/column_inference.tsv", table, raw_to_sample);
        writeInferredConfigDraft(options.output_dir + "/inferred_input_config.json", table, raw_to_sample, options.input_file, options.output_dir);
    }
    if (options.infer_columns && options.preview_schema &&
        (options.id_col.empty() || (!options.wide_to_long && options.time_col.empty()) || options.trait_cols.empty())) {
        PreparedInputResult res;
        res.output_dir = options.output_dir;
        res.raw_rows = static_cast<int>(table.rows.size());
        return res;
    }

    requireColumn(table, options.id_col, "--id-col");
    if (!options.wide_to_long) requireColumn(table, options.time_col, "--time-col");
    if (!options.iid_col.empty()) requireColumn(table, options.iid_col, "--iid-col");
    if (options.trait_cols.empty()) throw std::runtime_error("MISSING_REQUIRED_COLUMN: --trait-cols is required for --input-file MVP");
    for (const auto& c : options.trait_cols) requireColumn(table, c, "--trait-cols");
    for (const auto& c : options.record_qcovars) requireColumn(table, c, "--record-qcovars");
    for (const auto& c : options.record_ccovars) requireColumn(table, c, "--record-ccovars");
    for (const auto& c : options.subject_qcovars) requireColumn(table, c, "--subject-qcovars");
    for (const auto& c : options.subject_ccovars) requireColumn(table, c, "--subject-ccovars");

    std::vector<WideSpec> wide_specs;
    std::vector<std::string> trait_names;
    if (options.wide_to_long) {
        std::set<std::string> seen_traits;
        for (const auto& c : options.trait_cols) {
            WideSpec spec = parseWideSpec(c);
            wide_specs.push_back(spec);
            if (!seen_traits.count(spec.trait)) {
                seen_traits.insert(spec.trait);
                trait_names.push_back(spec.trait);
            }
        }
    } else {
        trait_names = options.trait_cols;
    }

    struct Candidate {
        size_t row_index = 0;
        std::string raw_id;
        std::string fid;
        std::string iid;
        std::string raw_time;
        std::string single_trait;
        std::string single_value;
    };
    std::vector<Candidate> candidates;
    std::vector<std::pair<std::string, std::string>> id_report_rows;
    int unmatched = 0;
    for (size_t r = 0; r < table.rows.size(); ++r) {
        const auto& row = table.rows[r];
        std::string raw_id = getCell(table, row, options.id_col);
        std::string fid, iid, match_status;
        if (!options.iid_col.empty()) {
            fid = raw_id;
            iid = getCell(table, row, options.iid_col);
            if (fid_iid_keys.count(key2(fid, iid))) match_status = "MATCHED_FID_IID";
        } else {
            auto it = raw_to_sample.find(raw_id);
            if (it != raw_to_sample.end()) {
                fid = it->second.first;
                iid = it->second.second;
                match_status = "MATCHED_RAW_ID";
            }
        }
        if (match_status.empty()) {
            ++unmatched;
            id_report_rows.push_back({raw_id, "ID_NOT_MATCHED"});
            continue;
        }
        id_report_rows.push_back({raw_id, match_status});
        if (options.wide_to_long) {
            for (const auto& spec : wide_specs) {
                std::string value = getCell(table, row, spec.column);
                if (isMissing(value)) continue;
                candidates.push_back({r, raw_id, fid, iid, spec.raw_time, spec.trait, value});
            }
        } else {
            candidates.push_back({r, raw_id, fid, iid, getCell(table, row, options.time_col), "", ""});
        }
    }

    std::vector<std::string> raw_times;
    raw_times.reserve(candidates.size());
    for (const auto& c : candidates) raw_times.push_back(c.raw_time);
    std::vector<double> bases = normalizeBaseTimes(raw_times, options.time_parser);

    std::unordered_map<std::string, int> duplicate_counts;
    int duplicate_keys = 0;
    int max_duplicate_count = 0;
    for (size_t i = 0; i < candidates.size(); ++i) {
        std::string base = numberToString(bases[i]);
        int n = ++duplicate_counts[key3(candidates[i].fid, candidates[i].iid, base)];
        max_duplicate_count = std::max(max_duplicate_count, n);
        if (n > 1) ++duplicate_keys;
    }
    int time_multiplier = resolveTimeMultiplier(options.time_index_multiplier, max_duplicate_count);

    std::vector<PreparedRecord> records;
    records.reserve(candidates.size());
    std::unordered_map<std::string, int> within_counts;
    std::map<std::string, int> trait_nonmissing;
    for (size_t r = 0; r < candidates.size(); ++r) {
        const auto& cand = candidates[r];
        const auto& row = table.rows[cand.row_index];
        std::string base = numberToString(bases[r]);
        std::string within_key = key3(cand.fid, cand.iid, base);
        int within = within_counts[within_key]++;
        if (within > 0 && options.time_unique_policy == "strict") {
            throw std::runtime_error("DUPLICATE_RECORD_KEY: " + within_key);
        }
        double time_value = bases[r];
        if (options.time_unique_policy == "within-time-index") time_value = bases[r] * time_multiplier + within;
        PreparedRecord rec;
        rec.fid = cand.fid;
        rec.iid = cand.iid;
        rec.raw_id = cand.raw_id;
        rec.raw_time = raw_times[r];
        rec.time_base = base;
        rec.time = numberToString(time_value);
        rec.within_time_index = within;
        rec.raw_row_index = static_cast<int>(cand.row_index) + 2;
        if (options.wide_to_long) {
            rec.values[cand.single_trait] = cand.single_value;
            if (!isMissing(cand.single_value)) trait_nonmissing[cand.single_trait]++;
        } else {
            for (const auto& c : options.trait_cols) {
                rec.values[c] = getCell(table, row, c);
                if (!isMissing(rec.values[c])) trait_nonmissing[c]++;
            }
        }
        for (const auto& c : options.record_qcovars) rec.values[c] = getCell(table, row, c);
        for (const auto& c : options.record_ccovars) rec.values[c] = getCell(table, row, c);
        for (const auto& c : options.subject_qcovars) rec.values[c] = getCell(table, row, c);
        for (const auto& c : options.subject_ccovars) rec.values[c] = getCell(table, row, c);
        records.push_back(std::move(rec));
    }

    if (options.preview_schema) {
        writeSchemaPreview(options.output_dir + "/schema_preview.tsv", table, options,
                           static_cast<int>(candidates.size()), unmatched, duplicate_keys, trait_nonmissing);
        PreparedInputResult res;
        res.output_dir = options.output_dir;
        res.raw_rows = static_cast<int>(table.rows.size());
        res.prepared_records = static_cast<int>(records.size());
        res.matched_records = static_cast<int>(candidates.size());
        res.unmatched_records = unmatched;
        res.selected_trait = options.input_trait.empty() && !trait_names.empty() ? trait_names.front() : options.input_trait;
        return res;
    }

    std::ofstream sample(options.output_dir + "/sample_manifest.tsv");
    sample << "FID\tIID\tRAW_ID\tMATCH_SOURCE\tMATCH_STATUS\n";
    for (const auto& rec : records) sample << rec.fid << "\t" << rec.iid << "\t" << rec.raw_id << "\t" << options.id_col << "\tMATCHED\n";

    std::ofstream record(options.output_dir + "/record_manifest.tsv");
    record << "FID\tIID\tTIME\tRECORD_ID\tRAW_TIME\tTIME_BASE\tWITHIN_TIME_INDEX\tRAW_ROW_INDEX\n";
    for (size_t i = 0; i < records.size(); ++i) {
        const auto& rec = records[i];
        record << rec.fid << "\t" << rec.iid << "\t" << rec.time << "\tR" << (i + 1)
               << "\t" << rec.raw_time << "\t" << rec.time_base << "\t" << rec.within_time_index
               << "\t" << rec.raw_row_index << "\n";
    }

    std::ofstream pheno_long(options.output_dir + "/pheno_long.tsv");
    pheno_long << "FID\tIID\tTIME\tTRAIT\tVALUE\n";
    for (const auto& rec : records) {
        for (const auto& trait : trait_names) {
            auto it = rec.values.find(trait);
            std::string v = it == rec.values.end() ? "NA" : it->second;
            pheno_long << rec.fid << "\t" << rec.iid << "\t" << rec.time << "\t" << trait << "\t" << (isMissing(v) ? "NA" : v) << "\n";
        }
    }

    for (const auto& trait : trait_names) {
        std::ofstream tp(options.output_dir + "/traits/" + trait + ".pheno.tsv");
        tp << "FID\tIID\tTIME\tVALUE\n";
        for (const auto& rec : records) {
            auto it = rec.values.find(trait);
            std::string v = it == rec.values.end() ? "NA" : it->second;
            if (!isMissing(v)) tp << rec.fid << "\t" << rec.iid << "\t" << rec.time << "\t" << v << "\n";
        }
    }

    for (const auto& trait : trait_names) {
        std::ofstream rp(options.output_dir + "/rrm/" + trait + ".pheno.tsv");
        rp << "FID\tIID\tTIME\tVALUE";
        for (const auto& c : options.record_qcovars) rp << "\t" << c;
        for (const auto& c : options.record_ccovars) rp << "\t" << c;
        rp << "\n";
        for (const auto& rec : records) {
            auto it = rec.values.find(trait);
            std::string v = it == rec.values.end() ? "NA" : it->second;
            if (!isMissing(v)) {
                rp << rec.fid << "\t" << rec.iid << "\t" << rec.time_base << "\t" << v;
                for (const auto& c : options.record_qcovars) {
                    auto cv = rec.values.find(c);
                    rp << "\t" << (cv == rec.values.end() || isMissing(cv->second) ? "NA" : cv->second);
                }
                for (const auto& c : options.record_ccovars) {
                    auto cv = rec.values.find(c);
                    rp << "\t" << (cv == rec.values.end() || isMissing(cv->second) ? "NA" : cv->second);
                }
                rp << "\n";
            }
        }
    }

    std::map<std::string, std::string> fid_by_iid;
    std::map<std::string, std::map<std::string, std::pair<double, int>>> subject_trait_stats;
    for (const auto& rec : records) {
        fid_by_iid[rec.iid] = rec.fid;
        for (const auto& trait : trait_names) {
            auto it = rec.values.find(trait);
            if (it == rec.values.end() || isMissing(it->second)) continue;
            double v = 0.0;
            if (!parseDouble(it->second, v)) continue;
            auto& stat = subject_trait_stats[rec.iid][trait];
            stat.first += v;
            stat.second += 1;
        }
    }
    for (const auto& trait : trait_names) {
        std::ofstream sp(options.output_dir + "/subject/" + trait + ".pheno.tsv");
        sp << "FID\tIID\tVALUE\n";
        for (const auto& kv : subject_trait_stats) {
            auto it = kv.second.find(trait);
            if (it == kv.second.end() || it->second.second == 0) continue;
            sp << fid_by_iid[kv.first] << "\t" << kv.first << "\t" << (it->second.first / it->second.second) << "\n";
        }
    }
    std::ofstream mv_subject(options.output_dir + "/subject/multivariate.pheno.tsv");
    mv_subject << "FID\tIID";
    for (const auto& trait : trait_names) mv_subject << "\t" << trait;
    mv_subject << "\n";
    for (const auto& kv : subject_trait_stats) {
        mv_subject << fid_by_iid[kv.first] << "\t" << kv.first;
        for (const auto& trait : trait_names) {
            auto it = kv.second.find(trait);
            if (it == kv.second.end() || it->second.second == 0) mv_subject << "\tNA";
            else mv_subject << "\t" << (it->second.first / it->second.second);
        }
        mv_subject << "\n";
    }

    std::map<std::string, std::set<std::string>> trait_ids;
    for (const auto& rec : records) {
        for (const auto& trait : trait_names) {
            auto it = rec.values.find(trait);
            if (it != rec.values.end() && !isMissing(it->second)) trait_ids[trait].insert(rec.iid);
        }
    }
    std::ofstream trait_manifest(options.output_dir + "/trait_manifest.tsv");
    trait_manifest << "trait\tfile\trows\tunique_ids\n";
    for (const auto& trait : trait_names) {
        trait_manifest << trait << "\t" << options.output_dir << "/traits/" << trait << ".pheno.tsv\t"
                       << trait_nonmissing[trait] << "\t" << trait_ids[trait].size() << "\n";
    }

    std::map<std::string, int> repeats_by_id;
    std::set<std::string> unique_time_keys;
    for (const auto& rec : records) {
        repeats_by_id[rec.iid]++;
        unique_time_keys.insert(key3(rec.fid, rec.iid, rec.time));
    }
    std::vector<int> repeats;
    for (const auto& kv : repeats_by_id) repeats.push_back(kv.second);
    std::sort(repeats.begin(), repeats.end());
    std::ofstream summary(options.output_dir + "/summary.tsv");
    summary << "metric\tvalue\n";
    summary << "matched_ids\t" << repeats_by_id.size() << "\n";
    summary << "matched_rows\t" << records.size() << "\n";
    summary << "min_repeats\t" << (repeats.empty() ? 0 : repeats.front()) << "\n";
    summary << "median_repeats\t" << (repeats.empty() ? 0 : repeats[repeats.size() / 2]) << "\n";
    summary << "max_repeats\t" << (repeats.empty() ? 0 : repeats.back()) << "\n";
    summary << "duplicate_time_keys\t" << (records.size() - unique_time_keys.size()) << "\n";
    summary << "max_same_day_records_per_subject\t" << max_duplicate_count << "\n";
    summary << "time_multiplier\t" << time_multiplier << "\n";
    for (const auto& trait : trait_names) summary << "nonmissing_" << trait << "\t" << trait_nonmissing[trait] << "\n";

    std::ofstream ids(options.output_dir + "/matched_ids.txt");
    for (const auto& kv : repeats_by_id) ids << kv.first << "\n";

    std::ofstream time_scale(options.output_dir + "/time_scale.tsv");
    time_scale << "item\tvalue\n";
    time_scale << "time_definition\tdays_or_time_base * time_multiplier + within_time_index\n";
    time_scale << "time_multiplier\t" << time_multiplier << "\n";
    time_scale << "max_within_time_records_per_subject\t" << max_duplicate_count << "\n";

    std::vector<std::string> all_covars = options.record_qcovars;
    all_covars.insert(all_covars.end(), options.record_ccovars.begin(), options.record_ccovars.end());
    all_covars.insert(all_covars.end(), options.subject_qcovars.begin(), options.subject_qcovars.end());
    all_covars.insert(all_covars.end(), options.subject_ccovars.begin(), options.subject_ccovars.end());
    std::ofstream cov_audit(options.output_dir + "/covariate_audit.tsv");
    cov_audit << "covariate\tlevel\ttype\tnonmissing\tmissing_rate\tunique_values\n";
    for (const auto& cov : all_covars) {
        int miss = 0, nonmiss = 0;
        std::set<std::string> vals;
        for (const auto& rec : records) {
            auto it = rec.values.find(cov);
            std::string v = it == rec.values.end() ? "" : it->second;
            if (isMissing(v)) ++miss;
            else { ++nonmiss; vals.insert(v); }
        }
        std::string level = (std::count(options.subject_qcovars.begin(), options.subject_qcovars.end(), cov) ||
                             std::count(options.subject_ccovars.begin(), options.subject_ccovars.end(), cov)) ? "subject" : "record";
        std::string type = (std::count(options.record_ccovars.begin(), options.record_ccovars.end(), cov) ||
                            std::count(options.subject_ccovars.begin(), options.subject_ccovars.end(), cov)) ? "categorical" : "quantitative";
        int total = miss + nonmiss;
        cov_audit << cov << "\t" << level << "\t" << type << "\t" << nonmiss << "\t"
                  << (total ? (double)miss / total : 0.0) << "\t" << vals.size() << "\n";
    }

    int subject_inconsistent_count = 0;
    std::ofstream subj_check(options.output_dir + "/subject_covar_check.tsv");
    subj_check << "covariate\tIID\tunique_values\tstatus\n";
    std::vector<std::string> subject_covars = options.subject_qcovars;
    subject_covars.insert(subject_covars.end(), options.subject_ccovars.begin(), options.subject_ccovars.end());
    for (const auto& cov : subject_covars) {
        std::map<std::string, std::set<std::string>> vals_by_id;
        for (const auto& rec : records) {
            auto it = rec.values.find(cov);
            if (it != rec.values.end() && !isMissing(it->second)) vals_by_id[rec.iid].insert(it->second);
        }
        for (const auto& kv : vals_by_id) {
            if (kv.second.size() > 1) {
                ++subject_inconsistent_count;
                subj_check << cov << "\t" << kv.first << "\t" << kv.second.size() << "\tMIXED_SUBJECT_COVARIATE\n";
            }
        }
    }

    int validation_warnings = 0;
    std::ofstream validation(options.output_dir + "/validation_report.tsv");
    validation << "code\tseverity\titem\tmessage\n";
    auto add_warning = [&](const std::string& code, const std::string& item, const std::string& msg) {
        ++validation_warnings;
        validation << code << "\twarning\t" << item << "\t" << msg << "\n";
    };
    if (!table.rows.empty()) {
        double unmatched_rate = (double)unmatched / (double)table.rows.size();
        if (unmatched_rate > 0.2) add_warning("ID_MATCH_RATE_LOW", options.id_col, "More than 20% raw rows did not match genotype IDs");
    }
    for (const auto& trait : trait_names) {
        if (trait_nonmissing[trait] == 0) add_warning("TRAIT_ALL_MISSING", trait, "Trait has no non-missing prepared records");
    }
    for (const auto& cov : all_covars) {
        int miss = 0, nonmiss = 0;
        std::set<std::string> vals;
        for (const auto& rec : records) {
            auto it = rec.values.find(cov);
            std::string v = it == rec.values.end() ? "" : it->second;
            if (isMissing(v)) ++miss;
            else { ++nonmiss; vals.insert(v); }
        }
        int total = miss + nonmiss;
        if (total > 0 && (double)miss / total > 0.5) add_warning("COVARIATE_HIGH_MISSING", cov, "Covariate missing rate exceeds 50%");
        bool categorical = std::count(options.record_ccovars.begin(), options.record_ccovars.end(), cov) ||
                           std::count(options.subject_ccovars.begin(), options.subject_ccovars.end(), cov);
        if (categorical && vals.size() > 50) add_warning("CATEGORICAL_TOO_SPARSE", cov, "Categorical covariate has more than 50 levels");
    }
    if (subject_inconsistent_count > 0) add_warning("MIXED_SUBJECT_COVARIATE", "subject_covariates", "At least one subject-level covariate varies within IID");

    writeUnresolvedColumns(options.output_dir + "/unresolved_columns.tsv", table, options, wide_specs);

    writeHeaderAndRows(options.output_dir + "/covar_record_q.tsv", {"FID", "IID", "TIME"}, records, options.record_qcovars, true);
    writeHeaderAndRows(options.output_dir + "/covar_record_c.tsv", {"FID", "IID", "TIME"}, records, options.record_ccovars, true);
    writeHeaderAndRows(options.output_dir + "/covar_subject_q.tsv", {"FID", "IID"}, records, options.subject_qcovars, false);
    writeHeaderAndRows(options.output_dir + "/covar_subject_c.tsv", {"FID", "IID"}, records, options.subject_ccovars, false);

    std::vector<std::string> model_q = options.record_qcovars;
    model_q.insert(model_q.end(), options.subject_qcovars.begin(), options.subject_qcovars.end());
    std::vector<std::string> model_c = options.record_ccovars;
    model_c.insert(model_c.end(), options.subject_ccovars.begin(), options.subject_ccovars.end());
    writeHeaderAndRows(options.output_dir + "/model_qcovar.tsv", {"FID", "IID", "TIME"}, records, model_q, true);
    writeHeaderAndRows(options.output_dir + "/model_covar.tsv", {"FID", "IID", "TIME"}, records, model_c, true);
    writeSubjectModelCovars(options.output_dir + "/model_subject_qcovar.tsv", records, options.subject_qcovars);
    writeSubjectModelCovars(options.output_dir + "/model_subject_covar.tsv", records, options.subject_ccovars);

    std::ofstream idr(options.output_dir + "/id_match_report.tsv");
    idr << "RAW_ID\tSTATUS\n";
    for (const auto& p : id_report_rows) idr << p.first << "\t" << p.second << "\n";

    std::ofstream audit(options.output_dir + "/input_audit.json");
    audit << "{\n";
    audit << "  \"input_file\": \"" << jsonEscape(options.input_file) << "\",\n";
    audit << "  \"input_model\": \"" << jsonEscape(options.input_model) << "\",\n";
    audit << "  \"id_col\": \"" << jsonEscape(options.id_col) << "\",\n";
    audit << "  \"iid_col\": \"" << jsonEscape(options.iid_col) << "\",\n";
    audit << "  \"time_col\": \"" << jsonEscape(options.time_col) << "\",\n";
    audit << "  \"time_parser\": \"" << jsonEscape(options.time_parser) << "\",\n";
    audit << "  \"time_unique_policy\": \"" << jsonEscape(options.time_unique_policy) << "\",\n";
    audit << "  \"time_index_multiplier\": " << time_multiplier << ",\n";
    audit << "  \"wide_to_long\": " << (options.wide_to_long ? "true" : "false") << ",\n";
    audit << "  \"infer_columns\": " << (options.infer_columns ? "true" : "false") << ",\n";
    audit << "  \"duplicate_id_time_keys\": " << duplicate_keys << ",\n";
    audit << "  \"mixed_subject_covariate_warnings\": " << subject_inconsistent_count << ",\n";
    audit << "  \"validation_warnings\": " << validation_warnings << ",\n";
    audit << "  \"raw_rows\": " << table.rows.size() << ",\n";
    audit << "  \"matched_records\": " << candidates.size() << ",\n";
    audit << "  \"unmatched_records\": " << unmatched << ",\n";
    audit << "  \"traits\": {";
    for (size_t i = 0; i < trait_names.size(); ++i) {
        const auto& t = trait_names[i];
        audit << (i ? ", " : "") << "\"" << jsonEscape(t) << "\": " << trait_nonmissing[t];
    }
    audit << "},\n";
    audit << "  \"reports\": [\"summary.tsv\", \"trait_manifest.tsv\", \"covariate_audit.tsv\", \"subject_covar_check.tsv\", \"validation_report.tsv\", \"unresolved_columns.tsv\", \"id_match_report.tsv\"";
    if (options.infer_columns) audit << ", \"column_inference.tsv\", \"inferred_input_config.json\"";
    audit << "]\n";
    audit << "}\n";

    PreparedInputResult res;
    res.output_dir = options.output_dir;
    res.raw_rows = static_cast<int>(table.rows.size());
    res.prepared_records = static_cast<int>(records.size());
    res.matched_records = static_cast<int>(candidates.size());
    res.unmatched_records = unmatched;
    res.selected_trait = options.input_trait.empty() ? trait_names.front() : options.input_trait;
    if (!std::count(trait_names.begin(), trait_names.end(), res.selected_trait)) {
        throw std::runtime_error("TRAIT_NOT_PREPARED: --input-trait is not one of the prepared trait names: " + res.selected_trait);
    }
    res.pheno_file = options.output_dir + "/traits/" + res.selected_trait + ".pheno.tsv";
    res.rrm_pheno_file = options.output_dir + "/rrm/" + res.selected_trait + ".pheno.tsv";
    res.subject_pheno_file = options.output_dir + "/subject/" + res.selected_trait + ".pheno.tsv";
    res.multivariate_pheno_file = options.output_dir + "/subject/multivariate.pheno.tsv";
    if (!model_q.empty()) res.qcovar_file = options.output_dir + "/model_qcovar.tsv";
    if (!model_c.empty()) res.covar_file = options.output_dir + "/model_covar.tsv";
    if (!options.subject_qcovars.empty()) res.subject_qcovar_file = options.output_dir + "/model_subject_qcovar.tsv";
    if (!options.subject_ccovars.empty()) res.subject_covar_file = options.output_dir + "/model_subject_covar.tsv";
    return res;
}

} // namespace cosmic
