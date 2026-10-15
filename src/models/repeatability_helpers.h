#pragma once

#include <set>
#include <string>
#include <vector>
#include <Eigen/Dense>
#include "lmm_engine.h"

namespace cosmic {

void writeAssocHeader(const std::string& filename, bool out_ssf, bool eb_lmm = false);
void appendAssocResult(const std::string& filename, const GwasResult& res, bool out_ssf, bool eb_lmm = false);
void appendAssocResults(const std::string& filename, const std::vector<GwasResult>& results, bool out_ssf, bool eb_lmm = false);
void accumulateEbResults(const std::vector<GwasResult>& results, double& effective_tests, long long& tested_markers);
void writeRepeatabilityEbSummary(const std::string& filename, double effective_tests, long long tested_markers);

std::vector<std::string> sortChroms(const std::set<std::string>& chroms_set);
bool uses_per_snp_scaling(int alg);
std::string sanitizeChromLabel(std::string chrom);
std::string getCacheFilePath(const std::string& prefix, const std::string& label);
std::string getCacheMetaPath(const std::string& prefix, const std::string& label);
bool readMatrixFloat32Mmap(const std::string& file, Eigen::MatrixXf& out);
void writeMatrixFloat32(const std::string& file, const Eigen::MatrixXf& mat);
bool readCacheMeta(const std::string& file, double* scaling_factor, std::set<std::string>* observed_chroms,
                   const std::string& expected_identity);
void writeCacheMeta(const std::string& file, double scaling_factor, const std::set<std::string>& observed_chroms,
                    const std::string& identity);

} // namespace cosmic
