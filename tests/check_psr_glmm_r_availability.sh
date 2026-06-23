#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RESULT_DIR="${ROOT}/tests/results/psr_gwas_smoke"
OUT="${RESULT_DIR}/psr_glmm_r_availability.tsv"
URL="${PSR_GLMM_R_URL:-http://106.13.144.206/psrglmm/index.html}"
REF_DIR="${ROOT}/external/reference/psr_glmm_r"

mkdir -p "${RESULT_DIR}"
timestamp="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
available="FALSE"
status="NA"
note="unavailable"

if command -v curl >/dev/null 2>&1; then
  tmp="${RESULT_DIR}/psr_glmm_r_index.html.tmp"
  status="$(curl -L --connect-timeout 5 --max-time 15 -o "${tmp}" -w '%{http_code}' -s "${URL}" || true)"
  if [[ "${status}" =~ ^2[0-9][0-9]$ && -s "${tmp}" ]]; then
    available="TRUE"
    note="index reachable; saved read-only copy under external/reference"
    mkdir -p "${REF_DIR}"
    cp "${tmp}" "${REF_DIR}/index.html"
  else
    note="index not reachable or empty; not a CI blocker"
  fi
  rm -f "${tmp}"
else
  note="curl not found; availability not checked"
fi

printf 'timestamp\turl\tavailable\thttp_status\tnote\n' > "${OUT}"
printf '%s\t%s\t%s\t%s\t%s\n' "${timestamp}" "${URL}" "${available}" "${status}" "${note}" >> "${OUT}"

cat > "${ROOT}/tests/PSR_GLMM_R_REFERENCE_STATUS.md" <<EOF
# PSR-GLMM/R Reference Status

Last checked: ${timestamp}

URL: ${URL}

Availability: ${available}

HTTP status: ${status}

Note: ${note}

The official PSR-GLMM/R resource is tracked as an external benchmark target only. It is not bundled with CosmicScanner, not required for CI, and not treated as a blocking dependency because availability and license status are not guaranteed. If the site becomes reachable, downloaded material is stored read-only under \`external/reference/psr_glmm_r/\` for adapter development and manual review.
EOF

echo "[PSR-GWAS] PSR-GLMM/R availability written to ${OUT}"
