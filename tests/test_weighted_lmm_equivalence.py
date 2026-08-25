#!/usr/bin/env python3
import math
from pathlib import Path

import numpy as np


def chisq1_p(stat):
    return math.erfc(math.sqrt(max(0.0, stat) / 2.0))


def gls_scan(y, X, g, V):
    C = np.column_stack([X, g])
    Vinv = np.linalg.inv(V)
    XtViX = C.T @ Vinv @ C
    beta = np.linalg.solve(XtViX, C.T @ Vinv @ y)
    cov = np.linalg.inv(XtViX)
    b = beta[-1]
    se = math.sqrt(cov[-1, -1])
    wald = (b / se) ** 2
    return b, se, wald, chisq1_p(wald)


def main():
    rng = np.random.default_rng(20260623)
    n = 24
    x1 = rng.normal(size=n)
    X = np.column_stack([np.ones(n), x1])
    Gs = rng.integers(0, 3, size=(n, 8)).astype(float)
    Gs -= Gs.mean(axis=0)
    g = Gs[:, 0]
    A = rng.normal(size=(n, n))
    K = A @ A.T
    K /= np.mean(np.diag(K))
    beta = np.array([0.4, -0.2])
    u = rng.multivariate_normal(np.zeros(n), 0.35 * K)
    y = X @ beta + 0.25 * g + u + rng.normal(scale=0.3, size=n)
    Rdiag = 0.2 + 1.5 * rng.random(n)
    Wsqrt = np.diag(1.0 / np.sqrt(Rdiag))

    # Layer 1: weighted GLS equivalence with explicit R.
    b1, se1, w1, p1 = gls_scan(y, X, g, np.diag(Rdiag))
    b2, se2, w2, p2 = gls_scan(Wsqrt @ y, Wsqrt @ X, Wsqrt @ g, np.eye(n))

    # Layer 2/3: weighted LMM equivalence with fixed V = sg*K + R.
    sg = 0.7
    V = sg * K + np.diag(Rdiag)
    Vw = sg * (Wsqrt @ K @ Wsqrt) + np.eye(n)
    b3, se3, w3, p3 = gls_scan(y, X, g, V)
    b4, se4, w4, p4 = gls_scan(Wsqrt @ y, Wsqrt @ X, Wsqrt @ g, Vw)

    max_scan_diff = 0.0
    for k in range(Gs.shape[1]):
        gw = Wsqrt @ Gs[:, k]
        vals_a = gls_scan(y, X, Gs[:, k], V)
        vals_b = gls_scan(Wsqrt @ y, Wsqrt @ X, gw, Vw)
        max_scan_diff = max(max_scan_diff, max(abs(a - b) for a, b in zip(vals_a, vals_b)))

    rows = [
        ("weighted_gls", abs(b1 - b2), abs(se1 - se2), abs(w1 - w2), abs(p1 - p2)),
        ("weighted_lmm_fixed_v", abs(b3 - b4), abs(se3 - se4), abs(w3 - w4), abs(p3 - p4)),
        ("weighted_lmm_scan_fixed_vc", max_scan_diff, max_scan_diff, max_scan_diff, max_scan_diff),
    ]
    out = Path("/disk228/huangjie/PanOmics/CosmicScanner/tests/results/psr_gwas_smoke/weighted_lmm_equivalence.tsv")
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8") as f:
        f.write("test\tbeta_abs_diff\tse_abs_diff\twald_abs_diff\tp_abs_diff\tpass\n")
        for name, db, ds, dw, dp in rows:
            ok = db < 1e-8 and ds < 1e-8 and dw < 1e-8 and dp < 1e-8
            f.write(f"{name}\t{db:.12g}\t{ds:.12g}\t{dw:.12g}\t{dp:.12g}\t{int(ok)}\n")
            if not ok:
                raise SystemExit(f"{name} failed weighted equivalence tolerance")
    print(f"[PSR-GWAS] Weighted LMM equivalence summary written to {out}")


if __name__ == "__main__":
    main()
