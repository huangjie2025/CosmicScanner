# Type-I Error 校准报告

**模型**: maternal, threshold, gxe
**重复次数**: 10

## 1. λ_GC 基因组膨胀因子

理想值 = 1.0。判定：PASS [0.95, 1.05] / WARN [0.90, 1.10] / FAIL 其他

| Model | Mean | SD | Min | Max | Status |
|-------|-----:|---:|----:|----:|:------:|
| maternal | 1.0014 | 0.0537 | 0.9162 | 1.0841 | PASS |
| threshold | 1.0030 | 0.0209 | 0.9666 | 1.0326 | PASS |
| gxe | 1.0212 | 0.0228 | 0.9780 | 1.0530 | PASS |

## 2. Type-I Error Rate

理想值 = α。判定：PASS ≤ 2α / WARN (2α, 3α] / FAIL > 3α

### α = 0.05

| Model | Mean | SE | Expected | Ratio | Status |
|-------|-----:|---:|---------:|------:|:------:|
| maternal | 0.048550 | 0.000862 | 5.00e-02 | 0.971 | PASS |
| threshold | 0.049500 | 0.001017 | 5.00e-02 | 0.990 | PASS |
| gxe | 0.052900 | 0.000908 | 5.00e-02 | 1.058 | PASS |

### α = 0.01

| Model | Mean | SE | Expected | Ratio | Status |
|-------|-----:|---:|---------:|------:|:------:|
| maternal | 0.009650 | 0.000495 | 1.00e-02 | 0.965 | PASS |
| threshold | 0.009500 | 0.000520 | 1.00e-02 | 0.950 | PASS |
| gxe | 0.011150 | 0.000579 | 1.00e-02 | 1.115 | PASS |

### α = 0.001

| Model | Mean | SE | Expected | Ratio | Status |
|-------|-----:|---:|---------:|------:|:------:|
| maternal | 0.001050 | 0.000217 | 1.00e-03 | 1.050 | PASS |
| threshold | 0.001250 | 0.000203 | 1.00e-03 | 1.250 | PASS |
| gxe | 0.000950 | 0.000304 | 1.00e-03 | 0.950 | PASS |

### α = 5e-08

| Model | Mean | SE | Expected | Ratio | Status |
|-------|-----:|---:|---------:|------:|:------:|
| maternal | 0.000000 | 0.000000 | 5.00e-08 | 0.000 | PASS |
| threshold | 0.000000 | 0.000000 | 5.00e-08 | 0.000 | PASS |
| gxe | 0.000000 | 0.000000 | 5.00e-08 | 0.000 | PASS |

## 3. KS 检验（p-value 均匀性）

H0: p-values ~ U(0,1)。P > 0.05 表示不拒绝均匀性

| Model | Mean KS-P | Min KS-P | % P>0.05 |
|-------|----------:|---------:|---------:|
| maternal | 0.6055 | 0.1144 | 100.00% |
| threshold | 0.6793 | 0.3222 | 100.00% |
| gxe | 0.5350 | 0.0587 | 100.00% |

## 4. 总结

- **maternal**: λ_GC=1.0014 ✓  Type-I ✓  → **PASS**
- **threshold**: λ_GC=1.0030 ✓  Type-I ✓  → **PASS**
- **gxe**: λ_GC=1.0212 ✓  Type-I ✓  → **PASS**

**总体判定**: ✓ ALL PASS

QQ 数据见: /tmp/type1_ext_final/qq_data.csv
原始汇总见: /tmp/type1_ext_final/type1_summary.csv
