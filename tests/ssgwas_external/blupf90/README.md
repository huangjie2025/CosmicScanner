# BLUPF90 / postGSf90 readiness kit (Phase 6.6, spec 13/22)

SKIP-safe: the BLUPF90 suite executables (`renf90`, `postGSf90`, `airemlf90`)
are **not** present in this environment, so nothing here is executed by CI.
This kit documents exactly how to run the Phase C comparison once the
reference executables are available.

Reference sources (official): https://nce.ads.uga.edu/html/projects/programs/

## Files

| File | Purpose |
|---|---|
| `renum_inputs.py` | convert CosmicScanner ssGBLUP inputs (PLINK bed/bim/fam + pedigree) to BLUPF90 renumbered input + pedigree files |
| `run_blupf90.sh` | end-to-end runner: renumber, write renf90.par (ssGBLUP, fixed VCs), run renf90 + postGSf90, emit `postGSf90.blp`-derived SNP effects |
| `README.md` (this file) | protocol + expected equivalence checks |

## Protocol (external comparison, Phase C)

1. `python3 renum_inputs.py --bed data.sim --ped data.ped --pheno data_y.phe --out work/`
   produces `work/input.ped`, `work/input.dat`, `work/input.pheno`.
2. `run_blupf90.sh` writes `renf90.par` with:

   ```
   DATAFILE: input.dat
   NUMBER_OF_TRAITS: 1
   PEDIGREE: input.ped
   PED_DEPTH: 3
   MODEL: y = mu + animal
   RANDOM: animal
   OPTION solution animal
   OPTION SNP_FILE 1 2000
   OPTION fix_eff_alpha 0.616532
   OPTION residual 1.85018
   ```

   (fixed VCs so the comparison is VCE-free; swap for `OPTION method VCE` only
   if a VCE comparison is intended.)
3. `renf90` solves the MME and writes `solutions`; `postGSf90` backsolves SNP
   effects to `postGSf90.blp` (columns: CHR SNP a1 a2 freq effect...).
4. Compare:
   - EBV: `solutions` animal effects vs CosmicScanner
     `*.ss_backsolve.ebv.txt` (Pearson / slope / max diff; both tools express
     EBV in `sigma_a` units, so scale must match at fixed VCs).
   - SNP effect: `postGSf90.blp` effect column vs `*.ss_backsolve.snp.txt`
     `EFFECT_ALLELE` (Pearson; absolute scale is not expected to match -
     postGSf90 outputs in residual-variance units, see spec 11/12).
   - Window variance: `postGSf90` window output (`OPTION window 10` /
     `OPTION blupf90_ped` setup) vs CosmicScanner `ss_backsolve.snp.txt`
     WINDOW_VAR / PCT_VAR columns (Definition A; semantic match only).

## Status

Not executed here (reference binaries absent). All equivalences above are
designed from the HIBLUP benchmark result: CosmicScanner EBV == HIBLUP EBV
(Pearson 1.0), so any BLUPF90-side drift is attributable to the BLUPF90
parameterization (e.g. base/genetic-group handling), not to CosmicScanner.
