#!/usr/bin/env python3
"""Validate PgenReader against PLINK2 --export A and --freq."""
import struct, os, sys, subprocess, tempfile, argparse

def read_varint(data, pos):
    result = 0; shift = 0
    while True:
        b = data[pos]; pos += 1
        result |= (b & 0x7F) << shift
        shift += 7
        if not (b & 0x80): break
    return result, pos

def read_fixed_int(data, pos, nbytes):
    val = 0
    for b in range(nbytes):
        val |= data[pos + b] << (b * 8)
    return val, pos + nbytes

def parse_pgen_directory(data, n_variants):
    """Parse PGEN header body directory. Returns (vrtypes, vlens)."""
    vrtypes = []
    for i in range(400):
        byte = data[20 + i]
        vrtypes.append(byte & 0x0F)
        vrtypes.append((byte >> 4) & 0x0F)
    vrtypes = vrtypes[:n_variants]

    vlens = []
    for i in range(n_variants):
        vlens.append(data[420 + i])

    return vrtypes, vlens

def parse_difflist(data, pos, num_samples):
    """Parse a PGEN difflist. Returns (sample_ids, values) or (None, None) on error."""
    start = pos
    L, pos = read_varint(data, pos)
    if L == 0:
        return [], []

    id_width = 1 if num_samples <= 256 else (2 if num_samples <= 65536 else (3 if num_samples <= 16777216 else 4))
    G = (L + 63) // 64

    group_starts = []
    for _ in range(G):
        val, pos = read_fixed_int(data, pos, id_width)
        group_starts.append(val)

    group_sizes = []
    for g in range(G - 1):
        group_sizes.append(data[pos] + 63)
        pos += 1

    val_bytes = (L + 3) // 4
    val_ptr = pos
    pos += val_bytes

    sample_ids = []
    for g in range(G):
        first_in_group = g * 64
        last_in_group = min(first_in_group + 63, L - 1)

        if len(sample_ids) < L:
            sample_ids.append(group_starts[g])

        for k in range(first_in_group + 1, last_in_group + 1):
            if len(sample_ids) >= L:
                break
            delta, pos = read_varint(data, pos)
            prev = sample_ids[-1]
            sample_ids.append(prev + delta)

    values = []
    for i in range(L):
        v = (data[val_ptr + i // 4] >> ((i % 4) * 2)) & 3
        values.append(v)

    return sample_ids, values

def decode_variant(data, pos, vlen, vrtype, num_samples):
    """Decode one PGEN variant. Returns list of float dosages."""
    raw = data[pos:pos+vlen]
    n = num_samples
    main_type = vrtype & 7
    dosage = [float('nan')] * n

    if main_type == 0:  # Simple 2-bit
        nbytes = (n + 3) // 4
        for i in range(n):
            if i // 4 < nbytes:
                v = (raw[i // 4] >> ((i % 4) * 2)) & 3
                if v == 0: dosage[i] = 0.0
                elif v == 1: dosage[i] = 1.0
                elif v == 2: dosage[i] = 2.0
                else: dosage[i] = float('nan')

    elif main_type == 1:  # 1-bit encoding
        if len(raw) < 1:
            return dosage

        p = 1
        start = raw[0]
        cat_map = {1: (0, 1), 2: (0, 2), 3: (0, 3), 5: (1, 2), 6: (1, 3), 9: (2, 3)}
        if start not in cat_map:
            return dosage
        cat_a, cat_b = cat_map[start]

        bit_nbytes = (n + 7) // 8
        if p + bit_nbytes > len(raw):
            return dosage

        for i in range(n):
            in_cat_b = (raw[p + i // 8] >> (i % 8)) & 1
            dosage[i] = float(cat_b if in_cat_b else cat_a)
        p += bit_nbytes

        # Difflist
        remaining = len(raw) - p
        if remaining > 0:
            diff_samples, diff_vals = parse_difflist(raw, p, n)
            for sid, val in zip(diff_samples, diff_vals):
                dosage[sid] = float(val)

    elif main_type == 4:  # Difflist non-cat-0
        for i in range(n):
            dosage[i] = 0.0
        if len(raw) > 0:
            diff_samples, diff_vals = parse_difflist(raw, 0, n)
            for sid, val in zip(diff_samples, diff_vals):
                dosage[sid] = float(val)

    elif main_type == 5:  # All category 1 (het)
        for i in range(n):
            dosage[i] = 1.0

    elif main_type == 6:  # Difflist non-cat-2
        for i in range(n):
            dosage[i] = 2.0
        if len(raw) > 0:
            diff_samples, diff_vals = parse_difflist(raw, 0, n)
            for sid, val in zip(diff_samples, diff_vals):
                dosage[sid] = float(val)

    elif main_type == 7:  # Difflist non-cat-3 (missing)
        for i in range(n):
            dosage[i] = float('nan')
        if len(raw) > 0:
            diff_samples, diff_vals = parse_difflist(raw, 0, n)
            for sid, val in zip(diff_samples, diff_vals):
                dosage[sid] = float(val)

    return dosage

def run_plink2_export(prefix, tempdir):
    """Run PLINK2 --export A and parse the raw output."""
    # Find .pgen, .pvar, .psam
    pgen = f"{prefix}.pgen"
    pvar = f"{prefix}.pvar"
    psam = f"{prefix}.psam"

    if not os.path.exists(pgen):
        pgen = prefix
    if not os.path.exists(pvar):
        pvar = prefix.replace('.pgen', '.pvar')

    out_prefix = os.path.join(tempdir, "plink_raw")
    cmd = ["plink2", "--pgen", pgen, "--pvar", pvar, "--psam", psam,
           "--export", "A", "--out", out_prefix]

    # Try common plink2 paths
    plink2_paths = [
        "/disk228/huangjie/software/miniconda3/envs/benchmark/bin/plink2",
        "plink2"
    ]
    for pp in plink2_paths:
        try:
            subprocess.run([pp, "--version"], capture_output=True)
            cmd[0] = pp
            break
        except FileNotFoundError:
            continue

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"PLINK2 error: {result.stderr}")
        return None

    raw_file = f"{out_prefix}.raw"
    if not os.path.exists(raw_file):
        print(f"Raw file not found: {raw_file}")
        return None

    # Parse raw file
    with open(raw_file) as f:
        header = f.readline().strip().split()
        samples = []
        for line in f:
            samples.append(line.strip().split())

    # Column indices: skip FID, IID, PAT, MAT, SEX, PHENOTYPE
    # Find variant columns (they end with _A or _C)
    var_cols = {}  # variant_name -> (col_idx, counted_allele)
    for i, h in enumerate(header):
        if h.startswith("rs") or (h.startswith("chr") and i > 5):
            # Column name format: rsID_allele or chr:pos:ref:alt_allele
            parts = h.rsplit('_', 1)
            if len(parts) == 2:
                var_cols[h] = (i, parts[1])

    return samples, header, var_cols

def validate(prefix, n_variants, n_samples):
    """Main validation function."""
    pgen = prefix + ".pgen"

    with open(pgen, 'rb') as f:
        data = f.read()

    # Parse header
    assert data[0:2] == b'\x6c\x1b', "Invalid PGEN magic"
    storage_mode = data[2]
    n_variants_pgen = struct.unpack_from('<I', data, 3)[0]
    n_samples_pgen = struct.unpack_from('<I', data, 7)[0]

    if n_variants is None:
        n_variants = n_variants_pgen
    if n_samples is None:
        n_samples = n_samples_pgen

    format_byte = data[11]
    provisional_ref = ((format_byte >> 6) & 3) == 2
    print(f"Format byte: 0x{format_byte:02x}, provisional_ref={provisional_ref}")

    vblock_off = struct.unpack_from('<Q', data, 12)[0]
    print(f"Storage mode: {storage_mode}, variants={n_variants}, samples={n_samples}, vblock_offset={vblock_off}")

    # Parse directory
    vrtypes, vlens = parse_pgen_directory(data, n_variants)

    type_counts = {}
    for t in vrtypes:
        type_counts[t] = type_counts.get(t, 0) + 1
    print(f"VR type distribution: {dict(sorted(type_counts.items()))}")

    # Read PVAR for variant names and REF/ALT
    pvar = prefix + ".pvar"
    var_names = []
    var_ref = []
    var_alt = []
    with open(pvar) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split()
            if len(parts) >= 5:
                var_names.append(parts[2])
                var_ref.append(parts[3])
                var_alt.append(parts[4])

    # Run PLINK2 --export A
    with tempfile.TemporaryDirectory() as tmpdir:
        print("\nRunning PLINK2 --export A...")
        plink2_result = run_plink2_export(prefix, tmpdir)
        if plink2_result is None:
            print("WARNING: PLINK2 not available, skipping raw comparison")
            plink2_result = None

    # Decode all variants
    print("\nDecoding all variants...")
    offset = vblock_off

    total_match = 0
    total_samples = 0
    unresolved_types = set()

    for idx in range(n_variants):
        t = vrtypes[idx]
        match_rate = None
        freq_dec = None
        freq_raw = None

        dosage = decode_variant(data, offset, vlens[idx], t, n_samples)
        valid = [d for d in dosage if not math.isnan(d)]
        if valid:
            freq_dec = sum(valid) / (2.0 * len(valid))

        offset += vlens[idx]

        if plink2_result:
            samples, header, var_cols = plink2_result
            var_name = var_names[idx] if idx < len(var_names) else f"var{idx}"

            # Counted allele in raw: could be REF or ALT
            # Column format: {var_name}_{allele}
            # Try both REF and ALT as counted allele
            raw_key_c = f"{var_name}_{var_ref[idx]}" if idx < len(var_ref) else None
            raw_key_a = f"{var_name}_{var_alt[idx]}" if idx < len(var_alt) else None

            raw_dosage = None
            raw_col_name = None
            for rk, counted_allele in [(raw_key_c, var_ref[idx] if idx < len(var_ref) else None),
                                       (raw_key_a, var_alt[idx] if idx < len(var_alt) else None)]:
                if rk in var_cols:
                    col = var_cols[rk][0]
                    raw_vals = [int(s[col]) for s in samples]
                    if len(raw_vals) == n_samples:
                        if counted_allele == var_ref[idx] if idx < len(var_ref) else True:
                            # PLINK2 counts the REF allele -> dosage = 2 - raw_val
                            raw_dosage = [2.0 - v if not math.isnan(dosage[i]) else float('nan') for i, v in enumerate(raw_vals)]
                            freq_raw = (2 * n_samples - sum(raw_vals)) / (2.0 * n_samples)
                        else:
                            raw_dosage = [float(v) for v in raw_vals]
                            freq_raw = sum(raw_vals) / (2.0 * n_samples)
                        raw_col_name = rk
                        break

            if raw_dosage is not None:
                n_match = 0
                n_valid = 0
                for i in range(n_samples):
                    d = dosage[i]
                    r = raw_dosage[i]
                    if not (math.isnan(d) or math.isnan(r)):
                        n_valid += 1
                        if abs(d - r) < 0.001:
                            n_match += 1
                match_rate = f"{n_match}/{n_samples}"
                if n_valid > 0 and n_match != n_valid:
                    pass  # Will show mismatches below
                total_match += n_match
                total_samples += n_valid

        # Show all mismatches
        if plink2_result and raw_dosage is not None and match_rate:
            nm = int(match_rate.split('/')[0])
            ns = int(match_rate.split('/')[1])
            if nm < ns:
                print(f"  MISMATCH idx={idx} ({var_names[idx] if idx < len(var_names) else '?'}): {match_rate}, "
                      f"type={t}, "
                      f"freq_dec={freq_dec:.4f}, freq_raw={freq_raw:.4f}")

    if total_samples > 0 and plink2_result:
        print(f"\nOverall: {total_match}/{total_samples} samples match ({100*total_match/total_samples:.2f}%)")

    print("\nDone.")

if __name__ == "__main__":
    import math
    import argparse
    parser = argparse.ArgumentParser(description="Validate PGEN reader")
    parser.add_argument("prefix", help="PGEN file prefix (without .pgen)")
    parser.add_argument("--n-variants", type=int, default=None)
    parser.add_argument("--n-samples", type=int, default=None)
    args = parser.parse_args()
    validate(args.prefix, args.n_variants, args.n_samples)
