#!/usr/bin/env python3
import array
import math
import sys


def read_floats(path):
    values = array.array("f")
    with open(path, "rb") as handle:
        values.frombytes(handle.read())
    if sys.byteorder != "little":
        values.byteswap()
    return values


expected = read_floats(sys.argv[1])
actual = read_floats(sys.argv[2])
tolerance = float(sys.argv[3])

if len(expected) != len(actual):
    raise SystemExit(f"float count differs: {len(expected)} != {len(actual)}")

max_difference = 0.0
for left, right in zip(expected, actual):
    if not math.isfinite(left) or not math.isfinite(right):
        if left != right:
            raise SystemExit("non-finite values differ")
        continue
    max_difference = max(max_difference, abs(left - right))

if max_difference > tolerance:
    raise SystemExit(f"maximum absolute difference {max_difference} exceeds {tolerance}")
