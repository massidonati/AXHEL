#!/usr/bin/env python3
# Copyright (C) 2026 University of Pisa - Dept. of Information Engineering
# SPDX-License-Identifier: Apache-2.0

import argparse
import csv
from pathlib import Path

KEYS = [
    "kernel", "poly_degree", "qi_bits", "modulus",
    "mod_factor", "input_mod_factor", "output_mod_factor",
]

METRICS = [
    "samples_total", "samples_kept", "outliers_removed",
    "raw_min_cycles", "raw_max_cycles", "q1_cycles", "q3_cycles", "iqr_cycles",
    "lower_fence_cycles", "upper_fence_cycles", "min_cycles", "max_cycles",
    "mean_cycles", "median_cycles", "stddev_cycles", "cv_percent",
    "mean_cycles_per_coeff", "median_cycles_per_coeff", "perf_scaled",
]


def load(path: Path):
    with path.open(newline="") as f:
        rows = list(csv.DictReader(f))
    return {tuple(row[k] for k in KEYS): row for row in rows}


def require_backend(rows, expected, label):
    wrong = sorted({row.get("backend", "") for row in rows.values() if row.get("backend") != expected})
    if wrong:
        raise SystemExit(f"{label} summary contains unexpected backend labels: {wrong}")


def as_float(row, name):
    return float(row[name])


def ratio(num, den):
    return f"{num / den:.6f}" if den else ""


def reduction(reference, optimized):
    return f"{100.0 * (1.0 - optimized / reference):.6f}" if reference else ""


def main():
    ap = argparse.ArgumentParser(
        description="Merge AXHEL ScalarStrict/Native/SVE summaries and compute vectorization speedups"
    )
    ap.add_argument("scalar_summary", type=Path)
    ap.add_argument("native_summary", type=Path)
    ap.add_argument("sve_summary", type=Path)
    ap.add_argument("output", type=Path)
    args = ap.parse_args()

    scalar = load(args.scalar_summary)
    native = load(args.native_summary)
    sve = load(args.sve_summary)
    require_backend(scalar, "ScalarStrict", "strict-scalar")
    require_backend(native, "Native", "native")
    require_backend(sve, "SVE", "SVE")

    keysets = {"scalar": set(scalar), "native": set(native), "sve": set(sve)}
    union = set().union(*keysets.values())
    if not all(keys == union for keys in keysets.values()):
        detail = ", ".join(f"missing_{name}={len(union - keys)}" for name, keys in keysets.items())
        raise SystemExit(f"ScalarStrict/Native/SVE scenario mismatch: {detail}")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fields = KEYS[:]
    for prefix in ("scalar", "native", "sve"):
        fields.extend(f"{prefix}_{m}" for m in METRICS)
    fields += [
        # Primary comparison: explicit SVE versus an environment with compiler
        # auto-vectorization disabled.
        "speedup_median_scalar_over_sve",
        "speedup_mean_scalar_over_sve",
        "median_cycle_reduction_sve_vs_scalar_percent",
        "mean_cycle_reduction_sve_vs_scalar_percent",
        # Secondary comparison: how much the compiler's normal Release policy
        # improves the public Native implementation by itself.
        "speedup_median_scalar_over_native",
        "speedup_mean_scalar_over_native",
        "median_cycle_reduction_native_vs_scalar_percent",
        "mean_cycle_reduction_native_vs_scalar_percent",
        # Existing comparison: explicit AXHEL SVE versus normal Native.
        "speedup_median_native_over_sve",
        "speedup_mean_native_over_sve",
        "median_cycle_reduction_sve_vs_native_percent",
        "mean_cycle_reduction_sve_vs_native_percent",
    ]

    def sort_key(key):
        row = scalar[key]
        return (row["kernel"], int(row["poly_degree"]), int(row["qi_bits"]))

    with args.output.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for key in sorted(union, key=sort_key):
            sc, n, sv = scalar[key], native[key], sve[key]
            row = {k: sc[k] for k in KEYS}
            for prefix, source in (("scalar", sc), ("native", n), ("sve", sv)):
                for metric in METRICS:
                    row[f"{prefix}_{metric}"] = source[metric]

            scmed, nmed, svmed = (as_float(r, "median_cycles") for r in (sc, n, sv))
            scmean, nmean, svmean = (as_float(r, "mean_cycles") for r in (sc, n, sv))

            row["speedup_median_scalar_over_sve"] = ratio(scmed, svmed)
            row["speedup_mean_scalar_over_sve"] = ratio(scmean, svmean)
            row["median_cycle_reduction_sve_vs_scalar_percent"] = reduction(scmed, svmed)
            row["mean_cycle_reduction_sve_vs_scalar_percent"] = reduction(scmean, svmean)

            row["speedup_median_scalar_over_native"] = ratio(scmed, nmed)
            row["speedup_mean_scalar_over_native"] = ratio(scmean, nmean)
            row["median_cycle_reduction_native_vs_scalar_percent"] = reduction(scmed, nmed)
            row["mean_cycle_reduction_native_vs_scalar_percent"] = reduction(scmean, nmean)

            row["speedup_median_native_over_sve"] = ratio(nmed, svmed)
            row["speedup_mean_native_over_sve"] = ratio(nmean, svmean)
            row["median_cycle_reduction_sve_vs_native_percent"] = reduction(nmed, svmed)
            row["mean_cycle_reduction_sve_vs_native_percent"] = reduction(nmean, svmean)
            writer.writerow(row)

    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
