#!/usr/bin/env python3
"""Plot timing results for gf2^8_mult from the detailed sketch summary CSV."""

import argparse
import csv
from pathlib import Path

try:
    import matplotlib.pyplot as plt
except ModuleNotFoundError as exc:
    raise SystemExit(
        "matplotlib is required. Install it in your Python environment, "
        "for example: python3 -m pip install matplotlib"
    ) from exc


DEFAULT_SUMMARY = (
    Path(__file__).resolve().parents[1]
    / "forme"
    / "sketch_experiments"
    / "aa_sketch_detail_20260618_032309"
    / "summary.csv"
)


def as_float(row, key, default=0.0):
    value = row.get(key, "")
    if value == "":
        return default
    return float(value)


def as_int(row, key):
    value = row.get(key, "")
    return int(float(value)) if value != "" else None


def load_rows(summary_csv, circuit):
    with summary_csv.open(newline="") as f:
        rows = list(csv.DictReader(f))

    selected = [
        row
        for row in rows
        if row.get("circuit") == circuit and row.get("variant") in {"reuse", "noreuse"}
    ]
    selected.sort(key=lambda row: (as_int(row, "p"), row["variant"]))
    if not selected:
        raise ValueError(f"No rows found for circuit={circuit!r} in {summary_csv}")
    return selected


def load_baseline(summary_csv, circuit):
    with summary_csv.open(newline="") as f:
        for row in csv.DictReader(f):
            if row.get("circuit") == circuit and row.get("variant") == "baseline":
                return row
    return None


def split_by_variant(rows):
    by_variant = {"reuse": [], "noreuse": []}
    for row in rows:
        by_variant[row["variant"]].append(row)
    return by_variant


def common_p_values(by_variant, variants):
    p_sets = [{as_int(row, "p") for row in by_variant[variant]} for variant in variants]
    common = sorted(set.intersection(*p_sets))
    if not common:
        raise ValueError(f"No common p values found for variants: {', '.join(variants)}")
    return common


def set_zero_based_ylim(ax, max_value, pad=0.08):
    upper = max_value * (1.0 + pad) if max_value > 0 else 1.0
    ax.set_ylim(bottom=0, top=upper)


def plot_total_time(rows, baseline, out_base):
    by_variant = split_by_variant(rows)
    max_time = 0.0

    fig, ax = plt.subplots(figsize=(7.0, 4.4), constrained_layout=True)
    styles = {
        "reuse": {"label": "with basis reuse", "marker": "o", "linewidth": 2.2},
        "noreuse": {
            "label": "without basis reuse",
            "marker": "s",
            "linewidth": 2.2,
            "linestyle": "--",
        },
    }

    for variant in ["reuse", "noreuse"]:
        data = by_variant[variant]
        if not data:
            raise ValueError(f"No rows found for variant={variant!r}")
        ps = [as_int(row, "p") for row in data]
        times = [as_float(row, "total_exec_s") for row in data]
        max_time = max(max_time, max(times))
        ax.plot(ps, times, **styles[variant])

    if baseline:
        baseline_time = as_float(baseline, "total_exec_s")
        max_time = max(max_time, baseline_time)
        ax.axhline(
            baseline_time,
            color="0.35",
            linewidth=1.5,
            linestyle=":",
            label=f"baseline ({baseline_time:.1f} s)",
        )

    ax.set_title("gf2^8_mult: total execution time")
    ax.set_xlabel("pruning parameter p")
    ax.set_ylabel("total execution time [s]")
    set_zero_based_ylim(ax, max_time)
    ax.grid(True, axis="y", alpha=0.3)
    ax.legend(frameon=False)

    fig.savefig(out_base.with_name(out_base.name + "_total_time.png"), dpi=300)
    fig.savefig(out_base.with_name(out_base.name + "_total_time.pdf"))
    plt.close(fig)


def plot_pruning_breakdown(rows, out_base, variant_mode):
    components = [
        ("sketch_build_ms", "build sketch"),
        ("sketch_basis_ms", "basis work"),
        ("sketch_mem_ms", "memory check"),
        ("sketch_other_ms", "other sketch"),
    ]
    colors = ["#4C78A8", "#72B7B2", "#F58518", "#B279A2"]
    by_variant = split_by_variant(rows)

    variants = ["reuse", "noreuse"] if variant_mode == "both" else [variant_mode]
    ps = common_p_values(by_variant, variants)

    fig, ax = plt.subplots(figsize=(8.0, 4.6), constrained_layout=True)

    if len(variants) == 1:
        width = 0.62
        x_positions = {p: i for i, p in enumerate(ps)}
        tick_positions = list(x_positions.values())
        tick_labels = [str(p) for p in ps]
    else:
        width = 0.36
        x_positions = {}
        tick_positions = []
        tick_labels = []
        for i, p in enumerate(ps):
            x_positions[(p, "reuse")] = i - width / 2
            x_positions[(p, "noreuse")] = i + width / 2
            tick_positions.append(i)
            tick_labels.append(str(p))

    for variant in variants:
        data_by_p = {as_int(row, "p"): row for row in by_variant[variant]}
        bottoms = [0.0] * len(ps)
        for (column, label), color in zip(components, colors):
            values = [as_float(data_by_p[p], column) / 1000.0 for p in ps]
            if len(variants) == 1:
                xs = [x_positions[p] for p in ps]
                legend_label = label
            else:
                xs = [x_positions[(p, variant)] for p in ps]
                legend_label = label if variant == variants[0] else None
            ax.bar(xs, values, width, bottom=bottoms, label=legend_label, color=color)
            bottoms = [bottom + value for bottom, value in zip(bottoms, values)]
    max_stacked_time = max(
        sum(as_float(row, column) / 1000.0 for column, _ in components)
        for variant in variants
        for row in by_variant[variant]
    )

    if len(variants) == 2:
        for p in ps:
            ax.text(
                x_positions[(p, "reuse")],
                -0.06,
                "R",
                ha="center",
                va="top",
                transform=ax.get_xaxis_transform(),
                fontsize=9,
            )
            ax.text(
                x_positions[(p, "noreuse")],
                -0.06,
                "N",
                ha="center",
                va="top",
                transform=ax.get_xaxis_transform(),
                fontsize=9,
            )
        subtitle = "R: with reuse, N: without reuse"
    else:
        subtitle = "with basis reuse" if variant_mode == "reuse" else "without basis reuse"

    ax.set_title(f"gf2^8_mult: pruning-time breakdown ({subtitle})")
    ax.set_xlabel("pruning parameter p")
    ax.set_ylabel("pruning time [s]")
    ax.set_xticks(tick_positions, tick_labels)
    set_zero_based_ylim(ax, max_stacked_time)
    ax.grid(True, axis="y", alpha=0.3)
    ax.legend(frameon=False, ncol=2)

    suffix = "_pruning_breakdown" if variant_mode == "both" else f"_pruning_breakdown_{variant_mode}"
    fig.savefig(out_base.with_name(out_base.name + suffix + ".png"), dpi=300)
    fig.savefig(out_base.with_name(out_base.name + suffix + ".pdf"))
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--summary", type=Path, default=DEFAULT_SUMMARY)
    parser.add_argument("--circuit", default="gf2^8_mult.tfc")
    parser.add_argument("--out-dir", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument(
        "--breakdown-variant",
        choices=["both", "reuse", "noreuse"],
        default="both",
        help="Which variant to show in the pruning-time stacked bar chart.",
    )
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    out_base = args.out_dir / "gf2_8_mult_timing"

    rows = load_rows(args.summary, args.circuit)
    baseline = load_baseline(args.summary, args.circuit)

    plot_total_time(rows, baseline, out_base)
    plot_pruning_breakdown(rows, out_base, args.breakdown_variant)

    print(f"Wrote figures to {args.out_dir}")


if __name__ == "__main__":
    main()
