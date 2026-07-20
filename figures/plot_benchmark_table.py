from pathlib import Path

import matplotlib.pyplot as plt


BENCHMARKS = [
    ("adder_8", 24, 399),
    ("barenco_tof_10", 19, 224),
    ("csla_mux_3", 16, 70),
    ("csum_mux_9", 30, 196),
    ("gf2^5_mult", 15, 175),
    ("gf2^6_mult", 18, 252),
    ("gf2^7_mult", 21, 343),
    ("gf2^8_mult", 24, 448),
    ("gf2^9_mult", 27, 567),
    ("gf2^10_mult", 30, 700),
    ("gf2^16_mult", 48, 1792),
    ("gf2^32_mult", 96, 7168),
    ("ham15-high", 20, 2457),
    ("ham15-med", 17, 574),
    ("ham15-low", 15, 161),
    ("mod_adder_1024", 28, 1995),
    ("qcla_adder_10", 36, 238),
    ("qcla_com_7", 24, 203),
    ("qcla_mod_7", 26, 413),
    ("tof_10", 19, 119),
]


def add_table(ax, rows):
    ax.axis("off")
    values = [[name, f"{n}", f"{t:,}"] for name, n, t in rows]
    table = ax.table(
        cellText=values,
        colLabels=["Benchmark", "n", "T gates"],
        colWidths=[0.62, 0.15, 0.23],
        cellLoc="left",
        colLoc="left",
        bbox=[0, 0, 1, 1],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(16)

    for (row, col), cell in table.get_celld().items():
        cell.set_edgecolor("#8b98a5")
        cell.set_linewidth(0.8)
        if row == 0:
            cell.set_facecolor("#234f70")
            cell.get_text().set_color("white")
            cell.get_text().set_weight("bold")
        else:
            cell.set_facecolor("#eef3f6" if row % 2 == 0 else "white")
        if col > 0:
            cell.get_text().set_ha("right")


def main():
    fig, axes = plt.subplots(1, 2, figsize=(16, 7.8))
    add_table(axes[0], BENCHMARKS[:10])
    add_table(axes[1], BENCHMARKS[10:])
    fig.suptitle("Benchmark circuits (n >= 15)", fontsize=24, fontweight="bold", y=0.97)
    fig.subplots_adjust(left=0.035, right=0.965, top=0.88, bottom=0.04, wspace=0.08)

    output = Path(__file__).resolve().parent / "benchmark_circuits_n15_table"
    fig.savefig(output.with_suffix(".png"), dpi=300, facecolor="white")
    fig.savefig(output.with_suffix(".pdf"), facecolor="white")


if __name__ == "__main__":
    main()
