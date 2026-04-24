"""Three-panel figure for proteomics abundance-missingness structure.

Panel A: Histogram of per-protein missingness rate for CPTAC ccRCC TMT
         mass spectrometry.
Panel B: Mean abundance binned by missingness rate for CPTAC TMT, with
         95 % CI error bars and a Kruskal-Wallis test.
Panel C: Histogram of per-protein missingness rate for TCGA BRCA/KIRC RPPA.
         Same axes as panel A so the platform contrast reads at a glance —
         a targeted antibody panel has almost no missingness, while bottom-up
         MS has missingness as a first-class feature of the assay.

Reads data/abundance_missingness.csv (combined, with `platform` column) and
writes figures/abundance_missingness.png.
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")  # headless backend; safe inside the container
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats


# Missingness-rate bins used in panel B.
BIN_EDGES = [-0.001, 0.0, 0.2, 0.5, 0.8, 1.01]
BIN_LABELS = ["0% (fully)", "0-20%", "20-50%", "50-80%", "80%+"]


def missingness_hist(ax, df, color, title, letter):
    """Draw a missingness-rate histogram on `ax` with log y-axis."""
    sns.histplot(
        data=df, x="missingness_rate", bins=40,
        color=color, edgecolor="white", ax=ax,
    )
    ax.set_yscale("log")
    ax.set_xlim(-0.02, 1.02)
    ax.set_xlabel("Missingness rate (fraction of samples)")
    ax.set_ylabel("Number of proteins (log scale)")
    ax.set_title(f"{letter}. {title}")
    n_fully = int((df["missingness_rate"] == 0).sum())
    n_partial = int((df["missingness_rate"] > 0).sum())
    ax.text(
        0.98, 0.95,
        (
            f"$n$ = {len(df):,} proteins\n"
            f"fully observed: {n_fully:,}\n"
            f"partially observed: {n_partial:,}"
        ),
        transform=ax.transAxes, ha="right", va="top",
        bbox=dict(facecolor="white", alpha=0.85, edgecolor="0.8"),
    )


def main():
    datadir = sys.argv[1] if len(sys.argv) > 1 else "data"
    figdir = sys.argv[2] if len(sys.argv) > 2 else "figures"
    os.makedirs(figdir, exist_ok=True)

    df_path = os.path.join(datadir, "abundance_missingness.csv")
    print(f"[plot] Reading {df_path} ...", flush=True)
    df = pd.read_csv(df_path)

    # Backwards compatibility: if there's no platform column (old CSV), treat
    # the whole file as CPTAC TMT.
    if "platform" not in df.columns:
        df["platform"] = "TMT"

    cptac = df[df["platform"] == "TMT"].reset_index(drop=True)
    rppa = df[df["platform"] == "RPPA"].reset_index(drop=True)

    cptac["missingness_bin"] = pd.cut(
        cptac["missingness_rate"], bins=BIN_EDGES, labels=BIN_LABELS,
    )

    # Kruskal-Wallis on CPTAC mean abundance across missingness bins.
    groups = [
        cptac.loc[cptac["missingness_bin"] == label, "mean_abundance"].dropna().values
        for label in BIN_LABELS
    ]
    H, p_kw = stats.kruskal(*groups)

    sns.set_theme(style="whitegrid", context="talk")
    fig, (ax_a, ax_b, ax_c) = plt.subplots(
        1, 3, figsize=(22, 7.5), constrained_layout=True,
    )

    # --- Panel A: CPTAC missingness histogram ---
    missingness_hist(
        ax_a, cptac, color="#4c78a8",
        title="CPTAC TMT missingness", letter="A",
    )

    # --- Panel B: CPTAC abundance by missingness bin ---
    palette = sns.color_palette("flare", n_colors=len(BIN_LABELS))
    sns.boxplot(
        data=cptac, x="missingness_bin", y="mean_abundance",
        order=BIN_LABELS, hue="missingness_bin", palette=palette, legend=False,
        fliersize=1.5, linewidth=1.2, ax=ax_b,
    )
    # Overlay bin means with 95 % CI error bars.
    grouped = cptac.groupby("missingness_bin", observed=True)["mean_abundance"]
    means = grouped.mean().reindex(BIN_LABELS)
    sems = grouped.sem().reindex(BIN_LABELS)
    ci95 = 1.96 * sems
    x_pos = list(range(len(BIN_LABELS)))
    ax_b.errorbar(
        x_pos, means.values, yerr=ci95.values,
        fmt="o-", color="black", linewidth=2.5, markersize=10,
        capsize=5, capthick=2, zorder=10, label="bin mean ± 95% CI",
    )
    ax_b.axhline(0, color="0.4", linestyle="--", linewidth=1, zorder=1)
    ax_b.set_xlabel("Missingness rate bin")
    ax_b.set_ylabel("Mean abundance (log-ratio, observed samples)")
    ax_b.set_title("B. CPTAC TMT abundance by bin")
    ax_b.set_ylim(-0.75, 0.75)
    counts = cptac["missingness_bin"].value_counts().reindex(BIN_LABELS).astype(int)
    ax_b.set_xticks(x_pos)
    ax_b.set_xticklabels(
        [f"{lbl}\n$n$={counts[lbl]:,}" for lbl in BIN_LABELS], fontsize=11,
    )
    ax_b.legend(loc="lower right", frameon=True, fontsize=11)
    ax_b.text(
        0.02, 0.98,
        (
            f"Kruskal-Wallis $H$ = {H:,.1f}\n"
            f"$p$ = {p_kw:.2e}"
        ),
        transform=ax_b.transAxes, ha="left", va="top",
        bbox=dict(facecolor="white", alpha=0.85, edgecolor="0.8"),
    )

    # --- Panel C: TCGA RPPA missingness histogram (contrast) ---
    if len(rppa) == 0:
        ax_c.axis("off")
        ax_c.text(
            0.5, 0.5, "No RPPA data in database",
            transform=ax_c.transAxes, ha="center", va="center", fontsize=14,
        )
    else:
        n_rppa_samples = None  # stored only in stats CSV; not needed for the plot
        missingness_hist(
            ax_c, rppa, color="#72b7b2",
            title="TCGA RPPA missingness (contrast)", letter="C",
        )

    fig.suptitle(
        "Platform-specific missingness: CPTAC TMT vs TCGA RPPA",
        fontsize=18,
    )

    out_path = os.path.join(figdir, "abundance_missingness.png")
    fig.savefig(out_path, dpi=150, bbox_inches="tight", pad_inches=0.4)
    plt.close(fig)
    print(f"[plot] Saved {out_path}", flush=True)


if __name__ == "__main__":
    main()
