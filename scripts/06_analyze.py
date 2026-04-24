"""Pre-analysis for Aim 1: proteomics abundance-missingness structure.

Tests the MNAR / detection-limit hypothesis that low-abundance proteins are
systematically more often missing than high-abundance proteins. Runs the same
per-protein summary on both platforms in the database:

  * CPTAC ccRCC TMT mass spectrometry  (platform='TMT')   — bottom-up MS
  * TCGA BRCA/KIRC RPPA antibody panel (platform='RPPA')  — targeted

Emits:
  - data/abundance_missingness.csv  with a `platform` column
  - data/stats_results.csv          with Spearman + Kruskal-Wallis per platform
"""

import sys
import os
import sqlite3
import numpy as np
import pandas as pd
from scipy import stats


def summarize_platform(con, platform):
    """Return (per-protein summary DataFrame, n_samples) for one platform."""
    df = pd.read_sql_query(
        "SELECT sample_id, protein_id, abundance_value "
        "FROM proteomics_matrix WHERE platform = ?",
        con, params=(platform,),
    )
    if df.empty:
        return None, 0

    n_samples = df["sample_id"].nunique()
    n_proteins = df["protein_id"].nunique()
    print(
        f"[analyze] {platform}: {len(df):,} observed rows across "
        f"{n_samples} samples x {n_proteins} proteins",
        flush=True,
    )

    g = df.groupby("protein_id")["abundance_value"]
    summary = pd.DataFrame({
        "n_observed": g.count(),
        "mean_abundance": g.mean(),
        "sd_abundance": g.std(),
    }).reset_index()
    summary["missingness_rate"] = 1.0 - summary["n_observed"] / n_samples
    summary = summary.dropna(subset=["mean_abundance"])
    summary = summary[summary["n_observed"] >= 3].reset_index(drop=True)
    summary["platform"] = platform
    # Tertiles by mean abundance (platform-local)
    summary["abundance_tertile"] = pd.qcut(
        summary["mean_abundance"], q=3,
        labels=["low", "medium", "high"],
    )
    return summary, n_samples


def stats_for(summary, platform, n_samples):
    """Run Spearman + Kruskal-Wallis on one platform's summary."""
    rho, p_rho = stats.spearmanr(
        summary["mean_abundance"], summary["missingness_rate"],
    )
    groups = [
        summary.loc[summary["abundance_tertile"] == t, "missingness_rate"].values
        for t in ["low", "medium", "high"]
    ]
    # Kruskal-Wallis requires at least some variance in each group; if a group
    # is degenerate (all identical — possible for nearly-complete platforms
    # where missingness is 0 for most proteins) scipy raises.
    try:
        H, p_kw = stats.kruskal(*groups)
    except ValueError:
        H, p_kw = float("nan"), float("nan")
    medians = {t: float(np.median(g)) if len(g) else float("nan")
               for t, g in zip(["low", "medium", "high"], groups)}

    return [
        {
            "platform": platform,
            "test": "spearman_correlation",
            "comparison": "mean_abundance vs missingness_rate",
            "statistic": float(rho),
            "p_value": float(p_rho),
            "n": int(len(summary)),
            "n_samples": int(n_samples),
            "notes": "negative rho => MNAR / detection-limit structure",
        },
        {
            "platform": platform,
            "test": "kruskal_wallis",
            "comparison": "missingness_rate across low/medium/high abundance tertiles",
            "statistic": float(H),
            "p_value": float(p_kw),
            "n": int(len(summary)),
            "n_samples": int(n_samples),
            "notes": (
                f"median missingness: low={medians['low']:.3f}, "
                f"medium={medians['medium']:.3f}, high={medians['high']:.3f}"
            ),
        },
    ]


def main():
    datadir = sys.argv[1] if len(sys.argv) > 1 else "data"
    db_path = sys.argv[2] if len(sys.argv) > 2 else os.path.join(datadir, "multiomics.db")

    print(f"[analyze] Reading proteomics from {db_path} ...", flush=True)
    con = sqlite3.connect(db_path)

    pieces = []
    stats_rows = []
    for platform in ("TMT", "RPPA"):
        summary, n_samples = summarize_platform(con, platform)
        if summary is None:
            print(f"[analyze] No rows for platform={platform}; skipping", flush=True)
            continue
        pieces.append(summary)
        stats_rows.extend(stats_for(summary, platform, n_samples))

    con.close()

    if not pieces:
        raise RuntimeError(
            "[analyze] No proteomics rows found for either TMT or RPPA. "
            "Did the fetch scripts run?"
        )

    combined = pd.concat(pieces, ignore_index=True)
    stats_df = pd.DataFrame(stats_rows)

    out_df = os.path.join(datadir, "abundance_missingness.csv")
    out_stats = os.path.join(datadir, "stats_results.csv")
    combined.to_csv(out_df, index=False)
    stats_df.to_csv(out_stats, index=False)

    print(f"[analyze] Saved plot-ready dataframe: {out_df} ({len(combined)} rows)")
    print(f"[analyze] Saved stats results:        {out_stats}")
    print()
    print(stats_df.to_string(index=False))


if __name__ == "__main__":
    main()
