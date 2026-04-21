"""Pre-analysis for Aim 1: CPTAC proteomics abundance-missingness structure.

Tests the MNAR / detection-limit hypothesis that low-abundance proteins are
systematically more often missing than high-abundance proteins. Produces:
  - a plot-ready dataframe of per-protein (mean_abundance, missingness_rate, ...)
  - a statistics CSV with Spearman correlation and Kruskal-Wallis results
"""

import sys
import os
import sqlite3
import numpy as np
import pandas as pd
from scipy import stats


def main():
    datadir = sys.argv[1] if len(sys.argv) > 1 else "data"
    db_path = sys.argv[2] if len(sys.argv) > 2 else os.path.join(datadir, "multiomics.db")

    print(f"[analyze] Reading CPTAC proteomics from {db_path} ...", flush=True)
    con = sqlite3.connect(db_path)
    cptac = pd.read_sql_query(
        "SELECT sample_id, protein_id, abundance_value "
        "FROM proteomics_matrix WHERE platform = 'TMT'",
        con,
    )
    con.close()

    if cptac.empty:
        raise RuntimeError(
            "[analyze] No CPTAC (platform='TMT') rows found in proteomics_matrix. "
            "Did script 04 run?"
        )

    n_samples = cptac["sample_id"].nunique()
    n_proteins = cptac["protein_id"].nunique()
    print(
        f"[analyze] Loaded {len(cptac):,} observed rows across "
        f"{n_samples} samples x {n_proteins} proteins",
        flush=True,
    )

    # Per-protein summary: abundance stats + missingness rate
    g = cptac.groupby("protein_id")["abundance_value"]
    summary = pd.DataFrame({
        "n_observed": g.count(),
        "mean_abundance": g.mean(),
        "sd_abundance": g.std(),
    }).reset_index()
    summary["missingness_rate"] = 1.0 - summary["n_observed"] / n_samples

    # Drop any proteins with zero variance or fewer than 3 observations
    # (Spearman needs variance; tertile split needs enough points)
    summary = summary.dropna(subset=["mean_abundance"])
    summary = summary[summary["n_observed"] >= 3].reset_index(drop=True)

    # Tertiles by mean abundance
    summary["abundance_tertile"] = pd.qcut(
        summary["mean_abundance"], q=3,
        labels=["low", "medium", "high"],
    )

    # --- Statistical tests ---
    # Spearman: is abundance negatively correlated with missingness?
    rho, p_rho = stats.spearmanr(summary["mean_abundance"], summary["missingness_rate"])

    # Kruskal-Wallis: do missingness rates differ across abundance tertiles?
    groups = [
        summary.loc[summary["abundance_tertile"] == t, "missingness_rate"].values
        for t in ["low", "medium", "high"]
    ]
    H, p_kw = stats.kruskal(*groups)

    medians = {t: float(np.median(g)) for t, g in zip(["low", "medium", "high"], groups)}

    stats_df = pd.DataFrame([
        {
            "test": "spearman_correlation",
            "comparison": "mean_abundance vs missingness_rate",
            "statistic": float(rho),
            "p_value": float(p_rho),
            "n": int(len(summary)),
            "notes": "negative rho => MNAR / detection-limit structure",
        },
        {
            "test": "kruskal_wallis",
            "comparison": "missingness_rate across low/medium/high abundance tertiles",
            "statistic": float(H),
            "p_value": float(p_kw),
            "n": int(len(summary)),
            "notes": (
                f"median missingness: low={medians['low']:.3f}, "
                f"medium={medians['medium']:.3f}, high={medians['high']:.3f}"
            ),
        },
    ])

    # --- Save outputs ---
    out_df = os.path.join(datadir, "abundance_missingness.csv")
    out_stats = os.path.join(datadir, "stats_results.csv")
    summary.to_csv(out_df, index=False)
    stats_df.to_csv(out_stats, index=False)

    print(f"[analyze] Saved plot-ready dataframe: {out_df} ({len(summary)} proteins)")
    print(f"[analyze] Saved stats results:        {out_stats}")
    print()
    print(stats_df.to_string(index=False))


if __name__ == "__main__":
    main()
