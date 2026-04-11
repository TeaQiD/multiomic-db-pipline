"""Build the multiomics SQLite database from all fetched CSVs."""

import sys
import os
import sqlite3
import datetime
import pandas as pd


SCHEMA = """
CREATE TABLE IF NOT EXISTS file_manifest (
    file_id TEXT PRIMARY KEY,
    cohort TEXT,
    assay TEXT,
    source_url TEXT,
    file_type TEXT,
    retrieved_at TEXT,
    version TEXT
);

CREATE TABLE IF NOT EXISTS samples (
    sample_id TEXT PRIMARY KEY,
    case_id TEXT,
    cohort TEXT,
    cancer_type TEXT,
    has_rnaseq INTEGER,
    has_proteomics INTEGER
);

CREATE TABLE IF NOT EXISTS clinical_annotations (
    sample_id TEXT,
    subtype TEXT,
    grade TEXT,
    stage TEXT,
    age REAL,
    sex TEXT
);

CREATE TABLE IF NOT EXISTS rnaseq_matrix (
    sample_id TEXT,
    gene_id TEXT,
    expression_value REAL
);

CREATE TABLE IF NOT EXISTS proteomics_matrix (
    sample_id TEXT,
    protein_id TEXT,
    abundance_value REAL,
    platform TEXT
);

CREATE TABLE IF NOT EXISTS masking_experiments (
    experiment_id TEXT PRIMARY KEY,
    cohort TEXT,
    missingness_type TEXT,
    masking_rate REAL,
    random_seed INTEGER,
    method TEXT
);

CREATE TABLE IF NOT EXISTS integrated_outputs (
    experiment_id TEXT,
    sample_id TEXT,
    latent_dim INTEGER,
    latent_value REAL,
    predicted_subtype TEXT,
    biomarker_rank INTEGER,
    biomarker_id TEXT
);
"""

NOW = datetime.datetime.now(datetime.timezone.utc).isoformat()


def load_csv(path, label):
    if not os.path.exists(path):
        print(f"[db] WARNING: {path} not found, skipping {label}")
        return None
    df = pd.read_csv(path, low_memory=False)
    print(f"[db] Loaded {len(df)} rows from {path}")
    return df


def build_samples(clinical_df, rnaseq_df, rppa_df, prot_df):
    """Combine all sources into the samples table."""
    rows = {}

    if clinical_df is not None:
        for _, row in clinical_df.iterrows():
            sid = _to_str(row.get("sample_submitter_id")) or _to_str(row.get("sample_id"))
            if not sid:
                continue
            cohort = _to_str(row.get("project_id"))
            cancer_map = {"TCGA-BRCA": "breast invasive carcinoma", "TCGA-KIRC": "kidney renal clear cell carcinoma"}
            rows[sid] = {
                "sample_id": sid,
                "case_id": _to_str(row.get("case_id")),
                "cohort": cohort,
                "cancer_type": cancer_map.get(cohort, cohort),
                "has_rnaseq": 0,
                "has_proteomics": 0,
            }

    if prot_df is not None:
        for sid in prot_df["sample_id"].unique():
            sid = str(sid)
            if sid not in rows:
                rows[sid] = {
                    "sample_id": sid,
                    "case_id": "",
                    "cohort": "CPTAC-CCRCC",
                    "cancer_type": "clear cell renal cell carcinoma",
                    "has_rnaseq": 0,
                    "has_proteomics": 0,
                }

    # Set flags
    rnaseq_samples = set(rnaseq_df["sample_id"].astype(str).unique()) if rnaseq_df is not None else set()
    rppa_samples = set(rppa_df["sample_id"].astype(str).unique()) if rppa_df is not None else set()
    prot_samples = set(prot_df["sample_id"].astype(str).unique()) if prot_df is not None else set()

    for sid, rec in rows.items():
        rec["has_rnaseq"] = 1 if sid in rnaseq_samples else 0
        rec["has_proteomics"] = 1 if (sid in rppa_samples or sid in prot_samples) else 0

    return pd.DataFrame(list(rows.values()))


def build_clinical_annotations(clinical_df, cptac_clin_df):
    rows = []

    if clinical_df is not None:
        for _, row in clinical_df.iterrows():
            sid = _to_str(row.get("sample_submitter_id")) or _to_str(row.get("sample_id"))
            if not sid:
                continue
            rows.append({
                "sample_id": sid,
                "subtype": _to_str(row.get("primary_diagnosis")),
                "grade": _to_str(row.get("tumor_grade")),
                "stage": _to_str(row.get("tumor_stage")),
                "age": _to_float(row.get("age_at_index")),
                "sex": _to_str(row.get("gender")),
            })

    if cptac_clin_df is not None:
        # Prefer exact columns when available; fall back to fuzzy match
        cols_lower = {c.lower(): c for c in cptac_clin_df.columns}
        age_col = cols_lower.get("age") or next((c for c in cptac_clin_df.columns if "age" in c.lower() and "smok" not in c.lower()), None)
        sex_col = cols_lower.get("sex") or cols_lower.get("gender")
        stage_col = cols_lower.get("tumor_stage_pathological") or next((c for c in cptac_clin_df.columns if "stage" in c.lower()), None)
        grade_col = cols_lower.get("histologic_grade") or next((c for c in cptac_clin_df.columns if "grade" in c.lower()), None)
        for _, row in cptac_clin_df.iterrows():
            sid = _to_str(row.get("sample_id"))
            if not sid:
                continue
            rows.append({
                "sample_id": sid,
                "subtype": "",
                "grade": _to_str(row[grade_col]) if grade_col else "",
                "stage": _to_str(row[stage_col]) if stage_col else "",
                "age": _to_float(row[age_col]) if age_col else None,
                "sex": _to_str(row[sex_col]) if sex_col else "",
            })

    return pd.DataFrame(rows)


def build_file_manifest(datadir, clinical_df, rnaseq_df, rppa_df, prot_df):
    entries = []
    sources = [
        ("tcga_clinical", "TCGA", "clinical", "https://api.gdc.cancer.gov/cases", "CSV"),
        ("tcga_rnaseq", "TCGA", "rnaseq", "https://api.gdc.cancer.gov/data", "TSV"),
        ("tcga_rppa", "TCGA", "rppa", "https://api.gdc.cancer.gov/data", "TXT"),
        ("cptac_ccrcc_proteomics", "CPTAC-CCRCC", "proteomics", "https://proteomics.cancer.gov/data-portal", "CSV"),
        ("cptac_ccrcc_clinical", "CPTAC-CCRCC", "clinical", "https://proteomics.cancer.gov/data-portal", "CSV"),
    ]
    for name, cohort, assay, url, ftype in sources:
        path = os.path.join(datadir, f"{name}.csv")
        if os.path.exists(path):
            entries.append({
                "file_id": name,
                "cohort": cohort,
                "assay": assay,
                "source_url": url,
                "file_type": ftype,
                "retrieved_at": NOW,
                "version": "latest",
            })
    return pd.DataFrame(entries)


def _to_float(val):
    try:
        f = float(val)
        if pd.isna(f):
            return None
        return f
    except (TypeError, ValueError):
        return None


def _to_str(val):
    """Convert to string, turning NaN / None / 'nan' into empty string."""
    if val is None:
        return ""
    try:
        if pd.isna(val):
            return ""
    except (TypeError, ValueError):
        pass
    s = str(val).strip()
    if s.lower() in ("nan", "none", "null"):
        return ""
    return s


def main():
    datadir = sys.argv[1] if len(sys.argv) > 1 else "data"
    db_path = sys.argv[2] if len(sys.argv) > 2 else os.path.join(datadir, "multiomics.db")
    os.makedirs(datadir, exist_ok=True)

    # Load CSVs
    clinical_df = load_csv(os.path.join(datadir, "tcga_clinical.csv"), "TCGA clinical")
    rnaseq_df = load_csv(os.path.join(datadir, "tcga_rnaseq.csv"), "TCGA RNA-seq")
    rppa_df = load_csv(os.path.join(datadir, "tcga_rppa.csv"), "TCGA RPPA")
    prot_df = load_csv(os.path.join(datadir, "cptac_ccrcc_proteomics.csv"), "CPTAC proteomics")
    cptac_clin_df = load_csv(os.path.join(datadir, "cptac_ccrcc_clinical.csv"), "CPTAC clinical")

    print(f"[db] Building database at {db_path} ...")
    con = sqlite3.connect(db_path)
    con.executescript(SCHEMA)

    # file_manifest
    manifest_df = build_file_manifest(datadir, clinical_df, rnaseq_df, rppa_df, prot_df)
    manifest_df.to_sql("file_manifest", con, if_exists="replace", index=False)
    print(f"[db] file_manifest: {len(manifest_df)} rows")

    # samples
    samples_df = build_samples(clinical_df, rnaseq_df, rppa_df, prot_df)
    samples_df.to_sql("samples", con, if_exists="replace", index=False)
    print(f"[db] samples: {len(samples_df)} rows")

    # clinical_annotations
    clin_ann_df = build_clinical_annotations(clinical_df, cptac_clin_df)
    clin_ann_df.to_sql("clinical_annotations", con, if_exists="replace", index=False)
    print(f"[db] clinical_annotations: {len(clin_ann_df)} rows")

    # rnaseq_matrix
    if rnaseq_df is not None:
        rnaseq_df[["sample_id", "gene_id", "expression_value"]].to_sql(
            "rnaseq_matrix", con, if_exists="replace", index=False
        )
        print(f"[db] rnaseq_matrix: {len(rnaseq_df)} rows")

    # proteomics_matrix (RPPA + CPTAC)
    prot_parts = []
    if rppa_df is not None:
        prot_parts.append(rppa_df[["sample_id", "protein_id", "abundance_value", "platform"]])
    if prot_df is not None:
        prot_parts.append(prot_df[["sample_id", "protein_id", "abundance_value", "platform"]])
    if prot_parts:
        prot_combined = pd.concat(prot_parts, ignore_index=True)
        prot_combined.to_sql("proteomics_matrix", con, if_exists="replace", index=False)
        print(f"[db] proteomics_matrix: {len(prot_combined)} rows")

    # Empty tables for future aims
    con.execute("DELETE FROM masking_experiments")
    con.execute("DELETE FROM integrated_outputs")
    con.commit()

    con.close()
    print(f"[db] Done. Database written to {db_path}")


if __name__ == "__main__":
    main()
