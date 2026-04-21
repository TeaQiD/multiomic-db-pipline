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


def build_samples(clinical_df, rnaseq_samples, proteomics_samples, cptac_sample_ids):
    """Combine all sources into the samples table.

    rnaseq_samples and proteomics_samples are sets collected while streaming
    the large CSVs into SQLite (so we never need to load them into pandas).
    cptac_sample_ids is the set of CPTAC sample_ids seen in the proteomics CSV
    that aren't in clinical_df (so we can add CPTAC-CCRCC rows to samples).
    """
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

    # Add CPTAC samples not already in rows (CPTAC has no TCGA-style clinical join)
    for sid in cptac_sample_ids:
        if sid not in rows:
            rows[sid] = {
                "sample_id": sid,
                "case_id": "",
                "cohort": "CPTAC-CCRCC",
                "cancer_type": "clear cell renal cell carcinoma",
                "has_rnaseq": 0,
                "has_proteomics": 0,
            }

    # Add TCGA RNA-seq samples that have no clinical row (common in test mode
    # when the sampled RNA-seq files don't overlap with the sampled clinical
    # cases; also possible in full mode for samples missing clinical metadata)
    for sid in rnaseq_samples:
        if sid not in rows:
            rows[sid] = {
                "sample_id": sid,
                "case_id": "",
                "cohort": "TCGA",
                "cancer_type": "",
                "has_rnaseq": 0,
                "has_proteomics": 0,
            }

    # Same for RPPA (TCGA proteomics) samples not in clinical/CPTAC
    tcga_prot_samples = proteomics_samples - cptac_sample_ids
    for sid in tcga_prot_samples:
        if sid not in rows:
            rows[sid] = {
                "sample_id": sid,
                "case_id": "",
                "cohort": "TCGA",
                "cancer_type": "",
                "has_rnaseq": 0,
                "has_proteomics": 0,
            }

    # Set has_* flags from the streaming sets
    for sid, rec in rows.items():
        rec["has_rnaseq"] = 1 if sid in rnaseq_samples else 0
        rec["has_proteomics"] = 1 if sid in proteomics_samples else 0

    return pd.DataFrame(list(rows.values()))


def stream_csv_to_table(con, csv_path, table, cols, chunksize=500_000):
    """Stream a CSV into a SQLite table in chunks; return the set of sample_ids seen.

    Uses append mode so the caller is responsible for clearing the table first.
    """
    samples = set()
    total = 0
    for chunk in pd.read_csv(csv_path, chunksize=chunksize, low_memory=False):
        chunk[cols].to_sql(table, con, if_exists="append", index=False)
        samples.update(chunk["sample_id"].astype(str).unique())
        total += len(chunk)
        print(f"[db]   {table}: {total:,} rows loaded ({len(samples)} unique samples so far)")
    return samples, total


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

    # Small CSVs: load fully (a few MB each)
    clinical_df = load_csv(os.path.join(datadir, "tcga_clinical.csv"), "TCGA clinical")
    cptac_clin_df = load_csv(os.path.join(datadir, "cptac_ccrcc_clinical.csv"), "CPTAC clinical")

    print(f"[db] Building database at {db_path} ...")
    # Remove any stale DB so we start clean (to_sql with append would otherwise grow)
    if os.path.exists(db_path):
        os.remove(db_path)
    con = sqlite3.connect(db_path)

    # Bulk-insert tuning — huge speedup for millions of rows
    con.execute("PRAGMA journal_mode = OFF")
    con.execute("PRAGMA synchronous = OFF")
    con.execute("PRAGMA temp_store = MEMORY")
    con.execute("PRAGMA cache_size = -200000")  # ~200 MB page cache

    con.executescript(SCHEMA)

    # --- Stream RNA-seq (the big one) ---
    rnaseq_csv = os.path.join(datadir, "tcga_rnaseq.csv")
    rnaseq_samples = set()
    rnaseq_total = 0
    if os.path.exists(rnaseq_csv):
        print("[db] Streaming rnaseq_matrix ...")
        rnaseq_samples, rnaseq_total = stream_csv_to_table(
            con, rnaseq_csv, "rnaseq_matrix",
            ["sample_id", "gene_id", "expression_value"],
        )
        con.commit()
        print(f"[db] rnaseq_matrix: {rnaseq_total:,} rows total")
    else:
        print(f"[db] WARNING: {rnaseq_csv} not found, skipping rnaseq_matrix")

    # --- Stream proteomics (RPPA + CPTAC into one table) ---
    proteomics_samples = set()
    cptac_sample_ids = set()
    prot_total = 0
    for label, csv_name, is_cptac in [
        ("RPPA", "tcga_rppa.csv", False),
        ("CPTAC proteomics", "cptac_ccrcc_proteomics.csv", True),
    ]:
        path = os.path.join(datadir, csv_name)
        if not os.path.exists(path):
            print(f"[db] WARNING: {path} not found, skipping {label}")
            continue
        print(f"[db] Streaming {label} into proteomics_matrix ...")
        samples, n = stream_csv_to_table(
            con, path, "proteomics_matrix",
            ["sample_id", "protein_id", "abundance_value", "platform"],
        )
        proteomics_samples.update(samples)
        if is_cptac:
            cptac_sample_ids.update(samples)
        prot_total += n
        con.commit()
    print(f"[db] proteomics_matrix: {prot_total:,} rows total")

    # --- Small tables built in memory ---
    manifest_df = build_file_manifest(datadir, clinical_df, None, None, None)
    manifest_df.to_sql("file_manifest", con, if_exists="replace", index=False)
    print(f"[db] file_manifest: {len(manifest_df)} rows")

    samples_df = build_samples(clinical_df, rnaseq_samples, proteomics_samples, cptac_sample_ids)
    samples_df.to_sql("samples", con, if_exists="replace", index=False)
    print(f"[db] samples: {len(samples_df)} rows "
          f"(has_rnaseq={int(samples_df['has_rnaseq'].sum())}, "
          f"has_proteomics={int(samples_df['has_proteomics'].sum())})")

    clin_ann_df = build_clinical_annotations(clinical_df, cptac_clin_df)
    clin_ann_df.to_sql("clinical_annotations", con, if_exists="replace", index=False)
    print(f"[db] clinical_annotations: {len(clin_ann_df)} rows")

    # Empty tables for future aims (created empty via SCHEMA, nothing to do)
    con.commit()

    # Create useful indices for downstream querying
    print("[db] Creating indices ...")
    con.executescript("""
        CREATE INDEX IF NOT EXISTS idx_rnaseq_sample ON rnaseq_matrix(sample_id);
        CREATE INDEX IF NOT EXISTS idx_rnaseq_gene   ON rnaseq_matrix(gene_id);
        CREATE INDEX IF NOT EXISTS idx_prot_sample   ON proteomics_matrix(sample_id);
        CREATE INDEX IF NOT EXISTS idx_prot_protein  ON proteomics_matrix(protein_id);
        CREATE INDEX IF NOT EXISTS idx_clin_sample   ON clinical_annotations(sample_id);
    """)
    con.commit()
    con.close()
    print(f"[db] Done. Database written to {db_path}")


if __name__ == "__main__":
    main()
