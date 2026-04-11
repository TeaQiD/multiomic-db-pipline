"""Fetch TCGA-BRCA and TCGA-KIRC RPPA protein expression data from GDC API."""

import sys
import os
import io
import time
import json
import requests
import pandas as pd

FILES_ENDPOINT = "https://api.gdc.cancer.gov/files"
DATA_ENDPOINT = "https://api.gdc.cancer.gov/data"
PROJECTS = ["TCGA-BRCA", "TCGA-KIRC"]
PAGE_SIZE = 500
SLEEP = 0.3


def query_rppa_manifest():
    """Return list of {file_id, file_name, sample_id} for RPPA files."""
    filters = {
        "op": "and",
        "content": [
            {
                "op": "in",
                "content": {"field": "cases.project.project_id", "value": PROJECTS},
            },
            {
                "op": "=",
                "content": {
                    "field": "data_type",
                    "value": "Protein Expression Quantification",
                },
            },
            {
                "op": "=",
                "content": {"field": "access", "value": "open"},
            },
        ],
    }
    fields = [
        "file_id",
        "file_name",
        "cases.project.project_id",
        "cases.samples.submitter_id",
    ]

    manifest = []
    from_offset = 0

    while True:
        params = {
            "filters": json.dumps(filters),
            "fields": ",".join(fields),
            "format": "json",
            "size": PAGE_SIZE,
            "from": from_offset,
        }
        resp = requests.get(FILES_ENDPOINT, params=params, timeout=120)
        if resp.status_code != 200:
            raise RuntimeError(f"[rppa] Files query HTTP {resp.status_code}: {resp.text[:200]}")

        data = resp.json()
        hits = data["data"]["hits"]
        if not hits:
            break

        for f in hits:
            sample_id = ""
            cases = f.get("cases", []) or []
            if cases:
                samples = cases[0].get("samples", []) or []
                if samples:
                    sample_id = samples[0].get("submitter_id", "")

            manifest.append({
                "file_id": f["file_id"],
                "file_name": f["file_name"],
                "sample_id": sample_id,
            })

        total = data["data"]["pagination"]["total"]
        from_offset += len(hits)
        print(f"[rppa] Manifest: {from_offset}/{total} file records")

        if from_offset >= total:
            break
        time.sleep(SLEEP)

    return manifest


def download_and_parse_rppa(file_id, sample_id):
    """Download a single RPPA file and return long-format rows.

    The GDC RPPA TSV has a header line followed by rows with columns:
    AGID, lab_id, catalog_number, set_id, peptide_target, protein_expression.
    We use peptide_target as protein_id and protein_expression as abundance.
    """
    resp = requests.get(
        f"{DATA_ENDPOINT}/{file_id}",
        timeout=120,
    )
    if resp.status_code != 200:
        print(f"[rppa]   WARNING: HTTP {resp.status_code} for file {file_id}, skipping")
        return []

    lines = resp.content.decode("utf-8").splitlines()
    if not lines:
        return []

    header = lines[0].split("\t")
    try:
        protein_col = header.index("peptide_target")
        value_col = header.index("protein_expression")
    except ValueError:
        print(f"[rppa]   WARNING: unexpected header in {file_id}: {header}")
        return []

    rows = []
    for line in lines[1:]:
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) <= max(protein_col, value_col):
            continue
        protein_id = parts[protein_col]
        try:
            abundance = float(parts[value_col])
        except ValueError:
            continue
        rows.append({
            "sample_id": sample_id,
            "protein_id": protein_id,
            "abundance_value": abundance,
            "platform": "RPPA",
        })
    return rows


def main():
    datadir = sys.argv[1] if len(sys.argv) > 1 else "data"
    os.makedirs(datadir, exist_ok=True)

    print("[rppa] Querying GDC RPPA file manifest ...")
    manifest = query_rppa_manifest()
    print(f"[rppa] Found {len(manifest)} RPPA files")

    all_rows = []
    for i, entry in enumerate(manifest):
        print(f"[rppa] Downloading file {i + 1}/{len(manifest)}: {entry['file_name']} ...")
        rows = download_and_parse_rppa(entry["file_id"], entry["sample_id"])
        all_rows.extend(rows)
        time.sleep(SLEEP)

    df = pd.DataFrame(all_rows)
    out_path = os.path.join(datadir, "tcga_rppa.csv")
    df.to_csv(out_path, index=False, encoding="utf-8")
    print(f"[rppa] Saved {len(df)} rows to {out_path}")


if __name__ == "__main__":
    main()
