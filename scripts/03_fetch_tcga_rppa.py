"""Fetch TCGA-BRCA and TCGA-KIRC RPPA protein expression data from GDC API."""

import sys
import os
import csv
import time
import json
import requests

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


def download_and_write_rppa(file_id, sample_id, writer, max_retries=5):
    """Download a single RPPA file and write rows directly to csv.writer.
    Returns number of rows written. Retries with exponential backoff on
    transient network errors (read timeouts, connection errors).
    """
    resp = None
    for attempt in range(1, max_retries + 1):
        try:
            resp = requests.get(f"{DATA_ENDPOINT}/{file_id}", timeout=120)
            break
        except (requests.exceptions.Timeout, requests.exceptions.ConnectionError) as e:
            if attempt == max_retries:
                print(
                    f"[rppa]   WARNING: {type(e).__name__} on {file_id} after "
                    f"{max_retries} attempts, skipping: {e}",
                    flush=True,
                )
                return 0
            backoff = 2 ** attempt  # 2, 4, 8, 16, 32 seconds
            print(
                f"[rppa]   {type(e).__name__} on {file_id} (attempt {attempt}/{max_retries}), "
                f"retrying in {backoff}s ...",
                flush=True,
            )
            time.sleep(backoff)

    if resp is None or resp.status_code != 200:
        code = resp.status_code if resp is not None else "no-response"
        print(f"[rppa]   WARNING: HTTP {code} for file {file_id}, skipping", flush=True)
        return 0

    lines = resp.content.decode("utf-8").splitlines()
    if not lines:
        return 0

    header = lines[0].split("\t")
    try:
        protein_col = header.index("peptide_target")
        value_col = header.index("protein_expression")
    except ValueError:
        print(f"[rppa]   WARNING: unexpected header in {file_id}: {header}", flush=True)
        return 0

    n = 0
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
        writer.writerow((sample_id, protein_id, abundance, "RPPA"))
        n += 1
    return n


def main():
    datadir = sys.argv[1] if len(sys.argv) > 1 else "data"
    os.makedirs(datadir, exist_ok=True)

    print("[rppa] Querying GDC RPPA file manifest ...", flush=True)
    manifest = query_rppa_manifest()
    print(f"[rppa] Found {len(manifest)} RPPA files", flush=True)

    out_path = os.path.join(datadir, "tcga_rppa.csv")
    total_rows = 0

    with open(out_path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(("sample_id", "protein_id", "abundance_value", "platform"))

        for i, entry in enumerate(manifest, 1):
            if i % 50 == 0 or i == 1 or i == len(manifest):
                print(
                    f"[rppa] File {i}/{len(manifest)}: {entry['file_name']} (total so far: {total_rows:,})",
                    flush=True,
                )
            total_rows += download_and_write_rppa(
                entry["file_id"], entry["sample_id"], writer
            )
            if i % 100 == 0:
                fh.flush()
            time.sleep(SLEEP)

    print(f"[rppa] Saved {total_rows:,} rows to {out_path}", flush=True)


if __name__ == "__main__":
    main()
