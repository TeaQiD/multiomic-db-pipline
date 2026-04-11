"""Fetch TCGA-BRCA and TCGA-KIRC RNA-seq (STAR Counts) from GDC API."""

import sys
import os
import io
import time
import json
import tarfile
import requests
import pandas as pd

FILES_ENDPOINT = "https://api.gdc.cancer.gov/files"
DATA_ENDPOINT = "https://api.gdc.cancer.gov/data"
PROJECTS = ["TCGA-BRCA", "TCGA-KIRC"]
PAGE_SIZE = 500
BATCH_SIZE = 100
SLEEP = 0.3


def query_file_manifest():
    """Return list of {file_id, file_name, sample_id, project_id} for STAR Counts TSVs."""
    filters = {
        "op": "and",
        "content": [
            {
                "op": "in",
                "content": {"field": "cases.project.project_id", "value": PROJECTS},
            },
            {
                "op": "=",
                "content": {"field": "analysis.workflow_type", "value": "STAR - Counts"},
            },
            {
                "op": "=",
                "content": {"field": "data_type", "value": "Gene Expression Quantification"},
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
        "cases.samples.sample_id",
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
            raise RuntimeError(f"[rnaseq] Files query HTTP {resp.status_code}: {resp.text[:200]}")

        data = resp.json()
        hits = data["data"]["hits"]
        if not hits:
            break

        for f in hits:
            project_id = ""
            sample_id = ""
            cases = f.get("cases", []) or []
            if cases:
                project_id = (cases[0].get("project") or {}).get("project_id", "")
                samples = cases[0].get("samples", []) or []
                if samples:
                    sample_id = samples[0].get("submitter_id", "")

            manifest.append({
                "file_id": f["file_id"],
                "file_name": f["file_name"],
                "project_id": project_id,
                "sample_id": sample_id,
            })

        total = data["data"]["pagination"]["total"]
        from_offset += len(hits)
        print(f"[rnaseq] Manifest: {from_offset}/{total} file records")

        if from_offset >= total:
            break
        time.sleep(SLEEP)

    return manifest


def parse_star_tsv(content: bytes, sample_id: str):
    """Parse a STAR Counts TSV file and return long-format rows."""
    lines = content.decode("utf-8").splitlines()
    rows = []
    for line in lines:
        if line.startswith("#") or line.startswith("N_"):
            continue
        parts = line.split("\t")
        if len(parts) < 8:
            continue
        gene_id = parts[0]
        try:
            fpkm = float(parts[7])  # fpkm_unstranded column
        except (ValueError, IndexError):
            continue
        rows.append({"sample_id": sample_id, "gene_id": gene_id, "expression_value": fpkm})
    return rows


def download_batch(file_ids, id_to_sample):
    """POST a batch of file IDs to GDC /data, untar, parse TSVs."""
    payload = json.dumps({"ids": file_ids})
    resp = requests.post(
        DATA_ENDPOINT,
        data=payload,
        headers={"Content-Type": "application/json"},
        timeout=600,
        stream=True,
    )
    if resp.status_code != 200:
        raise RuntimeError(f"[rnaseq] Data download HTTP {resp.status_code}")

    raw = io.BytesIO(resp.content)
    all_rows = []
    with tarfile.open(fileobj=raw, mode="r:gz") as tar:
        for member in tar.getmembers():
            if not member.name.endswith(".tsv"):
                continue
            # member.name is like "{file_uuid}/{filename}.tsv"
            parts = member.name.split("/")
            file_uuid = parts[0] if len(parts) > 1 else ""
            sample_id = id_to_sample.get(file_uuid, "")
            if not sample_id:
                continue
            fobj = tar.extractfile(member)
            if fobj is None:
                continue
            content = fobj.read()
            all_rows.extend(parse_star_tsv(content, sample_id))
    return all_rows


def main():
    datadir = sys.argv[1] if len(sys.argv) > 1 else "data"
    os.makedirs(datadir, exist_ok=True)

    print("[rnaseq] Querying GDC file manifest ...")
    manifest = query_file_manifest()
    print(f"[rnaseq] Found {len(manifest)} STAR Counts files")

    id_to_sample = {m["file_id"]: m["sample_id"] for m in manifest}
    file_ids = list(id_to_sample.keys())

    all_rows = []
    for i in range(0, len(file_ids), BATCH_SIZE):
        batch = file_ids[i : i + BATCH_SIZE]
        print(f"[rnaseq] Downloading batch {i // BATCH_SIZE + 1} ({len(batch)} files) ...")
        rows = download_batch(batch, id_to_sample)
        all_rows.extend(rows)
        print(f"[rnaseq]   Batch yielded {len(rows)} expression rows (total so far: {len(all_rows)})")
        time.sleep(SLEEP)

    df = pd.DataFrame(all_rows)
    out_path = os.path.join(datadir, "tcga_rnaseq.csv")
    df.to_csv(out_path, index=False, encoding="utf-8")
    print(f"[rnaseq] Saved {len(df)} rows to {out_path}")


if __name__ == "__main__":
    main()
