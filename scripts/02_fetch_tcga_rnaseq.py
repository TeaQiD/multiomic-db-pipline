"""Fetch TCGA-BRCA and TCGA-KIRC RNA-seq (STAR Counts) from GDC API."""

import sys
import os
import io
import csv
import time
import json
import tarfile
import requests

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


def parse_star_tsv_to_writer(content: bytes, sample_id: str, writer):
    """Parse a STAR Counts TSV and write rows directly to a csv.writer.

    Returns the number of rows written. Streams line-by-line so we never
    materialize the full gene set in memory.
    """
    n = 0
    text = content.decode("utf-8")
    for line in text.splitlines():
        if not line or line[0] == "#" or line.startswith("N_"):
            continue
        parts = line.split("\t")
        if len(parts) < 8:
            continue
        gene_id = parts[0]
        try:
            fpkm = float(parts[7])  # fpkm_unstranded column
        except (ValueError, IndexError):
            continue
        writer.writerow((sample_id, gene_id, fpkm))
        n += 1
    return n


def download_and_write_batch(file_ids, id_to_sample, writer):
    """POST a batch of file IDs, untar, parse TSVs, stream rows to the CSV writer.

    Returns the number of rows written for this batch.
    """
    payload = json.dumps({"ids": file_ids})
    resp = requests.post(
        DATA_ENDPOINT,
        data=payload,
        headers={"Content-Type": "application/json"},
        timeout=600,
    )
    if resp.status_code != 200:
        raise RuntimeError(
            f"[rnaseq] Data download HTTP {resp.status_code}: {resp.content[:200]!r}"
        )

    raw = io.BytesIO(resp.content)
    batch_rows = 0
    with tarfile.open(fileobj=raw, mode="r:gz") as tar:
        for member in tar:
            if not member.name.endswith(".tsv"):
                continue
            # member.name is "{file_uuid}/{filename}.tsv"
            uuid = member.name.split("/", 1)[0]
            sample_id = id_to_sample.get(uuid, "")
            if not sample_id:
                continue
            fobj = tar.extractfile(member)
            if fobj is None:
                continue
            batch_rows += parse_star_tsv_to_writer(fobj.read(), sample_id, writer)
    # Free the tar.gz bytes ASAP
    del raw
    return batch_rows


def main():
    datadir = sys.argv[1] if len(sys.argv) > 1 else "data"
    os.makedirs(datadir, exist_ok=True)

    print("[rnaseq] Querying GDC file manifest ...", flush=True)
    manifest = query_file_manifest()
    print(f"[rnaseq] Found {len(manifest)} STAR Counts files", flush=True)

    id_to_sample = {m["file_id"]: m["sample_id"] for m in manifest}
    file_ids = list(id_to_sample.keys())

    out_path = os.path.join(datadir, "tcga_rnaseq.csv")
    total_rows = 0
    n_batches = (len(file_ids) + BATCH_SIZE - 1) // BATCH_SIZE

    with open(out_path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(("sample_id", "gene_id", "expression_value"))

        for i in range(0, len(file_ids), BATCH_SIZE):
            batch_idx = i // BATCH_SIZE + 1
            batch = file_ids[i : i + BATCH_SIZE]
            print(
                f"[rnaseq] Batch {batch_idx}/{n_batches}: downloading {len(batch)} files ...",
                flush=True,
            )
            n = download_and_write_batch(batch, id_to_sample, writer)
            total_rows += n
            fh.flush()
            print(
                f"[rnaseq]   Batch {batch_idx}: +{n:,} rows (total {total_rows:,})",
                flush=True,
            )
            time.sleep(SLEEP)

    print(f"[rnaseq] Saved {total_rows:,} rows to {out_path}", flush=True)


if __name__ == "__main__":
    main()
