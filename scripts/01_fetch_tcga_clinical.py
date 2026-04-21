"""Fetch TCGA-BRCA and TCGA-KIRC clinical data from GDC API."""

import sys
import os
import time
import json
import requests
import pandas as pd

ENDPOINT = "https://api.gdc.cancer.gov/cases"
PROJECTS = ["TCGA-BRCA", "TCGA-KIRC"]
FIELDS = [
    "case_id",
    "submitter_id",
    "project.project_id",
    "demographic.gender",
    "demographic.age_at_index",
    "diagnoses.ajcc_pathologic_stage",
    "diagnoses.tumor_grade",
    "diagnoses.primary_diagnosis",
    "samples.sample_id",
    "samples.sample_type",
    "samples.submitter_id",
]
PAGE_SIZE = 500
SLEEP = 0.3

TEST_MODE = os.environ.get("PIPELINE_MODE", "test") == "test"
TEST_CASES_PER_PROJECT = 50  # tiny but covers both cohorts


def flatten_case(case):
    """Flatten one GDC case record into one row per sample."""
    rows = []
    project_id = case.get("project", {}).get("project_id", "")
    demographic = case.get("demographic", {}) or {}
    diagnoses = case.get("diagnoses", []) or []
    diag = diagnoses[0] if diagnoses else {}
    samples = case.get("samples", []) or []

    base = {
        "case_id": case.get("case_id", ""),
        "submitter_id": case.get("submitter_id", ""),
        "project_id": project_id,
        "gender": demographic.get("gender", ""),
        "age_at_index": demographic.get("age_at_index", ""),
        "tumor_stage": diag.get("ajcc_pathologic_stage", ""),
        "tumor_grade": diag.get("tumor_grade", ""),
        "primary_diagnosis": diag.get("primary_diagnosis", ""),
    }

    if not samples:
        rows.append({**base, "sample_id": "", "sample_type": "", "sample_submitter_id": ""})
    else:
        for s in samples:
            rows.append({
                **base,
                "sample_id": s.get("sample_id", ""),
                "sample_type": s.get("sample_type", ""),
                "sample_submitter_id": s.get("submitter_id", ""),
            })
    return rows


def fetch_project(project_id):
    mode_note = " [TEST MODE]" if TEST_MODE else ""
    print(f"[clinical] Querying {project_id}{mode_note} ...")
    filters = {
        "op": "=",
        "content": {"field": "project.project_id", "value": project_id},
    }
    all_rows = []
    from_offset = 0
    page_size = TEST_CASES_PER_PROJECT if TEST_MODE else PAGE_SIZE
    limit = TEST_CASES_PER_PROJECT if TEST_MODE else None

    while True:
        params = {
            "filters": json.dumps(filters),
            "fields": ",".join(FIELDS),
            "format": "json",
            "size": page_size,
            "from": from_offset,
        }
        resp = requests.get(ENDPOINT, params=params, timeout=120)
        if resp.status_code != 200:
            raise RuntimeError(f"[clinical] HTTP {resp.status_code} for {project_id}: {resp.text[:200]}")

        data = resp.json()
        hits = data["data"]["hits"]
        if not hits:
            break

        for case in hits:
            all_rows.extend(flatten_case(case))

        total = data["data"]["pagination"]["total"]
        from_offset += len(hits)
        print(f"[clinical]   {project_id}: {from_offset}/{total} cases fetched")

        if from_offset >= total:
            break
        if limit is not None and from_offset >= limit:
            break
        time.sleep(SLEEP)

    return all_rows


def main():
    datadir = sys.argv[1] if len(sys.argv) > 1 else "data"
    os.makedirs(datadir, exist_ok=True)

    all_rows = []
    for project in PROJECTS:
        all_rows.extend(fetch_project(project))

    df = pd.DataFrame(all_rows)
    out_path = os.path.join(datadir, "tcga_clinical.csv")
    df.to_csv(out_path, index=False, encoding="utf-8")
    print(f"[clinical] Saved {len(df)} rows to {out_path}")


if __name__ == "__main__":
    main()
