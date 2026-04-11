"""Fetch CPTAC ccRCC proteomics and clinical data via the cptac Python package.

The cptac package (>=1.5) exposes ccRCC data through several sources. We use:
  - proteomics: umich (TMT-labeled, protein-level abundance)
  - clinical:   mssm  (pan-cancer harmonized clinical metadata)

Both get_* calls require an explicit `source=` because the package's default-
selection code path is broken on Python 3.13 (generator vs len()).
"""

import sys
import os
import warnings

warnings.filterwarnings("ignore")

import pandas as pd  # noqa: E402
import cptac  # noqa: E402


def main():
    datadir = sys.argv[1] if len(sys.argv) > 1 else "data"
    os.makedirs(datadir, exist_ok=True)

    print("[cptac] Loading CPTAC ccRCC dataset (downloads on first use) ...")
    cc = cptac.Ccrcc()

    # --- Proteomics (umich TMT) ---
    print("[cptac] Extracting proteomics (source=umich, TMT) ...")
    prot = cc.get_proteomics(source="umich")

    # Columns are a MultiIndex ("Name", "Database_ID"); flatten to "NAME|ENSP..."
    if isinstance(prot.columns, pd.MultiIndex):
        prot.columns = [
            f"{name}|{db}" if db and str(db) != "nan" else str(name)
            for name, db in prot.columns
        ]

    prot.index.name = "sample_id"
    prot_long = (
        prot.reset_index()
        .melt(id_vars="sample_id", var_name="protein_id", value_name="abundance_value")
        .dropna(subset=["abundance_value"])
    )
    prot_long["platform"] = "TMT"

    prot_out = os.path.join(datadir, "cptac_ccrcc_proteomics.csv")
    prot_long.to_csv(prot_out, index=False, encoding="utf-8")
    print(f"[cptac] Saved {len(prot_long)} proteomics rows to {prot_out}")

    # --- Clinical (mssm harmonized) ---
    print("[cptac] Extracting clinical annotations (source=mssm) ...")
    clin = cc.get_clinical(source="mssm")
    clin.index.name = "sample_id"
    clin_out = os.path.join(datadir, "cptac_ccrcc_clinical.csv")
    clin.reset_index().to_csv(clin_out, index=False, encoding="utf-8")
    print(f"[cptac] Saved {len(clin)} clinical rows to {clin_out}")


if __name__ == "__main__":
    main()
