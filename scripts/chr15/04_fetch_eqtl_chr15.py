#!/usr/bin/env python3
"""
Phase 2 (working): fetch GTEx v8 eQTLs for the chr15 nicotinic locus
(CHRNA3/5, HYKK) from the EBI eQTL Catalogue REST API v2.

Key constraint discovered from the API:
  - the 'pos' window must be <= 1,000,000 bp. Larger windows are rejected
    with HTTP 400. So we split the target region into <=1 Mb sub-windows,
    fetch each, and concatenate.

Confirmed facts:
  - region param: pos=15:START-END  (chromosome 1-2 chars, NO 'chr')
  - datasets: whole blood QTD000356 (n=670), lung QTD000271 (n=510)
  - an/2 = sample size

Run on YOUR machine:
    python3 04_fetch_eqtl_chr15.py
"""
import os
import time
import requests
import pandas as pd

os.makedirs("results", exist_ok=True)   # create output dir if missing

API = "https://www.ebi.ac.uk/eqtl/api/v2"

CHROM = "15"
START = 78000000
END   = 79500000            # full target region (1.5 Mb -> will be split)
MAX_WINDOW = 1000000        # API hard limit

DATASETS = {
    "blood": "QTD000356",   # whole blood, n=670
    "lung":  "QTD000271",   # lung, n=510
}


def subwindows(start, end, step):
    """Yield (lo, hi) chunks of at most `step` bp covering [start, end]."""
    lo = start
    while lo <= end:
        hi = min(lo + step - 1, end)
        yield lo, hi
        lo = hi + 1


def fetch_window(dataset_id, chrom, lo, hi):
    """Fetch one <=1Mb window, paging with size+start if needed."""
    rows = []
    size = 1000
    start = 0
    region = f"{chrom}:{lo}-{hi}"
    while True:
        params = {"pos": region, "size": size, "start": start}
        r = requests.get(f"{API}/datasets/{dataset_id}/associations",
                         params=params, timeout=120)
        if r.status_code != 200:
            if r.status_code in (400, 404):
                # could be "no more pages" on later pages; stop quietly
                break
            print(f"    HTTP {r.status_code}: {r.text[:150]}")
            break
        batch = r.json()
        if not isinstance(batch, list) or not batch:
            break
        rows.extend(batch)
        if len(batch) < size:
            break
        start += size
        time.sleep(0.5)
    return rows


def main():
    print(f"== Fetch {CHROM}:{START}-{END} in <=1Mb chunks ==")
    for key, ds_id in DATASETS.items():
        print(f"\n  Tissue '{key}' (dataset {ds_id}):")
        all_rows = []
        for lo, hi in subwindows(START, END, MAX_WINDOW):
            print(f"    window {CHROM}:{lo}-{hi} ...", end=" ")
            rows = fetch_window(ds_id, CHROM, lo, hi)
            print(f"{len(rows)} rows")
            all_rows.extend(rows)
            time.sleep(0.5)
        if not all_rows:
            print("    no rows returned.")
            continue

        df = pd.DataFrame(all_rows)
        # safety: keep only rows truly inside the full region
        df["chromosome"] = df["chromosome"].astype(str)
        df["position"] = pd.to_numeric(df["position"], errors="coerce")
        df = df[(df["chromosome"] == CHROM) &
                (df["position"] >= START) & (df["position"] <= END)]

        before = len(df)
        if "rsid" in df.columns:
            df = df.drop(columns=["rsid"])
        df = df.drop_duplicates(subset=["variant", "molecular_trait_id"])
        after = len(df)

        out = f"results/eqtl_chr15_{key}.tsv"
        df.to_csv(out, sep="\t", index=False)
        print(f"    total {before} rows -> {after} after dedupe")
        print(f"    distinct genes: {df['gene_id'].nunique()}")
        print(f"    position range: {int(df['position'].min())}-{int(df['position'].max())}")
        print(f"    saved: {out}")

    print("\nDone. Next: harmonize each eQTL table against the GWAS chr15 slice.")


if __name__ == "__main__":
    main()
