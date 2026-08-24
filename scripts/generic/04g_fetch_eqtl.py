#!/usr/bin/env python3
"""
Fetch GTEx v8 eQTLs for ANY locus from the EBI eQTL Catalogue REST API v2.
Generic version: pass the locus name, chromosome, and window on the command line.

Datasets: whole blood QTD000356 (n=670), lung QTD000271 (n=510).
API constraint: window must be <= 1,000,000 bp (auto-split if larger).

Lancer depuis la racine Project/. Sortie : results/<locus>/eqtl_<locus>_<tissu>.tsv

Usage:
    python3 scripts/generic/04g_fetch_eqtl.py <locus_name> <chrom> <start> <end>
Example (CRP locus, chr1):
    python3 scripts/generic/04g_fetch_eqtl.py CRP 1 159200000 160200000
"""
import os
import sys
import time
import requests
import pandas as pd

API = "https://www.ebi.ac.uk/eqtl/api/v2"
MAX_WINDOW = 1000000

DATASETS = {
    "blood": "QTD000356",   # whole blood, n=670
    "lung":  "QTD000271",   # lung, n=510
    "liver": "QTD000266",   # liver, n=208 (added for CRP/HNF1A - D3)
}


def subwindows(start, end, step):
    lo = start
    while lo <= end:
        hi = min(lo + step - 1, end)
        yield lo, hi
        lo = hi + 1


def fetch_window(dataset_id, chrom, lo, hi):
    rows = []
    size = 1000
    start = 0
    region = f"{chrom}:{lo}-{hi}"
    while True:
        r = requests.get(f"{API}/datasets/{dataset_id}/associations",
                         params={"pos": region, "size": size, "start": start},
                         timeout=120)
        if r.status_code != 200:
            if r.status_code in (400, 404):
                break
            print(f"      HTTP {r.status_code}: {r.text[:150]}")
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


def main(locus, chrom, start, end):
    chrom = str(chrom)
    start = int(start)
    end = int(end)
    os.makedirs(f"results/{locus}", exist_ok=True)
    print(f"== Locus {locus}: fetch {chrom}:{start}-{end} in <=1Mb chunks ==")
    for key, ds_id in DATASETS.items():
        print(f"\n  Tissue '{key}' (dataset {ds_id}):")
        all_rows = []
        for lo, hi in subwindows(start, end, MAX_WINDOW):
            print(f"    window {chrom}:{lo}-{hi} ...", end=" ")
            rows = fetch_window(ds_id, chrom, lo, hi)
            print(f"{len(rows)} rows")
            all_rows.extend(rows)
            time.sleep(0.5)
        if not all_rows:
            print("    no rows returned.")
            continue

        df = pd.DataFrame(all_rows)
        df["chromosome"] = df["chromosome"].astype(str)
        df["position"] = pd.to_numeric(df["position"], errors="coerce")
        df = df[(df["chromosome"] == chrom) &
                (df["position"] >= start) & (df["position"] <= end)]

        before = len(df)
        if "rsid" in df.columns:
            df = df.drop(columns=["rsid"])
        df = df.drop_duplicates(subset=["variant", "molecular_trait_id"])
        after = len(df)

        out = f"results/{locus}/eqtl_{locus}_{key}.tsv"
        df.to_csv(out, sep="\t", index=False)
        print(f"    {before} rows -> {after} after dedupe")
        print(f"    distinct genes: {df['gene_id'].nunique()}")
        print(f"    position range: {int(df['position'].min())}-{int(df['position'].max())}")
        print(f"    saved: {out}")

    print(f"\nDone. Next: extract GWAS {locus} slice, then harmonize.")


if __name__ == "__main__":
    if len(sys.argv) < 5:
        print(__doc__); sys.exit(1)
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
