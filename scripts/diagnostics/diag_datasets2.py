#!/usr/bin/env python3
"""
Diagnostic v2: the /datasets endpoint may paginate. Fetch ALL datasets by
paging, then list every distinct study_label so we find how GTEx is named,
and print the blood/lung GTEx datasets.

Run on YOUR machine:
    python3 diag_datasets2.py
"""
import requests
import json
from collections import Counter

API = "https://www.ebi.ac.uk/eqtl/api/v2"


def get_all_datasets():
    """Page through the datasets endpoint until no more come back."""
    all_ds = []
    size = 1000
    start = 0
    while True:
        r = requests.get(f"{API}/datasets",
                         params={"quant_method": "ge", "size": size,
                                 "start": start},
                         timeout=60)
        if r.status_code != 200:
            print(f"HTTP {r.status_code} at start={start}")
            break
        batch = r.json()
        if not batch:
            break
        all_ds.extend(batch)
        if len(batch) < size:
            break
        start += size
    return all_ds


def main():
    ds = get_all_datasets()
    print(f"Total gene-expression datasets fetched: {len(ds)}\n")

    # All distinct study labels + counts
    studies = Counter(d.get("study_label") for d in ds)
    print("== Distinct study_label values ==")
    for s, c in sorted(studies.items()):
        print(f"  {c:>4}  {s}")

    # Find the GTEx-looking study label
    gtex_labels = [s for s in studies if s and "gtex" in s.lower()]
    print(f"\nGTEx-looking labels: {gtex_labels}")

    # Print GTEx blood + lung datasets
    print("\n== GTEx blood / lung datasets ==")
    for d in ds:
        sl = str(d.get("study_label", "")).lower()
        sg = str(d.get("sample_group", "")).lower()
        tl = str(d.get("tissue_label", "")).lower()
        if "gtex" in sl and ("blood" in sg or "blood" in tl
                             or "lung" in sg or "lung" in tl):
            print({k: d.get(k) for k in
                   ["dataset_id", "study_label", "sample_group",
                    "tissue_label", "condition_label", "sample_size"]})


if __name__ == "__main__":
    main()
