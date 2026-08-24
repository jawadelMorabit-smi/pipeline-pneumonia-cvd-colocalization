#!/usr/bin/env python3
"""
Diagnostic: inspect what the eQTL Catalogue API actually returns for datasets,
so we can fix the blood/lung filter in 04_fetch_eqtl_chr15.py.

Run on YOUR machine:
    python3 diag_datasets.py
"""
import requests
import json

API = "https://www.ebi.ac.uk/eqtl/api/v2"


def try_call(label, params):
    print(f"\n=== {label} : params={params} ===")
    r = requests.get(f"{API}/datasets", params=params, timeout=60)
    print(f"HTTP {r.status_code}")
    if r.status_code != 200:
        print(r.text[:300])
        return []
    data = r.json()
    print(f"returned {len(data)} datasets")
    return data


def main():
    # 1) What keys does a dataset record even have?
    all_ge = try_call("all gene-expression datasets", {"quant_method": "ge"})
    if all_ge:
        print("\nKeys in first record:")
        print(list(all_ge[0].keys()))
        print("\nFirst record in full:")
        print(json.dumps(all_ge[0], indent=2)[:800])

    # 2) Try to find GTEx records and print their study/tissue fields
    print("\n\n=== Scanning for GTEx + blood/lung across all ge datasets ===")
    hits = 0
    for d in all_ge:
        study = str(d.get("study_label", "")).lower()
        blob = json.dumps(d).lower()
        if "gtex" in study or "gtex" in blob:
            if any(k in blob for k in ["blood", "lung"]):
                hits += 1
                # print the identifying fields we care about
                print({k: d.get(k) for k in
                       ["dataset_id", "study_id", "study_label",
                        "sample_group", "tissue_label", "tissue_id",
                        "quant_method", "sample_size"]})
                if hits >= 12:
                    break
    if hits == 0:
        print("No GTEx blood/lung matches found by scanning. "
              "Print a few raw records to see labels:")
        for d in all_ge[:5]:
            print(json.dumps(d, indent=2)[:400])
            print("---")


if __name__ == "__main__":
    main()
