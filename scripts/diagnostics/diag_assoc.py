#!/usr/bin/env python3
"""
Probe the associations endpoint to see WHY it returns HTTP 400,
and find the region-parameter format the API accepts.

Run on YOUR machine:
    python3 diag_assoc.py
"""
import requests

API = "https://www.ebi.ac.uk/eqtl/api/v2"
DS = "QTD000356"   # GTEx whole blood

# Try several ways of specifying the region / paging.
attempts = [
    {"pos": "15:78500000-78600000", "size": 5},
    {"pos": "chr15:78500000-78600000", "size": 5},
    {"chromosome": "15", "bp_lower": 78500000, "bp_upper": 78600000, "size": 5},
    {"chromosome": 15, "position": 78520813, "size": 5},
    {"gene_id": "ENSG00000166669", "size": 5},   # HYKK, as a sanity fallback
    {"size": 5},                                  # no filter at all
]

for i, params in enumerate(attempts, 1):
    url = f"{API}/datasets/{DS}/associations"
    r = requests.get(url, params=params, timeout=60)
    print(f"\n=== attempt {i}: {params} ===")
    print("URL:", r.url)
    print("HTTP", r.status_code)
    if r.status_code == 200:
        data = r.json()
        n = len(data) if isinstance(data, list) else "?"
        print(f"OK — {n} rows")
        if isinstance(data, list) and data:
            print("keys:", list(data[0].keys()))
            print("first row:", data[0])
    else:
        # show the error body so we learn what the API wants
        print("BODY:", r.text[:400])
