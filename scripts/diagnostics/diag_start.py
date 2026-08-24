#!/usr/bin/env python3
"""
Isolate why the full fetch returns nothing while diag_assoc worked.
Hypothesis: the 'start' pagination param triggers a 400.

Run on YOUR machine:
    python3 diag_start.py
"""
import requests

API = "https://www.ebi.ac.uk/eqtl/api/v2"
DS = "QTD000356"
REGION = "15:78000000-79500000"

tests = [
    {"pos": REGION, "size": 5},                         # like diag (worked)
    {"pos": REGION, "size": 5, "start": 0},             # + start=0
    {"pos": REGION, "size": 1000},                      # bigger page, no start
    {"pos": REGION, "size": 1000, "start": 1000},       # second page
]

for i, params in enumerate(tests, 1):
    r = requests.get(f"{API}/datasets/{DS}/associations",
                     params=params, timeout=90)
    print(f"\n=== test {i}: {params} ===")
    print("URL:", r.url)
    print("HTTP", r.status_code)
    if r.status_code == 200:
        data = r.json()
        n = len(data) if isinstance(data, list) else "?"
        print(f"OK — {n} rows")
        if isinstance(data, list) and data:
            print("first pos:", data[0]["chromosome"], data[0]["position"])
            print("last  pos:", data[-1]["chromosome"], data[-1]["position"])
    else:
        print("BODY:", r.text[:300])
