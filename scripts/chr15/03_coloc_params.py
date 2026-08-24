#!/usr/bin/env python3
"""
Step 1 closeout: compute the coloc inputs N and s from the FinnGen
manifest case/control counts for J10_PNEUMONIA.

Read num_cases and num_controls from finngen_R13_manifest.tsv
(the J10_PNEUMONIA row), then run:

    python3 03_coloc_params.py <num_cases> <num_controls>

Example (NUMBERS BELOW ARE PLACEHOLDERS — use your real manifest values):
    python3 03_coloc_params.py 25000 475000
"""
import sys

def main(cases, controls):
    cases = int(cases)
    controls = int(controls)
    N = cases + controls
    s = cases / N
    print("== coloc inputs for J10_PNEUMONIA ==")
    print(f"  cases    = {cases:,}")
    print(f"  controls = {controls:,}")
    print(f"  N (total)= {N:,}        # coloc 'N'")
    print(f"  s (case proportion) = {s:.4f}   # coloc 's'")
    print(f"  type     = \"cc\"        # case-control")
    print()
    print("Plug N, s, type into dataset1 of the coloc.abf() call.")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(__doc__); sys.exit(1)
    main(sys.argv[1], sys.argv[2])
