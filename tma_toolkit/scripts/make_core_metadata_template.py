#!/usr/bin/env python
"""Write an empty core-metadata sheet for a TMA grid (e.g. 10 rows x 20 columns = 200 cores).

    python scripts/make_core_metadata_template.py --slide-id TMA2 --rows 10 --cols 20 \
        --extra grade stage msi_status --out config/core_metadata_TMA2.csv

Open the CSV in Excel, fill one row per core (patient_id, organ, tissue_type, ...), set
include=FALSE for cores to leave out, and keep it as CSV. Several slides: run once per
slide and concatenate, or pass --append to add to an existing sheet.
Then check it:  python scripts/make_core_metadata_template.py --check config/core_metadata_TMA2.csv
"""

import argparse
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from tma_toolkit.metadata import make_template, read_core_metadata  # noqa: E402


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--slide-id")
    p.add_argument("--tma-id", help="prefix of core ids (default = slide id)")
    p.add_argument("--rows", type=int)
    p.add_argument("--cols", type=int)
    p.add_argument("--extra", nargs="*", default=[], help="additional metadata columns")
    p.add_argument("--out")
    p.add_argument("--append", action="store_true")
    p.add_argument("--check", help="validate an existing sheet and print a summary")
    a = p.parse_args()

    if a.check:
        df = read_core_metadata(a.check)
        print(f"{len(df)} cores, {df.slide_id.nunique()} slide(s)")
        for c in ("organ", "tissue_type", "diagnosis", "include"):
            if c in df:
                print(df[c].value_counts(dropna=False).to_string(), "\n")
        if "patient_id" in df:
            print(f"{df.patient_id.nunique()} patients; cores per patient: "
                  f"{df.groupby('patient_id').size().describe()[['min', 'mean', 'max']].round(1).to_dict()}")
        return

    t = make_template(a.slide_id, a.rows, a.cols, a.tma_id, a.extra)
    if a.append and os.path.exists(a.out):
        t = pd.concat([pd.read_csv(a.out), t], ignore_index=True)
    t.to_csv(a.out, index=False)
    print(f"wrote {len(t)} rows -> {a.out}")


if __name__ == "__main__":
    main()
