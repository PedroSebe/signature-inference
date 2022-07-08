#!/usr/bin/env python
"""
Generate a BED file describing the genomic regions targeted by MSK-IMPACT based
on a spreadsheet from the supplemental data in the paper available at
<https://doi.org/10.1016/j.jmoldx.2014.12.006>
"""

import pandas as pd
import sys

df = pd.read_excel(sys.argv[1], skiprows=3)
df["Chr"] = "chr" + df["Chr"].astype(str)
df["start"] = df["start"] - 1
bed = df[["Chr","start","stop"]].assign(info=df["Type"]+":"+df["Gene"]+":"+df["Target"])
bed.to_csv(sys.stdout, sep="\t", header=False, index=False)