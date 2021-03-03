#!/usr/bin/env python3
# coding: utf-8

import pandas as pd
import sys, os
from pathlib import Path
from pysradb.search import SraSearch

genus = sys.argv[1]
if len(sys.argv) > 2:
    num_results = sys.argv[2]
else:
    num_results = 2

path = Path("/scratch/jho5ze/bionets/operons/SRA/genera/")

if not os.path.exists(path / genus):
    os.mkdir(path / genus)

path = path / genus

query = SraSearch(return_max=num_results, query=f"{genus} RNA-Seq")
query.search()
query.get_df().to_csv(path/f"{genus}_query_info.csv")

