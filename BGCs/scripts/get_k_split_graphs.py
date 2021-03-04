#!/usr/bin/env
import pandas as pd
import random
import subprocess
import os


smash = pd.read_csv("../misc_data/antismash_db_results.csv", sep="\t")
# smash[smash["NCBI accession"] == "NC_010612.1"] 
counts = smash["NCBI accession"].value_counts()
above_10_index = counts[counts > 5].index
smash = smash[smash["NCBI accession"].isin(above_10_index)]
smash["NCBI accession"] = smash["NCBI accession"].apply(lambda row: row.split(".")[0])
smash_dict = dict()
for accession, df in smash.groupby("NCBI accession"):
    smash_dict[accession] = df[["From", "To"]].to_numpy()
smash_nums = sorted([(k, len(value)) for k, value in smash_dict.items()], key = lambda tup: tup[1], reverse=True)
random.seed(42)
random.shuffle(smash_nums)
top = [i[0] for i in smash_nums[:10]]
