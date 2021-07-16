import pandas as pd
import numpy as np
from scipy.spatial.distance import pdist
import sys
sys.path.append('/scratch/jho5ze/bionets/covid/scripts')
import entropy_utils as eu
from screed import ScreedDB
msadb = ScreedDB("../data/msa_0625/usa_msa_0625.fasta")

def map_chars(line):
    char_map = {"A":0, "C":1, "T":2, "G":3}
    chars = []
    for c in line:
        if c in char_map:
            chars.append(char_map[c])
        else:
            chars.append(4)
            
    return chars

def position_difference(u, v):
    return sum(u == v) / u.shape[0]

def calculate_jaccard_from_accessions(accessions):
    sample_arr = []
    for line in eu.msa_from_screed_ids(accessions):
        if line[0] == ">":
            continue
        else:
            line_arr = map_chars(line)
            sample_arr.append(line_arr)
    sample_arr = np.array(sample_arr)
    x = pdist(sample_arr, position_difference)
    return sum(x) / x.shape[0]

location, start, end, uncert_file = sys.argv[1:5]
samples = sys.argv[5:]
jaccard = calculate_jaccard_from_accessions(samples)
with open(uncert_file, "w") as dest:
    dest.write(f"location,start,end,samples,jaccard\n")
    dest.write(f"{location},{start},{end},{len(samples)},{jaccard}\n")
