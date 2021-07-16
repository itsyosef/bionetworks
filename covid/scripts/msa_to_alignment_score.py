import pandas as pd
import numpy as np
from scipy.spatial.distance import pdist
from math import comb
import sys
sys.path.append('/scratch/jho5ze/bionets/covid/scripts')
import entropy_utils as eu
from screed import ScreedDB
msadb = ScreedDB("../data/msa_0625/usa_msa_0625.fasta")

def score_seqs(seq1, seq2):
    if len(seq1) != len(seq2):
        raise ValueError("Sequences must be the same length")
    score = 0
    gap_length = 0
    gap_seq = "1"
    shared_dash = 0
    for pos, char1 in enumerate(seq1):
        char2 = seq2[pos]
        if char1 == char2 == "-":
            shared_dash += 1
            continue
        elif char1 != "-" and char2 != "-":
            if char1 == char2:
                score +=1
            if gap_length > 0:
                score -= 5 + gap_length
                gap_length = 0
        elif char1 == "-":
            if gap_seq == "1":
                gap_length += 1
            elif gap_length > 0:
                score -= 5 + gap_length
                gap_length = 1
                gap_seq = "1"
        elif char2 == "-":
            if gap_seq == "2":
                gap_length += 1
            elif gap_length > 0:
                score -= 5 + gap_length
                gap_length = 1
                gap_seq = "2"

    return score / (len(seq1) - shared_dash)

def calculate_alignment_score_from_accessions(accessions):
    sample_arr = []
    for line in eu.msa_from_screed_ids(accessions):
        if line[0] == ">":
            continue
        else:
            sample_arr.append(line)
            
    total_score = 0
    for ix, seq1 in enumerate(sample_arr):
        for ij, seq2 in enumerate(sample_arr[ix+1:]):
            total_score += score_seqs(seq1, seq2)
            
    total_score /= comb(len(sample_arr), 2)
    return total_score

location, start, end, uncert_file = sys.argv[1:5]
samples = sys.argv[5:]
alignment_score = calculate_alignment_score_from_accessions(samples)
with open(uncert_file, "w") as dest:
    dest.write(f"location,start,end,samples,alignment_score\n")
    dest.write(f"{location},{start},{end},{len(samples)},{alignment_score}\n")
