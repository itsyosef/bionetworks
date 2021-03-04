#!/usr/bin/env python3
# coding: utf-8

import networkx as nx
import json
from Bio import SeqIO
import sys
import pandas as pd
import random
import numpy as np
import subprocess
import requests

def get_splits(k_split, num_splits=5, top_k = 100):

    smash = pd.read_csv("../misc_data/antismash_db_results.csv", sep="\t")
    top_k_accessions = smash["NCBI accession"].value_counts()[:top_k].index
    # above_10_index = counts[counts > 5].index
    smash = smash[smash["NCBI accession"].isin(top_k_accessions)]
    smash["NCBI accession"] = smash["NCBI accession"].apply(lambda row: row.split(".")[0])
    smash_dict = dict()
    for accession, df in smash.groupby("NCBI accession"):
        smash_dict[accession] = df[["From", "To"]].to_numpy()
    smash_nums = sorted([(k, len(value)) for k, value in smash_dict.items()], key = lambda tup: tup[1], reverse=True)
    random.seed(42)
    random.shuffle(smash_nums)
    splits = [[] for i in range(num_splits)]
    for ix, tup in enumerate(smash_nums):
        split = ix % num_splits
        splits[split].append(tup[0])

    test = splits.pop(k_split)
    train = [i for s in splits for i in s]
    
    test = subprocess.run(["../scripts/genome_accession_to_patric_genome_ids.sh"] + test, capture_output=True).stdout.decode("utf-8").split("\n")
    test = [i for i in test if i != ""]
    
    train = subprocess.run(["../scripts/genome_accession_to_patric_genome_ids.sh"] + train, capture_output=True).stdout.decode("utf-8").split("\n")
    train = [i for i in train if i != ""]
    
    return train, test

def chunker(data, chunk_size):
    for i in range(0, len(data), chunk_size):
        yield data[i:i+chunk_size]

def tsv_from_genome_ids(outfile, genome_ids, limit=25000000, stream=False):
    with open(outfile, "w") as dest:
        pass
    #Changed to 1 because the output was getting clipped and changing the limit didn't seem to fix it
    for gids in chunker(genome_ids, 1):
        selectors = ["ne(feature_type)","eq(annotation,PATRIC)",f"in(genome_id,({','.join(gids)}))"]
        genomes = f"and({','.join(selectors)})"    
        limit = f"limit({limit})"
        select = "select(genome_id,genome_name,accession,annotation,feature_type,patric_id,refseq_locus_tag,alt_locus_tag,uniprotkb_accession,start,end,strand,na_length,gene,product,figfam_id,plfam_id,pgfam_id,go,ec,pathway,aa_sequence_md5)&sort(+genome_id,+sequence_id,+start)"
        base = "https://www.patricbrc.org/api/genome_feature/"
        query = "&".join([genomes, limit, select])
        headers = {"accept":"text/tsv", "content-type": "application/rqlquery+x-www-form-urlencoded"}

        genome_ids = []

        r = requests.post(url=base, data=query, headers=headers, stream=stream)

        if r.encoding is None:
            r.encoding = "utf-8"
        with open(outfile, "a") as dest:
            dest.write(r.text)
            
def gen_pfam_graph_input(source, destination):
    with open(source) as src:
        with open(destination, "w") as dest:
            for ix, line in enumerate(src.readlines()):
                vals = line.strip().split("\t")
                if ix == 0:
                    names = vals
                    name_ids = {name:ix for ix, name in enumerate(names)}
                try:
                    vals = {name:vals[ix] for ix, name in enumerate(names)}
                except:
                    continue
                md5 = vals["aa_sequence_md5"].replace('"', '')
                pgfam = vals["pgfam_id"]

                #This means that the pgfam is NAN and should be skipped
                if len(pgfam) < 2:
                    continue

                dest.write(line)

                try:
                    pfams = seq_dict[md5]
                except Exception as e:
                    continue
                if len(pfams) > 0:
                    for pfam in pfams:
                        pfam_line = line.split("\t")
                        pfam_line[name_ids["start"]] = str(int(vals['start']) + int(pfam['start']))
                        pfam_line[name_ids["end"]] = str(int(vals['start']) + int(pfam['end']))
                        pfam_line[name_ids["na_length"]] = str(int(pfam['end']) - int(pfam['start']))
                        pfam_line[name_ids["pgfam_id"]] = pfam["sequence_domain"]
                        pfam_line[name_ids["figfam_id"]] = pfam["sequence_domain"]
                        pfam_line[name_ids["plfam_id"]] = pfam["sequence_domain"]
                        pfam_line[name_ids["annotation"]] = "CUSTOM"
                        dest.write(("\t".join(pfam_line)))


outprefix = sys.argv[1]
splits = sys.argv[2:]
                        
records = [r for r in SeqIO.parse("/sfs/lustre/bahamut/scratch/jho5ze/bionets/BGCs/genera/Mycobacterium/jsonhfasta/2021-02-24.jsonhfasta", "fasta")]
seq_dict = [json.loads(r.description.split(None, 1)[-1]) for r in records]
seq_dict = {x["md5"]:sorted([i for i in x.setdefault("sequence_domain", [{"source":"none", "start":"0"}]) if i["source"] == "pfam"], key = lambda x: int(x["start"])) for x in seq_dict}
# seq_dict = {"md5":[{"sequence_domain":"pfam", "source":"none", "start":"0", "end":"0"}]}

for split in splits:
    split = int(split)
    train, test = get_splits(split)
    tsv_from_genome_ids(outprefix + f"_split_{split}_train.tsv", train)
    tsv_from_genome_ids(outprefix + f"_split_{split}_test.tsv", test)
    gen_pfam_graph_input(outprefix + f"_split_{split}_train.tsv", outprefix + f"_split_{split}_train_with_pfams.tsv")
    gen_pfam_graph_input(outprefix + f"_split_{split}_test.tsv", outprefix + f"_split_{split}_test_with_pfams.tsv")

