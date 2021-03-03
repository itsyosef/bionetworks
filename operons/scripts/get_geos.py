import sqlite3
import pandas as pd
import numpy as np
import wget
import tarfile
from pathlib import Path
import sys, os

global con, cursor, gse_map

con = sqlite3.connect("GEOmetadb.sqlite")
cursor = con.cursor()
cols = "title,gse,summary,overall_design,supplementary_file"
gse_map = {i:ix for ix, i in enumerate(cols.split(","))}

def get_genus_ftps(genus):

    query = f"""SELECT * 
            FROM gse_ft WHERE gse_ft
            MATCH '{genus}'"""

    query = f"""SELECT *
            FROM gse_ft WHERE gse_ft
            MATCH 'RNA AND seq'
            AND gse_ft
            MATCH '{genus}'
            """
    results = cursor.execute(query).fetchall()

    positive_terms = "wild WT dRNA-".split()
    negative_terms = "plasmid".split()

    new_results = []
    for i in results:
        descriptive_text = i[gse_map["title"]] + i[gse_map["summary"]] + i[gse_map["overall_design"]]
        add = True
        for n in negative_terms:
            if n in descriptive_text:
                add = False
                break
        if add:
            new_results.append(i)

    results = new_results

    ftps = []

    for ix, i in enumerate(results):
        descriptive_text = i[gse_map["title"]] + i[gse_map["summary"]] + i[gse_map["overall_design"]]
        for p in positive_terms:
            if p in descriptive_text:
                if i[-1] is None:
                    continue
                if "RAW.tar" in i[-1]:
                    ftps.append([j for j in i[-1].split(";") if "RAW" in j][0].strip())
                break
    return ftps

def get_wigs_from_ftps(ftps, path="."):
    path = Path(path)
    for file in ftps:
        extracted = False
        print("Downloading: ", file)
        wget.download(file)
        file = file.split("/")[-1]
        with tarfile.open(file, "r:") as tar:
            names = tar.getmembers()
            wig_files = [n for n in names if ".wig" in n.name] # and ("wt" in n.name or "WT" in n.name) ]
            if len(wig_files) > 0:
                if not os.path.exists(path):
                    os.mkdir(path)

                tar.extractall(path=path, members=wig_files)
    
        os.remove(file)
        
genus = sys.argv[1].split("/")[-1].split("_")[0]
path = "/scratch/jho5ze/bionets/operons/geo_output/" + genus

ftps = get_genus_ftps(genus)
get_wigs_from_ftps(ftps, path=path)
