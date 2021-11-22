#!/usr/bin/env python3

import pandas as pd
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt
import networkx as nx
from glob import glob
from pathlib import Path
import json
from tqdm import tqdm
import sys
import os
from math import ceil

susc = False

with open("../data/edge_list_dict.json") as src:
    node_degrees = json.load(src)
    
base = Path("/scratch/jho5ze/bionets/epihiper/data/output/")
experiment = "exp5"
if not os.path.exists(f"../data/processed_{experiment}"):
    os.mkdir(f"../data/processed_{experiment}")
scenarios = sorted(glob(str(base / experiment / "*")))
replicates = sorted([i for i in glob(str(base / experiment / scenarios[0] / "*")) if "replicate" in i], key = lambda x: int(x.split("_")[-1]))

ix = int(sys.argv[1])
scenario = scenarios[ix // 40].split("/")[-1]
replicate = replicates[20 + (ix % 40)].split("/")[-1]

output = pd.read_csv(base / experiment / scenario / replicate / "output.csv.gz")
summary = pd.read_csv(base / experiment / scenario / replicate / "outputSummary.csv.gz")

if summary[[f"var2Isymp[current]", f"var2Iasymp[current]"]].iloc[-1].sum() < 1:
    sys.exit(0)

node_ids = [int(i) for i in node_degrees.keys()]
state_array = ["S" for i in range(7688058+1)]
state_dict = {state:[-1 for i in range(7688058+1)] for state in list(output.exit_state.unique()) + ["S"]}
for pid in node_ids:
    state_dict["S"][pid] = pid
dfs = []
susc_dfs = []
for tick in tqdm(range(400)):
    df = output[output.tick == tick]
    for pid, state in df.set_index("pid")["exit_state"].to_dict().items():
        prev_state = state_array[pid]
        state_array[pid] = state
        state_dict[prev_state][pid] = -1 
        state_dict[state][pid] = pid 
    if susc:
        bins = {i:0 for i in range(ceil(max(node_degrees.values()) / 50))}
        for pid in list(set(state_dict["S"]) - {-1}):
            bins[node_degrees[str(pid)]//50] += 1
        x = pd.DataFrame(bins, index=[tick])
        susc_dfs.append(x)
        if tick % 10 == 0:
            susc_data = pd.concat(susc_dfs)
            susc_data.to_csv(f"../data/processed_{experiment}/{scenario}_{replicate}_susc_data.csv", index=False)
    else:
        x = pd.DataFrame(list(set(state_dict["var1E"]) - {-1}), columns=["pid"])
        x["tick"] = tick
        x["variant"] = "var1"
        dfs.append(x)
        x = pd.DataFrame(list(set(state_dict["var2E"]) - {-1}), columns=["pid"])
        x["tick"] = tick
        x["variant"] = "var2"
        dfs.append(x)
        if tick % 10 == 0:
            variant_data = pd.concat(dfs)
            #variant_data["node_degree"] = variant_data.pid.apply(lambda row: node_degrees[str(row)])

            variant_data.to_csv(f"../data/processed_{experiment}/{scenario}_{replicate}_variant_data.csv", index=False)
#
    
#     x = pd.DataFrame(suscep_degrees.values(), columns=["node_degree"])
#     x["tick"] = tick
#     susc_dfs.append(x)
if susc:
    susc_data = pd.concat(susc_dfs)
    susc_data.to_csv(f"../data/processed_{experiment}/{scenario}_{replicate}_susc_data.csv", index=False)
else:
    variant_data = pd.concat(dfs)
    variant_data["node_degree"] = variant_data.pid.apply(lambda row: node_degrees[str(row)])

    # suscep_data = pd.concat(susc_dfs)

    variant_data.to_csv(f"../data/processed_{experiment}/{scenario}_{replicate}_variant_data.csv", index=False)
# suscep_data.to_csv(f"../data/processed_scenarios/{scenario}_{replicate}_suscep_data.csv", index=False)
