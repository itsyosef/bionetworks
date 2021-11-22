import pandas as pd
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt
import networkx as nx
from glob import glob
from pathlib import Path
import json
from tqdm import tqdm
import sys, os

with open("../data/edge_list_dict.json") as src:
    node_degrees = json.load(src)
    
    node_traits = pd.read_table("/project/biocomplexity/nssac/EpiHiperSynPop/v1.9.0/usa_va_2017_SynPop/va_persontrait_epihiper.txt", low_memory=True, sep=",", usecols=[0,2,3,5], header=1)
node_age = node_traits.set_index("pid")["age"].to_dict()

base = Path("/scratch/jho5ze/bionets/epihiper/data/output/")
experiment = "exp5"

if not os.path.exists(f"../plots/{experiment}"):
    os.mkdir(f"../plots/{experiment}")

scenarios = sorted(glob(str(base / experiment / "*")))
replicates = sorted([i for i in glob(str(base / experiment / scenarios[0] / "*")) if "replicate" in i], key = lambda x: int(x.split("_")[-1]))

ix = int(sys.argv[1])
scenario = scenarios[ix].split("/")[-1]

def scenario_to_title(scenario):
    scenario = scenario.split("/")[-1]
    scenario = scenario.replace("tau_", "Variant 1: ")
    scenario = scenario.replace("_var2_1.4_seed2p", "")
    scenario = scenario.replace("_5", ", Location: 5")
    scenario = scenario.replace("_seed2t_", ", Emergence Day: ")
#     scenario = scenario.replace("", "")
    return scenario

def plot_scenario(scenario, replicate_data, box=False, tick_spread=7, data_col = "node_degree", max_tick=280):
    fig, axs = plt.subplots(figsize=(120,50), ncols=5, nrows=4)
    box_str = "_box" if box else "_violin"
    for ix in range(20):
        y = ix // 4
        x = ix % 4

        ax = axs[x, y]
        data = replicate_data[ix]
        data = data[data.tick <= max_tick]
        data["plot_tick"] = data.tick.apply(lambda row: tick_spread * (row // tick_spread))
        split = True if len(data["variant"].unique()) == 2 else False
        if box:
            sns.boxplot(x="plot_tick", y=data_col, hue="variant", data=data, ax=ax)
        else:
            sns.violinplot(x="plot_tick", y=data_col, hue="variant", data=data, split=split, scale="count", inner="quartile", bw=0.2, ax=ax)
        ax.set_title(f"Replicate {ix+1}", fontsize=36)
        ax.legend().remove()
        ax.tick_params(axis="x", labelsize=20, labelrotation=45)
    plt.suptitle(f"{data_col} for Scenario: {scenario_to_title(scenario)}", y = 0.92, fontsize=58)
    plt.savefig(f"../plots/{experiment}/{data_col}{box_str}_{scenario}.png", dpi=100, facecolor="white", bbox_inches="tight")


    
    # for scenario in ["tau_0.02_var2_1.4_seed2p_51059_seed2t_120"]: #scenarios[2:]:
# for scenario in scenarios:
#     try:
#         scenario = scenario.split("/")[-1]
    # scenario = "tau_0.02_var2_1.4_seed2p_51059_seed2t_120"
replicate_data = []

print(scenario)

for replicate in replicates:
    variant_data = pd.read_csv(f"../data/processed_{experiment}/{scenario.split('/')[-1]}_{replicate.split('/')[-1]}_variant_data.csv")
    variant_data["node_degree"] = variant_data.pid.apply(lambda row: node_degrees[str(row)])
    variant_data["node_age"] = variant_data.pid.apply(lambda row: node_age[row])
#     variant_data["plot_tick"] = variant_data.tick.apply(lambda row: day_agg * (row // day_agg))
    replicate_data.append(variant_data)
plot_scenario(scenario, replicate_data)
plot_scenario(scenario, replicate_data, data_col="node_age")
plot_scenario(scenario, replicate_data, box=True)
plot_scenario(scenario, replicate_data, box=True, data_col="node_age")
# except:
# print(scenario)
# continue
#     break