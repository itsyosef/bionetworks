import pandas as pd
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt
from glob import glob
import pickle
import sys

def get_date_pop_prev(data_path, var="1", prev = 0.1, cumulative = True):
    if cumulative:
        state_mod = "[in]"
    else:
        state_mod = "[current]"
    num_nodes = 7605430
    data = pd.read_csv(data_path)
    data["var_pop_prev"] = data[[c for c in data.columns if c.startswith(f"var{var}") and state_mod in c and "I" in c]].sum(axis=1) / num_nodes
    data["v1"] = data[[c for c in data.columns if c.startswith(f"var1") and "[current]" in c and "I" in c]].sum(axis=1)
    data["v2"] = data[[c for c in data.columns if c.startswith(f"var2") and "[current]" in c and "I" in c]].sum(axis=1)
    data["var2_prev"] = data["v2"] / (data["v2"] + data["v1"])
    if cumulative:
        data["var_pop_prev"] = data["var_pop_prev"].cumsum()
        
    cutoff_tick = data[data["var_pop_prev"] <= prev].index.max()
    v2_prev = data.iloc[cutoff_tick]["var2_prev"]
    return cutoff_tick, v2_prev

base_data = "/gpfs/gpfs0/project/bii_nssac/COVID-19_USA_EpiHiper/rivanna/20210827-var_surv"
save_base = "/project/biocomplexity/biosurveillance/EpiHiper/postProcessing/nodeProfiles"

exp = "exp7"
scenarios = [i.split("/")[-1] for i in glob(f"{base_data}/{exp}/*")]

scenario = scenarios[int(sys.argv[1])]
sc_tick = int(scenario.split("_")[7])
cum_prev_cutoff = [0.1, 0.01, 0.001]

replicate_data_v1 = {offset:[] for offset in cum_prev_cutoff}
# v1_rep_date_prevs = []
replicate_data_v2 = {offset:[] for offset in cum_prev_cutoff}
replicate_data_v2_re = {offset:[] for offset in cum_prev_cutoff}
rep_date_prevs = []

for i in range(60):
    data_path = f"{base_data}/{exp}/{scenario}/replicate_{i+1}/output.csv.gz"
    all_data = pd.read_csv(data_path)
    for cutoff in cum_prev_cutoff:
        v1_tick, v1_prev = get_date_pop_prev(f"{base_data}/{exp}/{scenario}/replicate_{i+1}/outputSummary.csv.gz", var="1", prev=cutoff)
        v2_tick, v2_prev = get_date_pop_prev(f"{base_data}/{exp}/{scenario}/replicate_{i+1}/outputSummary.csv.gz", var="2", prev=cutoff)
        
        data_v1 = all_data[all_data.tick <= v1_tick]
        data_v2 = all_data[(all_data.tick <= v2_tick) & (all_data.tick >= sc_tick)]
        
        nodes_v1 = data_v1[data_v1.exit_state.isin(['var1Iasymp', 'var1Isymp'])].pid.unique()    
        nodes_v2 = data_v2[data_v2.exit_state.isin(["var2Isymp", "var2Iasymp", 'var2var1Isymp', 'var2var1Iasymp'])].pid.unique()    
        nodes_v2_re = data_v2[data_v2.exit_state.isin(['var2var1Isymp', 'var2var1Iasymp'])].pid.unique()    
        
        replicate_data_v1[cutoff].append([f"rep_{i+1}"] + list(nodes_v1))
        replicate_data_v2[cutoff].append([f"rep_{i+1}"] + list(nodes_v2))
        replicate_data_v2_re[cutoff].append([f"rep_{i+1}"] + list(nodes_v2_re))
        
        rep_date_prevs.append((1, scenario, i+1, cutoff, v1_tick, v1_prev))
        rep_date_prevs.append((2, scenario, i+1, cutoff, v2_tick, v2_prev))
        
#     break

with open(f"{save_base}/{exp}/{scenario}_node_V1_binary_infection_profile.pkl", "wb") as dest:
    pickle.dump(replicate_data_v1, dest)
    
with open(f"{save_base}/{exp}/{scenario}_node_V2_binary_infection_profile.pkl", "wb") as dest:
    pickle.dump(replicate_data_v2, dest)
    
with open(f"{save_base}/{exp}/{scenario}_node_V2_from_V1_reinfection_profile.pkl", "wb") as dest:
    pickle.dump(replicate_data_v2_re, dest)
    
pd.DataFrame.from_records(rep_date_prevs, columns = ["variant", "scenario", "replicate", "cutoff", "cutoff_tick", "v2_prev_at_cutoff"]).to_csv(f"{save_base}/{exp}/{scenario}_node_cumulative_cutoff_prevalence_info_for_profiles.csv", index=False)