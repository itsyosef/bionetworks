import pandas as pd
import numpy as np
import networkx as nx
from math import log
from pathlib import Path
from matplotlib import pyplot as plt
from copy import deepcopy
import sys
import os
import entropy_utils as eu

graph_out_path = Path(sys.argv[1])
uncert_file = sys.argv[2]
seq_paths = sys.argv[3:]

is_houston = "Houston" in uncert_file

num_genomes = len(seq_paths)

if not is_houston:
    graph, genomes = eu.generate_nx_from_msa(seq_paths)
    eu.add_root(graph, genomes)
    collapsed_graph = eu.collapse_graph(graph)

aa_graph, aa_genomes = eu.generate_nx_from_msa(seq_paths, aa_level=True)
eu.add_root(aa_graph, aa_genomes)
aa_collapsed_graph = eu.collapse_graph(aa_graph)

if not is_houston:
    eu.annotate_graph(graph)
    eu.annotate_graph(collapsed_graph)
eu.annotate_graph(aa_graph)
eu.annotate_graph(aa_collapsed_graph)

if not is_houston:
    graph_presence_info = eu.get_genome_node_presence(graph, genomes)
    collapsed_graph_presence_info = eu.get_genome_node_presence(collapsed_graph, genomes)
aa_graph_presence_info = eu.get_genome_node_presence(aa_graph, aa_genomes)
aa_collapsed_graph_presence_info = eu.get_genome_node_presence(aa_collapsed_graph, aa_genomes)

columns = []
values = [] 

if not is_houston:
    iter_info = zip([graph, collapsed_graph, aa_graph, aa_collapsed_graph], ["full", "collapsed", "aa_full", "aa_collapsed"], [graph_presence_info, collapsed_graph_presence_info, aa_graph_presence_info, aa_collapsed_graph_presence_info])
else:
    iter_info = zip([aa_graph, aa_collapsed_graph], ["aa_full", "aa_collapsed"], [aa_graph_presence_info, aa_collapsed_graph_presence_info])
    

for G, graph_name, presence_info in iter_info:
    for seqs in range(2):
        seq_name = f"_sequences_{seqs}"
        for genomes in range(2):
            if "full" in graph_name:
                if not (seqs ^ genomes):#full_both and full_none reduce to full_genomes and full_seq, respectively
                    continue
            genome_name = f"_genomes_{genomes}"
            columns.append(graph_name + seq_name + genome_name)
            uncert = str(eu.calc_uncertainty(G, presence_info, weigh_on_sequence=seqs, weigh_on_genomes=genomes))
            values.append(uncert)
            print()
            print("Done with: ", graph_name + seq_name + genome_name, ": ", uncert)
            print()
            
graph_name = str(graph_out_path).split("/")[-1]
graph_path = graph_out_path / (graph_name + ".gpkl")
if not is_houston:
    nx.write_gpickle(collapsed_graph, graph_path)
    
if not os.path.exists(uncert_file):
    with open(uncert_file, "w") as dest:
        dest.write(f"graph_name,{','.join(columns)}\n")
        
with open(uncert_file, "a") as dest:
    dest.write(f"{graph_name},{','.join(values)}\n")