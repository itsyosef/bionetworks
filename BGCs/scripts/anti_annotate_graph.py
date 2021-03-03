#!/usr/env/bin python

import pandas as pd
import numpy as np
import networkx as nx
from Bio import SeqIO
import graphviz as gv
import json
import sys

def is_feature_in_bgc(contig_name, start, stop):
    start = int(start)
    stop = int(stop)
    try:
        bgcs = smash_dict[contig_name]
    except:
        return False
    for b_start, b_stop in bgcs:
        if (start>= b_start) and (stop <= b_stop):
            return True
    return False

def annotate_nodes(pfams):
    for node_id in pfams.nodes:
        node = pfams.nodes[node_id]
        node_bgc = False
        if type(node["features"]) == str:
            features = json.loads(node["features"].replace('""', '"'))
        else:
            features = node["features"]
            
        if node["family"][:2] == "PF":
            node["resolution"] = "Pfam"
        else:
            node["resolution"] = "PGFam"

        for genome in features["info"]:
            for contig in features["info"][genome]:
                for ix, sequence in enumerate(features["info"][genome][contig]):
                    is_in_bgc = is_feature_in_bgc(contig, sequence["start"], sequence["end"])
                    features["info"][genome][contig][ix]["in_bgc"] = is_in_bgc
                    if is_in_bgc:
                        node_bgc = True

        node["features"] = features
        node["BGC"] = node_bgc
        
def annotate_edges(pfams):
    for u,v,a in pfams.edges(data=True):
        start = pfams.nodes[u]
        stop = pfams.nodes[v]
        edge_bgc = start["BGC"] and stop["BGC"]
        is_start_pfam = start["family"][:2] == "PF"
        is_stop_pfam = stop["family"][:2] == "PF"
        if is_start_pfam and is_stop_pfam:
            edge_type = "pfam"
        elif is_start_pfam ^ is_stop_pfam:
            edge_type = "mixed"
        else:
            edge_type = "pgfam"
        a["BGC"] = edge_bgc
        a["LinkType"] = edge_type

def process_graph(graph_path, save=True):
    pfams = nx.readwrite.read_gexf(graph_path)    
    annotate_nodes(pfams)
    annotate_edges(pfams)
    if save:
        save_path = graph_path.split(".")[0] + "_processed.gexf"
        nx.readwrite.write_gexf(pfams, save_path)        
        
graph = sys.argv[1]
process_graph(graph)