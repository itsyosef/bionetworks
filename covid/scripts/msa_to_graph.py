import pandas as pd
import numpy as np
import networkx as nx
from math import log
from pathlib import Path
from matplotlib import pyplot as plt
from copy import deepcopy
import sys
import os

def add_root(ggraph, num_genomes):
    no_in_nodes = []
    for node in ggraph.nodes():
        if ggraph.in_degree(node) == 0:
            no_in_nodes.append(node)
    if len(no_in_nodes) == 1:
        return no_in_nodes[0]
    else:
        ggraph.add_node("root", sequence="N", num_genomes=num_genomes)
        for node in no_in_nodes:
            ggraph.add_edge("root", node)
        return "root", no_in_nodes

def get_complement_nodes(ggraph, x_node): #, num_genomes):
    #Could benefit from some indexing
    #Could optimize by finding the consensus nodes and stopping 
    #     early there since all above or below are "implied"
    #     but would need to potentially iterate over the nodes 
    #     again since switching from negation to inclusion
    node_states = {node:True for node in ggraph.nodes}
        
    def follow_ancestors(node):
        if node_states[node] or node == x_node: #stop if node already visited
            node_states[node] = False
            for anc_node in ggraph.predecessors(node):
                follow_ancestors(anc_node)
            
    def follow_descendants(node):
        if node_states[node] or node == x_node: #stop if node already visited
            node_states[node] = False
            for dec_node in ggraph.successors(node):
                follow_descendants(dec_node)
    
    follow_ancestors(x_node)
    follow_descendants(x_node)
    
    return [node for node in node_states if node_states[node]]

def calc_node_chars(ggraph, node):
    return len(ggraph.nodes[node]["sequence"]) * ggraph.nodes[node]["num_genomes"]

def calc_uncertainty(ggraph, num_chars=None): #, num_genomes):
    
    char_per_node = False
    
    if num_chars is None: #Use each node as a character
        num_chars = len(graph.nodes())
        char_per_node = True
    
    uncert = 0
    
    for x_node in ggraph.nodes():
        x_node_info = ggraph.nodes[x_node]
        x_node_characters = 1 if char_per_node else calc_node_chars(ggraph, x_node)
        
        x_uncert = x_node_characters / num_chars
        
        #For the first term
        uncert -= x_uncert * log( x_uncert , 2)
        
        #For the second term
        y_nodes = get_complement_nodes(ggraph, x_node)
        y_set_chars = len(y_nodes) if char_per_node else sum([calc_node_chars(ggraph, y_node) for y_node in y_nodes])
        for y_node in y_nodes:
            
            y_node_characters = 1 if char_per_node else calc_node_chars(ggraph, y_node)
            y_uncert = ( y_node_characters / y_set_chars )
            
            uncert -= x_uncert * y_uncert * log( y_uncert , 2)
            
    return uncert

def generate_nx_from_msa(seq_paths):
    num_chars = 0
    genome_len = 34742

    graph = nx.DiGraph()
    for ix in range(genome_len):
        for letter in ["A","T","C","G","N"]:
            graph.add_node(f"{ix}_{letter}", genomes=[], sequence=letter, num_genomes=0)

    for seq in seq_paths:
        with open(seq) as src:
            for line in src.readlines():
                line = line.strip()
                if line[0] == ">":
                    genome = line[1:] #.split("/")[2]
                else:

                    prev_node = None
                    num_chars += len(line)
                    for ix, char in enumerate(line): #[34000:35000]):
                        if char == "-":
                            char = "N" #Preserves overall length of genome graph
                        if char not in ["A","T","C","G","N"]:
                            char = "N"
                            continue
                        node = f"{ix}_{char}"
                        graph.nodes[node]["genomes"].append(genome)
                        graph.nodes[node]["num_genomes"] += 1
                        if prev_node is not None and node not in graph.successors(prev_node):
                            graph.add_edge(prev_node, node)
                        prev_node = node
    empty_nodes = []
    for node in graph.nodes():
        if not graph.nodes[node]["genomes"]:
            empty_nodes.append(node)
            
    graph.remove_nodes_from(empty_nodes)
    
    return graph, num_chars

def collapse_graph(graph):
    new_graph = deepcopy(graph)
    redundant_nodes = []
    node_states = {node:False for node in new_graph.nodes} #tracks if nodes will be removed
    
    for main_node in new_graph:
        seq = new_graph.nodes[main_node]["sequence"]
        cur_node = main_node
        successors = list(new_graph.successors(cur_node))
        if len(successors) == 0 or len(successors) > 1 or node_states[successors[0]]:
            continue
        elif len(list(new_graph.predecessors(successors[0]))) > 1:
            continue
        while len(successors) <= 1:
            cur_node = successors[0]
            redundant_nodes.append(cur_node)
            node_states[cur_node] = True
            seq += new_graph.nodes[cur_node]["sequence"]
            successors = list(new_graph.successors(cur_node))
            if len(successors) == 0 or node_states[successors[0]]:
                break
            elif len(list(new_graph.predecessors(successors[0]))) > 1:
                break

        for node in new_graph.successors(cur_node):
            new_graph.add_edge(main_node, node)
                
        new_graph.nodes[main_node]["sequence"] = seq
    new_graph.remove_nodes_from(redundant_nodes)
        
    return new_graph 

graph_out_path = Path(sys.argv[1])
uncert_file = sys.argv[2]
seq_paths = sys.argv[3:]

num_genomes = len(seq_paths)
graph, num_chars = generate_nx_from_msa(seq_paths)
add_root(graph, num_genomes)
collapsed_graph = collapse_graph(graph)

full_uncert = calc_uncertainty(collapsed_graph, num_chars=num_chars)
collapsed_uncert = calc_uncertainty(collapsed_graph)

graph_name = str(graph_out_path).split("/")[-1]
graph_path = graph_out_path / (graph_name + ".gpkl")
nx.write_gpickle(collapsed_graph, graph_path)
if not os.path.exists(uncert_file):
    with open(uncert_file, "w") as dest:
        dest.write("graph_name,full_uncertainty,collapsed_uncertainty\n")
        
with open(uncert_file, "a") as dest:
    dest.write(f"{graph_name},{full_uncert},{collapsed_uncert}\n")