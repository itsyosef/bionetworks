import pandas as pd
import numpy as np
import networkx as nx
from math import log
from pathlib import Path
from matplotlib import pyplot as plt
from copy import deepcopy
import sys
import os

def add_root(ggraph, num_genomes, genomes):
    no_in_nodes = []
    for node in ggraph.nodes():
        if ggraph.in_degree(node) == 0:
            no_in_nodes.append(node)
    if len(no_in_nodes) == 1:
        return no_in_nodes[0]
    else:
        ggraph.add_node("root", sequence="N", num_genomes=num_genomes, genomes=genomes)
        for node in no_in_nodes:
            ggraph.add_edge("root", node)
        return "root", no_in_nodes

def get_complement_nodes(ggraph, x_node, induced_nodes_dict, num_genomes):
    x_genomes = ggraph.nodes[x_node]["genomes"]
    
    if len(x_genomes) == num_genomes: #All genomes included so all nodes induced
        return []
    
    induced_nodes = set().union(*[induced_nodes_dict[genome] for genome in x_genomes])
    all_nodes = set(ggraph.nodes())
    
    complement_nodes = all_nodes - induced_nodes
    return complement_nodes
    
# def get_complement_nodes(ggraph, x_node):
#     y_nodes = []
    
#     x_genomes = set(ggraph.nodes[x_node]["genomes"])
    
#     for node in ggraph.nodes():
#         node_genomes = set(ggraph.nodes[node]["genomes"])
#         if len(node_genomes & x_genomes) == 0:
#             y_nodes.append(node)
            
#     return y_nodes
    
# def get_complement_nodes(ggraph, x_node): #, num_genomes):
#     #Could benefit from some indexing
#     #Could optimize by finding the consensus nodes and stopping 
#     #     early there since all above or below are "implied"
#     #     but would need to potentially iterate over the nodes 
#     #     again since switching from negation to inclusion
#     node_in_y_states = {node:True for node in ggraph.nodes}
#     node_visited_states = {node:False for node in ggraph.nodes}
#     genomes = set(ggraph.nodes[x_node]["genomes"])
        
#     def follow_ancestors(node):
#         if not node_visited_states[node] or node == x_node: #stop if node already visited
#             node_visited_states[node] = True
#             if len(set(ggraph.nodes[node]["genomes"]) & genomes) > 0:
#                 node_in_y_states[node] = False
#                 for anc_node in ggraph.predecessors(node):
#                     follow_ancestors(anc_node)
            
#     def follow_descendants(node):
#         if not node_visited_states[node] or node == x_node: #stop if node already visited
#             node_visited_states[node] = True
#             if len(set(ggraph.nodes[node]["genomes"]) & genomes) > 0:
#                 node_in_y_states[node] = False
#                 for dec_node in ggraph.successors(node):
#                     follow_descendants(dec_node)
    
#     follow_ancestors(x_node)
#     follow_descendants(x_node)
    
#     return [node for node in node_in_y_states if node_in_y_states[node]]

def get_sequence_genome_weight(ggraph, node):
    return len(ggraph.nodes[node]["sequence"]) * ggraph.nodes[node]["num_genomes"]

def get_genome_weight(ggraph, node):
    return ggraph.nodes[node]["num_genomes"]

def get_sequence_weight(ggraph, node):
    return len(ggraph.nodes[node]["sequence"])

def get_weight(ggraph, node):
    return 1

def calc_uncertainty(ggraph, induced_nodes_dict, weigh_on_sequence=True, weigh_on_genomes=True): #, num_genomes):
    
    """
    Calculates the uncertainty/entropy of the genome graph using configurable weighing schemes
    
    ToDo:
        Determine if the genome weighing by itself should be, e.g. for the graph, the sum of the number of genomes in each node
    """
    
    num_genomes = len(induced_nodes_dict.keys())
    
    if weigh_on_sequence and num_chars is None:
        raise ValueError("If using sequence weighing, must provide the total number of characters in the graph")
    
    if weigh_on_genomes and weigh_on_sequence:
        func = get_sequence_genome_weight
    elif weigh_on_genomes:
        func = get_genome_weight
    elif weigh_on_sequence:
        func = get_sequence_weight
    else:
        func = get_weight
        
    graph_weight = sum([func(ggraph, node) for node in ggraph.nodes])
    uncert = 0
    
    terms = dict()
    
    for x_node in ggraph.nodes():
        
        x_weight = func(ggraph, x_node)
        
        x_uncert = x_weight / graph_weight
        
        #For the first term
        uncert -= x_uncert * log( x_uncert , 2)
        
        #For the second term
#         y_nodes = get_complement_nodes(ggraph, x_node)
        y_nodes = get_complement_nodes(ggraph, x_node, induced_nodes_dict, num_genomes)
        y_set_weight = sum([func(ggraph, y_node) for y_node in y_nodes])
        for y_node in y_nodes:
            
            y_node_weight = func(ggraph, y_node)
            y_uncert = ( y_node_weight / y_set_weight )
            
            uncert -= x_uncert * y_uncert * log( y_uncert , 2)
    return uncert

def generate_nx_from_msa(seq_paths):
    num_chars = 0
    genome_len = 34742
    genomes =[]
    
    graph = nx.DiGraph()
    for ix in range(genome_len):
        for letter in ["A","T","C","G"]: #,"N"]:
            graph.add_node(f"{ix}_{letter}", genomes=[], sequence=letter, num_genomes=0)

    for seq in seq_paths:
        with open(seq) as src:
            for line in src.readlines():
                line = line.strip()
                if line[0] == ">":
                    genome = line[1:].split("/")[2]
                    genomes.append(genome)
                else:

                    prev_node = None
                    num_chars += len(line)
                    for ix, char in enumerate(line): 
#                         if char == "-":
#                             char = "N" #Preserves overall length of genome graph
                        if char not in ["A","T","C","G"]: #,"N"]:
#                             char = "N"
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
    
    return graph, num_chars, genomes

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
graph, num_chars, genomes = generate_nx_from_msa(seq_paths)
add_root(graph, num_genomes, genomes)
collapsed_graph = collapse_graph(graph)

columns = []
values = []

graph_genome_induced_nodes = {genome:set() for genome in genomes}
collapsed_graph_genome_induced_nodes = {genome:set() for genome in genomes}
for node in graph.nodes():
    node_genomes = graph.nodes[node]["genomes"]
    for genome in genomes:
        if genome in node_genomes:
            graph_genome_induced_nodes[genome].add(node)
            
for node in collapsed_graph.nodes():
    node_genomes = collapsed_graph.nodes[node]["genomes"]
    for genome in genomes:
        if genome in node_genomes:
            collapsed_graph_genome_induced_nodes[genome].add(node) 

# for G, graph_name in zip([graph, collapsed_graph], ["full", "collapsed"]):
for G, graph_name, induced_nodes in zip([collapsed_graph], ["collapsed"], [collapsed_graph_genome_induced_nodes]):
    for seqs in range(2):
        seq_name = f"_sequences_{seqs}"
        for genomes in range(2):
            genome_name = f"_genomes_{genomes}"
            columns.append(graph_name + seq_name + genome_name)
            values.append(str(calc_uncertainty(G, induced_nodes, weigh_on_sequence=seqs, weigh_on_genomes=genomes)))

graph_name = str(graph_out_path).split("/")[-1]
graph_path = graph_out_path / (graph_name + ".gpkl")
nx.write_gpickle(collapsed_graph, graph_path)
if not os.path.exists(uncert_file):
    with open(uncert_file, "w") as dest:
        dest.write(f"graph_name,{','.join(columns)}\n")
        
with open(uncert_file, "a") as dest:
    dest.write(f"{graph_name},{','.join(values)}\n")