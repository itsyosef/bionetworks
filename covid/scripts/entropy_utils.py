import pandas as pd
import numpy as np
import networkx as nx
from math import log
from pathlib import Path
from matplotlib import pyplot as plt
from copy import deepcopy
from tqdm import tqdm
from screed import ScreedDB
import re
from Bio.Seq import Seq

msadb = ScreedDB("../data/msa_0408/usa_msa_0408.fasta")

def expand_ggraph(ggraph):
    """
    Expand the graph from normal run of GrnGraph so each node is a nucleotide
    """
    new_graph = nx.DiGraph()
    entry_points = {node:None for node in ggraph.nodes()}
    
    for node in ggraph.nodes():
        
        node_info = ggraph.nodes[node]
        
        if entry_points[node] is None:
            last_node = None
            entry_points[node] = node + "_0"
            start_idx = 0
        else:
            last_node = entry_points[node]
            start_idx = 1
            
        ids = set(ggraph.nodes[node]["ids"].split(","))
        num_genomes = len(ids)
        for ix, character in enumerate(node_info["sequence"][start_idx:]):
            ix += start_idx
            new_node = node + "_" + str(ix)
#             new_graph.add_node(new_node, attr_dict={"sequence":character, "num_genomes":num_genomes})
            new_graph.add_node(new_node, sequence = character, num_genomes = num_genomes)
            if last_node is not None:
                new_graph.add_edge(last_node, new_node, weight=num_genomes)
            last_node = new_node

        for child in ggraph.successors(node):
            child_info = ggraph.nodes[child]
            child_ids = set(child_info["ids"].split(","))
            
            entry_point = entry_points[child]
            if entry_point is not None:
                new_graph.add_edge(last_node, entry_point, weight=len(ids & child_ids))
            else:

                child_node = child + "_0"
#                 new_graph.add_node(child_node, attr_dict={"sequence":child_info["sequence"][0], "num_genomes":len(child_info["ids"].split(","))})
                new_graph.add_node(child_node, sequence=child_info["sequence"][0], num_genomes=len(child_info["ids"].split(",")))
                entry_points[child] = child_node
                new_graph.add_edge(last_node, child_node, weight=len(ids & child_ids))
                
    return new_graph

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
            new_graph.add_edge(main_node, node, weight=new_graph[cur_node][node]["weight"])
                
        new_graph.nodes[main_node]["sequence"] = seq
    new_graph.remove_nodes_from(redundant_nodes)
        
    return new_graph     

def annotate_graph(graph):
    for node in graph.nodes():
        node = graph.nodes[node]
        node["sequence_len"] = len(node["sequence"])
    
def add_root(ggraph, genomes):
    num_genomes = len(genomes)
    no_in_nodes = []
    for node in ggraph.nodes():
        if ggraph.in_degree(node) == 0:
            no_in_nodes.append(node)
    if len(no_in_nodes) == 1:
        return no_in_nodes[0]
    else:
        ggraph.add_node("root", sequence="A", num_genomes=num_genomes, genomes=genomes)
        for node in no_in_nodes:
            ggraph.add_edge("root", node, weight=ggraph.nodes[node]["num_genomes"])
        return "root", no_in_nodes

def get_sequence_genome_weight(ggraph, node):
    return ggraph.nodes[node]["sequence_len"] * ggraph.nodes[node]["num_genomes"]

def get_genome_weight(ggraph, node):
    return ggraph.nodes[node]["num_genomes"]

def get_sequence_weight(ggraph, node):
    return len(ggraph.nodes[node]["sequence"])

def get_weight(ggraph, node):
    return 1

def calc_uncertainty(ggraph, graph_presence_info, weigh_on_sequence=True, weigh_on_genomes=True):
    
    """
    Calculates the uncertainty/entropy of the genome graph using configurable weighing schemes
    
    ToDo:
        Determine if the genome weighing by itself should be, e.g. for the graph, the sum of the number of genomes in each node
    """
    
    sorted_genome_dict, sorted_nodes, genome_presences = graph_presence_info
    
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
    
#     for x_node in ggraph.nodes():
    for x_node in tqdm(ggraph.nodes()):
        
        #uninduced nodes
        genome_ixs = [sorted_genome_dict[g] for g in ggraph.nodes[x_node]["genomes"]]
        genome_x_presences = genome_presences[genome_ixs]
        mask = np.any(genome_x_presences.T, axis=1)
        y_nodes = sorted_nodes[~mask]
        
        x_weight = func(ggraph, x_node)
        
        x_uncert = x_weight / graph_weight
        
        #For the first term
        uncert -= x_uncert * log( x_uncert , 2)
        
        #For the second term
        y_set_weight = sum([func(ggraph, y_node) for y_node in y_nodes])
        for y_node in y_nodes:
            
            y_node_weight = func(ggraph, y_node)
            y_uncert = ( y_node_weight / y_set_weight )
            
            uncert -= x_uncert * y_uncert * log( y_uncert , 2)
    return uncert

def msa_from_seq_paths(seq_paths): 
    for seq in seq_paths:
        with open(seq) as src:
            for line in src.readlines():
                yield line
                
def msa_from_screed_ids(ids): 
    for seq_id in ids:
        if seq_id in msadb:
            record = msadb[seq_id]
            name = ">"+record["name"]
            sequence = str(record["sequence"])
            yield name
            yield sequence

def msa_from_string(msas):
    for line in msas.split("\n"):
        if line.strip() != "":
            yield line
            
def match_line_to_reference(line):
    with open("/scratch/jho5ze/bionets/covid/data/msa_0408/reference_msa.txt") as src:
        reference_line = src.readlines()[0]
    processed_line = ""
    for ix, char in enumerate(line):
        if reference_line[ix] != "-":
            processed_line += char
    return processed_line

def translate_line(line):
    get_dash = re.compile("-+")
    dashes = get_dash.findall(line)
    if dashes and min([len(i)%3 == 0 for i in dashes]) == False: #Some dash stretch is not divisible by 3
        return ""
    else:
        dashes.append("") #so we don't run out of next "-" stretches in the loop below
        aa_line = ""
        if line[0] == "-":
            aa_line += "-" * (len(dashes[0])//3)
            dashes.pop(0)
        for ix, nas in enumerate([i for i in get_dash.split(line) if i != ""]):
            aa_line += str(Seq(nas).translate())
            aa_line += "-" * (len(dashes[ix])//3)
        return aa_line
    
def translate_record(record):
    orfs = [(29558, 29674), \
            (28274, 29533), \
            (27894, 28259), \
            (27756, 27887), \
            (27394, 27759), \
            (27202, 27387), \
            (26523, 27191), \
            (26245, 26472), \
            (25393, 26220), \
            (21563, 25384), \
            (13468, 21555), \
            (266, 13483)]
    aa_seqs = []
    if type(record) == str:
        record_seq = match_line_to_reference(record)
    else:
        record_seq = match_line_to_reference(str(record["sequence"]))
    for start, stop in orfs:
        start -= 1 #1 based indexing from the preprocessing of orfs (from the GFF file)
        segment = record_seq[start:stop]
        aa_seqs.append(translate_line(segment))
    return "".join(aa_seqs)
        
def generate_nx_from_msa(msa, accepted_letters=None, screed_ids=True, aa_level=False):
    genome_len = 34742 #54
    genomes = []
    if accepted_letters is None:
        if aa_level:
            accepted_letters = [char for char in "ARNDCQEGHILKMFPSTWYVUOXBZJ*"]
        else:
            accepted_letters = ["A","T","C","G"]
    
    #Initialize the graph
    graph = nx.DiGraph()
    for ix in range(genome_len):
        for letter in accepted_letters: 
            graph.add_node(f"{ix}_{letter}", genomes=[], sequence=letter, num_genomes=0)
            
    if type(msa) == str:
        msa_yield = msa_from_string
    else: #Type is iterable
        if screed_ids:
            msa_yield = msa_from_screed_ids
        else:
            msa_yield = msa_from_seq_paths

    edges = dict()
    for line in msa_yield(msa):
        line = line.strip()
        if line[0] == ">":
            genome = line[1:]
            if len(genome.split("/")) >= 2:
                genome = genome.split("/")[2]
            genomes.append(genome)
        else:
            if aa_level:
                line = translate_record(line)
            prev_node = None
            for ix, char in enumerate(line): 
                if char not in accepted_letters: #Treats all other letters as a dash (and a dash like a dash)
                    continue
                node = f"{ix}_{char}"
                graph.nodes[node]["genomes"].append(genome)
                graph.nodes[node]["num_genomes"] += 1
                if prev_node is not None: 
                    edges.setdefault((prev_node, node), [0])[0] += 1
                prev_node = node
            
    edges = [(key[0], key[1], value[0]) for key, value in edges.items()]
    graph.add_weighted_edges_from(edges)
            
    empty_nodes = []
    for node in graph.nodes():
        if not graph.nodes[node]["genomes"]:
            empty_nodes.append(node)
            
    graph.remove_nodes_from(empty_nodes)

    return graph, genomes

def get_genome_node_presence(ggraph, genomes):
    sorted_genomes = np.array(sorted(genomes))
    sorted_nodes = np.array(sorted(ggraph.nodes()))
    
    gs = [[] for g in sorted_genomes]
    for node in sorted_nodes:
        node_genomes = ggraph.nodes[node]["genomes"]
        for gix, g in enumerate(sorted_genomes):
            presence = True if g in node_genomes else False
            gs[gix].append(presence)
    
    sorted_genome_dict = {genome:ix for ix, genome in enumerate(sorted_genomes)}
    
    return (sorted_genome_dict, sorted_nodes, np.array(gs))