import networkx as nx
import sys
import json
from networkx.readwrite import json_graph

in_graph = sys.argv[1]
out_graph = sys.argv[2]

g = nx.readwrite.read_gexf(in_graph)

json.dump(json_graph.node_link_data(g), open(out_graph, "w"))
