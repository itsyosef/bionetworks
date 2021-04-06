import numpy as np 
import pandas as pd
import networkx as nx
import math

from networkx.drawing.nx_agraph import graphviz_layout
import matplotlib.pyplot as plt
import functools

@functools.total_ordering
class Node():
    def __init__(self, value, index):
        self.left = None
        self.right = None
        self.parent = None
        self.value = value
        self.index = index

    def __eq__(self, other):
        return self.value == other.value and self.index == other.index
    
    def __lt__(self, other):
        return self.value < other.value or (self.value == other.value and self.index < other.index)
    
    def __str__(self):
        return f"Node: {self.value}"
    
class CartesianTree():
    
    """
    Need to make it totally ordered by index (ties broken by index so 
    if it comes first it is smaller)
    """
    
    def __init__(self, arr):
        self.root = Node(arr[0], 0)
        self.size = 1
        self.last_node = self.root
        
        for i in arr[1:]:
            self.insert(i)
        
    def insert(self, data):
        node_index = self.size
        new_node = Node(data, node_index)
        cur_node = self.last_node
        while (cur_node != self.root) and cur_node > new_node:
            cur_node = cur_node.parent
        if cur_node == self.root and cur_node > new_node:
            self.root = new_node
            new_node.left = cur_node
            cur_node.parent = new_node
        else:
            if not cur_node.right is None:
                new_node.left = cur_node.right
                new_node.left.parent = new_node
            cur_node.right = new_node
            new_node.parent = cur_node
        self.last_node = new_node
        self.size += 1
    
    def plot(self, name):
        G = nx.DiGraph()
        def add_nodes(node):
            G.add_node(node.value)
            if not node.parent is None:
                G.add_edge(node.parent.value, node.value)
            if not node.right is None:
                add_nodes(node.right)
            if not node.left is None:
                add_nodes(node.left)
        
        add_nodes(self.root)

        plt.title('draw_networkx')
        p=nx.drawing.nx_pydot.to_pydot(G)
        p.write_png(f"{name}.png")
        
    def get_cartesian_hash(self):
        cart_hash = []
        def update_hash(node):
            cart_hash.append("0")
            if node.left is None:
                cart_hash.append("1")
            else:
                update_hash(node.left)
            
            if node.right is None:
                cart_hash.append("1")
            else:
                update_hash(node.right)
                
        update_hash(self.root)
        return "".join(cart_hash)
    
    def get_in_order(self):
        nodes = []
        def add_node(node):
            if not node.left is None:
                add_node(node.left)
                
            nodes.append(node)
            if not node.right is None:
                add_node(node.right)

        add_node(self.root)
        return nodes    
        
    def get_pairwise_comps(self, comp=min):
        nodes = self.get_in_order()
        arr = np.zeros((len(nodes), len(nodes)))
        
        for i, first in enumerate(nodes):
            cur_extreme = first
            for j, second in enumerate(nodes[i:]):
                cur_extreme = comp(cur_extreme, second)
                arr[i][j+i] = cur_extreme.index
        return arr
                

def topo_sort(g, start_node):
    nodes = g.nodes()
    node_status = {i:False for i in nodes}
    sorted_nodes = []

    def topo_sort_helper(g, node):
        if not node_status[node]:
            node_status[node] = True
            for neighbor in g.neighbors(node):
                topo_sort_helper(g, neighbor)
            sorted_nodes.append(node)
    
    topo_sort_helper(g, start_node)
    return sorted_nodes[::-1]

def super_bubble(g, start_node):
    
    def validate_superbubble(start_node, end_node):
        start = start_node[2]
        end = end_node[2]
        
        outchild = outchild_rangearr.range_extreme(start, end-1)
        outparent = outparent_rangearr.range_extreme(start+1, end)
        
        if outchild != end:
            return -1
        if outparent == start:
            return start_node
        else:
            if outparent not in total_topo_order_to_node:
                return prev_entrances[outparent]
            elif total_topo_order_to_node[outparent][1] == "Entrance":
                return total_topo_order_to_node[outparent]
            else:
                return prev_entrances[outparent]

    
    def report_superbubble(start, exit):
        super_bubbles = []
        if start is None or exit is None or topo_order[start[0]] >= topo_order[exit[0]]:
            node_candidates.pop()
            return super_bubbles
        
        s = prev_entrances[exit[2]]
        
        valid = None
        
        while topo_order[s[0]] >= topo_order[start[0]]:
            valid = validate_superbubble(s, exit)
            if (valid == s) or (valid == alternative_entrances[s[2]]) or (valid == -1):
                break
            alternative_entrances[s[2]] = valid
            s = valid
            
        node_candidates.pop()
        
        if valid == s:
            super_bubbles.append((s, exit))
            while node_candidates and node_candidates[-1] != s:
                if node_candidates[-1][1] == "Exit":
                    next_s = order_to_node[topo_order[s[0]]+1]
                    super_bubbles.extend(report_superbubble(next_s, node_candidates[-1]))
                else:
                    node_candidates.pop()
        return super_bubbles
    
    sorted_nodes = topo_sort(g, start_node)
    node_to_topo = {node:ix for ix, node in enumerate(sorted_nodes)} 
    outparent = []
    outchild = []
    
    node_candidates = []
    prev_entrances = dict()
    
    prev_ent = None
    
    for ix, node in enumerate(sorted_nodes):
        prev_entrances[ix] = prev_ent
        is_entrance = False
        is_exit = False

        topo_orders = [np.inf]
        for parent in g.predecessors(node):
            topo_orders.append(node_to_topo[parent])
            if g.out_degree(parent) == 1 and not is_exit:
                is_exit = True
                node_candidates.append((node, "Exit", ix))
        outparent.append(min(topo_orders))
            
        topo_orders = [-1]
        for child in g.successors(node):
            topo_orders.append(node_to_topo[child])
            if g.in_degree(child) == 1 and not is_entrance:
                is_entrance = True
                node_candidates.append((node, "Entrance", ix))
                prev_ent = (node, "Entrance", ix)
        outchild.append(max(topo_orders))
#         sorted_nodes[ix] = (node, is_entrance, is_exit)

    outparent_rangearr = RangeMinMaxArray(outparent, comp=min)
    outchild_rangearr = RangeMinMaxArray(outchild, comp=max)
    
#     node_candidates, prev_entrances = identify_entrance_exit_nodes(g, start_node)
    alternative_entrances = {i:None for i in prev_entrances}
    topo_order = {i[0]:ix for ix, i in enumerate(node_candidates)}    
    order_to_node = {ix:i for ix, i in enumerate(node_candidates)}  
    #Works because we only ever need to use this to check if a node is an 
    #entrance. It overrides all nodes that are an exit and entrance with their
    #"entrance" entry because they always come after the exit entries for each
    #node...May need to reimplement
    total_topo_order_to_node = {i[2]:i for ix, i in enumerate(node_candidates)}    
#     node_candidates, prev_entrances = node_candidates[::-1], prev_entrances[::-1] #Need to operate on them in reverse order
    entrances = [node[0] for node in node_candidates if node[1] == "Entrance"]
    
    super_bubbles = []
    
    while node_candidates:
        if node_candidates[-1][1] == "Entrance":
            
            node_candidates.pop()
        else:
            
            super_bubbles.extend(report_superbubble(node_candidates[0], node_candidates[-1]))
    
    return super_bubbles, outparent

class RangeMinMaxArray():
    """
    Need to add an array to store the precomputed range extrema for
    all block index queries to make constant lookup of the extrema for
    all blocks fully contained in the query (and then only have three 
    extrema to consider)
    
    Adapted from 
    ftp://nozdr.ru/biblio/kolxoz/Cs/CsLn/Combinatorial%20Pattern%20Matching,%2017%20conf.,%20CPM%202006(LNCS4009,%20Springer,%202006)(ISBN%203540354557)(424s).pdf#page=46
        (Johannes Fischer and Volker Heun)
    https://en.wikipedia.org/wiki/Range_minimum_query
    https://en.wikipedia.org/wiki/Cartesian_tree
    """
    def __init__(self, data, comp=min):
        self.data = data
        self.comp=comp
        self.block_size = int(math.log(len(self.data), 2) // 1)
        self.blocks = [self.data[ix:ix + self.block_size] for ix in range(0,len(self.data),self.block_size)]
        self.block_extrema = [self.comp(arr) for arr in self.blocks]
        
        #Construct long query array
        array_len = len(self.blocks)
        log_len = math.floor(math.log(array_len, 2)) + 1
        self.long_query_array = np.zeros((array_len, log_len))
        for j in range(log_len):
            for ix, i in enumerate(self.block_extrema):
                if j == 0:
                    self.long_query_array[ix][0] = min(self.block_extrema[ix:ix+1])
                else:
                    second_index = ix + 2**(j-1)
                    if second_index >= array_len:
                        self.long_query_array[ix][j] = self.long_query_array[ix][j-1]
                    else:
                        self.long_query_array[ix][j] = min(self.long_query_array[ix][j-1], \
                                              self.long_query_array[second_index][j-1])

        self.build_cartesian_map()
        
    def __iter__(self):
        for i in self.data:
            yield i
            
    def __getitem__(self, key):
        return self.data[key]
    
    def get_pairwise_block_comps(self):
        
        for i, first in enumerate(self.block_extrema):
            cur_extreme = first
            for j, second in enumerate(nodes[i:]):
                cur_extreme = self.comp(cur_extreme, second)
                arr[i][j+i] = cur_extreme.index
        return arr
        
    def get_long_block_extrema(self, i, j):
        
        if i == j:
            return self.block_extrema[i]
        
        l = math.floor(math.log(j-i, 2))
        return self.comp(self.long_query_array[i][l], self.long_query_array[j - (2**l) +1][l])
            
    def build_cartesian_map(self):
        cartesian_map = dict() #I know this breaks the linear time algorithm (worst case is O(n) where the 
        #wikipedia article https://en.wikipedia.org/wiki/Range_minimum_query recommends an array of
        #size log n, but "key in dict" in constant on average by https://wiki.python.org/moin/TimeComplexity
        #So I'm saving some space and "readability")

        block_id_to_cart_hash = dict()
        for ix, block in enumerate(self.blocks):
            cart_tree = CartesianTree(block)
            cart_hash = cart_tree.get_cartesian_hash()
            if not cart_hash in cartesian_map:
                cartesian_map[cart_hash] = cart_tree.get_pairwise_comps(comp=self.comp)
            block_id_to_cart_hash[ix] = cart_hash 

        self.cartesian_map, self.block_id_to_cart_hash = cartesian_map, block_id_to_cart_hash
        
    def range_extreme(self, left, right):
        if left > right:
            raise ValueError("Left index must be less than or equal to the right index.")
        extrema = []
        left_block = math.floor(left/self.block_size)
        right_block = math.floor(right/self.block_size)
        
        if right_block - left_block > 1: #So there is at least one block completely between them
            extrema.append(self.get_long_block_extrema(left_block+1, right_block-1))
        
        #If the indices are in the same block, both left_extrema and right_extrema are the same
        
        left_block_offset = left_block * self.block_size
        left_block_left, left_block_right = (left - left_block_offset, \
                              min(self.block_size - 1, right - left_block_offset))
        
        left_block_cart_hash = self.block_id_to_cart_hash[left_block]
        left_extrema_index = int(self.cartesian_map[left_block_cart_hash][left_block_left][left_block_right])
        left_extrema = self.blocks[left_block][left_extrema_index]
        
        right_block_offset = right_block * self.block_size
        right_block_left, right_block_right = (max(0, left - right_block_offset), \
                                               right - right_block_offset)
        
        right_block_cart_hash = self.block_id_to_cart_hash[right_block]
        right_extrema_index = int(self.cartesian_map[right_block_cart_hash][right_block_left][right_block_right])
        right_extrema = self.blocks[right_block][right_extrema_index]
        
        extrema.append(left_extrema)
        extrema.append(right_extrema)
        
        return self.comp(extrema)
    

g = nx.DiGraph()
g.add_nodes_from(range(1,16))
starts = [1, 1, 2, 3, 3, 3, 4, 5, 5, 6, 6, 7, 8, 8, 9, 10, 11, 12, 13, 13, 15]
ends =  [2, 3, 3, 11, 5, 4, 8, 9, 6, 10, 7, 8, 13, 14, 10, 7, 12, 8, 15, 14, 14]
edges = list(zip(starts, ends))
g.add_edges_from(edges)

# import pdb; pdb.set_trace()
print(super_bubble(g, 1))