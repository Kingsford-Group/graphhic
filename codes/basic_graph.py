import sys
#from collections import defaultdict
import os
import json
import numpy as np
from time import time
#import time
import vg_pb2
import stream
from tqdm import tqdm
from collections import defaultdict, Counter
from Bio import SeqIO
import networkx as nx
import pickle, os
from copy import deepcopy
import pdb
from heapq import *
from gam_utils import *

def is_gz_file(filepath):
    with open(filepath, 'rb') as test_f:
        return test_f.read(2) == b'\x1f\x8b'

class Graph:
    """
        A general genome graph class.
        self.nodes --- dictionary that maps node idx to node label
        self.adj_list --- dictionary of dictionary that maps node to a list of neighbors with the weight of the edge
        self.in/outdeg --- dict that maps node to its in/out degrees
    """

    def __init__(self, nodes=None, edges=None, fname=None, get_degrees=True, chromosome="all"):

        
        assert chromosome in ["all","chrX","chrY"] + ["chr"+str(i) for i in range(1,23)]
        self.chromosome = chromosome
        if edges is None:
            edges = defaultdict(Counter)
        if nodes is None:
            nodes = dict()
        self.nodes = nodes
        self.node_weights = Counter()
        self.adj_list = edges
        self.adj_list_reverse = defaultdict(Counter)
        self.node_id_max = -1

        self.indeg = defaultdict(int)
        self.outdeg = defaultdict(int)

        self.indeg_weighted = defaultdict(int)
        self.outdeg_weighted = defaultdict(int)

        self.top_order = [] # nodes sorted in topological order
        self.top_dict = dict() # query topological order given node id

        #self.nodes_aln_states = defaultdict(int)

        if fname != None:
            suffix = fname.split(".")[-1]
            if suffix == "gfa":
                self.from_gfa(fname)
            elif suffix == "vg":
                self.from_vg(fname)
        else:
            self.get_degrees()
            self.get_adj_reverse()

    def __eq__(self, other):
        if isinstance(other, Graph):
            if self.nodes != other.nodes:
                print("nodes are different!")
                return False
            if self.node_weights != other.node_weights:
                print("node weights are different!")
                return False
            if self.adj_list != other.adj_list:
                print("adjacent lists are different!")
                return False
            if self.adj_list_reverse != other.adj_list_reverse:
                print("reverse adjacent lists are different!")
                return False
            if self.indeg != other.indeg:
                print("in degrees are different!")
                return False
            if self.outdeg != other.outdeg:
                print("out degrees are different!")
                return False
            if self.indeg_weighted != other.indeg_weighted:
                print("in degree weigths are different!")
                return False
            if self.outdeg_weighted != other.outdeg_weighted:
                print("out degree weigths are different!")
                return False
            return True
        return False

    def initialize_nodes_aln_states(self):
        self.nodes_aln_states = defaultdict(int)
        self.empty_nodes = []
        for node in self.nodes:
            self.nodes_aln_states[node] = 0

    def update_nodes_aln_states(self,path):
        for node in path:
            self.nodes_aln_states[node] = 1
    
    def update_empty_nodes(self):
        for node in self.nodes_aln_states:
            if self.nodes_aln_states[node]==0 and self.nodes[node]!="s" and self.nodes[node]!="t":
                self.empty_nodes.append(node)
    
    def clear(self):
        self.nodes = dict()
        self.adj_list = defaultdict(Counter)
        self.adj_list_reverse = defaultdict(Counter)

    def get_source_sink(self):
        sources = [v for v in self.nodes if self.indeg[v] == 0 and self.outdeg[v] != 0]
        sinks = [v for v in self.nodes if self.outdeg[v] == 0 and self.indeg[v] != 0]

        return sources, sinks

    def update_degree(self, node1, node2):
        self.indeg[node2] += 1
        self.outdeg[node1] += 1

        self.indeg_weighted[node2] += self.adj_list[node1][node2]
        self.outdeg_weighted[node1] += self.adj_list[node1][node2]

    def get_degrees(self):
        self.indeg = defaultdict(int)
        self.outdeg = defaultdict(int)
        for node1 in self.adj_list:
            for node2 in self.adj_list[node1]:
                self.update_degree(node1, node2)
        return self.indeg, self.outdeg

    def num_edges(self):
        count = 0
        for n in self.nodes:
            for m in self.adj_list[n]:
                count += 1
        return count

    def from_gfa(self, fname, skip=False):
        """Read from gfa file

        Args:
            fname (_type_): _description_
            skip (bool, optional): Whether skip non ATCGst nodes. Defaults to False.
        """
        self.clear()
        for line in open(fname):
            l = line.rstrip("\n").split("\t")
            if "S" in l[0]:
                if l[1][0]=="s":
                    l1 = l[1][1:]
                else:
                    l1 = l[1]
                # skip nodes with only N or other unknown characters in it
                if "s" ==l[2]:
                    self.nodes[int(l1)] = "$"
                elif "t"==l[2]:
                    self.nodes[int(l1)] = "#"
                else:
                    self.nodes[int(l1)] = l[2]
                '''
                if ("A" in l[2] or "C" in l[2] or "T" in l[2] or "G" in l[2]) or (not skip):
                    self.nodes[int(l[1])] = l[2]
                if ("s" in l[2]) or (not skip):
                    pdb.set_trace()
                    self.nodes[int(l[1])] = "$"
                if ("t" in l[2]) or (not skip):
                    pdb.set_trace()
                    self.nodes[int(l[1])] = "#"
                '''
            if "L" in l[0]:
                if l[1][0]=="s":
                    l1 = l[1][1:]
                    l3 = l[3][1:]
                else:
                    l1 = l[1]
                    l3 = l[3]
                n1 = int(l1)
                n2 = int(l3)
                if l[2]!="+" or l[4]!="+":
                    print(l[2],l[4])
                if n1 in self.nodes and n2 in self.nodes:
                    self.add_edge(n1, n2)

    def from_vg(self,fname):
        self.clear()
        '''
        if is_gz_file(fname):
            istream = stream.open(fname, "rb",gzip=)
        else:
            istream = stream.open(fname, "rb",gzip=False)
        '''
        with stream.open(fname, "rb",gzip=is_gz_file(fname)) as istream:
            for data in istream:
                try:
                    vg = vg_pb2.Graph()
                    vg.ParseFromString(data)
                    #pdb.set_trace()
                    '''
                    if len(vg.path)>=1:
                        print(len(vg.path),vg.path[0].name)
                    '''
                    if self.chromosome!="all" and len(vg.path)==1 and vg.path[0].name!=self.chromosome:
                        continue
                    #vg_list.append(vg)
                    for node in vg.node:
                        self.add_node(node.sequence,node.id)
                        #self.nodes[node.id] = node.sequence
                    for edge in vg.edge:
                        n1 = getattr(edge, "from")
                        n2 = int(edge.to)
                        '''
                        if n1 not in self.nodes or n2 not in self.nodes:
                            print("node currently not been readed")
                        '''
                        self.add_edge(n1,n2)
                        '''
                        if n1 in self.nodes and n2 in self.nodes:
                            self.add_edge(n1,n2)
                        '''
                    #pdb.set_trace()
                except Exception:
                    if data != b'VG':
                        print(data)
                        exit()
                    #pdb.set_trace()
                    continue
        #istream.close()

    '''
    def from_vg(self, fname):
        self.clear()

        # read protobuf format of vg
        vg_list = []
        with stream.open(fname, "rb") as istream:
            for data in istream:
                try:
                    vg = vg_pb2.Graph()
                    vg.ParseFromString(data)
                    vg_list.append(vg)
                except Exception:
                    continue

        for vg in vg_list:
            for node in vg.node:
                self.nodes[node.id] = node.sequence
                #if "A" in node.sequence or "C" in node.sequence or "T" in node.sequence or "G" in node.sequence or "s" in node.sequence or "t" in node.sequence:
                #   self.nodes[node.id] = node.sequence

        for vg in vg_list:
            for edge in vg.edge:
                n1 = getattr(edge, "from")
                n2 = int(edge.to)
                if n1 in self.nodes and n2 in self.nodes:
                    self.add_edge(n1,n2)
    '''

    def to_gfa(self, fname):
        with open(fname, "w") as ofile:
            ofile.write("\t".join(["H", "VN:Z:1.0"]) + "\n")

            for n in self.nodes:
                strings = self.nodes[n]
                if strings == "$":
                    strings = "s"
                if strings == "#":
                    strings = "t"
                ofile.write("\t".join(["S", str(n), strings]) + "\n")

            for u in self.adj_list:
                for v in self.adj_list[u]:
                    for _ in range(self.adj_list[u][v]):
                        ofile.write("\t".join(["L", str(u), "+", str(v), "+", "0M"]) + "\n")

    def get_nx_graph(self, flow=False):
        graph = nx.DiGraph()

        for u in self.adj_list:
            for v in self.adj_list[u]:

                if flow:
                    if self.nodes[u] == "#" and self.nodes[v] == "$":
                        continue
                    graph.add_edge(u, v, capacity=self.adj_list[u][v])
                else:
                    graph.add_edge(u, v)
        return graph

    def add_node(self, seq, nodeid=-100):

        # assign a node id if not provided
        if nodeid == -100:
            if len(self.nodes) == 0:
                node_id = 0
            else:
                for i in range(len(self.nodes)+1):
                    if i not in self.nodes:
                        nodeid = i
                        break
        self.nodes[nodeid] = seq
        self.update_max_node_id(nodeid)
        return nodeid

    def add_node_weights_from_aln(self, aln_dict):
        for aln_name in aln_dict:
            path = aln_dict[aln_name].path
            for n in path:
                self.node_weights[n] += 1

    def update_max_node_id(self,nodeid):
        self.node_id_max = max(self.node_id_max, nodeid + 1)

    def remove_node(self,nodeid):
        if nodeid not in self.nodes:
            print("Fail to remove node %i -- nodes do not exist" % nodeid)
            return
        seq = self.nodes.pop(nodeid)
        return seq

    def add_edge(self, node1, node2):
        '''
        if node1 not in self.nodes or node2 not in self.nodes:
            print("[Warning] (%i,%i)--- nodes do not exist" % (node1, node2))
            #print("Fail to add edge (%i,%i)--- nodes do not exist" % (node1, node2))
            #return
        '''

        if self.adj_list[node1][node2] == 0:
            self.indeg[node2] += 1
            self.outdeg[node1] += 1

        self.adj_list[node1][node2] = 1
        self.adj_list_reverse[node2][node1] = 1

        self.indeg_weighted[node2] += 1
        self.outdeg_weighted[node1] += 1

    def remove_edge(self,parent_node,node):
        #pdb.set_trace()
        if parent_node not in self.nodes or node not in self.nodes:
            if parent_node not in self.nodes:
                print("cannot find", parent_node)
            if node not in self.nodes:
                print("cannot find", node)
            print("Fail to remove edge (%i,%i)--- nodes do not exist" % (parent_node, node))
            assert 1==-1
            return 
        edge_weight = self.adj_list[parent_node].pop(node)
        edge_weight_2 = self.adj_list_reverse[node].pop(parent_node)
        assert edge_weight==edge_weight_2

        self.outdeg[parent_node] -=1
        self.outdeg_weighted[parent_node] -=edge_weight
        self.indeg[node] -=1
        self.indeg_weighted[node]-=edge_weight
        assert self.indeg[node] == len(self.adj_list_reverse[node]) 
        assert self.outdeg[parent_node] == len(self.adj_list[parent_node])
 
    def vg_graph_add_source_sink(self, print_msg = False):
        # only use for graphs that don't have source or sink nodes but are DAGs

        sources, sinks = self.get_source_sink()

        # add source
        source = self.add_node("$")
        for s in sources:
            self.add_edge(source, s)

        # add sink
        sink = self.add_node("#")
        for s in sinks:
            self.add_edge(s, sink)

        if print_msg:
            print("Old sources", sources)
            print("Old sinks", sinks)
            print("New source", source)
            print("New sink", sink)

        return self

    def get_adj_reverse(self):
        # record the father node
        self.adj_list_reverse = defaultdict(Counter)
        for u in self.adj_list:
            for v in self.adj_list[u]:
                self.adj_list_reverse[v][u]+=1

    def longest_path(self, node1, node2, score = "node_length"):
            """
                Find the longest path between node 1 and node 2 in the graph. Current implementation only counts the number of nodes in the path.

                Args:
                - node1, node2: two nodes in the graph. node1 is source and node2 is sink of the path
                - topological_sort (bool) --- whether redo topological sort
                - score (str) --- "node_length" or "node_weight"

                Return:
                - longest_path (list) --- a list of node ids
            """
            '''
            if top_order == None:
                top_order, top_dict = self.topological_sort()
            '''
            #pdb.set_trace()
            if self.top_order == []:
                _,_ = self.topological_sort()
            
            #pdb.set_trace()


            if score == "node_weight" and len(self.node_weights) == 0:
                assert False, "run self.add_node_weights_from_aln first."

            longest_path = [node2]

            curr_idx = self.top_dict[node1]
            end_idx = self.top_dict[node2]
            stack = [node1]

            dists = defaultdict(int)
            dists[node1] = 0
            parents = defaultdict(int)
            parents[node1] = -1
            for curr_node in self.top_order[self.top_dict[node1]: self.top_dict[node2]+1]:

                if curr_node != node1:
                    for parent in self.adj_list_reverse[curr_node]:

                        curr_score = len(self.nodes[parent])
                        if score == "node_weight":
                            curr_score = self.node_weights[parent]
                        
                        if score == "node_length":
                            if dists[parent] + curr_score > dists[curr_node]:
                                parents[curr_node] = parent
                                dists[curr_node] = dists[parent] + curr_score
                        if score == "node_weight":
                            if dists[parent] + curr_score >= dists[curr_node]:
                                parents[curr_node] = parent
                                dists[curr_node] = dists[parent] + curr_score
                curr_idx += 1
                if curr_idx <= end_idx:
                    stack.append(self.top_order[curr_idx])

            parent = parents[node2]
            while parent != -1:
                longest_path.append(parent)
                parent = parents[parent]
                
            return longest_path[::-1]

    def find_bubble(self, source):
        """
            Find the bubble that starts at a given node. If no bubble is found, return the source.
            TODO: currently this algorithm does not check for tips or cycles. It is not necessary right now but should be added for safety reasons.

            Args:
                - source: a node in the graph
            
            Return:
                - bubble --- a list of nodes sorted in topological order in the bubble
                - sink --- the end of the bubble
        """
        
        node_label = defaultdict(str)
        S = set([source])
        
        seen_nodes = set()
        visited_nodes = set()
        
        sink = -1
        
        while len(S) > 0:
            #pdb.set_trace()
            
            curr_node = S.pop()
            node_label[curr_node] = "visited"
            visited_nodes.add(curr_node)
            seen_nodes.discard(curr_node)
            
            for neighbor in self.adj_list[curr_node]:
                
                node_label[neighbor] = "seen"
                visited_nodes.discard(neighbor)
                seen_nodes.add(neighbor)
                
                all_visited_flag = True
                for parent in self.adj_list_reverse[neighbor]:
                    if node_label[parent] != "visited":
                        all_visited_flag = False
                        break
                if all_visited_flag:
                    S.add(neighbor)

            if len(S) == 1 and S == seen_nodes:
                sink = S.pop()
            
        if sink != -1:
            bubble = self.top_order[self.top_dict[source]:self.top_dict[sink]+1]
        else:
            bubble = [source]
        return bubble, sink

    def find_all_bubbles(self):
        """
            Find all bubble structures in the graph and order them into hierarchical structures. A bubble is another bubble's child if both 
            its start and end are included in the parent.
            
            Return:
            - all_bubbles --- a dictionary bubbles, key is the node id of bubble source
            - bubble_hierarchy --- a Graph object describing the tree structure of bubble hierarchy.
            - bubble_layers --- a list of list representing layers of bubbles as a result of doing BFS in bubble_hierarchy
        """

        self.topological_sort()
        
        bubble_hierarchy = Graph()
        bubble_hierarchy.add_node(nodeid=-1, seq=-1)

        all_bubbles = defaultdict()
        prev_bubble = (-1,-1)
        parent = -1
        #pdb.set_trace()

        for node in self.top_order:

            # treat the source and the sink of the graph as a bubble
            if self.top_dict[node] == 0 or self.top_dict[node] == len(self.top_order)-1:
                all_bubbles[node] = [node]
                bubble_hierarchy.add_node(nodeid=node, seq=node)
                bubble_hierarchy.add_edge(-1, node)
                    
            if len(self.adj_list[node]) ==1:
                continue

            bubble, sink = self.find_bubble(node)

            if sink != -1:
                all_bubbles[node] = bubble
                bubble_hierarchy.add_node(nodeid=node, seq=sink)

                while parent != -1 and self.top_dict[bubble_hierarchy.nodes[parent]] < self.top_dict[sink]:
                    parent = list(bubble_hierarchy.adj_list_reverse[parent].keys())[0]

                if parent == -1 or self.top_dict[bubble_hierarchy.nodes[parent]] >= self.top_dict[sink] :
                    bubble_hierarchy.add_edge(parent, node)
                    parent = node

        # put bubbles into layers
        bubble_layers = []

        source = -1
        Q = [(-1,source)]

        curr_layer = -1
        visited = set()
        #pdb.set_trace()
        while len(Q) > 0:

            curr_node = Q.pop()
            #pdb.set_trace()

            layer = curr_node[0]
            if layer != -1:

                if len(bubble_layers) == layer:
                    bubble_layers.append([curr_node[1]])
                else:
                    bubble_layers[layer].append(curr_node[1])
            for neighbor in bubble_hierarchy.adj_list[curr_node[1]]:
                if neighbor in visited:
                    continue
                Q.insert(0,(layer+1, neighbor))
                visited.add(neighbor)

        return all_bubbles, bubble_hierarchy, bubble_layers


    def topological_sort(self):
        """
            Perform topological sort and updates top_order and top_dict fields

            Return:
            - top_order
            - top_dict
        """
        sources,sinks = self.get_source_sink()

        self.top_order = []
        visited = set()
        pushed = set()
        stack = []
        for s in sources:
            stack.append(s)
            while len(stack) != 0:
                node = stack[-1]
                
                if node in pushed:
                    stack.pop()
                    continue
                
                visited.add(node)
                
                all_visited = True
                
                for neighbor in self.adj_list[node]:
                    if neighbor not in visited:
                        stack.append(neighbor)
                        all_visited = False
                        
                if all_visited:
                    new_node = stack.pop()
                    pushed.add(new_node)
                    self.top_order.append(new_node)
        self.top_order = self.top_order[::-1]
        self.top_dict = dict([(n,i) for i,n in enumerate(self.top_order)])
        return self.top_order, self.top_dict

    def get_total_characters(self):
        total_count = 0
        for n in self.nodes:
            try:
                total_count += len(self.nodes[n])
            except TypeError:
                print(n)
        return total_count

    def get_num_edges(self):
        count = 0
        for n1 in self.adj_list:
            for n2 in self.adj_list[n1]:
                count += 1
        return count

    def get_avg_node_size(self):
        return np.mean([len(self.nodes[n]) for n in self.nodes])

    def get_median_node_size(self):
        return np.median([len(self.nodes[n]) for n in self.nodes])

    def get_avg_node_size_gt_10(self):
        return np.mean([len(self.nodes[n]) for n in self.nodes if len(self.nodes[n]) > 10])

    def print_graph_stats(self, title=""):
        print("---- %s graph stats ----" % title)
        print("* Num of nodes", len(self.nodes))
        print("* Num of edges", self.get_num_edges())
        print("* Total length of nodes", self.get_total_characters())
        print("* Average length of nodes", self.get_avg_node_size())
        print("* Average length of nodes gt 10", self.get_avg_node_size_gt_10())
        print("* Median length of nodes", self.get_median_node_size())