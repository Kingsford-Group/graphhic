from basic_graph import *
from math import ceil

class CompactGraph:

    def __init__(self, old_graph = None):

        self.old_graph = old_graph
        self.new_graph = Graph()
        self.num_nodes = 0
        self.first_to_node = {}
        self.last_to_node = {}
        self.old_to_new = {}
        self.new_to_old = {}
        self.node_size = -1

    def get_compact_graph(self, **kwargs):
        """
        Wrapper for different compacting functions
        
        Input:
            **argument --- a keyword dictionary. The first keyword must be "mode = <str>" specifying which mode of compaction is done.
            Other keywords are specified in different modes below.

        Current supporting modes:
        1. mode = branching
            - merges all 1-in-1-out nodes. No other keywords or arguments need to be provided
        
        2. mode = empty_nodes
            - removes all nodes with zero read coverage.
            - Requires another argument, a dictionary of alignment objects.
            - example: {mode = "empty_node", aln_dict = {aln_name:Alignment}}

        3. mode = small_variations
            - merges nodes that represent small variations with size below a threshold
            - Requires another argument, a size threshold
            - example: {mode = "small_variations", threshold=10}
        """
        if "mode" not in kwargs:
            print("get_compact_graph: this function is not properly called. Must specify mode!")
            print("\t - branching")
            print("\t - empty_nodes, argument = aln_dict")
            print("\t - small_variations, argument = threshold")
            exit()

        mode = kwargs["mode"]

        if mode == "branching":
            print("Merging non-branching nodes", flush=True)
            self.get_compact_graph_branching()

        if mode == "empty_nodes":
            if "aln_dict" not in kwargs:
                print("get_compact_graph: this function is not properly called. Must specify aln_dict!")
                exit()

            print("Merging empty nodes", flush=True)
            self.get_compact_graph_empty_nodes(aln_dict = kwargs["aln_dict"])
        
        if mode == "small_variations":
            if "threshold" not in kwargs:
                print("get_compact_graph: this function is not properly called. Must specify threshold!")
                exit()

            print("Merging small variations", flush=True)
            self.get_compact_graph_small_variations(threshold = kwargs["threshold"])


    def path_projection(self,alignment):
        """
        Convert an aligned path in terms of nodes in compact_graph
        
        Args:
            compact_graph (CompactGraph) -- graph after compaction
            alignment (Alignment object) -- one alignment object 
        Return:
            list [new_path, offset, reverse] 
                - new_path (list) -- a list of nodes in compact_graph that represents the aligned path
                - offset (int) -- the start location of the alignment on the first aligned node in terms of bases
                - isReversed (int) -- whether the alignment orientation is reversed or not
        """
        #pdb.set_trace()
        if len(alignment.path)==0:
            return alignment
        #pdb.set_trace()
        new_path = []
        first_node = alignment.path[0]
        if first_node not in self.old_to_new:
            assert first_node in self.new_graph.nodes
            new_path.append(first_node)
            offset = alignment.offset
        else:
            first_new_node = self.old_to_new[first_node]
            assert first_new_node in self.new_graph.nodes
            new_path.append(first_new_node)
            offset = 0
            if alignment.isReversed==1:
                new_to_old_node_list = self.new_to_old[first_new_node][::-1]
                alignment.isReversed = 0
            elif alignment.isReversed==0:
                new_to_old_node_list = self.new_to_old[first_new_node]
            for new_to_old_node in new_to_old_node_list:
                if new_to_old_node==first_node:
                    break
                offset+=len(self.old_graph.nodes[new_to_old_node])
            offset+=alignment.offset
        for i in range(1,len(alignment.path)):
            old_node = alignment.path[i]
            if old_node in self.old_to_new:
                new_node = self.old_to_new[old_node]
            else:
                assert old_node in self.new_graph.nodes
                new_node = old_node
            if new_path[-1]!=new_node:
                new_path.append(new_node)

        alignment.path = new_path
        alignment.offset = offset
        return alignment


# ============================================================================================
#
#                            COMPACT MODE = REMOVE NON-BRANCHING NODES
#
# ============================================================================================

    def update_new_graph_nodes_aln_states_with_branching(self):
        self.new_graph.initialize_nodes_aln_states()
        for node in self.new_to_old:
            #aln_marker = 0
            for old_node in self.new_to_old[node]:
                self.new_graph.nodes_aln_states[node] = self.new_graph.nodes_aln_states[node] or self.old_graph.nodes_aln_states[old_node]


    def concat_node_label(self, node_list, new_node_id) -> str:
        """
            Merge a list of nodes into new node. The only function that updates old_to_new and new_to_old mapping.

            Args:
                node_list --- input list of nodes to be merged, node id in old_graph
                new_node_id --- purposed new id for the node in new_graph. It is required that the new node id has not appeared in new_graph yet
        
            Return:
                concatenated node label
        """
        s = ""
        assert new_node_id not in self.new_to_old.keys()
        self.new_to_old[new_node_id] = {}
        for i,n in enumerate(node_list):
            #s += self.old_graph.nodes[n]
            self.old_to_new[n] = new_node_id
            self.new_to_old[new_node_id][n] = [i,len(s)]
            s += self.old_graph.nodes[n]
        return s

    def add_compact_node(self, non_branch_path):
        """
            Merge a non-branching path as a new compact node in the new graph.
        """
        self.first_to_node[non_branch_path[0]] = self.num_nodes
        self.last_to_node[non_branch_path[-1]] = self.num_nodes
        self.new_graph.add_node(self.concat_node_label(non_branch_path,self.num_nodes), self.num_nodes)
        self.num_nodes += 1
    
    def get_compact_graph_branching(self):

        source, sink = self.old_graph.get_source_sink()
        if len(source) == 1 and len(sink) == 1:
            source = source[-1]
            sink = sink[-1]
        else:
            assert False
        self.old_graph.get_adj_reverse()

        for v in self.old_graph.nodes:
            # skip one-in-one-out nodes. They will be taken care of by creating non-branching paths
            if self.old_graph.indeg[v] != 1 or self.old_graph.outdeg[v] != 1 or v==source or v==sink:

                # keep source and sink as individual nodes
                if v == source or v == sink or (self.old_graph.indeg[v] > 1 and self.old_graph.outdeg[v] > 1):
                    self.add_compact_node([v])

                non_branch_path = []                
                if self.old_graph.outdeg[v] == 1 and v!=source and v!= sink:
                    non_branch_path.append(v)
                for u in self.old_graph.adj_list[v]:

                    while self.old_graph.indeg[u] == 1 and self.old_graph.outdeg[u] == 1 and u!=source and u!=sink:
                        non_branch_path.append(u)
                        u = list(self.old_graph.adj_list[u].keys())[-1]

                    if self.old_graph.indeg[u] == 1 and u!=sink:
                        non_branch_path.append(u)
                
                    if len(non_branch_path) > 0:
                        self.add_compact_node(non_branch_path)

                    non_branch_path = []

        for u in self.old_graph.adj_list:
            for v in self.old_graph.adj_list[u]:
                if u not in self.last_to_node or v not in self.first_to_node:
                    continue
                n1 = self.last_to_node[u]
                n2 = self.first_to_node[v]
                self.new_graph.adj_list[n1][n2] += self.old_graph.adj_list[u][v]
                self.new_graph.adj_list_reverse[n2][n1] += self.old_graph.adj_list_reverse[v][u]
                self.new_graph.update_degree(n1, n2)

# ============================================================================================
#
#                            COMPACT MODE = REMOVE EMPTY NODES
#
# ============================================================================================

    def get_compact_graph_empty_nodes(self, aln_dict):
        self.get_empty_nodes(aln_dict)
        self.remove_empty_nodes()
    
    #def update_nodes_aln_states(self,path):


    def get_empty_nodes(self,aln_dict):
        """
        Obtain nodes with zero read coverage given an alignment dictionary data structure
        """
        self.nodes_aln_states = defaultdict(int)
        for node in self.old_graph.nodes:
            self.nodes_aln_states[node] = 0

        for aln_name in aln_dict:
            path = aln_dict[aln_name].path
            if len(path)==0:
                continue
            #pdb.set_trace()
            for node in path:
                self.nodes_aln_states[node] = 1
        # nodes that no reads aligned
        self.empty_nodes = []
        for node in self.nodes_aln_states:
            if self.nodes_aln_states[node]==0 and self.old_graph.nodes[node]!="s" and self.old_graph.nodes[node]!="t":
                self.empty_nodes.append(node)
        # print(len(self.empty_nodes))
    
    def remove_empty_nodes_without_merge(self):
        sources, sinks = self.old_graph.get_source_sink()
        source = sources[0]
        sinks = sinks[0]

        self.new_graph = deepcopy(self.old_graph)
        self.empty_nodes = deepcopy(self.old_graph.empty_nodes)
        self.nodes_aln_states = deepcopy(self.old_graph.nodes_aln_states)

        while len(self.empty_nodes)>0:
            node = self.empty_nodes.pop(0)
            # if its parent nodes have more than one childs, and its child nodes have more than one parents, than this node is removable
            removable_marker = 1
            for parent_node in self.new_graph.adj_list_reverse[node]:
                if len(self.new_graph.adj_list[parent_node])<=1:
                    removable_marker = 0
                    break
            if removable_marker == 0:
                continue
            for child_node in self.new_graph.adj_list[node]:
                if len(self.new_graph.adj_list_reverse[child_node])<=1:
                    removable_marker = 0
                    break
            if removable_marker == 0:
                continue
            # begin remove
            # remove edges
            parent_nodes = list(self.new_graph.adj_list_reverse[node].keys())
            child_nodes = list(self.new_graph.adj_list[node].keys())
            for parent_node in parent_nodes:
                self.new_graph.remove_edge(parent_node,node)
            for child_node in child_nodes:
                self.new_graph.remove_edge(node,child_node)
            _ = self.new_graph.remove_node(node)
    
    def remove_empty_nodes_3(self):
        self.old_to_new = {}
        self.new_to_old = {}

        sources, sinks = self.old_graph.get_source_sink()
        source = sources[0]
        sink = sinks[0]

        self.new_graph = deepcopy(self.old_graph)
        self.empty_nodes = deepcopy(self.old_graph.empty_nodes)
        self.nodes_aln_states = deepcopy(self.old_graph.nodes_aln_states)

        while len(self.empty_nodes)>0:
            node = self.empty_nodes.pop(0)
            ### nodes may be merged
            '''
            if node >=59060 and node <=59067:
                pdb.set_trace()
            '''
            if node not in self.new_graph.nodes:
                continue
            # if its parent nodes have more than one childs, and its child nodes have more than one parents, than this node is removable
            removable_marker = 1
            for parent_node in self.new_graph.adj_list_reverse[node]:
                if len(self.new_graph.adj_list[parent_node])<=1:
                    removable_marker = 0
                    break
            if removable_marker == 0:
                continue
            for child_node in self.new_graph.adj_list[node]:
                if len(self.new_graph.adj_list_reverse[child_node])<=1:
                    removable_marker = 0
                    break
            if removable_marker == 0:
                continue
            # begin remove
            # remove edges
            parent_nodes = list(self.new_graph.adj_list_reverse[node].keys())
            child_nodes = list(self.new_graph.adj_list[node].keys())
            for parent_node in parent_nodes:
                self.new_graph.remove_edge(parent_node,node)
            for child_node in child_nodes:
                self.new_graph.remove_edge(node,child_node)
            _ = self.new_graph.remove_node(node)
            # merge nodes
            for anchor_node in parent_nodes+child_nodes:
                if anchor_node not in self.new_graph.nodes:
                    continue
                non_branch_path = []
                aln_marker = 0
                #first do the backward search
                v = anchor_node
                while self.new_graph.outdeg[v]==1 and v!=source:
                    non_branch_path.insert(0,v)
                    aln_marker = aln_marker or self.nodes_aln_states[v]
                    if self.new_graph.indeg[v]==1:
                        assert len(self.new_graph.adj_list_reverse[v])==1
                        v = list(self.new_graph.adj_list_reverse[v])[0]
                    else:
                        break
                if len(non_branch_path)==0:
                    continue
                #next do the forward search
                v = anchor_node
                while self.new_graph.outdeg[v]==1:
                    assert len(self.new_graph.adj_list[v])==1
                    v = list(self.new_graph.adj_list[v])[0]
                    if self.new_graph.indeg[v]==1 and v!=sink:
                        non_branch_path.append(v)
                        #pdb.set_trace()
                        aln_marker = aln_marker or self.nodes_aln_states[v]
                    else:
                        break
                if len(non_branch_path)>1:
                    # print("added")
                    #pdb.set_trace()
                    #trash_nodes = trash_nodes + non_branch_path
                    new_node = self.add_merge_node(non_branch_path)
                    prev_connect_nodes = list(self.new_graph.adj_list_reverse[non_branch_path[0]].keys())
                    for node1 in prev_connect_nodes:
                        #self.new_graph.adj_list[node1][new_node] += self.new_graph.adj_list[node1][non_branch_path[0]]
                        #self.new_graph.adj_list_reverse[new_node][node1] += self.new_graph.adj_list_reverse[non_branch_path[0]][node1]
                        #self.new_graph.update_degree(node1,new_node)
                        self.new_graph.add_edge(node1,new_node)
                        self.new_graph.remove_edge(node1,non_branch_path[0])
                    after_connect_nodes = list(self.new_graph.adj_list[non_branch_path[-1]].keys())
                    for node2 in after_connect_nodes:
                        #self.new_graph.adj_list[new_node][node2] +=self.new_graph.adj_list[non_branch_path[-1]][node2]
                        #self.new_graph.adj_list_reverse[node2][new_node] +=self.new_graph.adj_list_reverse[node2][non_branch_path[-1]]
                        #self.new_graph.update_degree(new_node,node2)
                        self.new_graph.add_edge(new_node,node2)
                        self.new_graph.remove_edge(non_branch_path[-1],node2)
                    #pdb.set_trace()
                    self.nodes_aln_states[new_node] = aln_marker
                    if aln_marker==0:
                        self.empty_nodes.append(new_node)
                    for trash_node in non_branch_path:
                        _ = self.new_graph.remove_node(trash_node)

    def remove_empty_nodes_2(self):
        #First step: check if the emoty node is removable
        #pdb.set_trace()
        #self.get_empty_nodes(aln_dict)
        #pdb.set_trace()
        self.old_to_new = {}
        self.new_to_old = {}

        sources, sinks = self.old_graph.get_source_sink()
        source = sources[0]
        sink = sinks[0]

        self.new_graph = deepcopy(self.old_graph)
        self.empty_nodes = deepcopy(self.old_graph.empty_nodes)
        self.nodes_aln_states = deepcopy(self.old_graph.nodes_aln_states)
        # print(len(self.new_graph.nodes))

        while len(self.empty_nodes)>0:
            #pdb.set_trace()
            node = self.empty_nodes.pop(0)
            # print(node)
            '''
            if node > 113151:
                pdb.set_trace()
            '''
            if node not in self.old_graph.nodes:
                continue
            # if its parent nodes have more than one childs, and its child nodes have more than one parents, than this node is removable
            removable_marker = 1
            for parent_node in self.new_graph.adj_list_reverse[node]:
                if len(self.new_graph.adj_list[parent_node])<=1:
                    removable_marker = 0
                    break
            if removable_marker==0:
                continue
            for child_node in self.old_graph.adj_list[node]:
                if len(self.old_graph.adj_list_reverse[child_node])<=1:
                    removable_marker = 0
                    break
            if removable_marker==0:
                continue
            # begin remove
            # remove edges
            parent_nodes = list(self.new_graph.adj_list_reverse[node].keys())
            child_nodes = list(self.new_graph.adj_list[node].keys())
            for parent_node in parent_nodes:
                self.new_graph.remove_edge(parent_node,node)
            for child_node in child_nodes:
                self.new_graph.remove_edge(node,child_node)
            _ = self.new_graph.remove_node(node)

            #merge nodes 
            trash_nodes = []
            for anchor_node in parent_nodes+child_nodes:
                non_branch_path = []
                aln_marker = 0
                #first do the backward search
                v = anchor_node
                if anchor_node in trash_nodes:
                    continue
                while self.old_graph.outdeg[v]==1 and v!=source:
                    non_branch_path.insert(0,v)
                    #pdb.set_trace()
                    aln_marker = aln_marker or self.nodes_aln_states[v]
                    if self.old_graph.indeg[v]==1:
                        assert len(self.new_graph.adj_list_reverse[v])==1
                        v = list(self.new_graph.adj_list_reverse[v])[0]
                    else:
                        break
                if len(non_branch_path)==0:
                    continue
                #next do the forward search
                v = anchor_node
                while self.new_graph.outdeg[v]==1:
                    assert len(self.new_graph.adj_list[v])==1
                    v = list(self.new_graph.adj_list[v])[0]
                    if self.old_graph.indeg[v]==1 and v!=sink:
                        non_branch_path.append(v)
                        #pdb.set_trace()
                        aln_marker = aln_marker or self.nodes_aln_states[v]
                    else:
                        break
                if len(non_branch_path)>1:
                    print("added")
                    trash_nodes = trash_nodes + non_branch_path
                    new_node = self.add_merge_node(non_branch_path)
                    prev_connect_nodes = list(self.new_graph.adj_list_reverse[non_branch_path[0]].keys())
                    for node1 in prev_connect_nodes:
                        self.new_graph.adj_list[node1][new_node] += self.new_graph.adj_list[node1][non_branch_path[0]]
                        self.new_graph.adj_list_reverse[new_node][node1] += self.new_graph.adj_list_reverse[non_branch_path[0]][node1]
                        self.new_graph.update_degree(node1,new_node)
                        self.new_graph.remove_edge(node1,non_branch_path[0])
                    after_connect_nodes = list(self.new_graph.adj_list[non_branch_path[-1]].keys())
                    for node2 in after_connect_nodes:
                        self.new_graph.adj_list[new_node][node2] +=self.new_graph.adj_list[non_branch_path[-1]][node2]
                        self.new_graph.adj_list_reverse[node2][new_node] +=self.new_graph.adj_list_reverse[node2][non_branch_path[-1]]
                        self.new_graph.update_degree(new_node,node2)
                        self.new_graph.remove_edge(non_branch_path[-1],node2)
                    #pdb.set_trace()
                    self.nodes_aln_states[new_node] = aln_marker
                    if aln_marker==0:
                        self.empty_nodes.append(new_node)
            for trash_node in trash_nodes:
                _ = self.new_graph.remove_node(trash_node)


    def remove_empty_nodes(self):
        #First step: check if the emoty node is removable
        #pdb.set_trace()
        #self.get_empty_nodes(aln_dict)
        #pdb.set_trace()
        self.old_to_new = {}
        self.new_to_old = {}

        sources, sinks = self.old_graph.get_source_sink()
        source = sources[0]
        sink = sinks[0]

        self.new_graph = deepcopy(self.old_graph)
        #self.empty_nodes = deepcopy(self.old_graph.empty_nodes)
        # print(len(self.new_graph.nodes))

        while len(self.empty_nodes)>0:
            #pdb.set_trace()
            node = self.empty_nodes.pop(0)
            # print(node)
            '''
            if node > 113151:
                pdb.set_trace()
            '''
            if node not in self.old_graph.nodes:
                continue
            # if its parent nodes have more than one childs, and its child nodes have more than one parents, than this node is removable
            removable_marker = 1
            for parent_node in self.new_graph.adj_list_reverse[node]:
                if len(self.new_graph.adj_list[parent_node])<=1:
                    removable_marker = 0
                    break
            if removable_marker==0:
                continue
            for child_node in self.old_graph.adj_list[node]:
                if len(self.old_graph.adj_list_reverse[child_node])<=1:
                    removable_marker = 0
                    break
            if removable_marker==0:
                continue
            # begin remove
            # remove edges
            parent_nodes = list(self.new_graph.adj_list_reverse[node].keys())
            child_nodes = list(self.new_graph.adj_list[node].keys())
            for parent_node in parent_nodes:
                self.new_graph.remove_edge(parent_node,node)
            for child_node in child_nodes:
                self.new_graph.remove_edge(node,child_node)
            _ = self.new_graph.remove_node(node)

            #merge nodes 
            trash_nodes = []
            for anchor_node in parent_nodes+child_nodes:
                non_branch_path = []
                aln_marker = 0
                #first do the backward search
                v = anchor_node
                if anchor_node in trash_nodes:
                    continue
                while self.old_graph.outdeg[v]==1 and v!=source:
                    non_branch_path.insert(0,v)
                    aln_marker = aln_marker or self.nodes_aln_states[v]
                    if self.old_graph.indeg[v]==1:
                        assert len(self.new_graph.adj_list_reverse[v])==1
                        v = list(self.new_graph.adj_list_reverse[v])[0]
                    else:
                        break
                if len(non_branch_path)==0:
                    continue
                #next do the forward search
                v = anchor_node
                while self.new_graph.outdeg[v]==1:
                    assert len(self.new_graph.adj_list[v])==1
                    v = list(self.new_graph.adj_list[v])[0]
                    if self.old_graph.indeg[v]==1 and v!=sink:
                        non_branch_path.append(v)
                        aln_marker = aln_marker or self.nodes_aln_states[v]
                    else:
                        break
                if len(non_branch_path)>1:
                    # print("added")
                    trash_nodes = trash_nodes + non_branch_path
                    new_node = self.add_merge_node(non_branch_path)
                    prev_connect_nodes = list(self.new_graph.adj_list_reverse[non_branch_path[0]].keys())
                    for node1 in prev_connect_nodes:
                        self.new_graph.adj_list[node1][new_node] += self.new_graph.adj_list[node1][non_branch_path[0]]
                        self.new_graph.adj_list_reverse[new_node][node1] += self.new_graph.adj_list_reverse[non_branch_path[0]][node1]
                        self.new_graph.update_degree(node1,new_node)
                        self.new_graph.remove_edge(node1,non_branch_path[0])
                    after_connect_nodes = list(self.new_graph.adj_list[non_branch_path[-1]].keys())
                    for node2 in after_connect_nodes:
                        self.new_graph.adj_list[new_node][node2] +=self.new_graph.adj_list[non_branch_path[-1]][node2]
                        self.new_graph.adj_list_reverse[node2][new_node] +=self.new_graph.adj_list_reverse[node2][non_branch_path[-1]]
                        self.new_graph.update_degree(new_node,node2)
                        self.new_graph.remove_edge(non_branch_path[-1],node2)
                    self.nodes_aln_states[new_node] = aln_marker
                    if aln_marker==0:
                        self.empty_nodes.append(new_node)
            for trash_node in trash_nodes:
                _ = self.new_graph.remove_node(trash_node)

    def concat_label(self,node_seq,new_node):
        s = ""
        assert new_node not in self.new_to_old.keys()
        self.new_to_old[new_node] = {}
        counter = 0
        for node in node_seq:
            #s+=self.new_graph.nodes[node]
            #self.old_to_new[node] = new_node

            # if part of the "old node seq" is already part of the merged nodes
            # need to update new_to_old and old_to_new dict for iterative merging
            if node in self.new_to_old:
                #pdb.set_trace()
                for v in self.new_to_old[node]:
                    self.old_to_new[v] = new_node
                    self.new_to_old[new_node][v] = [counter,len(s)+self.new_to_old[node][v][1]]
                    counter+=1
                del self.new_to_old[node]
            else:
                self.old_to_new[node] = new_node
                self.new_to_old[new_node][node] = [counter,len(s)]
                counter+=1
            s+=self.new_graph.nodes[node]
        return s

    def add_merge_node(self,node_seq):

        nodeid = self.new_graph.node_id_max
        seq = self.concat_label(node_seq, nodeid)
        _ = self.new_graph.add_node(seq, nodeid)
        #pdb.set_trace()
        # nodeid = self.node_id_max
        # while nodeid in self.nodes:
        #     nodeid+=1
        #self.nodes[nodeid] = seq
        # self.update_max_node_id
        # seq = self.concat_label(node_seq,nodeid)
        # self.nodes[nodeid] = seq
        return nodeid

# ============================================================================================
#
#                            COMPACT MODE = REMOVE SMALL VARIATIONS
#
# ============================================================================================
    
    def get_compact_graph_small_variations(self, threshold):

        self.merge_small_vars(threshold)
    
    def get_compact_graph_small_variations_2(self, threshold):

        self.merge_small_vars_2(threshold)
    
    def merge_small_vars_2(self, threshold=10):

        #self.new_graph = deepcopy(self.old_graph)

        #top_order, top_dict = self.old_graph.topological_sort()
        all_bubbles, bubble_hierarchy, bubble_layers = self.old_graph.find_all_bubbles()
        print("Number of bubbles:", len(all_bubbles), "Layers:", len(bubble_layers))

        self.old_to_new = defaultdict()
        self.new_to_old = defaultdict(list)

        processed = set()

        total_bubbles = len(all_bubbles)

        removed_chars = 0

        self.new_graph = deepcopy(self.old_graph)

        for layer_id, layer in enumerate(bubble_layers[::-1]):
            total_bubbles = len(layer)
            counter = 0
            for bubble_start in layer:
                #pdb.set_trace()
                bubble = all_bubbles[bubble_start]
                longest_path = self.new_graph.longest_path(bubble[0], bubble[-1])
                longest_path_index = {longest_path[i]:i for i in range(len(longest_path))}
                longest_path_edges = set([(longest_path[i], longest_path[i+1]) for i in range(len(longest_path)-1)])
                assert longest_path[0]==bubble[0]
                assert longest_path[-1]==bubble[-1]

                '''
                self.old_to_new[bubble[0]] = bubble[0]
                self.new_to_old[bubble[0]].append(bubble[0])
                self.old_to_new[bubble[-1]] = bubble[-1]
                self.new_to_old[bubble[-1]].append(bubble[-1])
                '''

                for i,node in enumerate(bubble):
                    if node in self.new_graph.nodes:
                        childs = list(self.new_graph.adj_list[node].keys())
                        if node in longest_path_index:
                            for child in childs:
                                if (node,child) not in longest_path_edges and child in longest_path_index:
                                    index_1 = longest_path_index[node]
                                    index_2 = longest_path_index[child]
                                    assert index_1<index_2-1
                                    char_inbetween = sum([len(self.new_graph.nodes[longest_path[i]]) for i in range(index_1+1,index_2)])
                                    if char_inbetween <= threshold:
                                        self.new_graph.remove_edge(node,child)
                        else:
                            parents = list(self.new_graph.adj_list_reverse[node].keys())
                            if len(parents)==1 and len(childs)==1 and parents[0] in longest_path_index and childs[0] in longest_path_index:
                                index_1 = longest_path_index[parents[0]]
                                index_2 = longest_path_index[childs[0]]
                                assert index_1<index_2
                                char_inbetween = sum([len(self.new_graph.nodes[longest_path[i]]) for i in range(index_1+1,index_2)])
                                if char_inbetween <= threshold and len(self.new_graph.nodes[node])<=threshold:
                                    removed_chars += len(self.new_graph.nodes[node])
                                    self.old_to_new[node] = childs[0]
                                    self.new_to_old[childs[0]].append(node)
                                    self.new_graph.remove_edge(parents[0],node)
                                    self.new_graph.remove_edge(node,childs[0])
                                    self.new_graph.remove_node(node)
                                    continue
                        if node not in self.old_to_new:
                            self.old_to_new[node] = node
                            self.new_to_old[node].append(node)
                    else:
                        assert node in self.old_to_new
                counter+=1
            print("processed %i of %i bubbles in layer %i" %(counter, total_bubbles, len(bubble_layers)-layer_id-1) )
        
        # remove invalid edges
        new_adj_list = defaultdict(dict)
        for n in self.new_graph.adj_list:
            if n not in self.new_graph.nodes: 
                continue
            for m in self.new_graph.adj_list[n]:
                if m not in self.new_graph.nodes: 
                    continue
                new_adj_list[n][m] = 1

        self.new_graph.adj_list = new_adj_list
        _ = self.new_graph.get_degrees()
        self.new_graph.get_adj_reverse()
        #print(removed_chars)

        # update new graph topological sort
        self.new_graph.topological_sort()
        
        return self.old_to_new




    def merge_small_vars(self, threshold=10):
        
        self.new_graph = deepcopy(self.old_graph)

        top_order, top_dict = self.old_graph.topological_sort()
        all_bubbles, bubble_hierarchy, bubble_layers = self.old_graph.find_all_bubbles()
        print("Number of bubbles:", len(all_bubbles), "Layers:", len(bubble_layers))

        self.old_to_new = defaultdict()
        self.new_to_old = defaultdict(list)

        processed = set()

        total_bubbles = len(all_bubbles)

        removed_chars = 0

        for layer_id, layer in enumerate(bubble_layers[::-1]):
            total_bubbles = len(layer)
            counter = 0
            #pdb.set_trace()
            for bubble_start in layer:
                bubble = all_bubbles[bubble_start]
                
                total_char = sum([len(self.old_graph.nodes[bubble[i]]) for i in range(1,len(bubble)-1)])

                # simplest case applies only to leaf 
                if total_char <= threshold :
                    processed.add(bubble_start)
                    new_adj_list = defaultdict(dict)
                    retained_nodes = self.old_graph.longest_path(bubble[0], bubble[-1])
                    retained_nodes_set = set(retained_nodes)

                    for i,node in enumerate(retained_nodes):
                        self.old_to_new[node] = node
                        self.new_to_old[node].append(node)

                    for i in range(len(retained_nodes) - 1):
                        new_adj_list[retained_nodes[i]][retained_nodes[i+1]] = 1
                    
                    retained_edges = set([(retained_nodes[i], retained_nodes[i+1]) for i in range(len(retained_nodes)-1)])
                    
                    for i,n in enumerate(bubble[:-1]):
                        neighbors = list(self.new_graph.adj_list[n].keys())
                        for m in neighbors:
                            if (n,m) not in retained_edges and n in self.new_graph.nodes:
                                self.new_graph.remove_edge(n,m)
                                if n != bubble_start and n!= bubble[-1] and n not in retained_nodes_set:
                                    self.new_graph.remove_node(n)
                                    removed_chars += len(self.old_graph.nodes[n])
                                self.old_to_new[n] = bubble[-1]
                                self.new_to_old[bubble[-1]].append(n)
                    counter += 1
                    continue

                #pdb.set_trace()
                curr_node = bubble[0]
                curr_idx = 0
                while curr_node != bubble[-1]:
                    if curr_node not in self.new_graph.nodes: 
                        curr_idx += 1
                        curr_node = bubble[curr_idx]
                        continue

                    # always keep the source of a bubble
                    self.old_to_new[curr_node] = curr_node
                    self.new_to_old[curr_node].append(curr_node)

                    if len(self.new_graph.adj_list[curr_node]) <= 1: 
                        curr_idx += 1
                        curr_node = bubble[curr_idx]
                        continue
                    
                    ## TODO: keep working on one node until it cannot be further simplified

                    indel_thresh = False
                    to_retain = []
                    neighbor = -1

                    for neighbor in self.new_graph.adj_list[curr_node]:
                        if neighbor in self.new_graph.nodes and len(self.new_graph.nodes[neighbor]) > threshold:
                            to_retain.append(neighbor)
                            indel_thresh = True
                    if len(to_retain) == 0:
                        to_retain.append(neighbor)
                    to_retain = set(to_retain)

                    neighbors = list(self.new_graph.adj_list[curr_node].keys())
                    
                    for n_idx, neighbor in enumerate(neighbors):
                        if neighbor not in to_retain and neighbor in self.new_graph.nodes and neighbor != bubble[-1]:

                            to_replace = curr_node
                            if len(to_retain) == 1:
                                to_replace = list(to_retain)[-1]
                            
                            self.new_graph.remove_edge(curr_node, neighbor)

                            # reorient edges from removed node
                            neighbors_children = list(self.new_graph.adj_list[neighbor].keys())
                            for child in neighbors_children:
                                if child in self.new_graph.nodes:
                                    self.new_graph.add_edge(curr_node, child)
                                    self.new_graph.remove_edge(neighbor, child)

                            # reorient edges to removed node
                            parents = list(self.new_graph.adj_list_reverse[neighbor].keys())
                            for parent in parents:
                                if self.old_graph.top_dict[parent] < self.old_graph.top_dict[curr_node]:
                                    if parent != curr_node:
                                        self.new_graph.add_edge(parent, curr_node)
                                    self.new_graph.remove_edge(parent, neighbor)
                                elif self.old_graph.top_dict[parent] > self.old_graph.top_dict[curr_node]:
                                    for child in neighbors_children:
                                        if parent != child:
                                            self.new_graph.add_edge(parent, child)
                                    self.new_graph.remove_edge(parent, neighbor)
                                else:
                                    continue

                            self.old_to_new[neighbor] = to_replace
                            self.new_to_old[to_replace].append(neighbor)
                            if neighbor in self.new_to_old:
                                for old_node in self.new_to_old[neighbor]:
                                    self.old_to_new[old_node] = to_replace
                                    self.new_to_old[to_replace].append(old_node)
                            self.new_graph.remove_node(neighbor)
                            removed_chars += len(self.old_graph.nodes[neighbor])
                        elif neighbor in to_retain or neighbor == bubble[-1]:
                            self.old_to_new[neighbor] = neighbor
                            self.new_to_old[neighbor].append(neighbor)
                        
                    curr_idx += 1
                    curr_node = bubble[curr_idx]

                counter += 1
            print("processed %i of %i bubbles in layer %i" %(counter, total_bubbles, len(bubble_layers)-layer_id-1) )

        # remove invalid edges
        new_adj_list = defaultdict(dict)
        for n in self.new_graph.adj_list:
            if n not in self.new_graph.nodes: 
                continue
            for m in self.new_graph.adj_list[n]:
                if m not in self.new_graph.nodes: 
                    continue
                new_adj_list[n][m] = 1

        self.new_graph.adj_list = new_adj_list
        _ = self.new_graph.get_degrees()
        self.new_graph.get_adj_reverse()
        print(removed_chars)
        
        return self.old_to_new

# ============================================================================================
#
#                                        SPLIT NODES 
#
# ============================================================================================

    def split_one_node(self, node_str):
        assert(self.node_size != -1)
        strings = []
        i = 0
        while i < len(node_str):
            strings.append(node_str[i:i + self.node_size])
            assert(len(node_str) >= len(node_str[i:i + self.node_size]))
            i += self.node_size
        return strings

    def split_nodes(self, node_size = 1000):
        '''
        Split nodes into sizes <= node_size
        - Do BFS and split when necessary
        - keep track of old_to_new mapping along the way
        '''
        self.old_to_new = {}
        self.new_to_old = {}
        self.old_graph.get_adj_reverse()
        self.node_size = node_size
        top_order, top_dict = self.old_graph.topological_sort()
        sources, sinks = self.old_graph.get_source_sink()
        source = sources[0]
        sink = sinks[0]
        queue = [source]
        prev_node = -1
        node_id = 0

        #visited = set()
        prev_nodes = []

        for curr_node in top_order:
            # curr_node = queue.pop(0)

            # if curr_node not in visited:
            assert curr_node not in self.old_to_new
            self.old_to_new[curr_node] = {}

            if curr_node in self.old_graph.adj_list_reverse:
                prev_nodes = list(self.old_graph.adj_list_reverse[curr_node].keys())
            # split node
            node_strings = self.split_one_node(self.old_graph.nodes[curr_node])
            added_nodes = []

            for i,s in enumerate(node_strings):
                assert node_id not in self.new_to_old
                self.new_graph.add_node(s, node_id)
                added_nodes.append(node_id)
                self.old_to_new[curr_node][i] = [node_id,i*self.node_size]
                self.new_to_old[node_id] = curr_node

                # add edge between nodes
                if prev_nodes != []:
                    for prev_node in prev_nodes:
                        self.new_graph.add_edge(self.old_to_new[prev_node][max(self.old_to_new[prev_node].keys())][0], node_id)
                    prev_nodes = []
                # add edge within nodes
                elif prev_node != -1:
                    self.new_graph.add_edge(prev_node, node_id)

                prev_node = node_id
                node_id += 1
            #self.old_to_new[curr_node] = added_nodes
            #visited.add(curr_node)
                
                # add neighbors to queue
                # for neighbor in self.old_graph.adj_list[curr_node]:
                    # queue.append(neighbor)

    def path_projection_split(self, alignment):
        """
        If aligned path has length one, 
            look at the offset to determine which one or two nodes a node is projected to
            node_idx = offset // self.node_size
        Otherwise, for the first node, always use self.old_to_new[node][-1]
                   for the last node, always use self.old_to_new[node][0]
                   for the other nodes, use all of self.old_to_new[node]

        Args:
            alignment --- Alignment object 
        Return:
            alignment --- changed alignment object
        """
        assert(isinstance(self.old_to_new[0], list))
        assert(self.node_size != -1)

        old_path = alignment.path
        new_path = []
        
        if len(old_path) == 1:
            node_offset = alignment.offset // self.node_size        
            node_start_offset = alignment.offset % self.node_size   # offset within the new node
            path_length = ceil((len(alignment.seq) + node_start_offset) / self.node_size)
            new_path = self.old_to_new[old_path[0]][node_offset:min(node_offset+path_length, len(self.old_to_new[old_path[0]]))]
            alignment.offset = node_start_offset
        
        else:
            new_path.append(self.old_to_new[old_path[0]][-1])
            alignment.offset = alignment.offset % self.node_size
            for idx in range(1, len(old_path) - 1):
                for new_node in self.old_to_new[old_path[idx]]:
                    new_path.append(new_node)
            new_path.append(self.old_to_new[old_path[-1]][0])
        
        alignment.path = new_path
        return alignment

def summarize_old_to_new_mappings_pair(mapping1, mapping2):
    """
    summarizes old_to_new mappings for a pair of mappings in the order of condensing

    Args:
        mapping1 -- old_to_new mapping of a former step
        mapping2 -- old_to_new mapping of a later step
    """
    summarized_old_to_new = defaultdict()
    summarized_new_to_old = defaultdict(list)

    for old_node in mapping1:
        try:
            summarized_old_to_new[old_node] = mapping2[mapping1[old_node]]
            summarized_new_to_old[mapping2[mapping1[old_node]]].append(old_node)
        except KeyError:
            print("KeyError: mapping1 node not found in mapping2 ", old_node, mapping1[old_node])

    for old_node in mapping2:
        to_map = [old_node]
        if old_node in summarized_old_to_new:
            assert mapping2[old_node] == summarized_old_to_new[old_node] , " ".join(["mismatched mapping", str(old_node), str(mapping2[old_node]), str(summarized_old_to_new[old_node])])
        else:
            summarized_old_to_new[old_node] = mapping2[old_node]
            summarized_new_to_old[mapping2[old_node]].append(old_node)

    return summarized_old_to_new, summarized_new_to_old
    
