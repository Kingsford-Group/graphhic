"""
    Utilities to view and manipulate alignment-to-graph records
"""

from collections import defaultdict, Counter, namedtuple
import matplotlib.pyplot as plt
from basic_graph import *
import re
import gzip
import pdb

class Alignment:
    """
    A simplified alignment class

    self.name --- read identifier
    self.path -- aligned path, a sequence of nodes
    self.offset -- start location on the first aligned node in terms of bases
    self.length -- read length
    self.isReversed -- whether the alignment orientation is reversed or not
    self.seq -- the sequence of the alignment
    """

    def __init__(self, name="", path=[], offset=0, length=0, isReversed=0, seq=""):
        self.name = name
        self.path = path
        self.offset = offset
        self.length = length
        self.isReversed = isReversed
        self.seq = seq

class Alignment_GAF:
    """
    A simplified alignment class from .gaf files

    self.name --- read identifier
    self.path -- aligned path, a sequence of nodes
    self.offset -- start location on the first aligned node in terms of bases
    self.offset_end -- end location on the aligned path in terms of bases
    self.length -- read length
    self.isReversed -- whether the alignment orientation is reversed or not, if it is 2, the read contains an inversion
    self.seq -- the sequence of the alignment
    """

    def __init__(self, name="", path=[], offset=-1, offset_end=-1, length=0, isReversed=-1, seq=""):
        self.name = name
        self.path = path
        self.offset = offset
        self.offset_end = offset_end
        self.length = length
        self.isReversed = isReversed
        self.seq = seq
    def print(self):
        print("name:",self.name)
        print("path:",self.path)
        print("offset:",self.offset)
        print("offset_end:",self.offset_end)
        print("length:",self.length)
        print("isReversed:",self.isReversed)

def get_aln_obj_gaf(record) -> Alignment_GAF:
    record = record.split('\t')
    name = record[0]
    length = int(record[1])
    if record[5]=="*":
        aln_obj = Alignment_GAF(name=name, length=length)
        return aln_obj
    path_info = list(filter(None, re.split('(>)|(<)', record[5])))
    path = []
    orient_list = []
    for n in path_info:
        if n in [">","<"]:
            orient_list.append(n)
        else:
            path.append(int(n))
    if orient_list.count('>') !=0 and orient_list.count('<') != 0:
        isReversed = 2
    elif orient_list.count('>') !=0:
        isReversed = 0
    else:
        isReversed = 1
    offset = int(record[7])
    offset_end = int(record[8])
    aln_obj = Alignment_GAF(name=name, 
                        path=path, 
                        offset=offset,
                        offset_end=offset_end,
                        length=length, 
                        isReversed=isReversed)
    return aln_obj


def write_aln_dict_pkl(aln_dict, fname):
    with open(fname, "wb") as ofile:
        pickle.dump(aln_dict, ofile)

def get_aln_obj_gam(record) -> Alignment:
    """
    Get Alignment object from one line of gam record
    """
    mapping = record.path.mapping
    path = []
    counter = 0
    reverse_true = 0
    first_offset = -1
    for pos in mapping:
        isReversed = int(pos.position.is_reverse)
        offset = pos.position.offset
        if counter == 0:
            reverse_true = isReversed
            first_offset = offset
        else:
            assert offset==0
            assert isReversed == reverse_true
        counter+=1
        path.append((pos.position.node_id, pos.rank))
    path = sorted(path, key=lambda x:x[-1])
    aln_obj = Alignment(name=record.name, 
                        path=[int(n[0]) for n in path], 
                        offset=int(first_offset), 
                        isReversed=reverse_true,
                        seq = record.sequence)
    return aln_obj

def get_aln_obj_list(aln_name, aln_list) -> Alignment:
    """
    Get Alignment object from a list
    """
    return Alignment(name=aln_name,
                     path=aln_list[0],
                     offset=aln_list[1],
                     isReversed=aln_list[2])

def read_aln_from_gam(gam_file, pair=1) -> tuple:
    """
    Read a dictionary of alignments from gam file

    Args:
        gam_file (str) --- path to gam file
        pair (int) -- rank in read pair. (1 if _1 in filename, 2 if _2 in filename)

    Returns:
        tuple(dict, int)
        alignments --- a dictionary of alignment objects. Key = read name, value = Alignment objects
        size --- number of records in gam file
    """
    cannot_read = []
    aln_dict = defaultdict()
    counter = 0
    with stream.open(gam_file,"rb",header=b'GAM') as istream:
        for data in tqdm(istream):
            aln = vg_pb2.Alignment()
            try:
                aln.ParseFromString(data)
                aln_obj = get_aln_obj_gam(aln)
                if len(aln_obj.path) == 0:
                    continue
                aln_dict[aln.name + "_" + str(pair)] = aln_obj
                counter += 1
            except Exception:
                if data != b'GAM':
                    exit()
                # print("cannot read")
                cannot_read.append(data)
    return aln_dict, counter
    
def read_aln_from_pkl(pkl_file) -> dict:
    """
    Read a dictionary of Alignments from pkl file
    This function assumes that the pkl file contains a dictionary of alignment objects
    """
    with open(pkl_file, "rb") as infile:
        aln_dict = pickle.load(infile)
    return aln_dict

def read_aln_from_pkl_prev(pkl_file) -> dict:
    """
    Read a dictionary of Alignments from pkl file
    This function assumes that the pkl file contains a dictionary of paths. Need to convert to Alignment objects
    """
    new_aln_dict = dict()
    with open(pkl_file, "rb") as infile:
        aln_dict = pickle.load(infile)
        counter = 0
        counter2 = 0

        for aln_name in aln_dict:
            aln_obj = get_aln_obj_list(aln_name, aln_dict[aln_name])
            if len(aln_obj.path) > 0:
                new_aln_dict[aln_name] = aln_obj
            if len(aln_obj.path) == 0 and counter < 10:
                print(aln_name, aln_obj.path, aln_obj.offset)
                counter += 1
            if len(aln_obj.path) > 10 and counter2 < 10:
                print(aln_name, aln_obj.path)
                counter2 += 1
    return new_aln_dict

def aln_path_length(aln_dict) -> tuple:

    length_list = []
    for name in aln_dict:
        length_list.append(len(aln_dict[name].path))
    
    length_counter = Counter(length_list)
    return length_list, length_counter

def print_counter(counter):

    keys = sorted(counter.keys())
    for key in keys:
        print(key, counter[key])

def check_order(path, top_dict): 
    """
    Check if the ordering of nodes on an alignment path follows topological order
    """ 
    prev_order = -1
    for n in path:
        curr_order = top_dict[n]
        if prev_order > curr_order:
            print("no!", n, path)
            assert(False)
            exit()

def aln_to_matrix(aln_dict, graph):
    """
        Convert alignments in a dictionary to an interaction matrix.
        If len(path) > 1:
            - for each adjacent pair of nodes, add 1 interaction to the pair
            - for node in path1:
                - for node in path2:
                    - M[node1][node2] += 1
    """
    top_order, top_dict = graph.topological_sort()
    for aln_name in aln_dict:
        check_order(aln_dict[aln_name].path, top_dict)

    ## select reads with both ends are length 1
    aln_names = [n for n in aln_dict if "_1" in n]
    aln_names2 = [n for n in aln_dict if "_2" in n]
    print(len(aln_names), len(aln_names2))

    num_nodes = len(graph.nodes)
    M = np.zeros((num_nodes, num_nodes))

    added = 0
    notfound = 0
    for name in tqdm(aln_names):
        try:
            aln1 = aln_dict[name]
            aln2 = aln_dict[name.split("_")[0] + "_2"]  # find mate pair

            if len(aln1.path) > 1:
                for idx in range(len(aln1.path)-1):
                    M[aln1.path[idx]][aln1.path[idx + 1]] += 1
                    M[aln1.path[idx + 1]][aln1.path[idx]] += 1

            if len(aln2.path) > 1:
                for idx in range(len(aln2.path)-1):
                    M[aln2.path[idx]][aln2.path[idx + 1]] += 1
                    M[aln2.path[idx + 1]][aln2.path[idx]] += 1

            for node1 in aln1.path:
                for node2 in aln2.path:
                    M[node1][node2] += 1
                    M[node2][node1] += 1

                added += 1
        except KeyError:
            notfound += 1

    print("added:", added, "pair of reads")
    print("not found:", notfound)
    return M

def get_final_location(record,compact_graphs,compact_modules):
    #pdb.set_trace()
    ### here forward_offset_start, forward_offset_end are the first and last base location indexed from zero
    assert len(compact_graphs)>=1
    assert len(compact_graphs)==len(compact_modules)
    start_graph = compact_graphs[0].old_graph
    start_node = record.path[0]
    end_node = record.path[-1]
    if record.isReversed == 0:
        forward_offset_start = record.offset
        forward_offset_end = record.offset_end - sum([len(start_graph.nodes[record.path[i]]) for i in range(len(record.path)-1)])
        assert forward_offset_end>=0
    elif record.isReversed == 1:
        forward_offset_start = len(start_graph.nodes[record.path[0]]) - record.offset -1
        forward_offset_end = sum([len(start_graph.nodes[record.path[i]]) for i in range(len(record.path))]) - record.offset_end - 1
        assert forward_offset_start>=0
        assert forward_offset_end>=0
    else:
        print("cannot handle this orientation!")
        exit(0)
    for i in range(len(compact_modules)):
        if compact_modules[i]=="branching":
            new_start_node = compact_graphs[i].old_to_new[start_node]
            new_end_node = compact_graphs[i].old_to_new[end_node]
            forward_offset_start += compact_graphs[i].new_to_old[new_start_node][start_node][1]
            forward_offset_end += compact_graphs[i].new_to_old[new_end_node][end_node][1]
            assert forward_offset_start < len(compact_graphs[i].new_graph.nodes[new_start_node])
            assert forward_offset_end < len(compact_graphs[i].new_graph.nodes[new_end_node])
            start_node = new_start_node
            end_node = new_end_node
        elif compact_modules[i]=="empty_nodes":
            if start_node not in compact_graphs[i].old_to_new:
                assert start_node in compact_graphs[i].new_graph.nodes
            else:
                new_start_node = compact_graphs[i].old_to_new[start_node]
                assert new_start_node in compact_graphs[i].new_graph.nodes
                forward_offset_start += compact_graphs[i].new_to_old[new_start_node][start_node][1]
                assert forward_offset_start < len(compact_graphs[i].new_graph.nodes[new_start_node])
                start_node = new_start_node
            if end_node not in compact_graphs[i].old_to_new:
                assert end_node in compact_graphs[i].new_graph.nodes
            else:
                new_end_node = compact_graphs[i].old_to_new[end_node]
                assert new_end_node in compact_graphs[i].new_graph.nodes
                forward_offset_end += compact_graphs[i].new_to_old[new_end_node][end_node][1]
                assert forward_offset_end < len(compact_graphs[i].new_graph.nodes[new_end_node])
                end_node = new_end_node
        elif compact_modules[i] == "small_variations":
            assert start_node in compact_graphs[i].old_to_new
            assert end_node in compact_graphs[i].old_to_new
        
            new_start_node = compact_graphs[i].old_to_new[start_node]
            new_end_node = compact_graphs[i].old_to_new[end_node]
            assert new_start_node in compact_graphs[i].new_graph.nodes
            assert new_end_node in compact_graphs[i].new_graph.nodes

            if start_node != new_start_node:
                forward_offset_start = 0
            if end_node != new_end_node:
                forward_offset_end = 0

            start_node = compact_graphs[i].old_to_new[start_node]
            end_node = compact_graphs[i].old_to_new[end_node]

        elif compact_modules[i]=="split":
            start_index = int(forward_offset_start/compact_graphs[i].node_size)
            new_start_node = compact_graphs[i].old_to_new[start_node][start_index][0]
            forward_offset_start -= compact_graphs[i].old_to_new[start_node][start_index][1]
            assert forward_offset_start>=0
            assert forward_offset_start < len(compact_graphs[i].new_graph.nodes[new_start_node])
            start_node = new_start_node

            end_index = int(forward_offset_end/compact_graphs[i].node_size)
            new_end_node = compact_graphs[i].old_to_new[end_node][end_index][0]
            forward_offset_end -= compact_graphs[i].old_to_new[end_node][end_index][1]
            assert forward_offset_end>=0
            assert forward_offset_end < len(compact_graphs[i].new_graph.nodes[new_end_node])
            end_node = new_end_node
    return start_node,end_node

def main():
    out_dir = "/mnt/disk26/user/yutongq/HiCGenomeGraph_data/graph_pruning/"

    orig_aln_dict = read_aln_from_pkl_prev(out_dir + "698_1_chr22_aln.pkl")
    print("Number of records (original):", len(orig_aln_dict))

    compact_aln_dict = read_aln_from_pkl_prev(out_dir + "698_1_chr22_aln_compact.pkl")
    print("Number of records (compact):", len(compact_aln_dict))

    new_compact_aln_dict =  read_aln_from_pkl_prev(out_dir + "698_1_chr22_aln_compact_new.pkl")
    print("Number of records (new compact):", len(new_compact_aln_dict))

    orig_length = aln_path_length(orig_aln_dict)
    compact_length = aln_path_length(compact_aln_dict)
    new_compact_length = aln_path_length(new_compact_aln_dict)

    print_counter(orig_length[-1])
    print_counter(compact_length[-1])
    print_counter(new_compact_length[-1])

    fig, ax = plt.subplots(1)
    ax.hist(orig_length[0], label="original")
    ax.hist(compact_length[0], alpha=0.5, label="condense_1")
    ax.hist(new_compact_length[0], alpha=0.5, label="condense_2")
    ax.legend()
    ax.set_xlabel("Path Length")
    ax.set_ylabel("Number of alignments")
    fig.savefig(out_dir + "path_length_distribution.png")

if __name__ == "__main__":
    main()