from basic_graph import *
from compact_graph import *
# import matplotlib.pyplot as plt
from gam_utils import *
import sys
import pdb
from collections import defaultdict, namedtuple
import argparse


parser = argparse.ArgumentParser('Graph Pruning')
parser.add_argument('--graph', type=str)
parser.add_argument('--f1', type=str)
parser.add_argument('--f2',type=str)
parser.add_argument('--prefix',type=str)
parser.add_argument('--chr',type=str)
args = parser.parse_args()

vg_file = args.graph
chromosome = args.chr

print("=== Reading Graph")
print(vg_file)
start = time()
graph_vg = Graph(fname=vg_file,chromosome=chromosome)
assert set(graph_vg.adj_list.keys()).issubset(set(graph_vg.nodes.keys()))
assert set(graph_vg.adj_list_reverse.keys()).issubset(set(graph_vg.nodes.keys()))
print("number of nodes in original graph:", len(graph_vg.nodes))
print("time reading graph:", time() - start)
print()

graph_vg.vg_graph_add_source_sink(print_msg=True)
#graph_vg.to_gfa(out_dir + "chr22.gfa")
graph_vg.print_graph_stats("vg")
print("total characters", graph_vg.get_total_characters())

print("=== Node Merging")
start = time()
compact_graph_branching_1 = CompactGraph(graph_vg)
compact_graph_branching_1.get_compact_graph(mode = "branching")
print("number of nodes after node merging:",len(compact_graph_branching_1.new_graph.nodes))
print("total characters", compact_graph_branching_1.new_graph.get_total_characters())
print("time node merging:", time() - start)
print()


print("=== Reading alignments")

assert len(compact_graph_branching_1.old_to_new) == len(graph_vg.nodes), "size mismatch"

compact_graph_branching_1.new_graph.initialize_nodes_aln_states()

gzip_write_1 = gzip.open("paired_"+args.f1,"wt")
gzip_write_2 = gzip.open("paired_"+args.f2,"wt")
gaf_file1 = gzip.open(args.f1,"rt")
gaf_file2 = gzip.open(args.f2,"rt")

line1 = gaf_file1.readline()
index1 = int(line1.split("\t")[0].split(".")[2])
line2 = gaf_file2.readline()
index2 = int(line2.split("\t")[0].split(".")[2])

counter=0
#counter1=0
#counter2=0
while line1!='' and line2!='':
    while index1<index2:
        '''
        counter1+=1
        if counter1%100000==0:
            print("counter1:",counter1)
        '''
        line1 = gaf_file1.readline()
        if line1=='':
            break
        index1 = int(line1.split("\t")[0].split(".")[2])
    while index1>index2 and line1!='':
        '''
        counter2+=1
        if counter2%100000==0:
            print("counter2:",counter2)
        '''
        line2 = gaf_file2.readline()
        if line2=='':
            break
        index2 = int(line2.split("\t")[0].split(".")[2])
    if index1==index2 and line1!='' and line2!='':
        counter+=1
        if counter%100000==0:
            print("Paired reads:",counter)
        aln1 = get_aln_obj_gaf(line1)
        aln2 = get_aln_obj_gaf(line2)
        #print(aln1.name,aln2.name)
        assert aln1.name==aln2.name
        new_path1 = set()
        new_path2 = set()
        assert len(aln1.path)>0
        assert len(aln2.path)>0
        assert aln1.path[0]<graph_vg.node_id_max
        assert aln2.path[0]<graph_vg.node_id_max
        _ = gzip_write_1.write(line1)
        _ = gzip_write_2.write(line2)
        for pos1 in aln1.path:
            #if pos1 < graph_vg.node_id_max:
                #marker_1=1
            new_path1.add(compact_graph_branching_1.old_to_new[pos1])
        for pos2 in aln2.path:
            #if pos2 < graph_vg.node_id_max:
                #marker_2=1
            new_path2.add(compact_graph_branching_1.old_to_new[pos2])
        #assert marker_1==marker_2
        _ = compact_graph_branching_1.new_graph.update_nodes_aln_states(new_path1)
        _ = compact_graph_branching_1.new_graph.update_nodes_aln_states(new_path2)
        line1 = gaf_file1.readline()
        if line1!='':
            index1 = int(line1.split("\t")[0].split(".")[2])
        line2 = gaf_file2.readline()
        if line2!='':
            index2 = int(line2.split("\t")[0].split(".")[2])

gzip_write_1.close()
gzip_write_2.close()
gaf_file1.close()
gaf_file2.close()

compact_graph_branching_1.new_graph.update_empty_nodes()

print("=== Removing empty nodes")
start = time()
compact_graph_nonempty = CompactGraph(deepcopy(compact_graph_branching_1.new_graph))
compact_graph_nonempty.remove_empty_nodes_3()
print("number of nodes after empty_node_removal:", len(compact_graph_nonempty.new_graph.nodes))
print("time empty_node_removal:", time() - start)
print()

print("=== Node Merging again")
compact_graph_branching_2 = CompactGraph(deepcopy(compact_graph_nonempty.new_graph))
compact_graph_branching_2.get_compact_graph(mode="branching")
print("number of nodes after node merging:",len(compact_graph_branching_2.new_graph.nodes))
print("total characters", compact_graph_branching_2.new_graph.get_total_characters())
print("time node merging:", time() - start)
print()

print("=== Removing small bubbles")
compact_graph_bubble = CompactGraph(deepcopy(compact_graph_branching_2.new_graph))
compact_graph_bubble.get_compact_graph_small_variations_2(threshold = 10)
print("number of nodes after small_sv_removal:", len(compact_graph_bubble.new_graph.nodes))
compact_graph_bubble.new_graph.print_graph_stats("small_variations")
print("time small_sv_removal:", time() - start)
print()

print("=== Node Merging again")
compact_graph_branching_3 = CompactGraph(deepcopy(compact_graph_bubble.new_graph))
compact_graph_branching_3.get_compact_graph(mode="branching")
#snp_merged_nonbranching_graph = compact_graph_branching_2.new_graph
print("number of nodes after node merging:",len(compact_graph_branching_3.new_graph.nodes))
print("total characters", compact_graph_branching_3.new_graph.get_total_characters())
print("time node merging:", time() - start)
print()

print("=== Node binning")
graph_split = CompactGraph(deepcopy(compact_graph_branching_3.new_graph))
graph_split.split_nodes(node_size = 10000)
graph_split.new_graph.print_graph_stats("graph_split")
graph_split.new_graph.to_gfa(args.prefix+"_"+chromosome+"_final_graph.gfa")


print("=== Generating graph-based contact matrix")

compact_graphs = [compact_graph_branching_1,compact_graph_nonempty,compact_graph_branching_2,compact_graph_bubble,compact_graph_branching_3,graph_split]
compact_modules = ["branching","empty_nodes","branching","small_variations","branching","split"]


num_nodes = len(graph_split.new_graph.nodes)
#print(num_nodes)
M = np.zeros((num_nodes, num_nodes))
inter_node_pairs = 0
total_pair = 0

final_out = open(args.prefix+"_"+chromosome+"_pair_out.txt", "w")

with gzip.open("paired_"+args.f1,"rt") as gaf_file1, gzip.open("paired_"+args.f2,"rt") as gaf_file2:
    for line1, line2 in tqdm(zip(gaf_file1,gaf_file2)):
        aln1 = get_aln_obj_gaf(line1)
        aln2 = get_aln_obj_gaf(line2)
        assert aln1.name==aln2.name
        assert len(aln1.path)>0
        assert len(aln2.path)>0
        assert aln1.path[0]<graph_vg.node_id_max
        assert aln2.path[0]<graph_vg.node_id_max
        total_pair +=1
        start1,end1 = get_final_location(aln1,compact_graphs,compact_modules)
        start2,end2 = get_final_location(aln2,compact_graphs,compact_modules)
        _ = final_out.write(aln1.name + "\t" + str(start1) + "\t" + str(end1) + "\t" + str(start2) + "\t" + str(end2) + "\n")
        M[start1,start2]+=1
        if start1!=start2:
            M[start2,start1]+=1
        if start1!=end1 or start2!=end2:
            inter_node_pairs+=1

final_out.close()
print("Inter-node pairs", inter_node_pairs)

np.save(args.prefix+"_"+chromosome+"_M.npy",M)