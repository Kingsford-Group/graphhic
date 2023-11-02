from graphtadfinder_utils_parallel import *
import scipy
import scipy.interpolate
import Bio
from Bio import SeqIO
import pdb
import argparse

parser = argparse.ArgumentParser('Heuristic DP')
parser.add_argument('--reference', type=str)
#parser.add_argument('--matrix', type=str)
parser.add_argument('--prefix',type=str)
parser.add_argument('--chromosomes',type=str)
#parser.add_argument('--threads', type=int,default=10)
args = parser.parse_args()

#pdb.set_trace()

chromosomes = args.chromosomes.split(",")
output_handle = open(args.prefix+"_reconstruction.fa","w")

with open(args.reference)as handle:
    for record in SeqIO.parse(handle, "fasta"):
        #print(record.id)
        if record.id in chromosomes:
            chromosome = record.id
            print("reconstructed:"+chromosome)
            reconstruct_seq = ""
            graph_path = args.prefix+"_"+chromosome+"_final_graph.gfa"
            graph = Graph(fname=graph_path)
            matrix_path = args.prefix+"_"+chromosome+"_M.npy"
            M = np.load(matrix_path)
            M_sum = (M.sum() + np.diag(M).sum())/2
            farmatusMU_gamma1 = np.load("./mu0_info/GM12878_"+chromosome+"_MAPQ0_g1.0_10kb.npy")
            M_healthy_sum = float(np.load("./mu0_info/GM12878_"+chromosome+"_MAPQ0_totalinteration.npy"))
            alpha = 1.0
            mu_gamma1 = scipy.interpolate.InterpolatedUnivariateSpline(farmatusMU_gamma1[:,0], farmatusMU_gamma1[:,1]*(M_sum/M_healthy_sum)**alpha)
            tad_score_config = {
                "base_reso":10000,
                "new_reso":10000,
                "mu":mu_gamma1,
                "base_gamma":1,
                "new_gamma":0.1,
                "offdiag_bins":300
            }
            nx_di_graph = graph.get_nx_graph()
            node_weight_dict = dict()
            for v in nx_di_graph.nodes():
                node_weight_dict[v] = len(graph.nodes[v])
            nx.set_node_attributes(nx_di_graph, node_weight_dict, "weight")
            V_topo_sort = list(nx.topological_sort(nx_di_graph))
            Boundaries_raw = np.load(args.prefix+"_"+chromosome+"_boundaries_raw.npy")
            exadj_mat = extended_adj_mat(nx_di_graph,V_topo_sort)
            genome_reconstructed = genome_reconstruction_2(Boundaries_raw,V_topo_sort[-1],V_topo_sort,M,nx_di_graph,exadj_mat,tad_score_config)
            assert nx.is_path(nx_di_graph,genome_reconstructed)
            np.save(args.prefix+"_"+chromosome+"_reconstructed_nodes.npy",np.array(genome_reconstructed))
            for node in genome_reconstructed:
                reconstruct_seq +=graph.nodes[node]
            if reconstruct_seq[0]=="$":
                reconstruct_seq = reconstruct_seq[1:]
            if reconstruct_seq[-1]=="#":
                reconstruct_seq = reconstruct_seq[:-1]
            record.seq = Bio.Seq.Seq(reconstruct_seq)
            SeqIO.write(record,output_handle,"fasta")
        else:
            SeqIO.write(record,output_handle,"fasta")

output_handle.close()
