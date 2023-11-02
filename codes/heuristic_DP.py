from graphtadfinder_utils_parallel import *
import scipy
import scipy.interpolate
import argparse

parser = argparse.ArgumentParser('Heuristic DP')
parser.add_argument('--graph', type=str)
parser.add_argument('--matrix', type=str)
parser.add_argument('--prefix',type=str)
parser.add_argument('--chr',type=str)
parser.add_argument('--threads', type=int,default=10)
args = parser.parse_args()

chromosome = args.chr
graph_path = args.graph
graph = Graph(fname=graph_path)

matrix_path = args.matrix
M = np.load(matrix_path)
#print(M.shape)
M_sum = (M.sum() + np.diag(M).sum())/2


farmatusMU_gamma1 = np.load("./mu0_info/GM12878_"+chromosome+"_MAPQ0_g1.0_10kb.npy")
#farmatusMU_gamma1[:,1] = farmatusMU_gamma1[:,0]*farmatusMU_gamma1[:,1]
#farmatusMU_gamma1[:,0] = farmatusMU_gamma1[:,0]+1
M_healthy_sum = float(np.load("./mu0_info/GM12878_"+chromosome+"_MAPQ0_totalinteration.npy"))
### estimate mu function
alpha=1.0
mu_gamma1 = scipy.interpolate.InterpolatedUnivariateSpline(farmatusMU_gamma1[:,0], farmatusMU_gamma1[:,1]*(M_sum/M_healthy_sum)**alpha)

tad_score_config = {
    "base_reso":10000,
    "new_reso":10000,
    "mu":mu_gamma1,
    "base_gamma":1,
    "new_gamma":0.1,
    "offdiag_bins":300
}

parallel_config = {"threads":args.threads,"num_threshold":200}

nx_di_graph = graph.get_nx_graph()
node_weight_dict = dict()
for v in nx_di_graph.nodes():
    node_weight_dict[v] = len(graph.nodes[v])
nx.set_node_attributes(nx_di_graph, node_weight_dict, "weight")

TAD_boundaries, OPT, Boundaries_raw = DPTAD_4_parallel(nx_di_graph,M,tad_score_config,parallel_config)

np.save(args.prefix+"_"+chromosome+"_graph_TAD_boundaries.npy",np.array(TAD_boundaries))
np.save(args.prefix+"_"+chromosome+"_OPT.npy",OPT)
#np.save("./results/round_2/chr22_693_smallvg_greedy2_OPT_D.npy",OPT_D)
np.save(args.prefix+"_"+chromosome+"_boundaries_raw.npy",Boundaries_raw)