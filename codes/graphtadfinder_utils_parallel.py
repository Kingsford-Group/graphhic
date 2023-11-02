import sys
import numpy as np
from basic_graph import *
from multiprocessing import Manager
import multiprocessing


def extended_adj_mat(graph,T):
    '''
    construct an extended adjacency matrix for graph
    graph: a networkx DiGraph instance
    returns a n*n matrix (n = number of nodes in graph) where 
    M[i, j] = 1 if node i can reach node j; otherwise M[i, j] = 0
    '''
    #T = list(nx.topological_sort(graph))
    N = len(T)
    M = np.zeros((N, N), dtype=bool)

    for i in range(len(T)-1, -1, -1):
        v = T[i]
        M[v, v] = 1
        for u in graph._succ[v]:
            M[v, u] = 1
            # wherever M[u, w] = 1, M[v, w] = 1
            M[v, np.where(M[u]==1)[0]] = 1
    return M


def boundaries_update(boundaries,end_node):
    TADs = []
    l = end_node
    while l!=-1:
        k = boundaries[l][1]
        if boundaries[l][2]==1:
            TADs.append([k,l])
        l = boundaries[l][0]
    return TADs

def compute_score_2(path,path_length,M,tad_score_config):
    #pdb.set_trace()
    bin_length = float(path_length/tad_score_config["base_reso"])
    if bin_length<3:
        mu = np.inf
    else:
        mu = tad_score_config["mu"](bin_length)*bin_length**(tad_score_config["base_gamma"]-tad_score_config["new_gamma"])
    M_sub = M[path][:,path]
    interations = (M_sub.sum() + np.diag(M_sub).sum())/2
    q_score = interations/(bin_length**tad_score_config["new_gamma"]) - mu
    return q_score

def genome_reconstruction_2(boundaries,end_node,T,M,graph,exadj_mat,tad_score_config):
    V_num = len(T)
    T_index = {T[j]:j for j in range(V_num)}
    node_l = end_node
    #pdb.set_trace()
    genome = []
    while node_l!=-1:
        if node_l==0:
            genome = [0]+genome
            break
        node_k = boundaries[node_l][1]
        index_k = T_index[node_k]
        index_l = T_index[node_l]
        #pdb.set_trace()
        _,kl_path = q_greedy_3(index1=index_k,index2=index_l,T=T,M=M,graph=graph,exadj_mat=exadj_mat,tad_score_config=tad_score_config)
        genome = kl_path + genome
        node_l = boundaries[node_l][0]
    return genome


def q_greedy_3(index1,index2,T,M,graph,exadj_mat,tad_score_config):
    def node_weighted_shortest_path(u,v,d):
        return graph.nodes[u]["weight"]
    #pdb.set_trace()
    T = np.array(T)
    nodes_inbetween = T[index1:index2+1]
    node_start = T[index1]
    node_end = T[index2]
    shortest_path = nx.dijkstra_path(graph, node_start, node_end, weight=node_weighted_shortest_path)
    shortest_path_length = sum([graph.nodes[node]["weight"] for node in shortest_path])
    shortest_path_score = compute_score_2(shortest_path,shortest_path_length,M,tad_score_config)
    best_path_save = {"score":shortest_path_score,"path_list":shortest_path}
    max_path_top = np.array([index1,index2])
    max_path = T[max_path_top]
    #best_path_save = {"score":-np.inf,"path_list":[]}
    max_path_length = graph.nodes[node_start]["weight"]+graph.nodes[node_end]["weight"]
    if node_start in graph._pred[node_end]:
        path_edges = 1
    else:
        path_edges = 0
    max_score = -np.inf
    sub_M = M[nodes_inbetween,:][:,nodes_inbetween]
    row_sum = np.sum(sub_M, axis=1)
    inds_row = np.argsort(row_sum)[::-1]
    ends_marker = 0
    for index in inds_row:
        node = nodes_inbetween[index]
        node_top = index+index1
        if node == node_start or node == node_end:
            ends_marker+=1
            continue
        insert_ind = np.searchsorted(max_path_top, node_top)
        pred_top_id = max_path_top[insert_ind - 1]
        succ_top_id = max_path_top[insert_ind]
        pred = T[pred_top_id]
        succ = T[succ_top_id]
        if exadj_mat[pred, node]==0 or exadj_mat[node, succ]==0:
            continue
        ### insert node
        path_edges = path_edges +int(pred in graph._pred[node])+int(node in graph._pred[succ])-int(pred in graph._pred[succ])
        max_path_top = np.insert(max_path_top, insert_ind, node_top)
        max_path = T[max_path_top]
        max_path_length = max_path_length + graph.nodes[node]["weight"]
        if path_edges == len(max_path)-1 and ends_marker==2:
            break
    assert nx.is_path(graph,max_path)==1
    max_score = compute_score_2(max_path,max_path_length,M,tad_score_config)
    if max_score > best_path_save["score"]:
        best_path_save = {"score":max_score,"path_list":list(max_path)}
    return best_path_save["score"],best_path_save["path_list"]

def DPTAD_4_chunk(prev_nodes,index2,tad_score_config):
    #print("yes1")
    #t0 = time()
    node_l = V_topo_sort[index2]
    OPT_tmp = -np.inf
    Boundaries_tmp = -np.ones(3, dtype=int)
    '''
    chunk_size = len(prev_nodes)//parallel_config["threads"]
    if id_parallel==parallel_config["threads"]-1:
        prev_nodes_chunk = prev_nodes[id_parallel*chunk_size:]
    else:
        prev_nodes_chunk = prev_nodes[id_parallel*chunk_size:(id_parallel+1)*chunk_size]
    '''
    #print("yes1")
    #t0 = time()
    for node_k in prev_nodes:
        j = V_topo_sort_index[node_k]
        if exadj_mat[node_k,node_l] == 0:
            continue
        qkl,_ = q_greedy_3(index1=j,index2=index2,T=V_topo_sort,M=M_1,graph=G_1,exadj_mat=exadj_mat,tad_score_config=tad_score_config)
        node_k_parents = list(G_1._pred[node_k])
        if len(node_k_parents) == 0:
            OPT_v = 0
            node_v = -1
        else:
            OPT_v = np.max(OPT[node_k_parents])
            node_v = node_k_parents[np.argmax(OPT[node_k_parents])]
        if OPT_tmp < OPT_v+qkl or Boundaries_tmp[2]==-1:
            OPT_tmp = OPT_v+qkl
            Boundaries_tmp[0] = node_v
            Boundaries_tmp[1] = node_k
            if qkl<0:
                Boundaries_tmp[2] = 0
            else:
                Boundaries_tmp[2] = 1
    #print("yes2")
    #print(time()-t0)
    return [OPT_tmp]+list(Boundaries_tmp)

def DPTAD_4_parallel(nx_di_graph,M,tad_score_config,parallel_config):
    #pdb.set_trace()
    global G_1
    G_1 = deepcopy(nx_di_graph)
    global M_1
    M_1 = deepcopy(M)
    global V_topo_sort
    V_topo_sort = list(nx.topological_sort(G_1))
    V_num = len(V_topo_sort)
    global V_topo_sort_index
    V_topo_sort_index = {V_topo_sort[j]:j for j in range(V_num)}

    global exadj_mat
    exadj_mat = extended_adj_mat(G_1,V_topo_sort)

    global OPT
    OPT = -np.ones(V_num)*np.inf
    OPT[0] = 0
    Boundaries = -np.ones((V_num, 3), dtype=int)

    for i in tqdm(range(1,V_num)):
        #pdb.set_trace()
        node_l = V_topo_sort[i]
        #OPT_tmp = -np.inf
        #OPT_D_tmp = -np.inf
        #node_v = -1
        #node_k = -1
        parent_one_layer = set(G_1._pred[node_l])
        parents_all = set(G_1._pred[node_l])
        #pdb.set_trace()
        '''
        if i>=1629:
            pdb.set_trace()
        '''
        for _ in range(tad_score_config["offdiag_bins"]):
            new_parent_one_layer = set()
            if len(parent_one_layer)==0:
                break
            for node_k in parent_one_layer:
                new_parent_one_layer = new_parent_one_layer.union(set(G_1._pred[node_k]))
            parent_one_layer = new_parent_one_layer.difference(parents_all)
            parents_all = parents_all.union(parent_one_layer)
        parents_all = list(parents_all)
        #pdb.set_trace()
        if len(parents_all)>=parallel_config["num_threshold"]:
            ### do parallelism
            #pdb.set_trace()
            chunk_size = len(parents_all)//parallel_config["threads"]
            parallel_args = [(parents_all[id_parallel*chunk_size:(id_parallel+1)*chunk_size],i,tad_score_config) for id_parallel in range(parallel_config["threads"]-1)]
            parallel_args.append((parents_all[(parallel_config["threads"]-1)*chunk_size:],i,tad_score_config))
            with multiprocessing.Pool(parallel_config["threads"]) as pool:
                chunk_OPT = np.array(pool.starmap(DPTAD_4_chunk,parallel_args))
        else:
            chunk_OPT = np.array(DPTAD_4_chunk(parents_all,i,tad_score_config))
        if chunk_OPT.ndim==1:
            OPT[node_l] = chunk_OPT[0]
            Boundaries[node_l][0] = chunk_OPT[1]
            Boundaries[node_l][1] = chunk_OPT[2]
            Boundaries[node_l][2] = chunk_OPT[3]
        else:
            max_index = chunk_OPT[:,0].argmax()
            OPT[node_l] = chunk_OPT[max_index][0]
            Boundaries[node_l][0] = chunk_OPT[max_index][1]
            Boundaries[node_l][1] = chunk_OPT[max_index][2]
            Boundaries[node_l][2] = chunk_OPT[max_index][3]
        print(len(parents_all))
        print(Boundaries[node_l],OPT[node_l])
        '''
        if Boundaries[node_l][2]==-1:
            pdb.set_trace()
        '''
    TAD_list = boundaries_update(Boundaries,V_topo_sort[-1])
    return TAD_list, OPT, Boundaries
