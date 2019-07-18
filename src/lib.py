import argparse
import numpy as np
import sys
import networkx as nx

# TODO: saranno funzioni davvero utili, queste?
def conductance(s, adj, m=None, m_instances=None):
    """
        py:function:: conductance(s, adj, m=None, m_instances=None)

        Computing the motif conductance of a subgraph s,
        with respect to a motif m, when given;
        when it's not, edge conductance is computed.

        :param s: a list, i.e. subgraph we're analyzing
        :param adj: numpy NxN matrix, i.e. the adjacency matrix of the graph
        :param m: (optional) a set, i.e. motif considered.
        :param m_instances: (optional) a list of sets, i.e. the instances of M inside the graph.
                            Great performance improvment when given.
    """

    if m!=None and m_instances==None:  # if not given, we should compute it.
        m_instances = motif_instances(adj, m)

    cut_size = len(cut(s, adj, m, m_instances))
    set_volume = min(vol(s, adj, m, m_instances), vol(complement_set(s), adj, m, m_instances))
    return cut_size /set_volume

def cut(s, adj, m=None, m_instances=None):
    """
        py:funcion:: cut(s,adj,m)

        Computing the motif cut of a subgraph,
        with respect to a motif m, when given;
        when it's not, edge cut is computed.

        :param s: a list or a set, i.e. the subgraph we're analyzing
        :param adj: numpy NxN matrix, i.e. the adjacency matrix of the graph
        :param m: a set, i.e. motif considered.

        :return: a set of edges, i.e. the desired cut
    """

    if m!=None and m_instances==None:  # if not given, we should compute it.
        m_instances = motif_instances(adj, m)

    s_set = s if type(s) is set else set(s)
    return ([
    	(i,j+i)
    	for i, row in enumerate(adj)  # Iterating through all the rows
    	for j, el in enumerate(row[i:])  # Bypassing the left-lower triangle bc of symmetry
    	if el!=0 and (i in s_set) and (j+i not in s_set)  # Returning edges on the cut
    ])

def vol(s, adj, m=None, m_instances=None):
    """
        py:funcion:: vol(s, adj, m)

        Computing the motif volume of a subgraph,
        with respect to a motif m, when given;
        when it's not, edge volume is computed.

        :param s: a list, i.e. the subgraph we're analyzing
        :param adj: numpy NxN matrix, i.e. the adjacency matrix of the graph
        :param m: a set, i.e. motif considered.

        :return: a set of edges, i.e. the desired cut
    """
    '''
        saving this code for later
    edges_in_s = [
    	(s[i], s[j+i])
    	for i, row in enumerate(adj[np.ix_(s,s)])
    	for j, el in enumerate(row[i:])
    	if el!=0
    ]
    '''

    if m!=None and m_instances==None:  # if not given, we should compute it.
        m_instances = motif_instances(adj, m)

def u_triangle(adj):
    '''
        py:funcion:: motif_count(adj, m)

        Computing the undirected triangle occurences inside the graph; this
        method is well explained by Chiba and Nishizeki (1985); but, instead of
        enumerating all the triangles, our objective is to create the new W
        matrix, with Wij = #triangles on the edge. In other words, we want to
        count the edges participating in occurences of m (undirected triangle).

        :param adj: numpy NxN matrix, i.e. the adjacency matrix of the graph

        :return: numpy NxN matrix, i.e. the new adj matrix representing the undirected triangle count
    '''
    a = np.copy(adj)
    n = len(a)
    w_mat = np.zeros((n,n), dtype=np.uint16)

    sorted_degrees = sorted(
        [(i,sum(row)) for i,row in enumerate(a)],
        reverse=True,
        key=lambda t:t[1]
    )

    for i, row in enumerate(a[:-2,:]):
        # Iterating through all the rows, i.e. every node once exactly
        v = sorted_degrees[i][0]  # Iterating from the highest degree to lowest
        colouring = [  # Colouring every neighbour of v
            True if w!=v and elem!=0 else False
            for w, elem in enumerate(a[v,:])
        ]
        for j, el in enumerate(row[v+1:]):  # for each coloured neighbour u
            if el!=0:  # the "coloured" condition is here
                u = v+j+1  # We state who "u" is
                for w, elem in enumerate(a[u,:]):
                    if elem!=0 and colouring[w]:
                        # If u and v share a neighbour, this is a triangle.
                        w_mat[u,v] += 1
                        w_mat[v,w] += 1
                        w_mat[w,u] += 1
                        # Because of undirected graphs symmetry
                        w_mat[v,u] += 1
                        w_mat[w,v] += 1
                        w_mat[u,w] += 1
                        # triangles.add((v,u,w))
                # Removing the colour from u
                colouring[u] = False
        # Erasing v from the graph
        a[v,:] = 0
        a[:,v] = 0

    return w_mat

def no(adj):
    return adj

def filter_disconnected_components(motif_mat, inputs):
    '''
        py:funcion:: filter_disconnected_components(motif_mat, inputs)

        This function filters out disconnected components from the graph.
        Disconnected components may appear in the hard version of the problem,
        when some edges don't occur in motifs, and therefore their weight is 0.
        This function returns a new filtered matrix and changes the input
        values, directly.


        :param adj: numpy NxN matrix, i.e. the motif matrix

        :return: numpy NxN matrix, i.e. the filtered motif matrix and other re-mapped inputs
    '''

    n = len(motif_mat)
    disconnected_vertices = [
        i if sum(motif_mat[i,:])==0 else 0
        for i in range(n)
    ]
    disconnected_vertices = [
        el for el in disconnected_vertices
        if el>0
    ]
    motif_mat_filtered = np.delete(
        np.delete(
            motif_mat,
            disconnected_vertices,
            0
        ),
        disconnected_vertices,
        1
    )
    k = len(motif_mat_filtered)
    print("Filtered out "+str(n - k)+" disconnected nodes")

    inputs['v_labels'] = np.delete(
        inputs['v_labels'],
        disconnected_vertices
    )
    inputs['heat'] = np.delete(
        inputs['heat'],
        disconnected_vertices
    )

    return motif_mat_filtered

def transition_matrix(motif_mat, soft=False, adj=None):
    n = len(motif_mat)

    w = np.zeros((n,n), dtype=np.float64)
    counting_motifs = motif_mat+adj if soft else motif_mat

    for i in range(n):
        w[i,:] = counting_motifs[i,:]/sum(counting_motifs[i,:])

    return w

def approx_check(r, d_w, epsilon):
    """
        py:function:: approx_check(r, d_w, epsilon)

        Given the two vectors needed for the 'push' operation, this function
        returns, if exists, the first possible vertex v (its number)
        that doesn't satisfy the approximation;
        insteda, if every vertex v satisfies such condition, '-1' is returned.

        :param r: vector r. See APPR paper.
        :param d_w: vector d_w. See APPR paper.
        :param epsilon: tolerance parameter (for the approximation).

        :return: a vertex that doesn't satisfy the approx., '-1' instead.
    """

    list_of_vertices = np.arange(len(r))
    np.random.shuffle(list_of_vertices)

    for v in list_of_vertices:
        if float(r[v])/float(d_w[v]) >= epsilon:
            return v

    return -1

def awppr(w, u, alpha, epsilon):
    """
        py:function:: awppr(adj, w, u, alpha, epsilon)

        This function returns the approximate weighted personalized
        page rank vector of node u, using the weighted matrix w, on the
        graph represented by the adjacency matrix adj.

        :param w: weighted adjacency matrix obtained with the motif analysis.
        :param u: vertex on which the AWPPR is computed.
        :param alpha: teleportation/random step parameter (see papers).
        :param epsilon: tolerance parameter (for the approximation).

        :return: the APPR vector on the weighted graph represented by W.
    """

    n = len(w)

    # Initializing
    p_tilde = np.zeros(n, dtype=np.float64)
    r = np.zeros(n, dtype=np.float64)
    r[u] = 1

    # np.set_printoptions(threshold=sys.maxsize)
    d_w = np.sum(w, axis=0)

    v = approx_check(r, d_w, epsilon)
    while (v != -1):  # repeat until the approx. condition is met
        # push operation
        ro = r[v] - (epsilon/2) * d_w[v]
        p_tilde[v] += (1-alpha)*ro
        r[v] = (epsilon/2) * d_w[v]
        # r update
        r += (w[v,:].astype(np.float64) / d_w.astype(np.float64)) * alpha * ro
        # should we iterate again?
        v = approx_check(r, d_w, epsilon)

    return p_tilde

def mappr():
    # NOTE: probabilmente non è necessario. Sfrutteremo AWPPR per
    # ri-creare quello che fa HotNet2. Ci interessa infatti ottenere la matrice
    # E per ottenere il grafo H come descritto nel paper di Vandin et. al.
    pass

def diffusion_matrix(w, heat, alpha, epsilon, delta, method='AWPPR'):
    n = len(w)

    # Creating F, first.
    f = np.zeros((n,n), dtype=np.float64)

    if method == 'AWPPR':  # approximate PPR method
        for u in range(n):
            p_u_tilde = awppr(w, u, alpha, epsilon)
            f[:,u] = p_u_tilde
    elif method == 'HOTNET2':  # inverting the matrix directly
        f = (1-alpha) * np.linalg.inv(np.eye(n) - alpha*w)
    else:
        raise ValueError(method + " is an unrecognized method. Aborting.")

    # Then, applying the diffusion process, thus obtaining E
    e = np.dot(f,np.diag(heat))
    # Now, obtaining the H matrix, via pruning.
    h = np.zeros((n,n), dtype=np.float64)
    h[e>=delta] = 1

    return h

def extract_strong_cc(h):
    H = nx.from_numpy_matrix(h, create_using=nx.DiGraph())
    ccs = [H.subgraph(c) for c in nx.strongly_connected_components(H)]
    subgraphs = [
        set(cc.nodes)
        for cc in ccs
        if len(set(cc.nodes))>1
    ]
    print(subgraphs)
    return subgraphs

"""
    TODO: chiedere al professore:
    W ora ha elementi >=1, è pesata così.
    HOTNET2, invece, usa una "transition matrix" che ha pesi <1.
    Come interpretare quella matrice di transizione, ora?

    In più, perché ignorare _completamente_ i lati che non partecipano ai motif?
    Perché semplicemente non lasciare peso unitario ai lati in generale,
    aumentando invece il peso dei lati che partecipano ai motif?
"""
