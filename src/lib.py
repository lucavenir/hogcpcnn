import argparse
import numpy as np
import sys

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

    n = len(adj)
    w_mat = np.zeros((n,n), dtype=np.uint16)
    # triangles = set()

    sorted_degrees = sorted(
        [(i,sum(row)) for i,row in enumerate(adj)],
        reverse=True,
        key=lambda t:t[1]
    )

    for i, row in enumerate(adj[:-2,:]):
        # Iterating through all the rows, i.e. every node once exactly
        v = sorted_degrees[i][0]  # Iterating from the highest degree to lowest
        colouring = [  # Colouring every neighbour of v
            True if w!=v and elem!=0 else False
            for w, elem in enumerate(adj[v,:])
        ]
        for j, el in enumerate(row[v+1:]):  # for each coloured neighbour u
            if el!=0:  # the "coloured" condition is here
                u = v+j+1  # We state who "u" is
                for w, elem in enumerate(adj[u,:]):
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
        adj[v,:] = 0
        adj[:,v] = 0

    return w_mat

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
    # NOTE: probabilmente non è neanche necessario. Sfrutteremo AWPPR per
    # ri-creare quello che fa HotNet2. Ci interessa infatti ottenere la matrice
    # E per ottenere il grafo H come descritto nel paper di Vandin et. al.
    pass

def diffusion_matrix(w, heat, alpha, epsilon, delta):
    n = len(w)

    # Creating F, first.
    f = np.zeros((n,n), dtype=np.float64)

    for u in range(n):
        p_u_tilde = awppr(w, u, alpha, epsilon)
        f[:,u] = p_u_tilde

    # Then, applying the diffusion process, thus obtaining E
    e = np.dot(f,np.diag(heat))
    # Now, obtaining the H matrix, via pruning.
    h = np.zeros((n,n), dtype=np.float64)
    h[e>=delta] = 1

    return h

# TODO: Debugging. Anche con i piccoli grafi. E con quelli random.
# TODO: Lanciare sul cluster con gli altri dataset.
# TODO: processo di diffusione non con AWPPR, ma con semplice inversione.
# TODO: chiedere al professore: W ora ha elementi >1, è pesata così. È un problema?
