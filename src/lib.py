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

def sort_by_degree(adj):
    '''
        py:funcion:: sort_by_degree(adj)

        Computing a vector representing the non-increasing degree order of the vertices.
        This function has been wrote because we're looking for an by their degree
        (from the highest to the lowest).

        :param adj: numpy NxN matrix, i.e. the adjacency matrix of the graph

        :return: numpy 1xN vector, i.e. the non-increasing order by degree from the adj matrix
    '''

    n = len(adj)

    vertices = [i for i in range(n)]  # Yet to be ordered
    nodes_degree = [sum(row) for row in adj]  # Computing the nodes' degree

    ordered_vertices = sorted(
        vertices,
        key=lambda i:nodes_degree[i],
        reverse=True
    )

    # We return it as an Numpy array; we just love Numpy
    return np.array(ordered_vertices)

def counting_quadrangles(adj):
    '''
        py:funcion:: counting_quadrangles(adj)

        This function takes in the adjacency matrix adj and returns a list of
        quadrangles inside the graph.

        As pointed out by the 1985 paper we're citing and using, one quadrangle
        can be expressed by two vertices v and w, and a set of vertices U, such that
        they're adjacent to both v and w. Such set U is called "bipartition" and
        to form a quadrangle, at least two vertices must be in there.
        The returned list will look like this:

        set([(w,v, (u1, u2, ..., u_l)), ...])

        :param a: numpy NxN matrix, i.e. the adjacency matrix

        :return: a list of tuples, representing the quadrangles in the graph.
    '''

    # Creating the to-be-returned set, first.
    quadrangles_list = list()
    n = len(adj)  # The graph 'size'
    a = np.copy(adj)  # Deepcopying the adj matrix, so that we won't modify the original one

    # Implementing, again, the Chiba-Nishizeki algorithm for 'quadrangles'.

    # As demanded by the algorithm, the vertices exploration must go by their
    # degree, with a non-increasing order
    ordered_vertices = sort_by_degree(adj)

    for v in ordered_vertices:
        # This list represents the (u1, ..., u_l) set cited in the paper ("U")
        bipartitions = [[] for i in range(n)]

        for u, el1 in enumerate(adj[v,:]):  # For each u adj to v
            if el1!=0:  # if u is neighbour of v (i.e. adj to v)
                # Then, the paper says that every neighbour of u is a
                # "neighbour of distance 2 from v". There we go:
                for w, el2 in enumerate(adj[u,:]):  # For each w adj to u
                    # if w is neighbour of u (i.e. adj) and such that w!=v
                    if el2!=0 and w!=v:
                        # We save every u "between" v and w.
                        bipartitions[w].append(u)

        # For each element in the bipartition list:
        for w, b in enumerate(bipartitions):
            # if the bipartition has at least two elements,
            # we have the data structure ready, representing the quadrangle(s).
            if len(b)>1:
                quadrangles_list.append((v,w,b))

        # Erasing v from the graph
        a[v,:] = 0
        a[:,v] = 0

    return quadrangles_list

def k_bron_kerbosch(a, vertices_labels, r, p, x, cliques_list=[], k=4):
    '''
        py:funcion:: k_bron_kerbosch(a, vertices_labels, r, p, x, cliques_list=set(), k=4)

        This function takes in the matrix 'a' and the 'vertices_labels' list obtained
        in the preprocessing phase (see sort_by_degree function), and returns a
        list of k-sized cliques inside the graph.
        This method uses an adaptation of the Bron-Kerbosch algorithm, that instead of
        searching for MAXIMAL cliques inside our graph, it searches for cliques of
        size k. This algorithm runs in O(nk), and it is efficient as long as
        k = 0(log(n)). This algorithm is recursive. Pseudocode:

        k_BronKerbosch(R, P, X, k):
            if |R| = k:
                report R as a k-clique
            else
                choose a pivot u in P ⋃ X
                for each vertex v in P-N(u):
                    BronKerbosch1(R ⋃ {v}, P ⋂ N(v), X ⋂ N(v))
                    P := P - {v}
                    X := X ⋃ {v}

        The returned list will look like this:

        [(v_1,v_2,...,v_k), ...]

        # TODO: fare la descrizione delle variabili e il loro ruolo
        :param a: numpy NxN+1 matrix, i.e. the new sorted_by_degree adj matrix
        :param vertices_labels: list of size N, i.e. the map between old and new indexes.
        :param k: positive integer, i.e. the size of the clique we're investigating
        :return: a list of tuples, representing the quadrangle in the graph.
    '''

    # If the length of the clique-set is k, we've found the k-sized clique
    # we've been looking for.
    if len(r)==k:
        return [r]  # We return it as a list, so that we can merge it later

    if len(p)==0:
        return [r]

    # Next, we want to implement:
    # "choose a pivot u in P ⋃ X"
    x_U_p = list(x|p)  # This is the union between p and x.
    if len(x_U_p)==0:
        return []

    # Pivoting phase.
    # We want to find the vertex which has the highest degree among the xUp set
    u_min = np.argmin(vertices_labels[x_U_p])
    u = x_U_p[u_min]  # Vertex u s.t. degree is maximum in the set P⋃X

    # Then, we extract every neighbour of u. We want it to be a set so that
    # we can do the p-N(u) described in the algorithm.
    neighbours_u = set([
    	i
    	for i, el in enumerate(a[vertices_labels[u],1:])
    	if el==1
    ])

    for v in p-neighbours_u:
        # Recursive step of the algorithm
        cliques_list += k_bron_kerbosch(
            a,
            vertices_labels,
            r|{v},
            p&neighbours_u,
            x&neighbours_u,
            cliques_list=cliques_list,
            k=k
        )
        p.remove(v)
        x.add(u)

    return cliques_list

def counting_triangles(adj):
    '''
        py:funcion:: counting_triangles(a, vertices_labels)

        This function takes in the adjacency matrix adj and returns a list of
        triangles inside the graph.

        As pointed out by Chiba and Nishizeki (1985), one triangle can be
        expressed by two vertices v and w that share a common neighbour u.
        The returned list will look like this:

        set([(w, v, u), ...])

        :param adj: numpy NxN matrix, i.e. the adjacency matrix

        :return: a list of tuples, representing the triangles in the graph.
    '''

    triangles_list = list()  # Used to save triangles.
    a = np.copy(adj)  # Deepcopying the adj matrix, so that we won't modify the original one

    # The following algorithm is a straightforward implementation of the methods
    # used in Chiba, Nishizeki (1985); you can check their paper for more info.

    # As demanded by the algorithm, the vertices exploration must go by their
    # degree, with a non-increasing order
    ordered_vertices = sort_by_degree(adj)

    for v in ordered_vertices[:-2]:
        colouring = [  # Colouring every neighbour of v
            True if w!=v and el!=0 else False
            for w, el in enumerate(a[v,:])
        ]

        for u, el1 in enumerate(a[v,:]):  # "for each 'coloured' neighbour of v, say u"
            if el1!=0 and u!=v:  # the 'coloured' condition is here
                for w, el2 in enumerate(a[u,:]):  # Then, search for common neighbours
                    if el2!=0 and colouring[w]:  # is w a common neighbour?
                        # triangle found
                        triangles_list.append((u,v,w))

                # Removing the colour from u
                colouring[u] = False

        # 'Erasing' v from the graph to prevent duplicates
        a[v,:] = 0
        a[:,v] = 0

    return triangles_list

def clique(adj, k=4):
    '''
        ... # TODO: descrizione
    '''

    n = len(adj)
    w_mat = np.zeros((n,n), dtype=np.uint32)

    # A little verification on the clique size given.
    int_k = int(k)  # We also want to make sure the size is an int number.
    if k<4:
        print("WARNING: Clique size is weird (given: "+str(k)+")")
        if k==2:
            print("Computing the classic HOTNET2 algorithm")
            print("You can obtain the same results with by typing the option:")
            print("-m no")
            print("when launching this script.")
            return no(adj)
        elif k==3:
            print("Computing the triangle counting algorithm")
            print("You can obtain the same results with by typing the option:")
            print("-m triangle")
            print("when launching this script.")
            return triangle(adj)
        elif k!=int_k:
            raise ValueError("A clique with floating point size isn't acceptable")
        else:
            raise ValueError("The size given is too small. We reccomend 4<=k<=7")

    
    for c in cliques_list:
        # c is a set containing a clique of k elements.
        # Therefore, every possible edge of the weighted graph receives a +1
        c_list = list(c)  # Because we want to index it in order to fully enumerate
        for index, i in enumerate(c_list):
            for j in c_list[index+1:]:
                w_mat[i,j] += 1

                # Because of symmetry
                w_mat[j,i] += 1

    print(w_mat)
    input()
    return w_mat

def tailed_triangle(adj):
    '''
        py:funcion:: tailed_triangle(adj)

        Computing the triangle (with a "tail") occurences inside the graph; this
        method exploits triangle counting (see triangle); for each triangle found,
        this method investigates if there's a neighbour for each of the triangle vertices
        which isn't inside the triangle itself. If the answer is yes, that is a
        tailed triangle.
        As always, our objective is to create the new W matrix, with
        Wij = #triangles on the edge.
        In other words, we want to count the edges participating in occurences
        of the motif considered (i.e. undirected tailed triangle).

        :param adj: numpy NxN matrix, i.e. the adjacency matrix of the graph

        :return: numpy NxN matrix, i.e. the new adj matrix representing the motif count
    '''

    n = len(adj)
    w_mat = np.zeros((n,n), dtype=np.uint32)

    # First, we call the function that actually counts triangles in the graph
    triangles_list = counting_triangles(adj)

    # Now, for each triangle found, we want to investigate if they've got a tail
    # A tail could be found in any of the three vertices of the triangle.
    for triangle in triangles_list:
        u = triangle[0]
        v = triangle[1]
        w = triangle[2]

        # For each vertex in the triangle, we're investigating the presence of
        # another neighbour which is not already participating in the triangle
        # (i.e. the tail itself)

        # Investigating u's neighbours
        for j, el in enumerate(adj[u,:]):
            if j!=v and j!=w and el==1:
                # tailed triangle found.
                w_mat[u,v] += 1
                w_mat[v,w] += 1
                w_mat[w,u] += 1
                w_mat[u,j] += 1
                # Because of undirected graphs symmetry
                w_mat[v,u] += 1
                w_mat[w,v] += 1
                w_mat[u,w] += 1
                w_mat[j,u] += 1

        # Investigating v's neighbours
        for j, el in enumerate(adj[v,:]):
            if j!=u and j!=w and el==1:
                # tailed triangle found.
                w_mat[u,v] += 1
                w_mat[v,w] += 1
                w_mat[w,u] += 1
                w_mat[v,j] += 1
                # Because of undirected graphs symmetry
                w_mat[v,u] += 1
                w_mat[w,v] += 1
                w_mat[u,w] += 1
                w_mat[j,v] += 1

        # Investigating w's neighbours
        for j, el in enumerate(adj[w,:]):
            if j!=v and j!=u and el==1:
                # tailed triangle found.
                w_mat[u,v] += 1
                w_mat[v,w] += 1
                w_mat[w,u] += 1
                w_mat[w,j] += 1
                # Because of undirected graphs symmetry
                w_mat[v,u] += 1
                w_mat[w,v] += 1
                w_mat[u,w] += 1
                w_mat[j,w] += 1

    return w_mat

def triangle(adj):
    '''
        py:funcion:: triangle(adj)

        Computing the triangle occurences inside the graph; this
        method is well explained by Chiba and Nishizeki (1985).
        As always, our objective is to create the new W matrix, with
        Wij = #triangles on the edge.
        In other words, we want to count the edges participating in occurences
        of the motif considered (i.e. triangle).

        :param adj: numpy NxN matrix, i.e. the adjacency matrix of the graph

        :return: numpy NxN matrix, i.e. the new adj matrix representing the motif count
    '''
    n = len(adj)
    w_mat = np.zeros((n,n), dtype=np.uint32)

    # First, we call the function that actually counts triangles in the graph
    triangles_list = counting_triangles(adj)

    # Now, for each triangle found, we transition this into our motif_count matrix.
    for triangle in triangles_list:
        u = triangle[0]
        v = triangle[1]
        w = triangle[2]

        # But we actually want to count the edges participating in this motif.
        w_mat[u,v] += 1
        w_mat[u,w] += 1
        w_mat[v,w] += 1
        # Because of symmetry, we repeat the process.
        w_mat[w,v] += 1
        w_mat[w,u] += 1
        w_mat[v,u] += 1

    return w_mat

def tailed_quadrangle(adj):
    '''
        py:funcion:: tailed_quadrangle(adj)

        Computing the tailed quadrangle occurences inside the graph; the quadrangle
        method is well explained by Chiba and Nishizeki (1985). For each quadrangle,
        it is easy to find a tailed motif iff a vertex in the quadrangle has a
        neighbour which is not inside the quadrangle itself.
        As always, our objective is to create the new W matrix, with
        Wij = #triangles on the edge.
        In other words, we want to count the edges participating in occurences
        of the motif considered (i.e. tailed quadrangles).

        :param adj: numpy NxN matrix, i.e. the adjacency matrix of the graph

        :return: numpy NxN matrix, i.e. the new adj matrix representing the motif count
    '''

    n = len(adj)
    w_mat = np.zeros((n,n), dtype=np.uint32)

    # First, we call the function that actually counts quadrangles in the graph
    quadrangles_list = counting_quadrangles(a, vertices_labels)

    # Now, for each quadrangle found, we want to investigate if they've got
    # a tail; a tail could be found in any of the four vertices of the quadrangle
    for q in quadrangles_list:
        v = q[0]
        w = q[1]
        u_set = list(q[2])

        # Just like in the quadrangle method, we're explorating every possible
        # combination of u_i, u_j (with 1<=i<j<=N)
        for i, u_i in enumerate(u_set):
            for u_j in u_set[i+1:]:
                # Here, we have our quadrangle:
                # (v,w,u_i,u_j)

                # For each quadrangle, we're now investigating the presence of a
                # "tail", similarly to what we've done in 'tailed_triangle'

                # Now, investigating v's neighbours
                for j, el in enumerate(adj[v,:]):
                    if j!=v and j!=u_i and j!=u_j and j!=w and el==1:
                        # tailed quadrangle found.
                        w_mat[v,u_i] += 1
                        w_mat[v,u_j] += 1
                        w_mat[u_j,w] += 1
                        w_mat[u_i,w] += 1
                        w_mat[v,j] += 1
                        # Because of undirected graphs symmetry
                        w_mat[u_i,v] += 1
                        w_mat[u_j,v] += 1
                        w_mat[w,u_j] += 1
                        w_mat[w,u_i] += 1
                        w_mat[j,v] += 1

                # Now, investigating w's neighbours
                for j, el in enumerate(adj[w,:]):
                    if j!=v and j!=u_i and j!=u_j and j!=w and el==1:
                        # tailed quadrangle found.
                        w_mat[v,u_i] += 1
                        w_mat[v,u_j] += 1
                        w_mat[u_j,w] += 1
                        w_mat[u_i,w] += 1
                        w_mat[w,j] += 1
                        # Because of undirected graphs symmetry
                        w_mat[u_i,v] += 1
                        w_mat[u_j,v] += 1
                        w_mat[w,u_j] += 1
                        w_mat[w,u_i] += 1
                        w_mat[j,w] += 1

                # Now, investigating u_j's neighbours
                for j, el in enumerate(adj[u_j,:]):
                    if j!=v and j!=u_i and j!=u_j and j!=w and el==1:
                        # tailed quadrangle found.
                        w_mat[v,u_i] += 1
                        w_mat[v,u_j] += 1
                        w_mat[u_j,w] += 1
                        w_mat[u_i,w] += 1
                        w_mat[u_j,j] += 1
                        # Because of undirected graphs symmetry
                        w_mat[u_i,v] += 1
                        w_mat[u_j,v] += 1
                        w_mat[w,u_j] += 1
                        w_mat[w,u_i] += 1
                        w_mat[j,u_j] += 1

                # Now, investigating u_i's neighbours
                for j, el in enumerate(adj[u_i,:]):
                    if j!=v and j!=u_i and j!=u_j and j!=w and el==1:
                        # tailed quadrangle found.
                        w_mat[v,u_i] += 1
                        w_mat[v,u_j] += 1
                        w_mat[u_j,w] += 1
                        w_mat[u_i,w] += 1
                        w_mat[u_i,j] += 1
                        # Because of undirected graphs symmetry
                        w_mat[u_i,v] += 1
                        w_mat[u_j,v] += 1
                        w_mat[w,u_j] += 1
                        w_mat[w,u_i] += 1
                        w_mat[j,u_i] += 1

    return w_mat

def quadrangle(adj):
    '''
        py:funcion:: quadrangle(adj)

        Computing the quadrangle occurences inside the graph; this
        method is well explained by Chiba and Nishizeki (1985).
        As always, our objective is to create the new W matrix, with
        Wij = #triangles on the edge.
        In other words, we want to count the edges participating in occurences
        of the motif considered (i.e. quadrangles).

        :param adj: numpy NxN matrix, i.e. the adjacency matrix of the graph

        :return: numpy NxN matrix, i.e. the new adj matrix representing the motif count
    '''

    n = len(adj)
    w_mat = np.zeros((n,n), dtype=np.uint32)

    # First, we call the function that actually counts quadrangles_list in the graph
    quadrangles_list = counting_quadrangles(adj)

    # The quadrangles found inside the graph are actually stored in this list
    # in a way well explained in Chiba and Nishizeki (1985). Check the function
    # called "counting_quadrangles" to check the data structure's format.
    for q in quadrangles_list:
        v = q[0]
        w = q[1]
        # u_set is the set of vertices that form a complete bipartition.
        # (i.e. nodes in u are adjacent to both v and w)
        u_set = list(q[2])

        # For this reason, to completely enumerate every quadrangle inside this
        # set, we should consider every possible combination in the u_set.
        for i, u_i in enumerate(u_set):
            for u_j in u_set[i+1:]:
                # For every quadrangle found, we save the occurences.
                w_mat[v,u_i] += 1
                w_mat[v,u_j] += 1
                w_mat[u_j,w] += 1
                w_mat[u_i,w] += 1

                # And, because of symmetry:
                w_mat[u_i,v] += 1
                w_mat[u_j,v] += 1
                w_mat[w,u_j] += 1
                w_mat[w,u_i] += 1

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
        i if sum(motif_mat[i,:])==0 else -1
        for i in range(n)
    ]
    disconnected_vertices = [
        el for el in disconnected_vertices
        if el>=0
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
    print("Filtered out "+str(n-k)+" disconnected nodes")

    # Re-write them by reference, in place, so that we don't need to return
    # the new indexes.
    # TODO: We've de-bugged this already, but I want to be sure it works. # DEBUG: further debugging
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
    h = np.zeros((n,n), dtype=np.uint32)
    h[e>=delta] = 1

    return h

def extract_strong_cc(h):
    H = nx.from_numpy_matrix(h, create_using=nx.DiGraph())

    subgraphs = [
        set(c)
        for c in nx.strongly_connected_components(H)
        if len(c)>1
    ]

    return subgraphs

"""
    TODO: chiedere al professore:
    - VA BENE LA SEGUENTE VERSIONE "SOFT"?
        Perché ignorare _completamente_ i lati che non partecipano ai motif?
        Perché semplicemente non lasciare peso unitario ai lati in generale,
        aumentando invece il peso dei lati che partecipano ai motif?
"""
