import argparse
import networkx as nx
import numpy as np


# TODO: la funzione che, data adj e il motif, torna W. (vedi fogli)



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

def motif_count(adj, m='u_tri'):
    '''
        py:funcion:: motif_count(adj, m)

        Computing the motif occurences inside the graph; in our application,
        we want to count the _edges_ participating in occurences of m.
        Many papers I've read, instead, focus on counting the vertices of such
        occurences.

        :param adj: numpy NxN matrix, i.e. the adjacency matrix of the graph
        :param m: a string, i.e. motif considered.

        :return: numpy NxN matrix, i.e. the new adj matrix representing the motif count
    '''
    w = None
    if m=='u_tri':
        w = u_triangle_count(adj)
    else:
        raise ValueError("Unknown motif!")

    # write(w) # TODO: scrivere la matrice su file per evitare di ricalcolarla ogni volta

def u_triangle_count(adj):
    '''
        py:funcion:: motif_count(adj, m)

        Computing the undirected triangle occurences inside the graph; this method
        has been inspired by Marcus, Shavitt in their 2010 paper, counting vertices
        that participate in instances of m. Instead, in our application,
        we want to count the _edges_ participating in occurences of m.

        :param adj: numpy NxN matrix, i.e. the adjacency matrix of the graph

        :return: numpy NxN matrix, i.e. the new adj matrix representing the undirected triangle count
    '''
    n = len(adj)
    w_mat = np.zeros((n,n), dtype=np.uint16)
    triangles = set()

    for i, row in enumerate(adj):
        # Iterating through all the rows, i.e. every node once exactly
        for j, el in enumerate(row[i+1:]):
            if el!=0:
                v = i
                u = i+j+1
                colouring = [
            		True if w!=v and elem!=0 else False
            		for w, elem in enumerate(adj[u,:])
            	]
                for w, elem in enumerate(adj[v,:]):
                    if w!=u and elem!=0 and colouring[w]:
                        triangles.add(frozenset({v,u,w}))
    print(triangles)
    input()
    for triple in triangles:
        nodes = list(triple)
        v = nodes[0]
        u = nodes[1]
        t = nodes[2]
        w_mat[u,v] += 1
        w_mat[v,u] += 1
        w_mat[v,t] += 1
        w_mat[t,v] += 1
        w_mat[t,u] += 1
        w_mat[u,t] += 1

    return w_mat
