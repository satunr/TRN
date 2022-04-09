import networkx as nx
import pickle
import numpy as np

def tier(G):

    t1, t2, t3 = [], [], []
    for u in G.nodes():

        if G.in_degree(u) == 0:
            t1.append(u)

        elif G.out_degree(u) == 0:
            t3.append(u)

        else:
            t2.append(u)

    return t1, t2, t3

Organisms = ['Ecoli', 'Yeast', 'Mouse', 'Human']

for o in Organisms:

    G = nx.read_gml(o + '_Ordered.gml')
    P = pickle.load(open('NMC_' + str(o[0]) + '.p', "rb"))

    t1, t2, t3 = tier(G)

    print (np.mean([P[u] for u in t1]))
    print (np.mean([P[u] for u in t2]))
    print (np.mean([P[u] for u in t3]))
    print ()



