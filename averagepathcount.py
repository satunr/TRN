import networkx as nx

from multiprocessing import Pool

def pathc(L):

    G = L[0]
    N = L[1]
    indices = L[2]
    i = L[3]
    limit = L[4]

    den = 0.0
    num = 0.0

    for u in N[indices[i]:indices[i + 1]]:
        if i == 0 and u % 10 == 0:
            print u

        for v in G.nodes():
            if nx.has_path(G,u,v):
                num += len(list(nx.all_simple_paths(G,source = u,target = v,cutoff = limit)))
                den += 1

    return (num,den)

G = nx.read_gml('Mouse_Ordered.gml')
G = nx.convert_node_labels_to_integers(G,first_label = 0)
#G_COM = nx.read_gml('ScaleM.gml')

cores = 8
p = Pool(processes = cores)

nG = len(G)
N = sorted(G.nodes())

#Limit of cutoff
limit = 6

indices = [int(nG/cores) * i for i in range(cores)]
indices.append(nG)

L = p.map(pathc,[[G,N,indices,i,limit] for i in range(cores)])

NUM = 0.0
DEN = 0.0

for each in L:

    NUM += each[0]
    DEN += each[1]

print "Mouse: ",float(NUM/DEN)
