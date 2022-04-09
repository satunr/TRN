import pickle
import networkx as nx
from multiprocessing import Pool




def centrality(L):

    G = L[0]
    N = L[1]
    indices = L[2]
    i = L[3]
    M = []

    nmc = {}
    for u in G.nodes():
        nmc[u] = 0.0

    nsm = [[] for u in G.nodes()]

    for u in N[indices[i]:indices[i + 1]]:
        if i == 0 and u % 10 == 0:
            print (u, 'Human')

        for v in G.nodes():
            if v == u:
                continue
            for w in G.nodes():
                if w == u or w == v:
                    continue

                if G.has_edge(u,v) and G.has_edge(v,w) and G.has_edge(u,w):

                    M.append((u,v,w))
                    nmc[u] += 1
                    nmc[v] += 1
                    nmc[w] += 1

                    nsm[u].extend([v,w])
                    nsm[v].extend([u,w])
                    nsm[w].extend([u,v])



    return [M,nsm,nmc]

#G = nx.erdos_renyi_graph(n = 20,p = 0.1,directed = True)

G = nx.read_gml('Human_Ordered.gml')
#G = nx.convert_node_labels_to_integers(G,first_label = 0)

nG = len(G)
N = sorted(G.nodes())
cores = 8
p = Pool(processes = cores)

#Limit of cutoff
limit = 6

indices = [int(nG/cores) * i for i in range(cores)]
indices.append(nG)

L = p.map(centrality,[[G,N,indices,i] for i in range(cores)])

NMC = {}
for u in G.nodes():
    NMC[u] = 0.0

NSM = [[] for u in G.nodes()]
M = []

for each in L:

    m = each[0]
    nsm = each[1]
    nmc = each[2]
    M.extend(m)

    for i in range(len(nsm)):

        NMC[i] += nmc[i]
        NSM[i].extend(nsm[i])

NSM = [list(set(each)) for each in NSM]

print (NSM)
print (NMC)
print (len(G),len(G.edges()))

pickle.dump(M,open('Motif_Human.p','wb'))
pickle.dump(NSM,open('NSM_H.p','wb'))
pickle.dump(NMC,open('NMC_H.p','wb'))
#nx.write_gml(G,'Yeast_Ordered.gml')
