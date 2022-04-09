import networkx as nx

from multiprocessing import Pool

def motifcount(L):

    G = L[0]
    N = L[1]
    i = L[2]
    indices = L[3]

    c = 0.0
    for u in N[indices[i]:indices[i + 1]]:
        if i == 0 and u % 100 == 0:
            print u

        for v in G.nodes():
            for w in G.nodes():

                if G.has_edge(u,v) and G.has_edge(v,w) and G.has_edge(u,w):
                    c = c + 1

    return c

print "E. COLI"
G = nx.read_gml('ScaleE.gml')
N = sorted(G.nodes())
nG = len(G)
cores = 8

indices = [int(nG/cores) * i  for i in range(cores)]
indices.append(nG)

p = Pool(processes = cores)
I = p.map(motifcount,[[G,N,i,indices] for i in range(cores)])

C = 0
for each in I:
    C += each

print C