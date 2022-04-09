import networkx as nx
import itertools

from multiprocessing import Pool
from copy import deepcopy

def pathcount(I):

    G = I[0]
    N = I[1]
    indices = I[2]
    i = I[3]
    limit = I[4]
    EMC = I[5]

    #Total number of paths
    l_UV = 0.0

    #Total number of paths created by FFL motifs
    c_UV = 0.0

    occurrences = 0

    for u in N[indices[i]:indices[i + 1]]:
        if i == 0 and u % 10 == 0:
            print (u)
        for v in G.nodes():

            if u == v or not nx.has_path(G,u,v):
                continue

            P = list(nx.all_simple_paths(G,source = u,target = v,cutoff = limit))
            #P = P[:100]

            Counted = [False for z in range(len(P))]

            #print "All paths:",P
            l_UV += len(P)
            iterate = 0

            while iterate < len(P):
                if Counted[iterate] == True:
                    iterate += 1
                    continue

                p = P[iterate]

                #print "Current path:",p
                m, flag = creatematrix(EMC,p)

                if flag == True:

                    pc = pathgenerate(m,p,P)
                    #print "pc",pc

                    for i in range(len(pc)):

                        try:
                            if not Counted[P.index(pc[i])]:
                                Counted[P.index(pc[i])] = True
                                c_UV += 1

                        except Exception as e:
                            continue

                iterate += 1

    #print c_UV,l_UV
    return (c_UV,l_UV)

def pathgenerate(m,p,P):

    pc = []

    a = [range(len(w)) for w in m]
    permutations = list(itertools.product(*a))

    for i in range(len(permutations)):

        #new path
        p_duplicate = deepcopy(p)

        for j in range(len(permutations[i])):
            p_duplicate.insert(j * 2 + 1,m[j][permutations[i][j]])

        #Remove all "None" terms from the new path
        p_duplicate = [each for each in p_duplicate if each != None]

        pc.append(p_duplicate)

    return pc

def creatematrix(EMC,p):
    m = []
    flag = False
    for i in range(len(p) - 1):

        v = EMC[(p[i], p[i + 1])]
        m.append(v)

        if len(v) > 1:
            flag = True

    return m,flag

def directEdgeCentrality(I):

    G = I[0]
    N = I[1]
    indices = I[2]
    i = I[3]

    emc = {}
    for e in G.edges():
        emc[(e[0],e[1])] = []

    for u in N[indices[i]:indices[i + 1]]:
        if i == 0 and u % 10 == 0:
            print (u)

        for v in G.nodes():
            if v == u:
                continue
            for w in G.nodes():
                if w == u or w == v:
                    continue
                if G.has_edge(u,v) and G.has_edge(v,w) and G.has_edge(u,w):
                    emc[(u,w)].append(v)

    return emc

cores = 4
p = Pool(processes = cores)

G = nx.DiGraph()
G.add_edges_from([(1,3),(1,4),(3,4),(3,2),(4,2)])
#G.add_edges_from([(0,1),(1,2),(0,3),(3,2),(1,3)])

#G = nx.read_gml('Ecoli_Ordered.gml')
#G = nx.convert_node_labels_to_integers(G,first_label = 0)

C = list(nx.simple_cycles(G))
print (len(C))

for c in C:

    print (c)
    if len(c) == 1 and G.has_edge(c[0],c[0]):
        G.remove_edge(c[0],c[0])

    elif G.has_edge(c[0],c[1]):
        G.remove_edge(c[0],c[1])

nG = len(G)
N = sorted(G.nodes())

#Limit of cutoff
limit = 6

indices = [int(nG/cores) * i for i in range(cores)]
indices.append(nG)

L = p.map(directEdgeCentrality,[[G,N,indices,i] for i in range(cores)])
L_UV = 0.0
C_UV = 0.0

EMC = {}

for e in G.edges():
    EMC[(e[0],e[1])] = [None]

for each in L:
    for e in each.keys():
        EMC[(e[0],e[1])].extend(each[e])

print (EMC)

L = p.map(pathcount,[[G,N,indices,i,limit,EMC] for i in range(cores)])

for i in range(len(L)):

    C_UV += L[i][0]
    L_UV += L[i][1]

print (C_UV,L_UV)
print ("E. coli: ",float(C_UV)/float(L_UV))

