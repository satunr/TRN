import networkx as nx
import pickle
import numpy as np
import matplotlib.pyplot as plt



def fit(x,y):

    print (x)

    for i in range(len(x) - 1):
        for j in range(i + 1,len(x)):

            if x[i] > x[j]:
                temp = x[i]
                x[i] = x[j]
                x[j] = temp

                temp = y[i]
                y[i] = y[j]
                y[j] = temp


    z = np.polyfit(x,y,2)
    p = np.poly1d(z)
    return x,[p(pt) for pt in x]

def tiers(G):

    t1 = []
    t2 = []
    t3 = []

    for u in G.nodes():
        if G.in_degree(u) == 0:
            t1.append(u)
        elif G.out_degree(u) == 0:
            t3.append(u)
        else:
            t2.append(u)

    return G,t1,t2,t3

def max_path_length(G):

    m = 0.0

    for u in G.nodes():
        for v in G.nodes():
            if nx.has_path(G,u,v):
                s = nx.shortest_path_length(G,u,v)
                if s > m:
                    m = s

    return m

'''
def edge_disjoint_paths():

    M = pickle.load(open('Motif_' + s + '.p','rb'))
    G = nx.read_gml(s + '_Ordered.gml')
    NMC = pickle.load(open('NMC_' +s[0] + '.p','rb'))
    t1, t2, t3 = tiers(G)

    print list(nx.edge_disjoint_paths(G,0,5))
'''


def correlate():

    s = 'Human'

    M = pickle.load(open('Motif_' + s + '.p','rb'))
    G = nx.read_gml(s + '_Ordered.gml')
    NMC = pickle.load(open('NMC_' +s[0] + '.p','rb'))
    t1, t2, t3 = tiers(G)

    '''
    TypeA = [0.0 for u in G.nodes()]
    TypeB = [0.0 for u in G.nodes()]

    for m in M:
        TypeA[m[0]] += 1
        TypeB[m[1]] += 1

    N = [u for u in G.nodes() if u in t2]

    plt.scatter([TypeA[u] for u in N],[TypeB[u] for u in N])
    x,y = fit([TypeA[u] for u in N],[TypeB[u] for u in N])
    plt.plot(x,y)

    plt.title('Ecoli')
    plt.xlabel('Type A')
    plt.ylabel('Type B')

    plt.show()
    '''

    N = [u for u in G.nodes() if NMC[u] >= 100]
    t1_count = 0
    t2_count = 0
    t3_count = 0

    for u in N:
        if u in t1:
            t1_count += 1

        elif u in t2:
            t2_count += 1

        else:
            t3_count += 1

    print t1_count,t2_count,t3_count



def path_length():

    M = pickle.load(open('Motif_Yeast.p','rb'))
    G = nx.read_gml('Yeast_Ordered.gml')
    NMC = pickle.load(open('NMC_Y.p','rb'))
    t1, t2, t3 = tiers(G)

    N = [u for u in G.nodes() if u in t2]
    P = [[0.0 for u in range(len(N))] for v in range(len(N))]

    R = {}
    Y = [NMC[u] for u in sorted(N)]
    X = [u for u in sorted(N)]

    Z = [x for _,x in sorted(zip(Y,X))]
    print (Z)

    for u in sorted(N):
        R[u] = Z.index(u)

    m = max_path_length(G) + 10.0

    for u in N:
        for v in N:


            if nx.has_path(G,u,v):
                s = nx.shortest_path_length(G,u,v)
                P[R[u]][R[v]] = s
            else:
                P[R[u]][R[v]] = m


            '''
            if G.has_edge(u,v):
                P[R[u]][R[v]] = 10
            else:
                P[R[u]][R[v]] = 0
            '''
    plt.imshow(P,cmap = 'summer',interpolation = 'nearest')
    plt.show()

def motif_lost():

    M = pickle.load(open('Motif_Yeast.p','rb'))
    G = nx.read_gml('Yeast_Ordered.gml')
    NMC = pickle.load(open('NMC_Y.p','rb'))

    t1,t2,t3 = tiers(G)

    N = [u for u in G.nodes() if u in t2 and NMC[u] > 0]
    B = nx.betweenness_centrality(G)

    #Centrality of type A and B
    TypeA_C = {}
    TypeB_C = {}
    for u in G.nodes():
        TypeA_C[u] = 0.0
        TypeB_C[u] = 0.0

    for m in M:
        TypeA_C[m[0]] += 1
        TypeB_C[m[1]] += 1


    for u in N:

        plt.scatter(TypeA_C[u],G.out_degree(u),c = 'blue')
        plt.scatter(TypeB_C[u],G.out_degree(u),c = 'green')

    x,y = fit([TypeA_C[u] for u in N],[G.out_degree(u) for u in N])
    plt.plot(x,y,label = 'Type A',c = 'blue')

    x,y = fit([TypeB_C[u] for u in N],[G.out_degree(u) for u in N])
    plt.plot(x,y,label = 'Type B',c = 'green')
    plt.title('Yeast')
    plt.xlabel('Node Motif Centrality')
    plt.ylabel('Out-degree')
    plt.legend()
    plt.show()

def tier_roles():


    G = nx.read_gml('Mouse_Ordered.gml')
    M = pickle.load(open('Motif_Mouse.p','rb'))
    NMC = pickle.load(open('NMC_M.p','rb'))
    G,t1,t2,t3 = tiers(G)


    print (len(NMC.keys()),len(t1),len(t2),len(t3))

    for m in M:

        u = m[0]
        v = m[1]
        w = m[2]

        if NMC[u] < 100 and u in t3:
            plt.scatter(NMC[u],1,s = 1)

        if NMC[v] < 100 and v in t3:
            plt.scatter(NMC[v],2,s = 1)

        if NMC[w] < 100 and w in t3:
            plt.scatter(NMC[w],3,s = 1)

    plt.title('Mouse')
    plt.ylim([0,4])
    plt.show()


    '''
    g = nx.DiGraph()
    g.add_nodes_from([u for u in G.nodes() if u in t2 and NMC[u] < 100])

    for u in g.nodes():
        for v in g.nodes():
            if G.has_edge(u,v):
                g.add_edge(u,v)

    print len(g)
    print len(g.edges())

    g = nx.convert_node_labels_to_integers(g,first_label = 0)

    cnt = 0
    for e in g.edges():

        u = e[0]
        v = e[1]

        if v >= u:
            continue
        if g.has_edge(v,u):
            cnt += 1

    print ("Number of cycles of length 2:",cnt)
    '''
    '''
    d = nx.greedy_color(g,strategy = nx.coloring.strategy_largest_first)
    print (d)

    color = [0.0 for u in g.nodes()]
    for u in g.nodes():

        if d[u] == 0:
            color[u] = 'red'

        elif d[u] == 1:
            color[u] = 'green'

        elif d[u] == 2:
            color[u] = 'brown'

        elif d[u] == 3:
            color[u] = 'black'

        else:
            color[u] = 'blue'

    pos = nx.circular_layout(g)
    nx.draw(g,pos,node_color = color)
    plt.show()


    M = pickle.load(open('Motif_Mouse.p','rb'))

    for m in M:

        u = m[0]
        v = m[1]
        w = m[2]

        if NMC[u] >= 100:
            plt.scatter(NMC[u],1)

        if NMC[v] >= 100:
            plt.scatter(NMC[v],2)

        if NMC[w] >= 100:
            plt.scatter(NMC[w],3)

    plt.title('Mouse')
    plt.show()
    '''

def motif_centrality_correlation(G,t1,t2,t3,NMC):

    X = []
    Y = []
    M = []

    for u in G.nodes():

        if u % 50 == 0:
            print (u,'Human')

        for v in G.nodes():
            if v == u:
                continue

            for w in G.nodes():
                if w == u or w == v:
                    continue

                if G.has_edge(u,v) and G.has_edge(v,w) and G.has_edge(u,w):

                    M.append((u,v,w))

                    if u in t2:
                        X.append(NMC[u])
                        Y.append(1)

                    if v in t2:
                        X.append(NMC[v])
                        Y.append(2)

                    if w in t2:
                        X.append(NMC[w])
                        Y.append(3)

    plt.title('Ecoli')
    pickle.dump(M,open('Motif_Yeast.p','wb'))

    plt.scatter(X,Y,s = 1)
    plt.show()

def dist(N,NMC):

    max_m = max([NMC[u] for u in N])

    Y = [0.0 for i in range(int(max_m) + 1)]
    X = [i for i in range(int(max_m) + 1)]

    for u in N:
        Y[int(NMC[u])] += 1

    plt.plot(X,Y)
    plt.xlim([-10,max(X)])
    plt.ylim([-10,max(Y)])

    plt.show()

def rich_club(G,N,NMC):

    H = G.to_undirected()
    max_m = max([NMC[u] for u in N])

    Y = [0.0 for i in range(int(max_m) + 1)]
    X = [i for i in range(int(max_m) + 1)]

    for k in range(int(max_m) + 1):
        n = [u for u in N if NMC[u] > k]
        e = [e for e in H.edges() if e[0] in n and e[1] in n]

        if len(n) > 1:
            Y[k] = float(2 * len(e))/float(len(n) * (len(n) - 1))

    plt.plot(X,Y)
    plt.show()


#G = nx.read_gml('Mouse_Ordered.gml')

#self_loops = 0.0


#print ("self loops:",self_loops)
#NMC = pickle.load(open('NMC_M.p','rb'))
#print (NMC)


#t1,t2,t3 = tiers(G)
'''
for u in G.nodes():
    if u in t1:
        plt.scatter(NMC[u],1,s = 5)
    elif u in t2:
        plt.scatter(NMC[u],2,s = 5)
    else:
        plt.scatter(NMC[u],3,s = 5)

plt.show()
'''

def tier_roles2(s):

    M = pickle.load(open('Motif_' + s + '.p','rb'))
    G = nx.read_gml(s + '_Ordered.gml')
    NMC = pickle.load(open('NMC_' + s[0] + '.p','rb'))
    G, t1, t2, t3 = tiers(G)

    print (len([u for u in G.nodes() if NMC[u] >= 100 and u in t1]))
    print (len([u for u in G.nodes() if NMC[u] >= 100 and u in t2]))
    print (len([u for u in G.nodes() if NMC[u] >= 100 and u in t3]))


tier_roles2('Human')

#C = nx.closeness_centrality(G,normalized = True)
#motif_centrality_correlation(G,t1,t2,t3,NMC)
#tier_roles()

#rich_club(G,t2,NMC)
#motif_lost()
#path_length()
#correlate()

