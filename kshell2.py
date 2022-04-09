import networkx as nx
import pickle
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['xtick.labelsize'] = 15
mpl.rcParams['ytick.labelsize'] = 15

font = {'family':'normal','size':20}

mpl.rc('font',**font)

def normalize(L):

    s = max(L.values())

    for u in L.keys():
        L[u] = float(L[u])/float(s)

    return L

def mean_centrality(x,y,c):

    maxk = max([x[u] for u in x.keys()])

    xavg = [i for i in range(maxk + 1)]
    yavg = [(0.0,0.0) for i in range(maxk + 1)]

    for u in x.keys():
        yavg[x[u]] = (yavg[x[u]][0] + y[u],yavg[x[u]][1] + 1.0)

    for i in range(len(yavg)):
        if yavg[i] == (0.0,0.0):
            yavg[i] = 0.0
        else:
            yavg[i] = float(yavg[i][0])/float(yavg[i][1])

    print (len(xavg),len(yavg))

    plt.plot(yavg,xavg,color = c,marker = 'o')

def correlate(x,y,G):

    xaxis = [x[u] for u in sorted(G.nodes())]

    yaxis = [y[u] for u in sorted(G.nodes())]

    for i in range(len(xaxis)):
        for j in range(i + 1,len(xaxis)):

            if xaxis[i] > xaxis[j]:
                temp = xaxis[i]
                xaxis[i] = xaxis[j]
                xaxis[j] = temp

                temp = yaxis[i]
                yaxis[i] = yaxis[j]
                yaxis[j] = temp

    plt.scatter(xaxis,yaxis,s = 3)


    z = np.polyfit(xaxis,yaxis,2)
    p = np.poly1d(z)

    plt.plot(xaxis,[p(pt) for pt in xaxis])
    #plt.ylim([-10,200])
    plt.show()

def checkdeg(G,i):

    return [u for u in G.nodes() if G.degree(u) <= i]


def kshelld(G):

    kshell = {}
    for u in G.nodes():
        kshell[u] = 0.0

    k = 0

    while len(G) > 0:

        print (k)
        while True:

            nlist = checkdeg(G,k)

            if len(nlist) <= 0:
                break

            for u in nlist:
                kshell[u] = k

            G.remove_nodes_from(nlist)

        k = k + 1

    return kshell

s = 'Yeast'
G = nx.read_gml(s + '_Ordered.gml')
G = G.to_undirected()

#NMC = pickle.load(open('NMC_' + s +'.p','r'))
#NMC = normalize(NMC)

M = pickle.load(open('Motif_' + s +'.p','r'))

TypeA = {}
TypeB = {}
TypeC = {}

for u in G.nodes():
    TypeA[u] = 0.0
    TypeB[u] = 0.0
    TypeC[u] = 0.0

for m in M:
    if not (m[0] == m[1] or m[1] == m[2] or m[0] == m[2]):
        TypeA[m[0]] += 1
        TypeB[m[1]] += 1
        TypeC[m[2]] += 1


kshell = kshelld(G.copy())
correlate(TypeA,kshell,G)
correlate(TypeB,kshell,G)
correlate(TypeC,kshell,G)

#print (kshell)



