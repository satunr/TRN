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


    #z = np.polyfit(xaxis[:1000],yaxis[:1000],2)
    #p = np.poly1d(z)

    #plt.plot(xaxis[:1000],[p(pt) for pt in xaxis[:1000]])
    #plt.ylim([-10,200])

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

G = nx.read_gml('Yeast_Ordered.gml')
G = G.to_undirected()
NMC = pickle.load(open('NMC_Y.p','r'))
NMC = normalize(NMC)
kshell = kshelld(G.copy())
mean_centrality(kshell,NMC,'red')

G = nx.read_gml('Ecoli_Ordered.gml')
G = G.to_undirected()
NMC = pickle.load(open('NMC_E.p','r'))
NMC = normalize(NMC)
kshell = kshelld(G.copy())
mean_centrality(kshell,NMC,'black')

G = nx.read_gml('Mouse_Ordered.gml')
G = G.to_undirected()
NMC = pickle.load(open('NMC_M.p','r'))
NMC = normalize(NMC)
kshell = kshelld(G.copy())
mean_centrality(kshell,NMC,'green')

G = nx.read_gml('Human_Ordered.gml')
G = G.to_undirected()
NMC = pickle.load(open('NMC_H.p','r'))
NMC = normalize(NMC)
kshell = kshelld(G.copy())
mean_centrality(kshell,NMC,'blue')

plt.savefig('kshell_correlate',dpi = 300)
plt.show()
