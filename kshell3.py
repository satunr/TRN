import networkx as nx
import pickle
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import math

mpl.rcParams['xtick.labelsize'] = 15
mpl.rcParams['ytick.labelsize'] = 15

font = {'family':'normal','size':20}

mpl.rc('font',**font)

def kshelld(G):

    H = G.to_undirected()
    ks = {u:0 for u in H.nodes()}

    #H = largest_component(H)

    d = 0
    while len(H) > 0:

        while True:

            nlist = [u for u in H.nodes() if H.degree(u) <= d]

            if len(nlist) == 0:
                d = d + 1
                break

            for u in nlist:
                ks[u] = d

            H.remove_nodes_from(nlist)

    ks = {k:float(ks[k]) for k in ks.keys()}

    return ks,d

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


def checkdeg(G,i):

    return [u for u in G.nodes() if G.degree(u) <= i]


def find_roles(G,M):

    r_a_y = {u:0 for u in G.nodes()}

    for m in M:
        if m[0] in r_a_y.keys():
            r_a_y[m[0]] += 1

    return r_a_y

path = '/Users/satyakiroy/PycharmProjects/Main/whyGRNLinux/GRN ORIGINALs/'

'''
G = nx.read_gml(path + 'Mapped_Ecoli.gml')
G = G.to_undirected()
motifs_e = pickle.load(open(path + 'motif_Ecoli.p','rb'))
r_a_e = find_roles(G,motifs_e)
r_a_e = normalize(r_a_e)
kshell,d = kshelld(G.copy())
mean_centrality(kshell,r_a_e,'red')

G = nx.read_gml(path + 'Mapped_Yeast.gml')
G = G.to_undirected()
motifs_y = pickle.load(open(path + 'motif_Yeast.p','rb'))
r_a_y = find_roles(G,motifs_y)
r_a_y = normalize(r_a_y)
kshell,d = kshelld(G.copy())
mean_centrality(kshell,r_a_y,'green')

G = nx.read_gml(path + 'Mapped_Mouse.gml')
G = G.to_undirected()
motifs_m = pickle.load(open(path + 'motif_Mouse.p','rb'))
r_a_m = find_roles(G,motifs_m)
r_a_m = normalize(r_a_m)
kshell,d = kshelld(G.copy())
mean_centrality(kshell,r_a_m,'blue')

'''

TRNs = ['Ecoli','Yeast','Mouse','Human']
colors = ['red','green','blue','black']

for s in TRNs:

    X = []
    Y = []

    G = nx.read_gml(path + 'Mapped_' + s + '.gml')
    G = G.to_undirected()
    motifs_h = pickle.load(open(path + 'motif_' + s + '.p','rb'))
    r_a_h = find_roles(G,motifs_h)
    max_r_a_h = max(r_a_h.values())
    #r_a_h = normalize(r_a_h)
    kshell,d = kshelld(G.copy())

    x = list(range(11))
    y = [(0.0,0.0) for i in range(len(x))]

    for u in G.nodes():

        v = int(math.ceil(float(r_a_h[u])/float(max_r_a_h) * 10.0))
        y[v] = (y[v][0] + kshell[u],y[v][1] + 1)

    z = []

    for i in range(len(y)):

        z.append(None)
        if y[i][1] > 0:
            z[-1] = float(y[i][0]) / float(y[i][1])

            X.append(float(i)/10.0)
            Y.append(z[-1])

    plt.plot(X[1:],Y[1:],color = colors[TRNs.index(s)],marker = 'o')

plt.savefig('/Users/satyakiroy/PycharmProjects/Main/whyGRNLinux/GRN ORIGINALs/kshell_correlate.png',dpi = 300)
plt.show()
#mean_centrality(kshell,r_a_h,'black')

#plt.savefig(path + 'kshell_correlate',dpi = 300)
#plt.xlim([-0.01,0.05])
#plt.show()
