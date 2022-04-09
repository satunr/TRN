import networkx as nx
import pickle

import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['xtick.labelsize'] = 15
mpl.rcParams['ytick.labelsize'] = 15

font = {'family': 'normal',
        'size' : 20}

mpl.rc('font',**font)


#G = nx.DiGraph()
#G.add_edges_from([(0,1), (0,2), (1,2), (1,3), (2,3)])

def dist(G):

    maxd = max([G.out_degree(u) for u in G.nodes()])

    xaxis = [i for i in range(maxd + 1)]
    yaxis = [0 for i in range(maxd + 1)]

    for u in G.nodes():
        yaxis[G.out_degree(u)] += 1

    plt.plot(xaxis,yaxis)
    plt.show()

def decay(NMC):

    Frequency = [len(NMC[u]) for u in range(len(NMC))]

    Total = sum(Frequency)/3
    print Total

    print Frequency

    yaxis = []
    yaxis.append(Total)

    i = 0
    while Total > 0:

        vertex = Frequency.index(max(Frequency))
        print Total

        Total = Total - Frequency[vertex]

        for m in NMC[vertex]:

            u = m[0]
            v = m[1]
            w = m[2]

            if u != vertex:
                NMC[u].remove(m)

            if v != vertex:
                NMC[v].remove(m)

            if w != vertex:
                NMC[w].remove(m)

        NMC[vertex] = []
        Frequency = [len(NMC[z]) for z in range(len(NMC))]

        i = i + 1
        yaxis.append(Total)

    print yaxis
    plt.plot([i for i in range(len(yaxis))],yaxis,marker = 'o',markersize = 2,linewidth = 2)


NMC = pickle.load(open('NMC_HN.p','r'))
decay(NMC)

#NMC = pickle.load(open('NMC_MLN.p','r'))
#decay(NMC)

NMC = pickle.load(open('R_NMC_HLN.p','r'))
decay(NMC)

plt.xlim([-5,120])
plt.ylim([-100,6000])

plt.savefig('motif_decay_human1.png',dpi = 300)
plt.show()
