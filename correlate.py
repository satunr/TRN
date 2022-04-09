import pickle
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import networkx as nx

mpl.rcParams['xtick.labelsize'] = 15
mpl.rcParams['ytick.labelsize'] = 15

font = {'family': 'normal',
        'size' : 20}

mpl.rc('font',**font)

def between(G):

    BC = nx.betweenness_centrality(G)
    return BC

def closeness(G):

    CC = nx.closeness_centrality(G)
    return CC

def normalize(D):

    S = max(D.values())
    for u in D.keys():
        D[u] = float(D[u])/float(S)


    return D

def plotter(G,NMC,PC,colors):

    BC = between(G)
    CC = closeness(G)
    DC = {u:G.out_degree(u) for u in G.nodes()}

    BC = normalize(BC)
    CC = normalize(CC)
    NMC = normalize(NMC)
    DC = normalize(DC)

    y = [CC[u] for u in sorted(NMC.keys())]
    x = [NMC[u] for u in sorted(NMC.keys())]

    correlation(x,y,colors)
    print (np.corrcoef(x,y))

def correlation(x,y,colors):

    for i in range(len(x)):
        for j in range(i + 1,len(x)):

            if x[i] > x[j]:
                temp = x[i]
                x[i] = x[j]
                x[j] = temp

                temp = y[i]
                y[i] = y[j]
                y[j] = temp

    plt.scatter(x,y,s = 2,color = colors)

    z = np.polyfit(x,y,2)
    p = np.poly1d(z)

    plt.xlim([0,0.7])
    #plt.ylim([0,1.0])

    plt.plot(x,[p(pt) for pt in x],linewidth = 2,color = colors)

    #plt.show()

color = ['r','g','b','black']
GY = nx.read_gml('Yeast_Ordered.gml')
PC_Y = [GY.out_degree(u) for u in GY.nodes()]

NMC_Y = pickle.load(open('NMC_Y.p','rb'))
#PC_Y = pickle.load(open('PC_Y.p','r'))
NMC_Y = normalize(NMC_Y)

GM = nx.read_gml('Mouse_Ordered.gml')
PC_M = [GM.out_degree(u) for u in GM.nodes()]
NMC_M = pickle.load(open('NMC_M.p','rb'))
#PC_M = pickle.load(open('PC_M.p','r'))
NMC_M = normalize(NMC_M)

GH = nx.read_gml('Human_Ordered.gml')
PC_H = [GH.out_degree(u) for u in GH.nodes()]

NMC_H = pickle.load(open('NMC_H.p','rb'))
#PC_H = pickle.load(open('PC_H.p','r'))
NMC_H = normalize(NMC_H)

GE = nx.read_gml('Ecoli_Ordered.gml')
PC_E = [GE.out_degree(u) for u in GE.nodes()]
NMC_E = pickle.load(open('NMC_E.p','rb'))
#PC_H = pickle.load(open('PC_H.p','r'))
NMC_E = normalize(NMC_E)


#G = nx.read_gml('Yeast_Ordered.gml')
#NMC = pickle.load(open('NMC_Y.p','r'))
#NMC = normalize(NMC)

plotter(GY,NMC_Y,PC_Y,color[0])
plotter(GH,NMC_H,PC_H,color[1])
plotter(GM,NMC_M,PC_M,color[2])
plotter(GE,NMC_E,PC_E,color[3])
plt.savefig('closeness_correlate.png',dpi = 300)
plt.show()


plt.show()
