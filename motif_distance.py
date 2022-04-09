import networkx as nx
import pickle
import numpy as np
import matplotlib.pyplot as plt


def fit(X,Y,order):

    z = np.polyfit(X,Y,order)
    p = np.poly1d(z)

    Z = [p(u) for u in X]
    print (len(Z))

    plt.plot(X[:2450],Z[:2450])


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


    return t1, t2, t3

s = 'Human'

G = nx.read_gml(s + '_Ordered.gml')
NMC = pickle.load(open('NMC_' + s[0] + '.p','rb'))
C = nx.closeness_centrality(G,normalized = True)

t1,t2,t3 = tiers(G)

print (np.mean([NMC[u] for u in t1]),np.std([NMC[u] for u in t1]))
print (np.mean([NMC[u] for u in t2]),np.std([NMC[u] for u in t2]))
print (np.mean([NMC[u] for u in t3]),np.std([NMC[u] for u in t3]))

'''
#----------------------------------------------------------------------
p = np.percentile(NMC.values(),100)
print (p)
N = []
for u in G.nodes():
    if NMC[u] <= p:
        N.append(u)


X = [NMC[u] for u in sorted(N)]
Y = [C[u] for u in sorted(N)]

plt.scatter(X,Y)


for i in range(len(X) - 1):
    for j in range(i + 1,len(X)):

        if X[j] < X[i]:
            temp = X[i]
            X[i] = X[j]
            X[j] = temp

            temp = Y[i]
            Y[i] = Y[j]
            Y[j] = temp


'''

#------------------------------------------

print (NMC)
dim = 25

SP = [[(0.0,0.0) for _ in range(dim)] for _ in range(dim)]
F = [[0.0 for _ in range(dim)] for _ in range(dim)]
#print SP

for u in G.nodes():

    if NMC[u] > dim - 1:
        continue

    for v in G.nodes():
        if NMC[v] > dim - 1:
            continue

        if nx.has_path(G,u,v):
            s = nx.shortest_path_length(G,source = u,target = v)

            val = SP[int(NMC[u])][int(NMC[v])]

            num = val[0]
            den = val[1]

            num = num + s
            den = den + 1

            SP[int(NMC[u])][int(NMC[v])] = (num,den)

for i in range(len(SP)):
    for j in range(len(SP[0])):

        num = SP[i][j][0]
        den = SP[i][j][1]

        if den == 0.0:
            continue

        F[i][j] = float(num)/float(den)

F = np.array(F)
print (F)

plt.imshow(F, cmap='hot', interpolation='nearest')
plt.savefig('/Users/satyakiroy/PycharmProjects/Main/whyGRNLinux/Results/3/Human/4-res.png',dpi = 300)
plt.show()

'''
max_mc = max(NMC.values())
V = [[] for i in range(int(max_mc) + 1)]


for u in G.nodes():
    V[int(NMC[u])].append(C[u])

X = []
Y = []
for i in range(len(V)):

    if len(V[i]) > 0:
        X.append(i)
        Y.append(np.average(V[i]))


print (len(X))

plt.scatter(X,Y)
#fit(X[:77],Y[:77],4)

plt.xlabel('Node Motif centrality Value')
plt.ylabel('Average Closeness Centrality')

plt.show()
'''