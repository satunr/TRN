from networkx.algorithms.connectivity import build_auxiliary_edge_connectivity
from networkx.algorithms.connectivity import build_auxiliary_node_connectivity

import networkx as nx

from networkx.exception import NetworkXNoPath
# Define the default maximum flow function to use for the undelying
# maximum flow computations
from networkx.algorithms.flow import edmonds_karp
from networkx.algorithms.flow import preflow_push
from networkx.algorithms.flow import shortest_augmenting_path
default_flow_func = edmonds_karp
# Functions to build auxiliary data structures.
from networkx.algorithms.flow import build_residual_network

try:
    from itertools import filterfalse as _filterfalse
except ImportError:  # Python 2
    def _filterfalse(predicate, iterable):
        # https://docs.python.org/3/library/itertools.html
        # filterfalse(lambda x: x%2, range(10)) --> 0 2 4 6 8
        if predicate is None:
            predicate = bool
        for x in iterable:
            if not predicate(x):
                yield x

__all__ = [
    'edge_disjoint_paths',
    'node_disjoint_paths',
]

def _unique_everseen(iterable):
    # Adapted from https://docs.python.org/3/library/itertools.html examples
    "List unique elements, preserving order. Remember all elements ever seen."
    # unique_everseen('AAAABBBCCDAABBB') --> A B C D
    seen = set()
    seen_add = seen.add
    for element in _filterfalse(seen.__contains__, iterable):
        seen_add(element)
        yield element

def node_disjoint_paths(G, s, t, flow_func=None, cutoff=None, auxiliary=None,
                        residual=None):

    if s not in G:
        raise nx.NetworkXError('node %s not in graph' % s)
    if t not in G:
        raise nx.NetworkXError('node %s not in graph' % t)

    if auxiliary is None:
        H = build_auxiliary_node_connectivity(G)
    else:
        H = auxiliary

    mapping = H.graph.get('mapping', None)
    if mapping is None:
        raise nx.NetworkXError('Invalid auxiliary digraph.')

    # Maximum possible edge disjoint paths
    possible = min(H.out_degree('%sB' % mapping[s]),
                   H.in_degree('%sA' % mapping[t]))
    if not possible:
        raise NetworkXNoPath

    if cutoff is None:
        cutoff = possible
    else:
        cutoff = min(cutoff, possible)

    kwargs = dict(flow_func=flow_func, residual=residual, auxiliary=H,
                  cutoff=cutoff)

    # The edge disjoint paths in the auxiliary digraph correspond to the node
    # disjoint paths in the original graph.
    paths_edges = edge_disjoint_paths(H, '%sB' % mapping[s], '%sA' % mapping[t],
                                      **kwargs)
    for path in paths_edges:
        # Each node in the original graph maps to two nodes in auxiliary graph
        yield list(_unique_everseen(H.node[node]['id'] for node in path))



def edge_disjoint_paths(G, s, t, flow_func=None, cutoff=None, auxiliary=None,
                        residual=None):
    if s not in G:
        raise nx.NetworkXError('node %s not in graph' % s)
    if t not in G:
        raise nx.NetworkXError('node %s not in graph' % t)

    if flow_func is None:
        flow_func = default_flow_func

    if auxiliary is None:
        H = build_auxiliary_edge_connectivity(G)
    else:
        H = auxiliary

    # Maximum possible edge disjoint paths
    possible = min(H.out_degree(s), H.in_degree(t))
    if not possible:
        raise NetworkXNoPath

    if cutoff is None:
        cutoff = possible
    else:
        cutoff = min(cutoff, possible)

    # Compute maximum flow between source and target. Flow functions in
    # NetworkX return a residual network.
    kwargs = dict(capacity='capacity', residual=residual, cutoff=cutoff,
                  value_only=True)
    if flow_func is preflow_push:
        del kwargs['cutoff']
    if flow_func is shortest_augmenting_path:
        kwargs['two_phase'] = True
    R = flow_func(H, s, t, **kwargs)

    if R.graph['flow_value'] == 0:
        raise NetworkXNoPath

    # Saturated edges in the residual network form the edge disjoint paths
    # between source and target
    cutset = [(u, v) for u, v, d in R.edges(data=True)
              if d['capacity'] == d['flow'] and d['flow'] > 0]
    # This is equivalent of what flow.utils.build_flow_dict returns, but
    # only for the nodes with saturated edges and without reporting 0 flows.
    flow_dict = {n: {} for edge in cutset for n in edge}
    for u, v in cutset:
        flow_dict[u][v] = 1

    # Rebuild the edge disjoint paths from the flow dictionary.
    paths_found = 0
    for v in list(flow_dict[s]):
        if paths_found >= cutoff:
            # preflow_push does not support cutoff: we have to
            # keep track of the paths founds and stop at cutoff.
            break
        path = [s]
        if v == t:
            path.append(v)
            yield path
            continue
        u = v
        while u != t:
            path.append(u)
            try:
                u, _ = flow_dict[u].popitem()
            except KeyError:
                break
        else:
            path.append(t)
            yield path
            paths_found += 1

#G = nx.DiGraph()
#G.add_edges_from([(1,2), (1,3), (2,4), (3,4), (4,5), (4,6), (5,7), (6,7)])
#print (list(node_disjoint_paths(G, 1, 7)))
