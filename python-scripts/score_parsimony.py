#this implementation uses the MAT Cython API implementation of the Fitch algorithm, which is substantially
#more efficient than biopython. https://github.com/jmcbroome/matreePy
import mat

def get_parsimony_scoring(pbf, regionf):
    tree = mat.MATree(pbf)
    nrd = {}
    with open(regionf) as inf:
        for entry in inf:
            node, region = entry.strip().split()
            nrd[node] = region
    parsimony_assigned = tree.simple_parsimony(nrd)
    return parsimony_assigned, tree

def assign_parsimony_clusters(pbf = 'sim.mat.collapsed.pb', regionf = 'simulated_regions.txt'):
    print("Computing parsimony scores...")
    parsimony_assigned, tree = get_parsimony_scoring(pbf,regionf)
    pclusters = {}
    print("Calculating parsimony cluster assignments...")
    for leaf in tree.get_leaves_ids():
        #get its state.
        #iterate up towards the root until a parent without that state is encountered.
        lstate = parsimony_assigned[leaf]
        for n in tree.rsearch(leaf):
            if not n.is_leaf():
                #disregard ambiguous assignments for cluster construction
                if parsimony_assigned[n.id] != lstate: #and parsimony_assigned[n.name] != '-':
                    pclusters[leaf] = n.id
                    break
    return pclusters,parsimony_assigned