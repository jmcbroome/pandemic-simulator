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
    print("Computing parsimony based on {} leaf assignments.".format(len(nrd)))
    parsimony_assigned = tree.simple_parsimony(nrd)
    return parsimony_assigned, tree, nrd

def collapse_polytomy_assignment(passd):
    '''
    Takes sets of resolved polytomy nodes (formatted as "node_XX_polytomy_YYY") 
    and returns a dictionary of combined nodes states by majority among their polytomy
    '''
    groups = {}
    fixed = {}
    for k,v in passd.items():
        nspl = k.split("_")
        if "polytomy" not in nspl:
            fixed[k] = v
            continue
        node = nspl[0] + "_" + nspl[1]
        if node not in groups:
            groups[node] = {}
        if v not in groups[node]:
            groups[node][v] = 0
        groups[node][v] += 1
    cd = {k:max(v, key=v.get) for k,v in groups.items()}
    cd.update(fixed)
    return cd

def assign_parsimony_clusters(pbf = 'sim.mat.collapsed.pb', regionf = 'simulated_regions.txt'):
    print("Computing parsimony scores...")
    bifurc_parsimony_assignments, colltree, original_labeling = get_parsimony_scoring(pbf,regionf)
    parsimony_assignments = collapse_polytomy_assignment(bifurc_parsimony_assignments)
    print(parsimony_assignments)
    tree = mat.MATree(pbf)

    pclusters = {}
    print("Calculating parsimony cluster assignments...")
    for leaf in tree.get_leaves_ids():
        #get its state.
        #iterate up towards the root until a parent without that state is encountered.
        lstate = original_labeling[leaf]
        print(leaf,lstate)
        for n in tree.rsearch(leaf):
            if not n.is_leaf():
                #disregard ambiguous assignments for cluster construction
                if parsimony_assignments.get(n.id,'0') != lstate: #and parsimony_assigned[n.name] != '-':
                    pclusters[leaf] = n.id
                    print('assigned to cluster ' + n.id)
                    break
    return pclusters,parsimony_assignments