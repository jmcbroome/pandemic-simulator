#use parsimony to attempt to recapture true cluster states.
from Bio import Phylo
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from tqdm import tqdm
import subprocess

def make_msa(regiondict):
    #convert regions to a length one sequence msa for scoring.
    conv = {'0':"A",'1':"G"}
    msav = []
    for n,r in regiondict.items():
        sr = SeqRecord(Seq(conv[r]),id=n)
        msav.append(sr)
    return MultipleSeqAlignment(msav)

def get_clade_scores(tree,alignment):
    #literally just the scoring function from biopython but instead of adding everything at the end to a useless int it actually returns assignments.
    #biopython...
    if not tree.is_bifurcating():
        raise ValueError("The tree provided should be bifurcating.")
    if not tree.rooted:
        tree.root_at_midpoint()
    # sort tree terminals and alignment
    terms = tree.get_terminals()
    terms.sort(key=lambda term: term.name)
    alignment.sort()
    if not all(t.name == a.id for t, a in zip(terms, alignment)):
        raise ValueError(
            "Taxon names of the input tree should be the same with the alignment."
        )
    # term_align = dict(zip(terms, alignment))
    score = 0
    for i in range(len(alignment[0])):
        # parsimony score for column_i
        score_i = 0
        # get column
        column_i = alignment[:, i]
        # skip non-informative column
        if column_i == len(column_i) * column_i[0]:
            continue

        # start calculating score_i using the tree and column_i

        # Fitch algorithm without the penalty matrix
        # init by mapping terminal clades and states in column_i
        clade_states = dict(zip(terms, [{c} for c in column_i]))
        for clade in tqdm(tree.get_nonterminals(order="postorder")):
            clade_childs = clade.clades
            left_state = clade_states[clade_childs[0]]
            right_state = clade_states[clade_childs[1]]
            state = left_state & right_state
            if not state:
                state = left_state | right_state
                score_i = score_i + 1
            clade_states[clade] = state
    return clade_states

def get_parsimony_scoring(newickf = 'sim.mat.nwk', regionf = 'simulated_regions.txt'):
    tree = Phylo.read(newickf,'newick')
    nrd = {}
    with open(regionf) as inf:
        for entry in inf:
            node, region = entry.strip().split()
            nrd[node] = region
    msa = make_msa(nrd)
    csd = get_clade_scores(tree,msa) 
    rconv = {"A":'0',"G":'1'}
    parsimony_assigned = {}
    for k,v in csd.items():
        if len(v) == 1:
            parsimony_assigned[k.name] = rconv[list(v)[0]]
        else:
            #ambiguous. Not sure if this is actually the right way to interpet a fitch-parsimony set of states but whatever.
            parsimony_assigned[k.name] = '-'
    return parsimony_assigned, tree

def do_resolution(pbf="sim.mat.collapsed.pb"):
    print("Resolving polytomies...")
    subprocess.check_call('matUtils extract -R -i '+pbf+' -o '+pbf[:-3]+'.resolved.pb',shell=True)
    subprocess.check_call('matUtils extract -i '+pbf[:-3]+'.resolved.pb -t '+pbf[:-3]+'.resolved.nwk',shell=True)

def assign_parsimony_clusters(newickf = 'sim.mat.nwk', regionf = 'simulated_regions.txt'):
    print("Computing parsimony scores...")
    parsimony_assigned, tree = get_parsimony_scoring(newickf,regionf)
    #to get the clusters, for each leaf, iterate back from the tips until an internal node is encountered with a different state
    #then group them by the internal nodes
    pclusters = {}
    print("Calculating parsimony cluster assignments...")
    for leaf in tqdm(tree.get_terminals()):
        #get its state.
        #iterate up towards the root until a parent without that state is encountered.
        l = leaf
        lstate = parsimony_assigned[leaf.name]
        path = tree.get_path(leaf.name)
        for n in path[::-1]:
            if not n.is_terminal():
                #disregard ambiguous assignments for cluster construction
                if parsimony_assigned[n.name] != lstate: #and parsimony_assigned[n.name] != '-':
                    pclusters[leaf.name] = n.name
                    break

        # while l:
            # print(dir(l))
            # l = l.parent
            # if parsimony_assigned[l.name] != lstate:
            #     # if l not in pclusters:
            #     #     pclusters[l] = set()
            #     # pclusters[l].add(leaf)
            #     pclusters[leaf.name] = l.name
            #     break
    return pclusters,parsimony_assigned