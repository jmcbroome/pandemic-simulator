import pandas as pd
import numpy as np
import ete3
import glob
from sklearn.metrics import adjusted_rand_score

events = {}
with open("migrations.txt") as inf:
    for entry in inf:
        spent = entry.strip().split()
        if spent[0] == "Node":
            continue
        events[spent[0]] = spent[-1]
nodestates = {}
trueclusters = {}

tree = ete3.Tree("newick_output.nwk",format=1)
stree = ete3.Tree("sim.mat.nwk",format=1)
#the trees are equivalent, but matUtils doesn't retain internal node IDs, so we need to traverse these equivalent trees and 
#match up their names in traversal order
tnames = [n.name for n in tree.traverse()]
snames = [n.name for n in stree.traverse()]
assert len(tnames) == len(snames)

#for each node, traverse upwards from it until it encounters the latest event
#and take that state. Not a particularly efficient algorithm, but sufficient here.
#additionally, record clusters of samples descended from distinct migration events.
nodestates = {}
counter = {}
for l in tree.traverse():
    key = l.name
    checking_leaf = l.is_leaf()
    assigned = False
    while l:
        if l.name in events:
            nodestates[key] = events[l.name]
            if checking_leaf:
                if events[l.name] not in counter:
                    counter[events[l.name]] = 0
                counter[events[l.name]] += 1

                if l.name not in trueclusters:
                    trueclusters[l.name] = []
                trueclusters[l.name].append(key)
            assigned = True
            break
        else:
            l = l.up
    if not assigned:
        nodestates[key] = '0'
        if checking_leaf:
            if '0' not in counter:
                counter['0'] = 0
            counter['0'] += 1
#collect the confidences of internal nodes inferred by matUtils.
print(counter)
dfvs = []
for rf in glob.glob("assignments_out/*"):
    region = rf.split("/")[-1].split("_")[0]
    tdf = pd.read_csv(rf,sep='\t')
    tdf.rename(columns={'confidence_continuous':region + "_confidence"},inplace=True)
    dfvs.append(tdf)
assdf = dfvs[0]
for sd in dfvs[1:]:
    assdf = assdf.merge(sd,on="sample",how='outer')
assdf.replace(np.nan,0,inplace=True)
# print("DEBUG:", assdf.shape[0], len(snames))
# print(trueclusters)    
nd = {sn:nodestates[tnames[i]] for i,sn in enumerate(snames)}
counter = {}
for k,v in nd.items():
    if k in assdf['sample']:
        if v not in counter:
            counter[v] = 0
        counter[v] += 1
print(counter)
assdf['TrueState'] = assdf['sample'].apply(lambda x:nd[x])
#remove leaves, as they have guaranteed correctness and are noninformative.
assdf = assdf[assdf['sample'].apply(lambda x:"node_" in x)]
states = assdf.TrueState.value_counts().index

def get_predicted(d):
    for s in states:
        pv = d[s + "_confidence"]
        if pv > 0.5:
            return s
assdf['PredState'] = assdf.apply(get_predicted,axis=1)
cvc = (assdf.PredState == assdf.TrueState).value_counts(normalize=True)
print(assdf)

#now, we calculate the Adjusted Rand Index of cluster assignments for all simulated samples which are descended from at least one transition
#in to or out of a region in the data.
idf = pd.read_csv("simulated_fullout.collapsed.tsv",sep="\t").dropna()
iv = {}
for k,v in trueclusters.items():
    for sv in v:
        iv[sv]= k 
if len(iv) == 0:
    ari = 1 #not particularly meaningful. No migrations = no samples from the other region = no introductions = nothing to detect = can't be wrong??
else:
    idf['TC'] = idf['sample'].apply(lambda x:iv.get(str(x),np.nan))
    idf = idf.dropna()
    assert (idf.shape[0] > 0)
    ari = adjusted_rand_score(idf.introduction_node, idf.TC)

#save useful results.
with open("results.txt","w+") as outf:
    print("Percentage of internal nodes correctly assigned overall to full tree: " + str(cvc[True]), file=outf)
    print("Adjusted Rand Index of cluster labels on collapsed tree: " + str(ari), file = outf)
    print("Confusion Matrix of internal nodes", file=outf)
    print("________________", file=outf)
    print("Matrix\t" + "\t".join(["PredictedState=" + s for s in states]), file = outf)
    for ts in states:
        pred_states = assdf[(assdf.TrueState == ts)].PredState.value_counts()
        prow = "TrueState=" + ts + "\t" + "\t".join([str(pred_states[s]) for s in states])
        print(prow, file = outf)
