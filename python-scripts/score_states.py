import pandas as pd
import numpy as np
import ete3
import glob

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
stree = ete3.Tree("sim.mat.nwk",1)
#the trees are equivalent, but matUtils doesn't retain internal node IDs, so we need to traverse these equivalent trees and 
#match up their names in traversal order
tnames = [n.name for n in tree.traverse()]
snames = [n.name for n in stree.traverse()]
assert len(tnames) == len(snames)

#for each node, traverse upwards from it until it encounters the latest event
#and take that state. Not a particularly efficient algorithm, but sufficient here.
nodestates = {}
for l in tree.traverse():
    key = l.name
    assigned = False
    while l:
        if l.name in events:
            nodestates[key] = events[l.name]
            assigned = True
            break
        else:
            l = l.up
    if not assigned:
        nodestates[key] = '0'

#collect the confidences of internal nodes inferred by matUtils.
dfvs = []
for rf in glob.glob("assignments_out/*"):
    region = rf.split("/")[-1].split("_")[0]
    tdf = pd.read_csv(rf,sep='\t')
    tdf[region + "_confidence"] = tdf['confidence_continuous']
    dfvs.append(tdf)
assdf = dfvs[0]
for sd in dfvs[1:]:
    assdf = assdf.merge(sd,on="sample",how='outer')
assdf.replace(np.nan,0,inplace=True)

#apply confidence scores to the state data.
nd = {sn:nodestates[tnames[i]] for i,sn in enumerate(snames)}
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
#save useful results.
with open("results.txt","w+") as outf:
    print("Percentage of internal nodes correctly assigned overall: " + str(cvc[True]), file=outf)
    print("Confusion Matrix", file=outf)
    print("________________", file=outf)
    print("Matrix\t" + "\t".join(["PredictedState=" + s for s in states]), file = outf)
    for ts in states:
        pred_states = assdf[(assdf.TrueState == ts)].PredState.value_counts()
        prow = "TrueState=" + ts + "\t" + "\t".join([str(pred_states[s]) for s in states])
        print(prow, file = outf)
