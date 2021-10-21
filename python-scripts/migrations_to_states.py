import ete3
events = {}
with open("migrations.txt") as inf:
    for entry in inf:
        spent = entry.strip().split()
        if spent[0] == "Node":
            continue
        events[spent[0]] = spent[-1]
tree = ete3.Tree("newick_output.nwk",format=1)

with open("simulated_regions.txt","w+") as outf:
    leafstates = {}
    for l in tree.get_leaves():
        key = l.name
        while l:
            if l.name in events:
                print(key + "\t" + events[l.name], file = outf)
                break
            else:
                l = l.up