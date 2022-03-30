import pandas as pd
import subprocess
import time
from yaml import load, dump
import yaml
try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper

cyml = open('config.yaml')
config = yaml.load(cyml,Loader=Loader)
cyml.close()

def adjust_migration(nr, pref = 'doublepop'):
    #this will only set a two-population even migration rate matrix.
    with open("parameter-files/" + pref + ".mg",'w+') as outf:
        print("#Migration_format_version 0.0.1",file=outf)
        print(str(1-nr) + " " + str(nr),file=outf)
        print(str(nr) + " " + str(1-nr),file=outf)
        
def update_parameters(scale = 0.0009, migration = 0.01):
    adjust_migration(migration)
    config['phastsim-params']['scale'] = str(scale)
    nf = open('config.yaml','w')
    print(yaml.dump(config),file=nf)
    nf.close()

def read_summary_out(f = "summary_output.txt"):
    sd = {}
    with open(f) as inf:
        for entry in inf:
            k,v = entry.strip().split(":")
            v = float(v)
            sd[k] = v
    return sd

def read_results(f = "results.txt"):
    rd = {}
    tm = []
    ptm = []
    ln = 0
    with open(f) as inf:
        for entry in inf:
            if ln == 0:
                rd["Internal Assignment"] = float(entry.strip().split(":")[-1])
            elif ln == 1:
                rd["ARI"] = float(entry.strip().split(":")[-1])
            elif ln == 5 or ln == 6:
                s = entry.strip().split()
                tm.append([int(s[1]),int(s[2])])
            elif ln == 8:
                rd['Parsimony Internal Assignment'] = float(entry.strip().split(":")[-1])
            elif ln == 10:
                rd["ParsimonyARI"] = float(entry.strip().split(":")[-1])
            elif ln == 13 or ln == 14:
                s = entry.strip().split()
                ptm.append([int(s[1]),int(s[2])])
            ln += 1
    rd['ConfusionMatrix'] = tm
    rd['ParsimonyConfusionMatrix'] = ptm
    return rd

def collect_results():
    subprocess.check_call("matUtils summary -i sim.mat.collapsed.pb > summary_output.txt",shell=True)
    sd = read_summary_out()
    rd = read_results()
    sd.update(rd)
    return sd

def do_simulation_run(scale, migration):
    update_parameters(scale, migration)
    subprocess.check_call("snakemake --forceall -s simulate_introductions.smk -c1",shell=True)
    results = collect_results()
    return results

simmat = {k:[] for k in ["Scale","MigRate","NodesCollapsed","MutationsPerNode","MutationParsimony","ARI","TreeDepth","IAC", "ParsimonyIAC", "ParsimonyARI"]}
for s in [1e-4,5e-4,1e-3,5e-3]:
    try:
        for mr in [0.01,0.05,0.1,0.25]:
            print(s,mr)
            stime = time.time()
            try:
                rd = do_simulation_run(s,mr)
            except:
                print("Failed?")
                try:
                    rd = do_simulation_run(s,mr)
                except:
                    print("Failed twice.")
                    try:
                        rd = do_simulation_run(s,mr)
                    except:
                        print("Third failure- giving up on this parameter set.")
                        print(s,mr)
                        continue
            nc = 1999999 - rd['Total Nodes in Tree']
            mpn = rd['Total Tree Parsimony']/1999999
            p = rd['Total Tree Parsimony']
            ari = rd['ARI']
            td = rd['Mean Tree Depth']
            iac = rd['Internal Assignment']
            piac = rd['Parsimony Internal Assignment']
            pari = rd['ParsimonyARI']

            simmat['Scale'].append(s)
            simmat['MigRate'].append(mr)
            simmat['NodesCollapsed'].append(nc)
            simmat['MutationsPerNode'].append(mpn)
            simmat['MutationParsimony'].append(p)
            simmat['ARI'].append(ari)
            simmat['TreeDepth'].append(td)
            simmat['IAC'].append(iac)
            simmat['ParsimonyIAC'].append(piac)
            simmat['ParsimonyARI'].append(pari)
            print("completed in ", time.time() - stime)
    except KeyboardInterrupt:
        print("Interrupt requested, saving current results")
        break
smdf = pd.DataFrame(simmat)
smdf.to_csv("simulations_r13.tsv",sep='\t')
