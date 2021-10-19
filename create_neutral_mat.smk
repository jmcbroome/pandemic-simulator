#this workflow generates a completely simulated, simple SARS-CoV-2 MAT based on the config file and summarizes the results
#it does NOT incorporate any kind of selection, but can incorporate population-level parameters.
configfile: "config.yaml"
rule collapse_mat:
    input: 
        "sim.mat.pb"
    output: 
        "sim.mat.collapsed.pb"
    shell: 
        "matUtils extract -i {input} -o {output} -O"
        "matUtils summary -i {output}"

rule vgsim:
    output:
        "newick_output.nwk"
    shell:
        "python3 {config[vgexec]} {config[rt]} -it {config[it]} -s {config[samples]} -pm {config[ppmg]}.pp {config[ppmg]}.mg -su {config[sust]}.su -st {config[sust]}.st --createNewick"

rule phastsim:
    input:
        "newick_output.nwk"
    output:
        "sim.mat.pb"
    shell:
        "phastSim --output sim --reference {config[ref]} --scale {config[scale]} --createMAT --treeFile {input} --noHierarchy"