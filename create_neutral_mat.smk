#this workflow generates a completely simulated, simple SARS-CoV-2 MAT based on the config file and summarizes the results
#it does NOT incorporate any kind of selection, but can incorporate population-level parameters.
configfile: "config.yaml"
rule all:
    input:
        "sim.mat.collapsed.pb"

rule collapse_mat:
    input: 
        "sim.mat.pb"
    output: 
        "sim.mat.collapsed.pb"
    shell: 
        "matUtils extract -i {input} -o {output} -O"

rule vgsim:
    output:
        "newick_output.nwk",
        "migrations.txt"
    shell:
        "python3 {config[executables][vgexec]} {config[vgsim-params][rt]} -it {config[vgsim-params][it]} -s {config[vgsim-params][samples]} -pm {config[vgsim-params][ppmg]}.pp {config[vgsim-params][ppmg]}.mg -su {config[vgsim-params][sust]}.su -st {config[vgsim-params][sust]}.st --createNewick --writeMigrations"

rule phastsim:
    input:
        "newick_output.nwk"
    output:
        "sim.mat.pb"
    shell:
        "phastSim --output sim --reference {config[phastsim-params][ref]} --scale {config[phastsim-params][scale]} --createMAT --treeFile {input}"