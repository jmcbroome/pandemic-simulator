configfile: "config.yaml"
include: "create_neutral_mat.smk"

rule score_states:
    input:
        "simulated_fullout.tsv"
        "sim.mat.nwk"
        "newick_output.nwk"
        "migrations.txt"
    output:
        "results.txt"
    script:
        "python-scripts/score_states.py"

rule define_states:
    input:
        "migrations.txt"
    output:
        "simulated_regions.txt"
    script:
        "python-scripts/migrations_to_states.py"

rule run_introduce:
    input:
        "simulated_regions.txt"
        "sim.mat.pb"
    output:
        "simulated_fullout.tsv"
    shell:
        "matUtils introduce -i {input[1]} -s {input[0]} -o {output} -D assignments_out"

rule get_sim_nwk:
    input:
        "sim.mat.pb"
    output:
        "sim.mat.nwk"
    shell:
        "matUtils extract -i {input} -t {output}"
