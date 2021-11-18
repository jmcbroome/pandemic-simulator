# pandemic-simulator

This repository contains snakemake workflows and parameter files for the simulation of pandemic-scale data, particularly of viruses which resemble SARS-CoV-2. These are intended to support future pandemic preparedness and provide simulated ground truth for the validation and testing of statistics and analyses on real SARS-CoV-2 data.

Dependencies include:

[VGsim](https://github.com/Genomics-HSE/VGsim)

[phastSim](https://github.com/NicolaDM/phastSim)

[matUtils](https://github.com/yatisht/usher)

[Snakemake](https://snakemake.readthedocs.io/en/stable/)

And the python packages:

pandas, numpy, ete3, scikit-learn

Ensure all dependencies are installed and available on your path.

## The Config File
Our configuration file contains all relevant parameters, executable path information, and paths to files containing parameters as needed. Each snakemake pipeline relies on this configuration file to determine its basic characteristics.

## Snakefiles

### create_neutral_mat.smk

This is the simplest pipeline. It applies VGSim and phastSim to generate a single simulated MAT output.

#### Simulation Parameters Used

VGsim: rates (rt), iterations (it) and samples, migration and populations (ppmg), suspectibility (sust)

phastSim: scale, reference (ref)

### simulate_introductions.smk

This pipeline is used to simulate and validate the results of the "matUtils introduce" heuristic for phylogeographic state identification. 

#### Simulation Parameters Used

VGsim: rates (rt), iterations (it) and samples, migration and populations (ppmg), suspectibility (sust)

phastSim: scale, reference (ref)
