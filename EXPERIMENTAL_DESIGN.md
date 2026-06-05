# Aims

1. To understand the interface between the W6/32 antibody and the HLA-B*27:05 molecule and beta-2-microglobulin in the experimentally determined structure

2. To understand how that interface might change between HLA-B*27:05 and HLA-A*02:01 in the PRAME (SLLQHLIGL) expreimental structure

3. To understand how similar the prediction of the HLA-A*02:01 / PRAME single chain trimer is to the experimental structure 

4. To map the W6/32 interface onto the HLA-A*02:01 / PRAME single chain trimer prediction for protein design

# Information

W6/32 is a pan HLA conformationally sensitive antibody. We're focussing on the heavy chain as the light chain is sterically hindered by the peptide -> beta2m linker in the single chain trimer format

PRAME is an HLA-A*02:01 restricted epitope. Its sequence is SLLQHLIGL

All the structures are aligned in the same co-ordinate frame.

The predictions are not yet aligned into the co-ordinate frame. The canonical structure to align to is in the 'structures' folder (it is the antigen binding domain of the 1HHK structure).

# Conventions

1. Each tool/skill must be self contained so that it can be split out into a separate Python library or skill file. We want to build components that can be composed, rather than a bespoke pipeline. 

2. Each tool should take an input filepath as an argument and an output folder as an argument

3. The approach taken should always be:

Plan first, with the plan appended to PLAN.md

Then execute

Then update the CLAUDE.md file and a CHANGELOG.md file which should act as a lab notebook

On approval of a piece of work, commit the work to the Git repository after I've approved the commit message


# Jobs to do

1. Build a tool/skill which takes a PDB file and given a set of chains for the epitope and a set of chains for the paratope, generates a JSON file detailing the interactions. 

2. Build a tool/skill which takes a PDB file and aligns it to a canonical structure. 

3. Build a tool/skill which analyses the fine scale differences between two structures/predictions

e.g. how similar is one experimental structure with another
e.g. how good is the prediction when compared to a ground truth structure

4. Use the tools/skills to understand design questions which are:

a. How much does W6/32 change the conformation of the HLA-B*27:05 molecule on binding? Compare the W6/32 bound structure to the HLA-B*27:05 structure with the RRFSRSPIRR peptide. Focus on the W6/32 interface area only (we're not interested in changes related to the different peptides). In your plan, state which skills you'll use for this. The experimental structures are already in the same co-ordinate frame.

b. Can we have confidence to use the Boltz2 prediction of the HLA-A*02:01/PRAME single chain trimer as a design input to protein design. Is it a prediction that has strong similarity to the HLA-A*02:01/PRAME experimental structure? Is it good in general? Is the peptide conformation similar? Is the W6/32 interface well modelled?

c. Extension task. If we were to introduce an additional disulphide bridge between the alpha3 domain of the MHC molecule and the linker that joins beta2m to the alpha1 domain, where would the best location for that be?


# Libraries to be used

BioPython - for structure alignment and manipulation
PDBE Arpeggio - for assessing contacts - https://github.com/PDBeurope/arpeggio
BioPandas - for generating RMSDs etc

# Environment 

Python 3.14
uv for package management



