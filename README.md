## About GRAFree

GRAFree2 (GRaphical footprint based Alignment-Free method version 2) is a python based program for deriving the phylogenetic tree from the genomic sequence data of different species based on the graphical representation of the sequences. The program consists of deriving phylogenetic tree from the set of genomic sequences, generating bootstrap sequences, and generate the bootstrap trees from the bootstrap sequences.


## Required packages

The GRAFree is developed on python 3.6.10 and Linux based system. So this program is only compatible for python 3.6.10. The other prerequisites are following,
- Numpy 1.19.0
- Biopython v1.77
- Json 0.1.1


## Input directory

As input all the sequences should kept as a single file in FASTA format. GRAFree2 considers each sequence file as separate dataset.


## The GRAFree options

GRAFree can only be executed through the command prompt. It is necessary to make all the python file executable before running the program (use `chmod +x *.py`). The command to execute the GRAFree from the command prompt is as `./GRAFree.py` with the command line arguments. The options are as follows,

`--dataset`	Location of the dataset. Each dataset is a file having multiple of sequences in FASTA format. GRAFree2 can take multiple of such datasets placed in a directory.
`--outdir`	Location of the output. Inside --outdir, there may have multiple of directories for separate datasets. Each directory contains three sub directories -- GFP, Drift, and FeatVect.


## Execution

`python grafree2.py --dataset [dataset path] --outdir [output location]`






