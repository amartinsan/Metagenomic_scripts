#!/bin/bash

Export PATH or use conda
#export PATH=/scratch/share/apps/miniconda2/bin/:$PATH
#source activate graftm
#conda activate graftm

# You need the taxnonomy and sequences to make a package (.gpkg)

graftM create --sequences sequences_create.fasta \
              --taxonomy taxonomy_create.tsv \
              --threads 1 \
              --output example.gpkg

#source deactivate
