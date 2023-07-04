#!/bin/bash


#export PATH=/scratch/share/apps/miniconda2/bin/:$PATH
#source activate graftm
#conda activate graftm

# Uses a package (.gpkg) to search homologs in  some sequences


graftM graft \
     --forward /$$$$$$$PATH/*****_R1.fastq \
     --reverse /$$$$$$$PATH/****_R2.fastq \
     --graftm_package PAQUETE_DESEADO.gpkg \
     #--min_orf_length  \
     #--filter_minimum  \
     #--evalue \
     #--threads $NSLOTS \
     #--filter_minimum  \

#source deactivate
