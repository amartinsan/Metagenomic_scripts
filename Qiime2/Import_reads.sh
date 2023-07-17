#!/bin/bash

# Using qiime2.2022.11, remember to check updates and/or changes in code commands

#Load your enviromento or activate qiime2

  #source /share/apps/External/Miniconda3-4.7.10/bin/activate
  #export PATH=/share/apps/External/Miniconda3-4.6.14/envs/qiime2-2022.11/bin:$PATH
  #conda activate qiime2-2022.11

#For single-end reads
  qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path manifesto.csv \
  --output-path demux.qza \
  --input-format SingleEndFastqManifestPhred33

#For pair-end reads

  qiime tools import \
  --type SampleData[PairedEndSequencesWithQuality] \
  --input-format PairedEndFastqManifestPhred33 \
  --input-path manifesto.csv\
  --output-path demux.qza


#Using qiime to cut adapters from de seqs, depends on your data

#qiime cutadapt trim-paired \
#--i-demultiplexed-sequences untrimmed_demux.qza \
#--p-cores $NSLOTS \
#--p-front-f CCTACGGGNGGCWGCAG \
#--p-front-r GACTACHVGGGTATCTAATCC \
#--o-trimmed-sequences demux.qza

