#!/bin/bash

# Using qiime2.2022.11, remember to check updates and/or changes in code commands
#Load your enviromento or activate qiime2
  #source /share/apps/External/Miniconda3-4.7.10/bin/activate
  #export PATH=/share/apps/External/Miniconda3-4.6.14/envs/qiime2-2022.11/bin:$PATH
  #conda activate qiime2-2022.11

#Qiime2 makes de .qza (a kinda zip) a visualization. remember to have a metadata wiht usefull information and the same sample names

qiime feature-table tabulate-seqs \
 --i-data rep-seqs.qza \
 --o-visualization rep-seqs.qzv

qiime metadata tabulate \
 --m-input-file denoising-stats.qza \
 --o-visualization denoising-stats.qzv
 
#You can also fitler ASVs based on frequence or more

qiime feature-table filter-features \
  --i-table table.qza \
  --p-min-frequency 10 \
  --o-filtered-table fil-table.qza

qiime feature-table summarize  \
 --i-table fil-table.qza  \
 --o-visualization fil-table.qzv \
 --m-sample-metadata-file metadata.tsv

qiime feature-table summarize  \
 --i-table table.qza  \
 --o-visualization table.qzv \
 --m-sample-metadata-file metadata.tsv
