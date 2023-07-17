#!/bin/bash
# Using qiime2.2022.11, remember to check updates and/or changes in code commands
#Load your enviromento or activate qiime2

  #source /share/apps/External/Miniconda3-4.7.10/bin/activate
  #export PATH=/share/apps/External/Miniconda3-4.6.14/envs/qiime2-2022.11/bin:$PATH
  #conda activate qiime2-2022.11
  
#ANCOM can be applied to identify features that are differentially abundant (i.e. present in different abundances) across sample groups. 
#--pwhere filters table for only selected data

qiime feature-table filter-samples \
  --i-table table.qza \
  --m-metadata-file metadata.tsv \
  --p-where "[place]='zone3'" #In this example varible name place, variable zone3 \
  --o-filtered-table filter-table.qza

qiime composition add-pseudocount \
  --i-table filter-table.qza \
  --o-composition-table comp-filter-table.qza

 metadata column tiene que ser categorico
qiime composition ancom \
  --i-table comp-filter-table.qza \
  --m-metadata-file metadata.tsv \
  --m-metadata-column lugar \
  --o-visualization ancom-subject.qzv
