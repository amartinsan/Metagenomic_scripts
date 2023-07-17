#!/bin/bash
# Using qiime2.2022.11, remember to check updates and/or changes in code commands
#Load your enviromento or activate qiime2

  #source /share/apps/External/Miniconda3-4.7.10/bin/activate
  #export PATH=/share/apps/External/Miniconda3-4.6.14/envs/qiime2-2022.11/bin:$PATH
  #conda activate qiime2-2022.11

# Differential abundance is used to determine which features are significantly more/less abundant in different groups of samples


#Option 1: Correlation-clustering
qiime gneiss correlation-clustering \
  --i-table table.qza \
  --o-clustering hierarchy.qza
# metadata-column has to be a categoric variable

#qiime gneiss dendrogram-heatmap \
  --i-table table.qza \
  --i-tree hierarchy.qza \
  --m-metadata-file metadata.tsv \
  --m-metadata-column lugar \
  --p-color-map seismic \
  --o-visualization heatmap.qzv

#Option 2: Gradient-clustering

# gradient column must be a categoric variable

qiime gneiss gradient-clustering \
  --i-table table.qza \
  --m-gradient-file metadata.tsv \
  --m-gradient-column lugar \
  --o-clustering gradient-hierarchy.qza

qiime gneiss dendrogram-heatmap \
  --i-table table.qza \
  --i-tree hierarchy.qza \
  --m-metadata-file metadata.tsv \
  --m-metadata-column lugar \
  --p-color-map seismic \
  --o-visualization Gheatmap.qzv
