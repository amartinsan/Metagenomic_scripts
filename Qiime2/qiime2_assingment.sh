#!/bin/bash
# Using qiime2.2022.11, remember to check updates and/or changes in code commands
#Load your enviromento or activate qiime2

  #source /share/apps/External/Miniconda3-4.7.10/bin/activate
  #export PATH=/share/apps/External/Miniconda3-4.6.14/envs/qiime2-2022.11/bin:$PATH
  #conda activate qiime2-2022.11

#In this script the asingment runs for Vsearchand SKlearn with all of the 16s. 
#Taxonomy asingment takes a long time so use a good amount of memory
#The really cool thing on qiime2 is you can train the sklearn algorithm with the primers or 16s rRNA regions.
  
  qiime feature-classifier classify-consensus-vsearch \
 --i-query rep-seqs.qza \
 --i-reference-reads $$PATH_TO/Silva_99OTUS_sequence.qza \
 --i-referencex-taxonomy $$PATH_TO/Silva_taxonomy.qza \
 --o-classification taxo_vsearch.qza \
 --p-n-threads $NSLOTS

qiime metadata tabulate \
  --m-input-file taxo_vsearch.qza \
  --o-visualization taxo_vsearch.qzv

qiime taxa barplot \
 --i-table table.qza \
 --i-taxonomy taxo_vsearch.qza \
 --m-metadata-file metadata.tsv \
 --o-visualization vsearch_barplots.qzv

#Sklearn 

qiime feature-classifier classify-sklearn \
 --i-classifier $$PAth_to/silva-132-99-nb-classifier.qza \
 --i-reads rep-seqs.qza \
 --o-classification taxo_sklearn99.qza \

qiime metadata tabulate \
  --m-input-file taxo_sklearn99.qza \
  --o-visualization taxo_sklearn99.qza

qiime taxa barplot \
 --i-table table.qza \
 --i-taxonomy taxo_sklearn99.qza \
 --m-metadata-file metadata.tsv \
 --o-visualization ske99_barplots.qzv
