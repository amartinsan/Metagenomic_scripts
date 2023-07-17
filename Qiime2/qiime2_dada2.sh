#!/bin/bash
# Using qiime2.2022.11, remember to check updates and/or changes in code commands
#Load your enviromento or activate qiime2

  #source /share/apps/External/Miniconda3-4.7.10/bin/activate
  #export PATH=/share/apps/External/Miniconda3-4.6.14/envs/qiime2-2022.11/bin:$PATH
  #conda activate qiime2-2022.11

#Normaly dada2 is kinda strict with the denoising and you can lose a lot of sequences if you did not cut adapters or have low overlap in your reads.

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs demux.qza \
  --p-trim-left-f $n \
  --p-trim-left-r $n \
  --p-trunc-len-f $n \
  --p-trunc-len-r $n \
  --o-table table.qza \
  --o-representative-sequences rep-seqs.qza \
  --o-denoising-stats denoising-stats.qza

  #Some usefull options for not losing sequences are:
    
    --p-min-overlap $n #Determine min overlap CAN INCREASE CHIMERAS
    --p-min-fold-parent-over-abundance $n minimal potential parents CAN INCREASE CHIMERAS and false negatives
    --p-max-ee-f $n or --p-max-ee-R $n #number of expected error in error rate normaly 2 but 3 is not a bad number
    #Remember only using these if you have a really low number of denoised reads.
