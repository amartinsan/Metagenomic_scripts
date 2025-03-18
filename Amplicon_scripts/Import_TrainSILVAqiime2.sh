#### Para poder entrenar al classificador es necesario importar primero las secuencias a qiime2
## SILVA ARB 138
# Import Silva reference sequences
qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path silva_nr99_v138.1_wSpecies_train_set.fa.gz \
  --output-path silva-138-99-seqs.qza
#########################################################################
#########################################
###########
# Extract taxonomy from Silva FASTA headers (run this in bash)
# zcat silva_nr99_v138.1_wSpecies_train_set.fa.gz | grep "^>" | \
#  sed 's/^>//' | sed 's/ /\t/' > silva_taxonomy.tsv
###################################################
###############################
####

# Import Silva taxonomy
qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path silva_taxonomy.tsv \
  --output-path silva-138-99-tax.qza

###############-------->ESTA PARTE REQUIERE TIEMPO CORRER EN NOHUP O SERVIDOR<--------------- #########################

# Extract reference reads for your V3-V4 region
qiime feature-classifier extract-reads \
  --i-sequences silva-138-99-seqs.qza \
  --p-f-primer CCTACGGGNGGCWGCAG \
  --p-r-primer GACTACHVGGGTATCTAATCC \
  --p-trunc-len 0 \
  --o-reads silva-138-99-v3v4-region.qza

# Train the classifier with extracted V3-V4 reads
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads silva-138-99-v3v4-region.qza \
  --i-reference-taxonomy silva-138-99-tax.qza \
