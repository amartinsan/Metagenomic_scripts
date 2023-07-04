##mapuniprot2tax.sh##

python /$$$PATH$$/pair_seq_to_accesion.py \
      --fasta sequences_cdhit100.fasta \
      --mapping map_uniprot2accession.tsv

#load local libraries
if [ -n $R_LIBS ]; then
   export R_LIBS=/$$PATH$$$//local_R_library:$R_LIBS
else
   export R_LIBS=/$$PATH$$$/local_R_library
fi

# with accesion2taxonomy we get the taxonomy of each sequence
Rscript $$$PATH$$$/accession2taxonomy.R \
       sequences.fasta \
       mapa.txt
# Change the name of  the sequence file
mv out_cl.fas sequences_create.fasta
#Merge both accesion and taxonomy into one file
paste out_seqid_cl.txt out_taxonomy_cl.txt > taxonomy_create.tsv
