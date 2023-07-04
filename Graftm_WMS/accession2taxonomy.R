## R script
## This script takes two files as input: the first file contains two columns, a sequence_id column
## and an accession_number column; the second file is a fasta sequence file, where the sequence_id
## corresponds to the first file.
## The output consists of two files: the first file contains a table composed of sequence_id
## and lineage columns in green_genes format (k_;p_;c_;l_;f_;g_;s_) only for those accession numbers
## that return a valid taxa_id. The second file is the fasta file that contains the subset of
## sequences where a lineage was found.

## The script uses the taxonomizr library
## https://cran.r-project.org/web/packages/taxonomizr/README.html


## Check if the "nodes.dmp", "names.dmp", and "accessionTaxa.sql" files exist
###################################


#read fasta sequence
fasta_seqs<-read.fasta(file = args[1], as.string=TRUE, seqtype="AA" )

#read nodes and names from NCBI taxonomy
#taxaNodes<-read.nodes.sql('/$$$RUTA DEL ARCHIVO$$$/nodes.dmp')
#taxaNames<-read.names.sql('/$$$RUTA DEL ARCHIVO$$$/names.dmp')

#get taxa_id for each accesion number, the result is a vector of integers
taxaId<-accessionToTaxa(accessions[[2]],"/$$$RUTA DEL ARCHIVO$$$/accessionTaxa.sql", version="base")
taxaId
#get linage from taxa_id, the resutl is a matrix with seven columns k p c o f g s
lineage<-getTaxonomy(taxaId,"/$$$RUTA DEL ARCHIVO$$$/accessionTaxa.sql")

#logical subsetting, retrieve only those records with a valid accesion_number
new_fasta_seqs<-fasta_seqs[!is.na(x=taxaId)]
write.fasta(sequences=new_fasta_seqs, names=names(new_fasta_seqs), file.out="out.fas")
new_seqid<-accessions[[1]][!is.na(x=taxaId)]
lineage.frame <- as.data.frame(x=lineage, optional=TRUE, make.names=NA)
final_lineage<-lineage.frame[!is.na(x=taxaId),]
write.table(x=new_seqid, file="out_seqid.txt", quote=F, row.names=F, col.names=F )
write.table(x=final_lineage, file="out_taxonomy.txt", sep="; ", quote=F, col.names=F, row.names=F, na="unclassified")

#retrieve a second set of files for those records where a phylum is reported
new_fasta_seqs<-fasta_seqs[!is.na(x=lineage.frame[[2]])]
write.fasta(sequences=new_fasta_seqs, names=names(new_fasta_seqs), file.out="out_cl.fas")
new_seqid<-accessions[[1]][!is.na(x=lineage.frame[[2]])]
lineage.frame <- as.data.frame(x=lineage, optional=TRUE, make.names=NA)
final_lineage<-lineage.frame[!is.na(x=lineage.frame[[2]]),]
write.table(x=new_seqid, file="out_seqid_cl.txt", quote=F, row.names=F, col.names=F )
write.table(x=final_lineage, file="out_taxonomy_cl.txt", sep="; ", quote=F, col.names=F
            , row.names=F, na="unclassified")
