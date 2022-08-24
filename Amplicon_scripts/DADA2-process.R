library("psadd") #plot_krona()
library("phyloseq")
library("dada2")
library("microbiome")
library("biomformat")
library("tidyverse")

############################

#DEFINE PATHS

#path="/opt/Alchemycode/HG-2.0/seqs/"
list.files(path)
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE)

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
#DADA2 cuts
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,maxN=0, maxEE=c(3,3), truncQ=2, rm.phix=TRUE,compress=TRUE, multithread=T)
#out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(285,150),truncQ=2, rm.phix=TRUE,compress=TRUE,multithread=F)
errF <- learnErrors(filtFs, multithread=F)
errR <- learnErrors(filtRs, multithread=F)
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
#Merge
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE,minOverlap =5)
seqtab <- makeSequenceTable(mergers)
#Quimeras
#minFoldParentOverAbundanceminFoldParentOverAbundance controla la definicion de quimeras, 2 o mas genera muchas cuando la calidad baja
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE, minFoldParentOverAbundance = 2)
getN <- function(x) sum(getUniques(x))
#track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))

# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)

track <- cbind(out, getN(dadaFs), getN(dadaRs), getN(mergers), rowSums(seqtab.nochim))

colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

#Taxonomy assing 
taxa <- assignTaxonomy(seqtab.nochim, "C:/Users/AMS/Downloads/silva_nr99_v138.1_wSpecies_train_set.fa.gz",multithread=F,tryRC=TRUE)
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
#Import to phyloseq object
ps <- phyloseq(otu_table(t(seqtab.nochim), taxa_are_rows=T),tax_table(taxa))
#change ASVs names#....
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
#Filtrar las ASVs sin phylum y sin dominio (artefactos)
ps <- subset_taxa(ps, !is.na(Phylum) & !Phylum %in% c(""))
ps <- subset_taxa(ps, !is.na(Kingdom) & !Kingdom %in% c(""))
#Add names as metadata
sample_names(ps) = sample.names
metadataMAB=data.frame(Nombres=colnames(otu_table(ps)))
rownames(metadataMAB)=metadataMAB$Nombres
sample_data(ps)=metadataMAB
#generete taxa table for 
mab2OtuT <- data.frame(otu_table(ps))
mab2Taxa <- data.frame(tax_table(ps))
mab2join <- merge(mab2Taxa,mab2OtuT,by = 0,all = T)
mab2join$Row.names = NULL
mab2join$Consensus = NULL
mab2PS=mab2join
#Compositional "abundanciar relativa"
ABrel_ps = transform(ps, "compositional") 
mab2OtuT <- data.frame(otu_table(ABrel_ps))
mab2Taxa <- data.frame(tax_table(ABrel_ps))
mab2join <- merge(mab2Taxa,mab2OtuT,by = 0,all = T)
#EXTRAS retirar columna de OTUid y Consensus
mab2join$Row.names = NULL
mab2join$Consensus = NULL

#Krona per sample

for (i in sample_names(ps)) {
  KronaPerSAmple = subset_samples(ps,Nombres==i)
  psadd::plot_krona(KronaPerSAmple,output = paste("/opt/Alchemycode/WORK",i,sep = "_"),variable = "Nombres")
}
##Barplots
