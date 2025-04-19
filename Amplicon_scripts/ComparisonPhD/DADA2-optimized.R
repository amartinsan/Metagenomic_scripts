
#Comparison script Coatza vs GoM

library("psadd")
library("phyloseq")
library("dada2")
library("microbiome")
library("biomformat")
library("tidyverse")
library("tictoc")

# Timer
tic.clearlog()
tic("Total workflow")

# Cores setup
#T <- parallel::detectCores() - 1

path <- getwd() # Use current directory

# Create QIIME2 export directory if it exists dont creat
qiime2_dir <- file.path(path, "qiime2_export")
if(!dir.exists(qiime2_dir)) {
  dir.create(qiime2_dir)
}

##################################BEGIN##################################
# List files in current directory
#list.files(path)

# Adjust file pattern for already merged reads
###########REVIZAR O DEJAR SIEMPRE EL MISMO PATRON (SUFIJO)############
fnFs <- sort(list.files(path, pattern=".fastq.gz", full.names = TRUE))
# Exclude non-sample files (like the analysis text file)
fnFs <- fnFs[!grepl("fastq_length|LenghtAnalysis", fnFs)]

# Extract sample names
sample.names <- sapply(strsplit(basename(fnFs), ".fastq", fixed=TRUE), `[`, 1)

# Check for duplicates
if(any(duplicated(sample.names))) {
  print("Warning: Duplicate sample names found!")
  print(sample.names[duplicated(sample.names)])
  # Use a more unique identifier if duplicates exist
  sample.names <- basename(fnFs) %>% sub(".fastq.gz", "", .)
}
# Check again for duplicates after fallback
if(any(duplicated(sample.names))) {
  stop("Still have duplicate sample names after fallback method!")
}

# Set up filtered file paths
filtFs <- file.path(path, "filtered", paste0(sample.names, "_filt.fastq"))

# Create filtered directory if it doesn't exist
if(!dir.exists(file.path(path, "filtered"))) {
  dir.create(file.path(path, "filtered"))
}

names(filtFs) <- sample.names

#################################FILTERING##################################
#####Modified for single-end reads#####
message("Filtering and trimming reads...")
message("Para este amplcion de 500 se dejaron 465 pb, aprox. el tamaño del amplicon V3-V4")
out <- filterAndTrim(fnFs, filtFs,
                   trimRight=37,
                   minLen=150,  # Adjust based on your read quality
                   maxEE=3,
                   truncQ=3,
                   rm.phix=TRUE,
                   compress=TRUE,
                   verbose=TRUE,
                   multithread=T)


head(out)
tic("Learn Errors")
errF <- learnErrors(filtFs, multithread=T)
toc(log = TRUE)

tic("Dereplicate-Denoise")
derepFs <- derepFastq(filtFs)
names(derepFs) <- sample.names
dadaFs <- dada(derepFs, err=errF, multithread=T)
toc(log = TRUE)

# Modified for single-end reads (skip merge step)
seqtab <- makeSequenceTable(dadaFs)

tic("Remove Chimeras")
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE, minFoldParentOverAbundance = 2)

toc(log = TRUE)

# Track reads through pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoised", "nonchim")
rownames(track) <- sample.names
head(track, n=25)
# Write tracking table
write.table(track, sep="\t", quote=F, file=file.path(path, "Track.txt"))
saveRDS(seqtab.nochim, file.path(path, "seqtab.nochim.rds"))
print("DADA2 pipeline completed successfully.")
#########################################################
# EXPORT FOR QIIME2
#########################################################
# 1. Export ASV table (feature table)
asv_table <- t(seqtab.nochim) # QIIME2 expects features as rows
write.table(asv_table, file.path(qiime2_dir, "asv_table.tsv"),
            sep="\t", quote=FALSE)

# 2. Export sequences as FASTA (representative sequences)
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- paste0("ASV", seq_along(asv_seqs))
asv_fasta <- c(rbind(paste0(">", asv_headers), asv_seqs))
writeLines(asv_fasta, file.path(qiime2_dir, "rep-seqs.fasta"))

# 3. Create a mapping file between sequence and ASV IDs
asv_map <- data.frame(
  ASV_ID = asv_headers,
  Sequence = asv_seqs
)
write.table(asv_map, file.path(qiime2_dir, "asv_map.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

# 4. Create a simple metadata file with sample IDs
metadata <- data.frame(
  SampleID = sample.names,
  row.names = sample.names
)
write.table(metadata, file.path(qiime2_dir, "metadata.tsv"),
            sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

# 5. Generate BIOM file (easier import method)
biom_data <- make_biom(data=asv_table)
write_biom(biom_data, file.path(qiime2_dir, "feature-table.biom"))
print("Exported files for QIIME2:")
#########################################################
# PART 1: Basic Taxonomy and Save Initial Phyloseq
#########################################################
tic("Assign Basic Taxonomy")
print("Beginning basic taxonomy assignment...")
taxa <- assignTaxonomy(seqtab.nochim,
                     "/opt/Alchemycode/SilvaDB/silva_nr99_v138.1_wSpecies_train_set.fa.gz",
                     multithread=T,
                     tryRC=TRUE)
print("Basic taxonomy assignment complete.")
toc(log = TRUE)

# Create basic phyloseq object
tic("Create Basic Phyloseq Object and save")
ps_basic <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), tax_table(taxa))


# Rename ASVs
dna <- Biostrings::DNAStringSet(taxa_names(ps_basic))
names(dna) <- taxa_names(ps_basic)
ps_basic <- merge_phyloseq(ps_basic, dna)
taxa_names(ps_basic) <- paste0("ASV", seq(ntaxa(ps_basic)))
ps_basic <- subset_taxa(ps_basic, !is.na(Phylum) & !Phylum %in% c("") & Kingdom != "Unknown")
print("Phyloseq object created and removed unknown phylums")
# Add sample metadata
metadataMAB <- data.frame(Nombres=sample_names(ps_basic))
rownames(metadataMAB) <- metadataMAB$Nombres
sample_data(ps_basic) <- metadataMAB
save(ps_basic, file=file.path(path, "bacfile.rda"))
print("Basic phyloseq object created.")
toc(log = TRUE)


# Export basic taxonomy for QIIME2
#tax_cols <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
#tax_table_basic <- as.data.frame(taxa)
#tax_table_basic$Feature.ID <- asv_headers
#tax_table_basic <- tax_table_basic[, c("Feature.ID", tax_cols)]

# Clean up taxonomy strings for QIIME2 compatibility
#for (col in tax_cols) {
#  tax_table_basic[[col]] <- gsub("^[a-z]__", "", tax_table_basic[[col]])
#}

# Write basic taxonomy table
#write.table(tax_table_basic, file.path(qiime2_dir, "taxonomy_basic.tsv"),
#            sep="\t", quote=FALSE, row.names=FALSE)


# Create taxonomy table
#tic("Create Basic Taxonomy Table")
#mab2OtuT <- data.frame(otu_table(ps_basic))
##mab2Taxa <- data.frame(tax_table(ps_basic))
#mab2join <- merge(mab2Taxa, mab2OtuT, by = 0, all = TRUE)
#mab2join$Row.names <- NULL
#mab2join$Consensus <- NULL
#mab2PS_basic <- mab2join
#write.table(mab2join, sep="\t", quote=FALSE, row.names=FALSE,file="TaxonomyTable.txt")

##generete taxa table for
#mab2OtuT <- data.frame(otu_table(ps))
#mab2Taxa <- data.frame(tax_table(ps))
#mab2join <- merge(mab2Taxa,mab2OtuT,by = 0,all = T)
#mab2join$Row.names = NULL
#mab2join$Consensus = NULL
#mab2PS=mab2join
#Compositional "abundanciar relativa"
#ABrel_ps = transform(ps, "compositional")
#mab2OtuT <- data.frame(otu_table(ABrel_ps))
#mab2Taxa <- data.frame(tax_table(ABrel_ps))
#mab2join <- merge(mab2Taxa,mab2OtuT,by = 0,all = T)
#EXTRAS retirar columna de OTUid y Consensus
#mab2join$Row.names = NULL
#mab2join$Consensus = NULL


#toc(log = TRUE)

#########################################################
# PART 2: Add Species and Save Complete Phyloseq
#########################################################

tic("Add Species Information")
taxa_with_species <- tryCatch({
  addSpecies(taxa, "/opt/Alchemycode/SilvaDB/silva_species_assignment_v138.1.fa.gz")
}, error = function(e) {
  message("Error in addSpecies: ", e$message)
  message("Saving current progress...")
  save(taxa, file=file.path(path, "taxa_before_species.rda"))
  stop("Species assignment failed. Saved progress to taxa_before_species.rda")
})
toc(log = TRUE)

# Export taxonomy for QIIME2
#tax_cols <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
#tax_table <- as.data.frame(taxa_with_species)
#tax_table$Feature.ID <- asv_headers
#tax_table <- tax_table[, c("Feature.ID", tax_cols)]

# Clean up taxonomy strings for QIIME2 compatibility
#for (col in tax_cols) {
#  tax_table[[col]] <- gsub("^[a-z]__", "", tax_table[[col]])
#}

# Write taxonomy table
#write.table(tax_table, file.path(qiime2_dir, "taxonomyspecies.tsv"),
#            sep="\t", quote=FALSE, row.names=FALSE)

# Create complete phyloseq object with species
tic("Create Complete Phyloseq Object")
ps_complete <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), tax_table(taxa_with_species))
# Rename ASVs
dna <- Biostrings::DNAStringSet(taxa_names(ps_complete))
names(dna) <- taxa_names(ps_complete)
ps_complete <- merge_phyloseq(ps_complete, dna)
taxa_names(ps_complete) <- paste0("ASV", seq(ntaxa(ps_complete)))
ps_complete <- subset_taxa(ps_complete, !is.na(Phylum) & !Phylum %in% c("") & Kingdom != "Unknown")

# Add sample metadata
metadataMAB <- data.frame(Nombres=sample_names(ps_complete))
rownames(metadataMAB) <- metadataMAB$Nombres
sample_data(ps_complete) <- metadataMAB
save(ps_complete, file=file.path(path, "bacfilespecies.rda"))
print("Complete phyloseq object created with species information.")
toc(log = TRUE)

# Create complete taxonomy table
#tic("Create Complete Taxonomy Table")
#mab2OtuT <- data.frame(otu_table(ps_complete))
#mab2Taxa <- data.frame(tax_table(ps_complete))
#mab2join <- merge(mab2Taxa, mab2OtuT, by = 0, all = TRUE)
#mab2join$Row.names <- NULL
#mab2join$Consensus <- NULL
#mab2PS_complete <- mab2join
#write.table(mab2join, sep="\t", quote=FALSE, row.names=FALSE,file="TaxonomyTableSPECIES.txt")

# Save complete phyloseq object
#save(ps_complete, mab2join, file=file.path(path, "bacfilespecies.rda"))
#toc(log = TRUE)

# End total timer
toc(log = TRUE)

# Display timing log
cat("\n===== Execution Time Summary =====\n")
tic.log(format = TRUE)

# Print information about exported files
cat("\n===== QIIME2 Export Summary =====\n")
cat("Files exported for QIIME2 in directory:", qiime2_dir, "\n")
cat("To import into QIIME2, create an import script in the qiime2_export directory\n")
