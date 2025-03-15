# Load required libraries
library("psadd") #plot_krona()
library("phyloseq")
library("dada2")
library("microbiome")
library("biomformat")
library("tidyverse")
library("reticulate") # For Python/sklearn integration
library("Biostrings")

# ---- Setup Python Environment ----
# Point to your Python environment - update this with your path from 'which python3'
# after activating your sklearn_bio environment
use_python("~/environments/sklearn_bio/bin/python")

# Import sklearn libraries
sklearn <- import("sklearn.naive_bayes")
sklearn_preprocessing <- import("sklearn.feature_extraction.text")
numpy <- import("numpy")

# ---- Define Paths and Parameters ----
path <- "/opt/Alchemycode/Adrian/BAC/seqs/Trimm/second"
silva_fasta_path <- "/opt/Alchemycode/SilvaDB/silva_nr99_v138.1_wSpecies_train_set.fa.gz"
silva_taxonomy_path <- "/opt/Alchemycode/SilvaDB/silva_nr99_v138.1_wSpecies_train_set.fa.gz" # Assuming tax is in same file
silva_species_path <- "/opt/Alchemycode/SilvaDB/silva_species_assignment_v132.fa.gz"

# Define your V3-V4 primers
forward_primer <- "CCTACGGGNGGCWGCAG"
reverse_primer <- "GACTACHVGGGTATCTAATCC"

# ---- Functions for sklearn Taxonomy Assignment ----
# Function to extract V3-V4 region from reference sequences
extract_v3v4_region <- function(ref_seqs_path, forward_primer, reverse_primer) {
  # Load reference sequences
  message("Loading reference sequences...")
  ref_seqs <- readDNAStringSet(ref_seqs_path)
  
  # Function to find primer sites and extract region
  extract_region <- function(seq, fwd_primer, rev_primer) {
    # Convert to character for easier handling
    seq_char <- as.character(seq)
    
    # Find primer binding sites (allowing for some mismatches)
    # Note: Using fixed=FALSE to allow for degenerate bases in primers
    fwd_match <- gregexpr(fwd_primer, seq_char, fixed = FALSE)[[1]]
    rev_match <- gregexpr(rev_primer, seq_char, fixed = FALSE)[[1]]
    
    # If both primers found, extract region
    if (fwd_match[1] > 0 && rev_match[1] > 0) {
      start_pos <- fwd_match[1] + nchar(fwd_primer)
      end_pos <- rev_match[1] - 1
      
      # Ensure valid coordinates
      if (start_pos < end_pos) {
        return(substr(seq_char, start_pos, end_pos))
      }
    }
    return(NA)
  }
  
  # Apply extraction to reference sequences (this may take time)
  message("Extracting V3-V4 regions from reference sequences...")
  extracted_regions <- character(length(ref_seqs))
  
  # Process in batches to reduce memory usage
  batch_size <- 1000
  num_batches <- ceiling(length(ref_seqs) / batch_size)
  
  for (i in 1:num_batches) {
    message(paste0("Processing batch ", i, " of ", num_batches))
    start_idx <- (i-1) * batch_size + 1
    end_idx <- min(i * batch_size, length(ref_seqs))
    
    batch_seqs <- ref_seqs[start_idx:end_idx]
    extracted_regions[start_idx:end_idx] <- sapply(
      batch_seqs, 
      extract_region, 
      fwd_primer = forward_primer, 
      rev_primer = reverse_primer
    )
  }
  
  # Filter out sequences where extraction failed
  valid_idx <- !is.na(extracted_regions)
  message(paste0("Successfully extracted ", sum(valid_idx), " regions out of ", length(ref_seqs), " reference sequences"))
  
  filtered_regions <- extracted_regions[valid_idx]
  filtered_taxa <- names(ref_seqs)[valid_idx]
  
  # Return as DNAStringSet with taxonomy labels
  result_seqs <- DNAStringSet(filtered_regions)
  names(result_seqs) <- filtered_taxa
  
  return(result_seqs)
}

# Function to parse Silva taxonomy strings into a structured format
parse_silva_taxonomy <- function(tax_strings) {
  # Create a dataframe to hold parsed taxonomy
  tax_df <- data.frame(
    Kingdom = character(length(tax_strings)),
    Phylum = character(length(tax_strings)),
    Class = character(length(tax_strings)),
    Order = character(length(tax_strings)),
    Family = character(length(tax_strings)),
    Genus = character(length(tax_strings)),
    Species = character(length(tax_strings)),
    stringsAsFactors = FALSE
  )
  
  # Process each taxonomy string
  for (i in seq_along(tax_strings)) {
    # Extract taxonomy from FASTA header
    # Example format: ">D0__Bacteria;D1__Proteobacteria;D2__Gammaproteobacteria..."
    tax_string <- tax_strings[i]
    
    # Extract taxonomy levels - handle different Silva formats
    if (grepl(";", tax_string)) {
      # Split by semicolon
      tax_parts <- strsplit(tax_string, ";")[[1]]
      
      # Different possible formats
      if (grepl("^>*D[0-9]__", tax_parts[1])) {
        # Silva 132 format: ">D0__Bacteria;D1__Proteobacteria..."
        for (j in seq_along(tax_parts)) {
          level_name <- sub("^>*D[0-9]__", "", tax_parts[j])
          if (j == 1) tax_df$Kingdom[i] <- level_name
          else if (j == 2) tax_df$Phylum[i] <- level_name
          else if (j == 3) tax_df$Class[i] <- level_name
          else if (j == 4) tax_df$Order[i] <- level_name
          else if (j == 5) tax_df$Family[i] <- level_name
          else if (j == 6) tax_df$Genus[i] <- level_name
          else if (j == 7) tax_df$Species[i] <- level_name
        }
      } else if (grepl("^>*[kpcofgs]__", tax_parts[1])) {
        # Standard format: ">k__Bacteria;p__Proteobacteria..."
        for (j in seq_along(tax_parts)) {
          level_prefix <- substr(tax_parts[j], 1, 3)
          level_name <- sub("^>*[kpcofgs]__", "", tax_parts[j])
          
          if (level_prefix == "k__" || level_prefix == ">k_") tax_df$Kingdom[i] <- level_name
          else if (level_prefix == "p__" || level_prefix == "p_") tax_df$Phylum[i] <- level_name
          else if (level_prefix == "c__" || level_prefix == "c_") tax_df$Class[i] <- level_name
          else if (level_prefix == "o__" || level_prefix == "o_") tax_df$Order[i] <- level_name
          else if (level_prefix == "f__" || level_prefix == "f_") tax_df$Family[i] <- level_name
          else if (level_prefix == "g__" || level_prefix == "g_") tax_df$Genus[i] <- level_name
          else if (level_prefix == "s__" || level_prefix == "s_") tax_df$Species[i] <- level_name
        }
      }
    } else {
      # Single entry - likely just a species name
      tax_df$Species[i] <- tax_string
    }
  }
  
  return(tax_df)
}

# Function to train sklearn classifier
train_sklearn_classifier <- function(training_seqs, training_taxa) {
  message("Setting up sklearn classifier...")
  
  # Set up k-mer function in Python
  py_run_string("
import numpy as np
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.naive_bayes import MultinomialNB

def kmerize_seq(seq, k=8):
    # Convert sequence to k-mers
    return ' '.join([seq[i:i+k] for i in range(len(seq)-k+1)])
  ")
  
  # Convert DNAStringSet to character vector
  seq_chars <- as.character(training_seqs)
  
  message("Creating k-mer representation of training sequences...")
  # Create k-mer representation in Python
  py$seq_chars <- seq_chars
  py_run_string("
# Convert sequences to k-mers
kmerized_seqs = [kmerize_seq(seq) for seq in seq_chars]

# Create CountVectorizer
vectorizer = CountVectorizer(analyzer='word', ngram_range=(1, 1))
X_train = vectorizer.fit_transform(kmerized_seqs)
  ")
  
  # Pass taxonomy dataframe to Python
  py$training_taxa <- training_taxa
  
  message("Training classifiers for each taxonomic level...")
  # Train classifiers for each taxonomic level
  py_run_string("
# Initialize classifiers
kingdom_clf = MultinomialNB()
phylum_clf = MultinomialNB()
class_clf = MultinomialNB()
order_clf = MultinomialNB()
family_clf = MultinomialNB()
genus_clf = MultinomialNB()
species_clf = MultinomialNB()

# Train classifiers
kingdom_clf.fit(X_train, training_taxa['Kingdom'])
phylum_clf.fit(X_train, training_taxa['Phylum'])
class_clf.fit(X_train, training_taxa['Class'])
order_clf.fit(X_train, training_taxa['Order'])
family_clf.fit(X_train, training_taxa['Family'])
genus_clf.fit(X_train, training_taxa['Genus'])
species_clf.fit(X_train, training_taxa['Species'])
  ")
  
  message("Classifier training complete.")
  
  # Return the trained model details
  return(list(
    trained = TRUE,
    vectorizer = py$vectorizer,
    classifiers = list(
      Kingdom = py$kingdom_clf,
      Phylum = py$phylum_clf,
      Class = py$class_clf,
      Order = py$order_clf,
      Family = py$family_clf,
      Genus = py$genus_clf,
      Species = py$species_clf
    )
  ))
}

# Function to classify ASVs using the trained sklearn classifier
classify_with_sklearn <- function(asv_seqs, trained_model) {
  message("Classifying ASVs with sklearn...")
  
  # Convert ASV sequences to character vector
  asv_chars <- as.character(asv_seqs)
  
  # Pass to Python
  py$asv_chars <- asv_chars
  py_run_string("
# Convert ASVs to k-mers
asv_kmerized = [kmerize_seq(seq) for seq in asv_chars]

# Transform to feature vectors
X_asv = vectorizer.transform(asv_kmerized)

# Predict for each level
kingdom_pred = kingdom_clf.predict(X_asv)
phylum_pred = phylum_clf.predict(X_asv)
class_pred = class_clf.predict(X_asv)
order_pred = order_clf.predict(X_asv)
family_pred = family_clf.predict(X_asv)
genus_pred = genus_clf.predict(X_asv)
species_pred = species_clf.predict(X_asv)

# Get probabilities
kingdom_prob = np.max(kingdom_clf.predict_proba(X_asv), axis=1)
phylum_prob = np.max(phylum_clf.predict_proba(X_asv), axis=1)
class_prob = np.max(class_clf.predict_proba(X_asv), axis=1)
order_prob = np.max(order_clf.predict_proba(X_asv), axis=1)
family_prob = np.max(family_clf.predict_proba(X_asv), axis=1)
genus_prob = np.max(genus_clf.predict_proba(X_asv), axis=1)
species_prob = np.max(species_clf.predict_proba(X_asv), axis=1)
  ")
  
  # Construct taxonomy table
  tax_table <- data.frame(
    Kingdom = py$kingdom_pred,
    Phylum = py$phylum_pred,
    Class = py$class_pred,
    Order = py$order_pred,
    Family = py$family_pred,
    Genus = py$genus_pred,
    Species = py$species_pred,
    stringsAsFactors = FALSE
  )
  
  # Add confidence values
  confidence <- data.frame(
    Kingdom_conf = py$kingdom_prob,
    Phylum_conf = py$phylum_prob,
    Class_conf = py$class_prob,
    Order_conf = py$order_prob,
    Family_conf = py$family_prob,
    Genus_conf = py$genus_prob,
    Species_conf = py$species_prob
  )
  
  # Apply confidence threshold (using 0.8 as in QIIME2)
  for (i in 1:nrow(tax_table)) {
    if (confidence$Phylum_conf[i] < 0.8) {
      tax_table$Phylum[i] <- NA
      tax_table$Class[i] <- NA
      tax_table$Order[i] <- NA
      tax_table$Family[i] <- NA
      tax_table$Genus[i] <- NA
      tax_table$Species[i] <- NA
    }
    else if (confidence$Class_conf[i] < 0.8) {
      tax_table$Class[i] <- NA
      tax_table$Order[i] <- NA
      tax_table$Family[i] <- NA
      tax_table$Genus[i] <- NA
      tax_table$Species[i] <- NA
    }
    else if (confidence$Order_conf[i] < 0.8) {
      tax_table$Order[i] <- NA
      tax_table$Family[i] <- NA
      tax_table$Genus[i] <- NA
      tax_table$Species[i] <- NA
    }
    else if (confidence$Family_conf[i] < 0.8) {
      tax_table$Family[i] <- NA
      tax_table$Genus[i] <- NA
      tax_table$Species[i] <- NA
    }
    else if (confidence$Genus_conf[i] < 0.8) {
      tax_table$Genus[i] <- NA
      tax_table$Species[i] <- NA
    }
    else if (confidence$Species_conf[i] < 0.8) {
      tax_table$Species[i] <- NA
    }
  }
  
  # Set row names
  rownames(tax_table) <- names(asv_seqs)
  
  return(tax_table)
}

# ---- Main DADA2 Workflow ----
# Setup
message("Starting DADA2 workflow...")
list.files(path)

# Find fastq files
fnFs <- sort(list.files(path, pattern="_1.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_1."), `[`, 1)

# Setup filtered files
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))

# Create filtered directory if it doesn't exist
if(!dir.exists(file.path(path, "filtered"))) {
  dir.create(file.path(path, "filtered"))
}

# Filter and trim
message("Filtering and trimming reads...")
out <- filterAndTrim(fnFs, filtFs, 
                   trimRight=45, 
                   maxEE=3, 
                   truncQ=3, 
                   rm.phix=TRUE, 
                   compress=TRUE, 
                   verbose=TRUE, 
                   multithread=TRUE)

head(out, n=20)

# Learn error rates
message("Learning error rates...")
errF <- learnErrors(filtFs, multithread=TRUE)

# Dereplicate
message("Dereplicating sequences...")
derepFs <- derepFastq(filtFs)

# Denoise
message("Denoising sequences...")
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)

# Make sequence table
message("Creating sequence table...")
seqtab <- makeSequenceTable(dadaFs)

# Save intermediate results after ASV generation
message("Saving intermediate ASV table...")
saveRDS(seqtab, "seqtab_before_chimera.rds")
saveRDS(dadaFs, "dadaFs_results.rds")
saveRDS(derepFs, "derepFs_results.rds")

# Remove chimeras
message("Removing chimeras...")
seqtab.nochim <- removeBimeraDenovo(seqtab, 
                                  method="consensus", 
                                  multithread=TRUE, 
                                  verbose=TRUE, 
                                  minFoldParentOverAbundance = 1.5)

# Save post-chimera ASV table
message("Saving post-chimera ASV table...")
saveRDS(seqtab.nochim, "seqtab_nochim.rds")

# Track reads through pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "nonchim")
rownames(track) <- sample.names
head(track)
write.table(track, sep="\t", quote=FALSE, file="Track.txt")

# ---- Perform BOTH Taxonomy Assignments ----
# Check if we can resume from a previous run
if (file.exists("taxa_standard.rds")) {
  message("Found existing taxa_standard.rds file, loading previous results...")
  taxa_standard <- readRDS("taxa_standard.rds")
} else {
  # 1. Standard DADA2 taxonomy assignment
  message("Performing standard DADA2 taxonomy assignment...")
  taxa_standard <- assignTaxonomy(seqtab.nochim, silva_fasta_path, multithread=TRUE, tryRC=TRUE)
  taxa_standard <- addSpecies(taxa_standard, silva_species_path)
  
  # Create a copy of the standard taxa for the standard phyloseq object
  saveRDS(taxa_standard, "taxa_standard.rds")
}

# 2. sklearn-based taxonomy assignment
message("Checking for previous sklearn taxonomy results...")

# First check if we have the trained model already
sklearn_model_exists <- file.exists("trained_sklearn_model.rds")
v3v4_refs_exists <- file.exists("v3v4_refs.rds")
parsed_taxa_exists <- file.exists("parsed_taxa.rds")
taxa_sklearn_exists <- file.exists("taxa_sklearn.rds")

if (taxa_sklearn_exists) {
  message("Found existing sklearn taxonomy results, loading...")
  taxa_sklearn_matrix <- readRDS("taxa_sklearn.rds")
} else {
  message("Starting sklearn-based taxonomy assignment...")
  
  # Check if we can resume from extracted regions
  if (v3v4_refs_exists && parsed_taxa_exists) {
    message("Loading previously extracted V3-V4 regions and taxonomy...")
    v3v4_refs <- readRDS("v3v4_refs.rds")
    parsed_taxa <- readRDS("parsed_taxa.rds")
  } else {
    # Extract and prepare reference sequences for training
    message("Extracting V3-V4 regions from reference sequences...")
    v3v4_refs <- extract_v3v4_region(silva_fasta_path, forward_primer, reverse_primer)
    saveRDS(v3v4_refs, "v3v4_refs.rds")
    
    # Extract taxonomy from sequence headers
    ref_taxa <- names(v3v4_refs)
    
    # Parse taxonomy
    message("Parsing taxonomy from reference sequences...")
    parsed_taxa <- parse_silva_taxonomy(ref_taxa)
    saveRDS(parsed_taxa, "parsed_taxa.rds")
  }
  
  # Check if we already have a trained model
  if (sklearn_model_exists) {
    message("Loading previously trained sklearn model...")
    trained_model <- readRDS("trained_sklearn_model.rds")
  } else {
    # Train classifier
    message("Training sklearn classifier...")
    trained_model <- train_sklearn_classifier(v3v4_refs, parsed_taxa)
    saveRDS(trained_model, "trained_sklearn_model.rds")
  }
  
  # Prepare ASV sequences for classification
  asv_seqs <- colnames(seqtab.nochim)
  asv_headers <- paste0("ASV", seq_along(asv_seqs))
  asv_objects <- DNAStringSet(asv_seqs)
  names(asv_objects) <- asv_headers
  
  # Classify with sklearn
  message("Classifying ASVs with sklearn...")
  taxa_sklearn <- classify_with_sklearn(asv_objects, trained_model)
  
  # Convert to matrix for phyloseq
  taxa_sklearn_matrix <- as.matrix(taxa_sklearn)
  
  # Save sklearn results
  saveRDS(taxa_sklearn_matrix, "taxa_sklearn.rds")
}

# ---- Create Phyloseq Objects for Both Methods ----
# Create standard phyloseq object
message("Creating standard phyloseq object...")
ps_standard <- phyloseq(otu_table(t(seqtab.nochim), taxa_are_rows=TRUE), tax_table(taxa_standard))

# Add DNA sequences
dna <- Biostrings::DNAStringSet(taxa_names(ps_standard))
names(dna) <- taxa_names(ps_standard)
ps_standard <- merge_phyloseq(ps_standard, dna)
taxa_names(ps_standard) <- paste0("ASV", seq(ntaxa(ps_standard)))

# Filter taxa
ps_standard <- subset_taxa(ps_standard, !is.na(Phylum) & !Phylum %in% c("") & Kingdom != "Unknown")

# Add sample names as metadata
sample_names(ps_standard) <- sample.names
metadataMAB <- data.frame(Nombres=colnames(otu_table(ps_standard)))
rownames(metadataMAB) <- metadataMAB$Nombres
sample_data(ps_standard) <- metadataMAB

# Create taxonomy table
mab2OtuT_standard <- data.frame(otu_table(ps_standard))
mab2Taxa_standard <- data.frame(tax_table(ps_standard))
mab2join_standard <- merge(mab2Taxa_standard, mab2OtuT_standard, by=0, all=TRUE)
mab2join_standard$Row.names <- NULL
if("Consensus" %in% colnames(mab2join_standard)) {
  mab2join_standard$Consensus <- NULL
}

# Save standard results
save(ps_standard, mab2join_standard, file="bacfile_standard.rda")

# Create sklearn phyloseq object
message("Creating sklearn phyloseq object...")
ps_sklearn <- phyloseq(otu_table(t(seqtab.nochim), taxa_are_rows=TRUE), tax_table(taxa_sklearn_matrix))

# Add DNA sequences
dna <- Biostrings::DNAStringSet(taxa_names(ps_sklearn))
names(dna) <- taxa_names(ps_sklearn)
ps_sklearn <- merge_phyloseq(ps_sklearn, dna)
taxa_names(ps_sklearn) <- paste0("ASV", seq(ntaxa(ps_sklearn)))

# Filter taxa
ps_sklearn <- subset_taxa(ps_sklearn, !is.na(Phylum) & !Phylum %in% c("") & Kingdom != "Unknown")

# Add sample names as metadata
sample_names(ps_sklearn) <- sample.names
metadataMAB <- data.frame(Nombres=colnames(otu_table(ps_sklearn)))
rownames(metadataMAB) <- metadataMAB$Nombres
sample_data(ps_sklearn) <- metadataMAB

# Create taxonomy table
mab2OtuT_sklearn <- data.frame(otu_table(ps_sklearn))
mab2Taxa_sklearn <- data.frame(tax_table(ps_sklearn))
mab2join_sklearn <- merge(mab2Taxa_sklearn, mab2OtuT_sklearn, by=0, all=TRUE)
mab2join_sklearn$Row.names <- NULL
if("Consensus" %in% colnames(mab2join_sklearn)) {
  mab2join_sklearn$Consensus <- NULL
}

# Save sklearn results
save(ps_sklearn, mab2join_sklearn, file="bacfile_sklearn.rda")

# Create a comparison metrics between the two methods
message("Creating taxonomy comparison metrics...")
asv_ids <- taxa_names(ps_standard)

# Create comparison table
tax_comparison <- data.frame(
  ASV_ID = asv_ids,
  Standard_Phylum = tax_table(ps_standard)[,"Phylum"],
  Sklearn_Phylum = tax_table(ps_sklearn)[,"Phylum"],
  Standard_Class = tax_table(ps_standard)[,"Class"],
  Sklearn_Class = tax_table(ps_sklearn)[,"Class"],
  Standard_Order = tax_table(ps_standard)[,"Order"],
  Sklearn_Order = tax_table(ps_sklearn)[,"Order"],
  Standard_Family = tax_table(ps_standard)[,"Family"],
  Sklearn_Family = tax_table(ps_sklearn)[,"Family"],
  Standard_Genus = tax_table(ps_standard)[,"Genus"],
  Sklearn_Genus = tax_table(ps_sklearn)[,"Genus"],
  Standard_Species = tax_table(ps_standard)[,"Species"],
  Sklearn_Species = tax_table(ps_sklearn)[,"Species"]
)

# Calculate agreement percentages
agreement <- data.frame(
  Level = c("Phylum", "Class", "Order", "Family", "Genus", "Species"),
  Agreement_Count = numeric(6),
  Total_Count = numeric(6),
  Agreement_Percent = numeric(6)
)

for (i in 1:nrow(agreement)) {
  level <- agreement$Level[i]
  std_col <- paste0("Standard_", level)
  skl_col <- paste0("Sklearn_", level)
  
  # Count where both methods assigned a taxonomy (neither is NA)
  both_assigned <- sum(!is.na(tax_comparison[[std_col]]) & !is.na(tax_comparison[[skl_col]]))
  
  # Count where both methods agree
  agree_count <- sum(tax_comparison[[std_col]] == tax_comparison[[skl_col]], na.rm = TRUE)
  
  agreement$Agreement_Count[i] <- agree_count
  agreement$Total_Count[i] <- both_assigned
  agreement$Agreement_Percent[i] <- ifelse(both_assigned > 0, 
                                           round(agree_count / both_assigned * 100, 2), 
                                           0)
}

# Save comparison results
write.csv(tax_comparison, "taxonomy_comparison.csv", row.names = FALSE)
write.csv(agreement, "taxonomy_agreement_metrics.csv", row.names = FALSE)

message("Processing complete. Results saved to bacfile_standard.rda and bacfile_sklearn.rda")
message("Taxonomy comparison saved to taxonomy_comparison.csv and taxonomy_agreement_metrics.csv")
