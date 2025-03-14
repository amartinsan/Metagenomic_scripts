#Load required libraries
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

# Remove chimeras
message("Removing chimeras...")
seqtab.nochim <- removeBimeraDenovo(seqtab,
                                  method="consensus",
                                  multithread=TRUE,
                                  verbose=TRUE,
                                  minFoldParentOverAbundance = 1.5)

# Track reads through pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "nonchim")
rownames(track) <- sample.names
head(track)
write.table(track, sep="\t", quote=FALSE, file="Track.txt")

# ---- sklearn Taxonomy Assignment ----
# Two options: either continue with regular assignTaxonomy or use sklearn
message("Would you like to use sklearn for taxonomy assignment? (1=Yes, 0=No)")
use_sklearn <- as.logical(as.integer(readline(prompt = "Enter choice: ")))

if (use_sklearn) {
  message("Starting sklearn-based taxonomy assignment...")

  # Extract and prepare reference sequences for training
  v3v4_refs <- extract_v3v4_region(silva_fasta_path, forward_primer, reverse_primer)

  # Extract taxonomy from sequence headers
  ref_taxa <- names(v3v4_refs)

  # Parse taxonomy
  parsed_taxa <- parse_silva_taxonomy(ref_taxa)

  # Train classifier
  trained_model <- train_sklearn_classifier(v3v4_refs, parsed_taxa)

  # Prepare ASV sequences for classification
  asv_seqs <- colnames(seqtab.nochim)
  asv_headers <- paste0("ASV", seq_along(asv_seqs))
  asv_objects <- DNAStringSet(asv_seqs)
  names(asv_objects) <- asv_headers

  # Classify with sklearn
  tax_sklearn <- classify_with_sklearn(asv_objects, trained_model)

  # Convert to matrix for phyloseq
  taxa <- as.matrix(tax_sklearn)

  # Save sklearn results
  saveRDS(taxa, file="sklearn_taxonomy.rds")

  message("sklearn taxonomy assignment complete.")
} else {
  # Use standard DADA2 taxonomy assignment as in your original script
  message("Using standard DADA2 taxonomy assignment...")
  taxa <- assignTaxonomy(seqtab.nochim, silva_fasta_path, multithread=TRUE, tryRC=TRUE)
  taxa <- addSpecies(taxa, "/opt/Alchemycode/SilvaDB/silva_species_assignment_v132.fa.gz")
}

# ---- Create Phyloseq Object ----
message("Creating phyloseq object...")
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
ps <- phyloseq(otu_table(t(seqtab.nochim), taxa_are_rows=TRUE), tax_table(taxa))

# Add DNA sequences
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))

# Filter taxa
ps <- subset_taxa(ps, !is.na(Phylum) & !Phylum %in% c("") & Kingdom != "Unknown")

# Add sample names as metadata
sample_names(ps) <- sample.names
metadataMAB <- data.frame(Nombres=colnames(otu_table(ps)))
rownames(metadataMAB) <- metadataMAB$Nombres
sample_data(ps) <- metadataMAB

# Create taxonomy table
mab2OtuT <- data.frame(otu_table(ps))
mab2Taxa <- data.frame(tax_table(ps))
mab2join <- merge(mab2Taxa, mab2OtuT, by=0, all=TRUE)
mab2join$Row.names <- NULL
mab2join$Consensus <- NULL
mab2PS <- mab2join

# Save results
save(ps, mab2join, file="bacfile.rda")

message("Processing complete. Results saved to bacfile.rda")
