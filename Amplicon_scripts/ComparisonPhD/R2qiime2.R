# Export to QIIME2
library(phyloseq)
library(biomformat)
library(Biostrings)

# First, check the orientation of your phyloseq object
#assuming the phyloseq is named ps_bacteria
#The rep-seqs.fasta have duplicates nothing a >>>>> sort | uniq <<<<<<<< in bash can't solve
print(ps_bacteria)

# We need to transpose the OTU table properly
# First, extract the raw OTU table
otu <- otu_table(ps_bacteria)

# Force the OTU table to have taxa as rows
if (!taxa_are_rows(otu)) {
  # Create a new OTU table with taxa as rows
  otu_transposed <- t(otu)
  # Set the correct orientation flag
  attr(otu_transposed, "taxa_are_rows") <- TRUE
  
  # Replace the OTU table in the phyloseq object
  otu_table(ps_bacteria) <- otu_transposed
  print("Properly transposed OTU table to have taxa as rows")
}

# Now check that the transposition worked
print(paste("Taxa are rows now:", taxa_are_rows(ps_bacteria)))
otu_dim <- dim(otu_table(ps_bacteria))
print(paste("New OTU table dimensions:", otu_dim[1], "taxa and", otu_dim[2], "samples"))

# Extract feature IDs
feature_ids <- taxa_names(ps_bacteria)

# Export OTU table
otu_table_mat <- as(otu_table(ps_bacteria), "matrix")
print(paste("OTU matrix dimensions:", nrow(otu_table_mat), "x", ncol(otu_table_mat)))

# Check if dimensions match
if (length(feature_ids) == nrow(otu_table_mat)) {
  rownames(otu_table_mat) <- feature_ids
  otu_biom <- make_biom(data=otu_table_mat)
  write_biom(otu_biom, "feature_table.biom")
  print("Successfully exported OTU table")
} else {
  print("ERROR: Dimensions still don't match")
}

# Export taxonomy
tax_table_mat <- as(tax_table(ps_bacteria), "matrix")
qiime2_tax <- data.frame(
  `Feature ID` = feature_ids,
  Taxon = apply(tax_table_mat, 1, function(row) {
    # Handle NA values in taxonomy
    row[is.na(row)] <- ""
    # Create the taxonomy string
    tax_string <- paste(paste0(
      c("k", "p", "c", "o", "f", "g", "s"),
      "__",
      row
    ), collapse=";")
    # Clean up empty levels
    tax_string <- gsub(".__", "", tax_string, fixed=FALSE)
    tax_string <- gsub(".__$", "", tax_string, fixed=FALSE)
    return(tax_string)
  }),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

write.table(qiime2_tax, "taxonomy.tsv", sep="\t", quote=FALSE, row.names=FALSE)
print("Successfully exported taxonomy file")

# Export representative sequences
if(!is.null(refseq(ps_bacteria))) {
  seq_data <- refseq(ps_bacteria)
  names(seq_data) <- feature_ids
  writeXStringSet(seq_data, "rep_seqs.fasta", format="fasta")
  print("Successfully exported representative sequences")
}
