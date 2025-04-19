#!/bin/bash

# Properly activate conda
source /opt/bin/miniconda3/etc/profile.d/conda.sh
conda activate

# Create output directory if it doesn't exist
mkdir -p trimmed_output

# Process each FASTQ file
for f in *.fastq.gz; do
    # Run trimmomatic with CROP to cut all reads to exactly 500 bp
    # No MINLEN parameter so shorter reads won't be discarded
    trimmomatic SE -threads 1 ${f} trimmed_output/${f} CROP:500
done
