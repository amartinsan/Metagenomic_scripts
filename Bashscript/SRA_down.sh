#!/bin/bash

#Download sequences from NCBI SRA database

#get the SRA list and metadata file (if they have)

#######SRR_Acc_List.txt##############

#Get the SRA toolkit to download 

    sudo apt-get update
    sudo apt-get install sra-toolkit

    mkdir NCBI_download
    cd NCBI_download

#Prefetch SRA

    prefetch --option-file SRR_Acc_List.txt

#This will create a folder for each SRR, move all .sra to one folder before convertion or make a loop

    mkdir All_SRA_Files
    mv */*.sra All_SRA_Files/

#Convert from SRA to FASTQ 

    fastq-dump --split-files *.sra

#Check if you have all the files and sequences
#Bonus check average lenght of sequences in a fastq file
    awk 'NR%4 == 2 { total += length($0) } END { printf "Average sequence length: %.2f\n", total/(NR/4) }' your_file.fastq

#Downloaded .fastq files have the SRR name and not the one i want or the one in metadata.
#Download the metadata from NCBI SRA viewver
#Make a file that has the name of the sequence and the equivalent name:

    SRR12345678_1.fastq new_name_1.fastq
    SRR23456789_1.fastq new_name_2.fastq
    SRR34567890_1.fastq new_name_3.fastq

#While loop to change the name using a metadata file as reference

# Path to the metadata file
metadata_file="metadatafast.txt"

# Loop through lines in metadata file
while IFS=$'\t' read -r original_name new_name; do
    # Rename the file
    mv "$original_name" "$new_name"
done < "$metadata_file"


# Rename corresponding reverse reads if applicable
   reverse_name="${original_name/_1.fastq/_2.fastq}"
   new_reverse_name="${new_name/_1.fastq/_2.fastq}"
   mv "$reverse_name" "$new_reverse_name"
done < "$metadata_file"
