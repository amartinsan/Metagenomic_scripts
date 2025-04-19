#!/bin/bash

analyze_fastq_length() {
    local file=$1
    local sample_name=$(basename "$file" | sed 's/.fastq.gz$//')

    gzip -cd "$file" | awk -v sample="$sample_name" '
        NR%4==2 {
            len = length($0)
            lengths[len]++
            sum += len
            count++
            if (NR==2 || len < min) min = len
            if (len > max) max = len
        }
        END {
            if (count > 0) {
                avg = sum / count
                printf "%-30s | %8.2f | %8d | %8d | %8d\n", sample, avg, count, min, max
            }
        }'
}

print_header() {
    echo "--------------------------------------------------------------------------------------------------------"
    printf "%-30s | %8s | %8s | %8s | %8s\n" "Sample" "Avg Len" "Reads" "Min Len" "Max Len"
    echo "--------------------------------------------------------------------------------------------------------"
}

analyze_all_files() {
    shopt -s nullglob
    for file in *.fastq.gz; do
        [ -f "$file" ] || continue
        analyze_fastq_length "$file"
    done
}

# Main process
echo "Processing FASTQ Files"
print_header
analyze_all_files
echo "--------------------------------------------------------------------------------------------------------"

# Save to file
output_file="fastq_length_analysis_$(date +%Y%m%d).txt"
{
    echo "FASTQ Length Analysis Results - $(date)"
    echo ""
    print_header
    analyze_all_files
    echo "--------------------------------------------------------------------------------------------------------"
} > "$output_file"

echo -e "\nResults have been saved to $output_file"
