Using fastqc and multiqc to analyze quality and make an easy report.

(fastp is another option)


	mkdir Forward ; mkdir Reverse

	#mv *R2.fastq Forward/ ; mv R1.fastq Reverse/

	for file in *.fastq.gz ; do    fastqc $file ; done  && multiqc .
