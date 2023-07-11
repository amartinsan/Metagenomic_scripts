Using fastqc and multiqc to analyze quality of seqs and make an easy report.

-needs fastqc and multiqc installed

(fastp is another option)


	#mkdir Forward ; mkdir Reverse

	#mv *R2.fastq Forward/ ; mv R1.fastq Reverse/

	for file in *.fastq.gz ; do    fastqc $file ; done  && multiqc .

AWK scripts for different stuff	

Get sequence lenght of each sequences in  .fasta 

        awk '/^>/ {if (seqlen){print seqlen}; print ;seqlen=0;next; } { seqlen += length($0)}END{print seqlen}' $$DESIRED.fasta$$ > fasta_lenght.txt
	
	#awk '/^>/{if (l!="") print l; print; l=0; next}{l+=length($0)}END{print l}' $FASTA of interest$$$ > fasta.lenght.txt

Calculate mean length of the sequences in a .fasta

	awk '{/>/&&++a||b+=length()}END{print b/a}' input.fasta > meanleneght_FASTANAME.fasta

Add a string after the header, maybe a suffix or another header

	awk '/^>/ {$0=$0 "WRITE-TEXT"}1' input.fasta > extraheder.fasta

Now with sed, prefix/suffix
	
	sed "s/>/>WRITE-TEXT/" input.fasta > output.fasta
	sed 's/^>/>SampleA/g' input.fasta > output.fasta


Replace the whole header with any string 

	awk '/^>/{print ">WRITE-TEXT" ++i; next}{print}' input.fasta > replaceHeader.fasta
	sed 's/>.*/>WRITE-text/' input.fasta > replaceHeader.fasta


Remove duplicated sequences, you can also use CD-HITo.


	sed -e '/^>/s/$/@/' -e 's/^>/#/' input.fasta | tr -d '\n' |
	 tr "#" "\n" | tr "@" "\t" | 
	sort -u -t $'\t' -f -k 2,2  | 
	sed -e 's/^/>/' -e 's/\t/\n/' > NOduplicate_output.fasta



Snipet to re-arrenge columns positions, wors on a tsv, substitute "\t" for any other separator.e

	awk 'BEGIN {OFS="\t"}; { print $1,$2,$3,$4,$15 " " $14}' input.txt  > rearrengeTEXT.txt

Count sequencs in .fasta/.fastq or .gz file

	grep -c "^>" input.fasta
	grep ">" input.fasta| wc -l 
	grep -c "^@" input.fasta

In a .gz file (fasta)

	zcat *.fastq.gz | echo $((`wc -l`/4))

Count the total number of sequences present in each different read files 
```
for i in *.fasta; do 
	grep -c "^>" $i; 
done

for i in *.fastq; do 
	grep -c "^+$" $i; 
done

for i in *.fastq.gz; do 
	echo $(zcat ${i} | wc -l)/4|bc;
done
```

Convert a .fastq file to .fasta 
__NOTE__ This assumes that each FASTQ entry spans only four lines as is customary)
```
sed -n '1~4s/^@/>/p;2~4p' input.fasta > output.fasta
sed '/^@/!d;s//>/;N' input.fasta > output.fasta
```

Extract sequences that are longer from 21 to 25 bp from a dataset
```
cat input.fasta | paste - - - - | awk 'length($2)  >= 21 && length($2) <= 25' | sed 's/\t/\n/g' > output.fasta
awk 'BEGIN {OFS = "\n"} {header = $0 ; getline seq ; getline qheader ; getline qseq ; if (length(seq) >= 21 && length(seq) <= 25) {print header, seq, qheader, qseq}}' input.fasta > output.fasta
```

Extract sequences the first e.g. 400 lines of a .fasta file
```
sed -n -e "1,400p" input.fasta > output.fasta
```

Count the total number of bases in a FASTA file
```
grep -v ">" input.fasta | wc | awk '{print $3 - $1}' > output.fasta

```

Do arithmetic operations on the bash command line
```
echo $((1 + 1))
echo $((1 - 1))
echo $((1 * 1))
echo $((1 / 1))
echo $(((1+3) / (1+1)))
echo "scale=1; 1/2" | bc 
```

Rename delimiters or headers
```
rename ' ' '_' *
rename 's/test-this/REPLACESTRING/g' *
```

Rename file names
```
for f in *.fastq ; do 
	NEW=${f%.fastq}R1.fastq; mv ${f} "${NEW}"; 
done

for i in *.fastq; do 
	mv "$i" "ITS_R1_$i"; # adding something
done

for i in *.fastq; do 
	mv $i ${i//ITS_R1_/}; # removing something
done

for file in *\ *; do
	mv "$file" "${file// /_}"
done

```

Modify seqeunce header by separating columns
```
cut -d '|' -f 1 < file.fastq 
```

Check identical sequences with bash script and collapse them
__NOTE__ make sure you remove the output.fasta before you run the script twice!
```
grep -v '>\|^$'  inputfile.fasta | sort -u > Sequence
while read SEQUENCE; do
	grep -w -m 1 -B1 $SEQUENCE input.fasta >> output.fasta
done < Sequence
```

Substitute first _ with . in Illumina fastq reads files
```for f in *fastq; do 
	sed '/^@/s/_/\./' "$f" > "$f.new"  &&  mv "$f.new" "$f"
done
