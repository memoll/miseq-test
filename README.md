# miseq-test

## Study through sequences
### Unzip files
Gunzip $file.fastq.gz

### View the size of the files
ls -lh

### Count number of fastq files in the folder
ls *.fastq | wc -l

### Count the number of sequences in each file
grep "@M" -c *.fastq

### Count primer sequences in each file 
#### use . instead of wild-card nucleotides
grep "$primer_sequence" -c *.fastq

* F-primer 799F (R1): AACMGGATTAGATACCCKG

grep "AAC.GGATTAGATACCC.G" -c *.fastq

* R-primer 1115R (R2): AGGGTTGCGCTCGTTG

grep "AGGGTTGCGCTCGTTG" -c *.fastq

## Trim adapters
#### Trim using bbduk (in case adapters need to be trimmed)
for i in *.fastq; do /data/apps/bbmap/bbduk.sh in=$i out=adapter-trimmed/$i-trimmed.fastq ref=/data/apps/bbmap/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 tbo tpe; done

## Make contigs using Mothur (for primers with different barcode lengths):
#### Add the oligo_16s.txt to the folder
* ffastq, rfastq : forward and reverse fastq
* oligos : takes a .txt file that contains the sequences of the paired primers and barcodes and their sample identifier
* checkorient : if =t and mothur cannot find the barcodes and primers, it will search the reverse compliment
* pdiffs : maximum number of differences to the primer
* bdiffs : maximum number of differences to the barcode

mothur/

make.contigs(ffastq=file_R1_001.fastq, rfastq=file_R2_001.fastq, oligos=oligo_16s.txt, checkorient=t, processors=8, pdiffs=1, bdiffs=1)