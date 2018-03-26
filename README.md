# miseq-test

## Study through sequences
### Unzip files
Gunzip $file.fastq.gz

### View the size of the files
ls -lh

### Count the number of fastq files in the folder
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
#### If there are several fastq files
for i in *.fastq; do /data/apps/bbmap/bbduk.sh in=$i out=adapter-trimmed/$i-trimmed.fastq ref=/data/apps/bbmap/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 tbo tpe; done

#### If there are only 2fastq files (R1 and R2)
/data/apps/bbmap/bbduk.sh -Xmx1g in=file_R#_001.fastq out=adapter-trimmed/file-trimmed_R#.fastq ref=/data/apps/bbmap/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 tbo tpe

