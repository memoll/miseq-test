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

### Count primer sequences in each file (use . instead of wild-card nucleotides)
grep "$primer_sequence" -c *.fastq

R-primer 799F (R1): AACMGGATTAGATACCCKG
grep "AAC.GGATTAGATACCC.G" -c *.fastq

F-primer 1115R (R2): AGGGTTGCGCTCGTTG
grep "AGGGTTGCGCTCGTTG" -c *.fastq

Trim using bbduk (in case adapters need to be trimmed)
for i in *.fastq; do /data/apps/bbmap/bbduk.sh in=$i out=adapter-trimmed/$i-trimmed.fastq ref=/data/apps/bbmap/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 tbo tpe; done