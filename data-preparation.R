## Mona Parizadeh - January 2018

## Data preparation for DADA2

### As mentioned in DADA2 tutorial, for the analysis we'll need paired-end sequences in fastq format that have been split (or “demultiplexed”) by sample and also the barcodes/adapters have already been removed.
### If you are using barcodes with the same length, we may want to use "extract_barcodes.py" and "split_libraries_fastq.py" or "split_libraries.py" from Qiime1
### "split_libraries_fastq.py" doesn't have an option to extract various length barcodes, in that case, we'll need to convert fastq files to fasta and qual files and then use "split_libraries.py".

## Split fastq files to fasta and qual
### Convert trimmed sequences format from fastq to fasta and qual
python /data/shared/scripts/fastq_to_fasta_qual.py $file-trimmed
