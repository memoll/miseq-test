## Mona Parizadeh - March 2018

## How to prepare mapping file?

#### 1. In excel: #SampleID - BarcodeSequence	- LinkerPrimerSequence -	Barcode	- Reverse	- Description; save as .tsv

#### 2. Import mapping_file.tsv in R:
mapping_file <- read.table("Google Drive/mapping-miseq-run-test-2016.tsv")
colnames(mapping_file) <- (c("#SampleID", "BarcodeSequence", "LinkerPrimerSequence", "Barcode", "Reverse", "Description"))
mapping_file

##### Verify the length of the longest barcode (we'll need it to make fake barcodes later)
max_BarcodeSequence <- max(nchar(as.character(mapping_file$BarcodeSequence)))
max_Barcode <- max(nchar(as.character(mapping_file$Barcode)))

#### 3. Add primers to the end of each barcode sequence (in order to make fake barcodes with the same length later)
BarcodeSequence_primer <-  paste0(mapping_file$BarcodeSequence,mapping_file$LinkerPrimerSequence)
mapping_file$BarcodeSequence <- BarcodeSequence_primer
Barcode_primer <- paste0(mapping_file$Barcode,mapping_file$Reverse)
mapping_file$Barcode <- Barcode_primer

#### 4. Cut all the new barcodes to the same length (length of the longest barcode)

##### length of the shortest barcode
#### BarcodeSequence column:
#nchar(mapping_file$BarcodeSequence)
BarcodeSequence_trimmed <- strtrim(mapping_file$BarcodeSequence,max_BarcodeSequence)
nchar(BarcodeSequence_trimmed)
mapping_file$BarcodeSequence <- BarcodeSequence_trimmed

##### Barcode column:
#nchar(mapping_file$Barcode)
Barcode_trimmed <- strtrim(mapping_file$Barcode,max_Barcode)
nchar(Barcode_trimmed)
mapping_file$Barcode <- Barcode_trimmed

#### 5. Add Barcode (reverse/r2) to BarcodeSequence (forward/r2)
BarcodeSequence_r1_r2 <- paste0(mapping_file$BarcodeSequence,mapping_file$Barcode)
mapping_file$BarcodeSequence <- BarcodeSequence_r1_r2

##### Remove the Barcode (reverse/r2) column
mapping_file$Barcode <- NULL
mapping_file

#### 6. Export the mapping_file as a .txt file
write.table(mapping_file, "R_projects/miseq-test/mapping-miseq-run-test-2016.txt", sep = "\t", row.names = F, quote = F)

#### 7. Validate the mapping file (qiime):
#In Terminal: validate_mapping_file.py -m mapping-miseq-run-test-2016.txt -o validate-mapping

#check if you need to replace M,K,etc in barcodes with "."

