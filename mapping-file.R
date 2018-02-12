
mapping_file <- read.table("../../Dropbox/PhD/Sequencing/Mona_run_test/mapping-miseq-test - Sheet1.tsv")
colnames(mapping_file) <- (c("#SampleID", "BarcodeSequence", "LinkerPrimerSequence", "ReversePrimer", "Description"))
write.table(mapping_file, "mapping-miseq-test-R1.txt", sep = "\t", row.names = F, quote = F)



