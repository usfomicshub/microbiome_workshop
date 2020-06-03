# The purpose is to obtain Amplicon Sequence variant table for all the samples with microbiome data (16S V4 region).
# The present analysis is based on dada2 package in R
# For more details, please follow the link https://benjjneb.github.io/dada2/tutorial.html
# We will start with demultiplexed fastq files for all samples and this analysis is for paired-end data
# Thus, for each sample, there will be two files (_R1_001.fastq which is the forward and _R2_001.fastq which is the reverse on Illumina platform

# Let's make a directory to save the results.
setwd("/home/sarkar/Documents/microbiome_workshop/demo_results_new")

# Let's locate a directory where all the fastq files are stored and give it a name.
# If the fastq files are in zipped format, they should be unzipped first

demo_microbiome_fasqfiles <- "/home/sarkar/Documents/microbiome_workshop/demo_fastq/demofastqsamples"

# Let's load the appropriate R packages (dada2 and few others) for the analysis. It should be installed in advance
library(dada2)
library(GUniFrac)
library(simEd)
library(tictoc)
# It is important to note down the package version of each of the required R package
packageVersion("dada2")
packageVersion("GUniFrac")
packageVersion("simEd")
packageVersion("tictoc")


# In order to make the results reproducible when carried out multiple times, set.seed function is used.
set.seed(56456)

# We should check if the fastq source directory contains all our fastq files
list.files(demo_microbiome_fasqfiles)

# Now we are going to make two sub-directories to separate the forward and the reverse reads
demo_F <- sort(list.files(demo_microbiome_fasqfiles, pattern="_R1_001.fastq", full.names = TRUE))
demo_R <- sort(list.files(demo_microbiome_fasqfiles, pattern="_R2_001.fastq", full.names = TRUE))

# Extract filenames of all samples for future
demo_samplenames <- sapply(strsplit(basename(demo_F), "_"), '[', 1)

# Let's check how is our data with quality plots (first let's check the forward)
pdf('demo_F_quality.pdf', width = 12, height = 8, pointsize = 8)
plotQualityProfile(demo_F[1:4], n=1e+06)
dev.off()

# Now let's check how is our data for reverse reads
pdf('demo_R_quality.pdf', width = 12, height = 8, pointsize = 8)
plotQualityProfile(demo_R[1:4], n=1e+06)
dev.off()


# The next step is to filter the sequences appropriately depending on the data. 
# Our target is to discard the bad sequences, trim the ends and save the good reads to a new directory

# The first step here is to make a directory to which the good sequences will be stored
demo_goodF <- file.path(demo_microbiome_fasqfiles, "demo_good_filtered", paste0(demo_samplenames, "F_good.fastq.gz"))
demo_goodR <- file.path(demo_microbiome_fasqfiles, "demo_good_filtered", paste0(demo_samplenames, "R_good.fastq.gz"))
names(demo_goodF) <- demo_samplenames
names(demo_goodR) <- demo_samplenames

# This is the very important step of filter and trimming each fastq. These parameters are flexible and should depend on your data
# We also want to know how much time the most important steps take for this task. 
# So, we are using tictoc to find it out
# Let's find out how much time it takes to clean the data

tic(msg = NULL, quiet = TRUE)
demo_good_proper <- filterAndTrim(demo_F, demo_goodF, demo_R, demo_goodR, trimLeft = c(17, 21), truncLen = c(145, 135), maxN = 0, truncQ = 2, minQ=1, maxEE = c(2, 4), rm.phix = TRUE, n = 1e+5, compress = TRUE, verbose = TRUE)

# save the output of previous step
write.table(demo_good_proper, "demo_filteredout.txt", sep = "\t")

datacleantime <- toc(log = FALSE, quiet = TRUE)

# So, the total time to clean the data is
time_clean <- datacleantime$toc - datacleantime$tic

# Let's find out how much time it takes to generate the error models from the data
tic(msg = NULL, quiet = TRUE)

# Now let us calculate the error rates for the forward and reverse sequences.
demo_error_F <- learnErrors(demo_goodF, nbases = 1e+07, randomize = TRUE, MAX_CONSIST = 12, multithread = TRUE, verbose = TRUE)

# Let us plot the forward error rate
pdf('demo_error_F_plot.pdf', width = 10, height = 10, pointsize = 8)
plotErrors(demo_error_F, obs = TRUE, nominalQ = TRUE)
dev.off()

# Let us calculate the reverse error rate and plot the error graph
demo_error_R <- learnErrors(demo_goodR, nbases = 1e+07, randomize = TRUE, MAX_CONSIST = 12, multithread = TRUE, verbose = TRUE)
pdf('demo_error_R_plot.pdf', width = 10, height = 10, pointsize = 8)
plotErrors(demo_error_R, obs = TRUE, nominalQ = TRUE)
dev.off()

errmodel <- toc(log = FALSE, quiet = TRUE)

# So, the total time to get the error models from the data is
time_errormodel <- errmodel$toc - errmodel$tic

# The next step is to dereplicate the sequences in each sample
derep_demo_F <- derepFastq(demo_goodF, n = 1e+06, verbose = TRUE)
derep_demo_R <- derepFastq(demo_goodR, n = 1e+06, verbose = TRUE)

names(derep_demo_F) <- demo_samplenames
names(derep_demo_R) <- demo_samplenames

# Now it is time to run the actual algorithm of dada2 to determine the ASVs in the dataset
# This step is run separately for forward and reverse sets

# We also want to know how much time this steps takes:
tic(msg = NULL, quiet = TRUE)
demo_dada_F <- dada(derep_demo_F, err=demo_error_F, pool = TRUE, multithread = TRUE)
demo_dada_R <- dada(derep_demo_R, err=demo_error_R, pool = TRUE, multithread = TRUE)

asvtime <- toc(log = FALSE, quiet = TRUE)

# So, the total time to generate the ASVs from the data is
time_asv_generation <- asvtime$toc - asvtime$tic

# Let's see how many sequence variants we have got in the forward set
demo_dada_F[[1]]

# Now we are going to merge the forward and the reverse sets (the paired-end reads)

# Again, we are going to calculate the time for this step
tic(msg = NULL, quiet = TRUE)
demo_merged <- mergePairs(demo_dada_F, derep_demo_F, demo_dada_R, derep_demo_R, minOverlap = 20, maxMismatch = 0, verbose = TRUE)

mergetime <- toc(log = FALSE, quiet = TRUE)

time_merge_reads <- mergetime$toc - mergetime$tic

# Let's make a sequence table of all the ASVs
demo_sequence_table <- makeSequenceTable(demo_merged, orderBy = "abundance")

# We can check the distribution of the ASVs by length
table(nchar(getSequences(demo_sequence_table)))

# Another important step: Remove the chimeric sequences

tic(msg = NULL, quiet = TRUE)
demo_nochim <- removeBimeraDenovo(demo_sequence_table, method = "consensus", minFoldParentOverAbundance = 1, verbose = TRUE, multithread = TRUE)

chimeratime <- toc(log = FALSE, quiet = TRUE)

time_chimera_removal <- chimeratime$toc - chimeratime$tic

# Let's see how many ASVs remains
dim(demo_nochim)

# Let's see what proportion of sequences we retained after filtering for chimera
sum(demo_nochim)/sum(demo_sequence_table)

# Now we have gone through all the filtering, trimming, cleanup etc. Here we can check how many sequences we were able to retain after each step.
# This test is important for trouble-shooting purposes. Let's work this out with a function
fetch_numbers <- function(a) sum(getUniques(a))
demo_track_steps <- cbind(demo_good_proper, sapply(demo_dada_F, fetch_numbers), sapply(demo_dada_R, fetch_numbers), sapply(demo_merged, fetch_numbers), rowSums(demo_nochim))
colnames(demo_track_steps) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nochim")
rownames(demo_track_steps) <- demo_samplenames

# Let's save the output to a new file
write.table(demo_track_steps, "demo_filtering_steps_track.txt", sep = "\t")

# The next step is to assign taxonomy to all the ASVs
# We will use the Silva database v.132 for this purpose

tic(msg = NULL, quiet = TRUE)
demo_taxonomy <- assignTaxonomy(demo_nochim, "/home/sarkar/Documents/microbiome_workshop/demo_results/silva_nr_v132_train_set.fa", minBoot = 80, verbose = TRUE, multithread = TRUE)
write.table(demo_taxonomy, "demo_taxaout.txt", sep = "\t")

taxotime <- toc(log = FALSE, quiet = TRUE)

time_taxonomy <- taxotime$toc - taxotime$tic

# Let's create a table by replacing the ASV sequences with ids (ASV_1, ASV_2 etc.) and their corressponding classifications
demo_taxa_summary <- demo_taxonomy
row.names(demo_taxa_summary) <- NULL
head(demo_taxa_summary)

# Let's make a file listing all the ASVs and their sequences in fasta format
demo_asv_seqs <- colnames(demo_nochim)
demo_asv_headers <- vector(dim(demo_nochim)[2], mode = "character")
for (i in 1:dim(demo_nochim)[2]) {demo_asv_headers[i] <- paste(">ASV", i, sep = "_")}
demo_asv.fasta <- c(rbind(demo_asv_headers, demo_asv_seqs))
write(demo_asv.fasta, "demo_out_asv.fasta")

# At this step, we need to make a table of ASV counts for each sample (which is going to be most important for all statistical analyses)
demo_asv_tab <- t(demo_nochim)
row.names(demo_asv_tab) <- sub(">", "", demo_asv_headers)
write.table(demo_asv_tab, "demo_asv_counts.tsv", sep = "\t", quote=F, col.names = NA)

# Finally, let's make a table with the taxonomy of all the ASVs
demo_asv_taxa <- demo_taxonomy
row.names(demo_asv_taxa) <- sub(">", "", demo_asv_headers)
write.table(demo_asv_taxa, "demo_asvs_taxonomy.tsv", sep = "\t", quote=F, col.names = NA)
dim(demo_asv_taxa)

# Perform rarefaction to subsample equal number of reads to reduce bias

demo_asv_transpose <- t(demo_asv_tab)

demo_rarefied <- Rarefy(demo_asv_transpose, depth = min(rowSums(demo_asv_transpose)))

write.table(demo_rarefied$otu.tab.rff, "demo_rar_table.txt", sep = "\t")

write.table(demo_rarefied$discard, "demo_rar_discard.txt", sep = "\t")

# Remove those ASVs whose sum is zero in the total dataset

demo_final <- demo_rarefied$otu.tab.rff[, which(colSums(demo_rarefied$otu.tab.rff) != 0)]

write.table(demo_final, "demo_final_count_table.txt", sep = "\t")

# Let's see how much time the more important steps take
cat("\nThe total time for data clean:", time_clean , "seconds\n")
cat("The total time for error models:", time_errormodel , "seconds\n")
cat("The total time for generating ASVs:", time_asv_generation , "seconds\n")
cat("The total time for merging paired-end reads:", time_merge_reads , "seconds\n")
cat("The total time to remove chimera:", time_chimera_removal , "seconds\n")
cat("The total time to assign taxonomies to ASVs:", time_taxonomy , "seconds\n")

# That's the end of this task
## ** THANK YOU EVERYONE ** ##



