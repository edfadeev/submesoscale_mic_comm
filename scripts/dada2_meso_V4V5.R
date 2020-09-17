#load libraries and set random seed
library(ggplot2); packageVersion("ggplot2")
library(dada2); packageVersion("dada2")
library(gridExtra); packageVersion("gridExtra")
library(DECIPHER); packageVersion("DECIPHER")
library(phangorn); packageVersion("phangorn")
set.seed(123)


# Sort ensures forward/reverse reads are in same order
fnFs <- sort(list.files("Clipped", pattern="_R1.fastq", full.names = TRUE))
fnRs <- sort(list.files("Clipped", pattern="_R2.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_clip_R1.fastq"), `[`, 1)

filt_path <- file.path("Report")
if(!file_test("-d", filt_path)) dir.create(filt_path)

# quality check
QualityProfileFs <- list()
for(i in 1:length(fnFs)) {
  QualityProfileFs[[i]] <- list()
  QualityProfileFs[[i]][[1]] <- plotQualityProfile(fnFs[i])
}
pdf(file.path("Report","RawProfileForward.pdf"))
for(i in 1:length(fnFs)) {
  do.call("grid.arrange", QualityProfileFs[[i]])  
}
dev.off()
rm(QualityProfileFs)

QualityProfileRs <- list()
for(i in 1:length(fnRs)) {
  QualityProfileRs[[i]] <- list()
  QualityProfileRs[[i]][[1]] <- plotQualityProfile(fnRs[i])
}
pdf(file.path("Report","RawProfileReverse.pdf"))
for(i in 1:length(fnRs)) {
  do.call("grid.arrange", QualityProfileRs[[i]])  
}
dev.off()
rm(QualityProfileRs)

# Make directory and filenames for the filtered fastqs
filt_path <- file.path("Filtered")
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path("Filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path("Filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

#separate the different runs
fnFs_hi<- sort(file.path("Clipped",paste(c(1:55), "_clip_R1.fastq", sep = "")))
fnRs_hi<- sort(file.path("Clipped",paste(c(1:55), "_clip_R2.fastq", sep = ""))) 
filtFs_hi<- sort(file.path("Filtered",paste(c(1:55), "_F_filt.fastq.gz", sep = "")))
filtRs_hi<- sort(file.path("Filtered",paste(c(1:55), "_R_filt.fastq.gz", sep = ""))) 
#Filter and trim
out_hi <- filterAndTrim(fnFs_hi, filtFs_hi, fnRs_hi, filtRs_hi, truncLen = c(230, 195),
                        maxN=0, maxEE=c(3,3), truncQ=2, rm.phix=TRUE,
                        compress=TRUE, multithread=TRUE)
# Learn errors 
errF_hi <- learnErrors(filtFs_hi, multithread = TRUE, randomize = TRUE, MAX_CONSIST = 30, verbose = TRUE)
errR_hi <- learnErrors(filtRs_hi, multithread = TRUE, randomize = TRUE, MAX_CONSIST = 30, verbose = TRUE)
# Sample Inference 
dadaFs_hi <- dada(filtFs_hi, err=errF_hi, multithread=TRUE, verbose = TRUE)
dadaRs_hi <- dada(filtRs_hi, err=errR_hi, multithread=TRUE, verbose = TRUE)
#Merge paired reads
mergers_hi <- mergePairs(dadaFs_hi, filtFs_hi, dadaRs_hi, filtRs_hi, verbose=TRUE, minOverlap = 10)

#493:432
fnFs_mi1 <- sort(file.path("Clipped",paste(c(56:60), "_clip_R1.fastq", sep = "")))
fnRs_mi1 <- sort(file.path("Clipped",paste(c(56:60), "_clip_R2.fastq", sep = "")))
filtFs_mi1 <- sort(file.path("Filtered",paste(c(56:60), "_F_filt.fastq.gz", sep = "")))
filtRs_mi1 <- sort(file.path("Filtered",paste(c(56:60), "_R_filt.fastq.gz", sep = "")))
#Filter and trim
out_mi1 <- filterAndTrim(fnFs_mi1, filtFs_mi1, fnRs_mi1, filtRs_mi1, truncLen = c(230, 195),
                         maxN=0, maxEE=c(3,3), truncQ=2, rm.phix=TRUE,
                         compress=TRUE, multithread=TRUE)
# Learn errors
errF_mi1 <- learnErrors(filtFs_mi1, multithread = TRUE, randomize = TRUE, MAX_CONSIST = 30, verbose = TRUE)
errR_mi1 <- learnErrors(filtRs_mi1, multithread = TRUE, randomize = TRUE, MAX_CONSIST = 30, verbose = TRUE)
# Sample Inference 
dadaFs_mi1 <- dada(filtFs_mi1, err=errF_mi1, multithread=TRUE, verbose = TRUE)
dadaRs_mi1 <- dada(filtRs_mi1, err=errR_mi1, multithread=TRUE, verbose = TRUE)
#Merge paired reads
mergers_mi1 <- mergePairs(dadaFs_mi1, filtFs_mi1, dadaRs_mi1, filtRs_mi1, verbose=TRUE, minOverlap = 10)

##summary of filtering

# quality check
QualityProfileFs <- list()
for(i in 1:length(filtFs)) {
  QualityProfileFs[[i]] <- list()
  QualityProfileFs[[i]][[1]] <- plotQualityProfile(filtFs[i])
}
pdf(file.path("Report","FiltProfileForward.pdf"))
for(i in 1:length(filtFs)) {
  do.call("grid.arrange", QualityProfileFs[[i]])  
}
dev.off()
rm(QualityProfileFs)

QualityProfileRs <- list()
for(i in 1:length(filtRs)) {
  QualityProfileRs[[i]] <- list()
  QualityProfileRs[[i]][[1]] <- plotQualityProfile(filtRs[i])
}
pdf(file.path("Report","FiltProfileReverse.pdf"))
for(i in 1:length(filtRs)) {
  do.call("grid.arrange", QualityProfileRs[[i]])  
}
dev.off()
rm(QualityProfileRs)

# Plot error profiles
pdf(file.path("Report","ErrorProfiles.pdf"))
plotErrors(errF_hi, nominalQ = TRUE)+ggtitle("FWD-MiSeq- 493:415")
plotErrors(errF_mi1, nominalQ = TRUE)+ggtitle("FWD-MiSeq- 493:432")
plotErrors(errR_hi, nominalQ = TRUE)+ggtitle("REV-MiSeq- 493:415")
plotErrors(errR_mi1, nominalQ = TRUE)+ggtitle("REV-MiSeq- 493:432")
dev.off()


#write out filtered read counts
write.csv(rbind(out_hi,out_mi1), 
          file= file.path("Report","dada2_filterAndTrim_output.csv"))

#merge the hiseq and miseq into a single sequence table
seqtab<- mergeSequenceTables(table1= makeSequenceTable(mergers_hi), 
                             table2 = makeSequenceTable(mergers_mi1))

save.image("MESO_V4V5_dada2_sep_runs.Rdata")

#Combine together sequences that are identical 
seqtab1 <- collapseNoMismatch(seqtab, verbose = TRUE)

dim(seqtab1)

save.image("MESO_V4V5_dada2_sep_runs.Rdata")

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab1)))

seqtab.nochim <- removeBimeraDenovo(seqtab1, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

#proportion of chimeras
sum(seqtab.nochim)/sum(seqtab1)

# inspect output: remove singletons and 'junk' sequences
# read lengths modified for V34 amplicons / based upon output table where majority of reads occurs
seqtab.nochim2 <- seqtab.nochim[, nchar(colnames(seqtab.nochim)) %in% c(360:400) & colSums(seqtab.nochim) > 1]
dim(seqtab.nochim2)
summary(rowSums(seqtab.nochim2)/rowSums(seqtab.nochim))

#assign taxonomy
taxa <- assignTaxonomy(seqtab.nochim2, "../tax/silva_nr_v138_train_set.fa.gz", multithread=TRUE, tryRC = TRUE, verbose = TRUE)
taxa <- addSpecies(taxa, "../tax/silva_species_assignment_v138.fa.gz", tryRC = TRUE, verbose = TRUE)

save.image("MESO_V4V5_dada2_sep_runs.Rdata")


sample.order <- names(c(dadaFs_mi1,dadaFs_mi2,dadaFs_mi3,dadaFs_mi4,dadaFs_mi5))

# get summary tables 
getN <- function(x) sum(getUniques(x))
track <- cbind(rbind(out_mi1,out_mi2,out_mi3,out_mi4,out_mi5), 
               sapply(c(dadaFs_mi1,dadaFs_mi2,dadaFs_mi3,dadaFs_mi4,dadaFs_mi5), getN),
               sapply(c(mergers_mi1,mergers_mi2,mergers_mi3,mergers_mi4,mergers_mi5), getN),
               rowSums(seqtab.nochim[sample.order,]),
               rowSums(seqtab.nochim2[sample.order,]))
colnames(track) <- c("input", "filtered", "denoised", "merged", "nochim", "tabled")
rownames(track) <- sample.names
track <- data.frame(track)

# write output
write.table(track, file.path("Report","libs_summary.csv"),sep = "\t", quote = F)
write.table(t(seqtab.nochim2), file.path("Report","seqtab.txt"), quote = F, sep = "\t")
write.table(taxa, file.path("Report","taxonomy_table.txt"), sep = "\t", quote = F)



