######################################################################
######################################################################
#This script will reproduce the Rdata files 
#that are already in the repository!
#Therefore, no need to run it.
######################################################################
######################################################################

#set working directory
wd <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd)

#install required packages
#.inst_packages <- c("ggplot2", "gridExtra", "phyloseq",
#                    "igraph","UpSetR","vegan","cowplot",
#                    "dplyr","ggpmisc","VennDiagram",
#                    "reshape2","DESeq2","tidyr", "iNEXT", "olsrr")

#.inst <- .inst_packages %in% installed.packages()
#if(any(!.inst)) {
#  BiocManager::install(.inst_packages, ask = F)
#}

#load libraries
library(phyloseq); packageVersion("phyloseq")
library(vegan); packageVersion("vegan")
library(tidyr); packageVersion("tidyr")
library(DESeq2); packageVersion("DESeq2")
library(ggplot2); packageVersion("ggplot2")

#load functions
source('./scripts/miseqR.R')
source('./Scripts/clr_trans.R')

#####################################
#Parse for Phyloseq
#####################################
OTU<- read.csv("./dada2_output/bac-arch_seqtab_nochim2.txt", h=T, sep="\t")
TAX<- as.matrix(read.csv("./dada2_output/bac-arch_taxonomy_table.txt", h=T,sep = "\t"))
ENV <- read.csv("./dada2_output/MESO_samples_meta.txt", sep = "\t" , h = T, row.names = 1, fill = T, na.strings=c("","NA"))

# Check order of samples
all.equal(rownames(OTU), rownames(TAX))

#creating Phyloseq dataset
OTU <- otu_table(OTU, taxa_are_rows = TRUE)
TAX <- tax_table(TAX)
meta <- sample_data(ENV)
PS107_merged <- phyloseq(OTU, TAX, meta)

#####################################
#Fix categories 
#####################################
sample_data(PS107_merged)$Fraction<- factor(
  sample_data(PS107_merged)$Fraction, 
  levels = c("FL", "PA"))

sample_data(PS107_merged)$StationName<- factor(
  sample_data(PS107_merged)$StationName, 
  levels = c("T1","T2","3","T4","T5"))

#####################################
#Plot total number of reads and OTUs per sample
#####################################
readsumsdf <- data.frame(nreads = sort(taxa_sums(PS107_merged), TRUE), sorted = 1:ntaxa(PS107_merged), 
                         type = "ASVs")
readsumsdf <- rbind(readsumsdf, data.frame(nreads = sort(sample_sums(PS107_merged), 
                                                         TRUE), sorted = 1:nsamples(PS107_merged), type = "Samples"))
title = "Total number of reads"
p <- ggplot(readsumsdf, aes(x = sorted, y = nreads)) + geom_bar(stat = "identity")
OTU_libs_overview <-  p + ggtitle(title) + scale_y_log10() + facet_wrap(~type, 1, scales = "free")

#####################################
# Compute prevalence of each feature
#####################################
prevdf <-  apply(X = otu_table(PS107_merged),
                 MARGIN = ifelse(taxa_are_rows(PS107_merged), yes = 1, no = 2),
                 FUN = function(x){sum(x > 0)})

# Add taxonomy and total read counts to this data.frame
prevdf.tax  <-  data.frame(Prevalence = prevdf,
                      TotalAbundance = taxa_sums(PS107_merged),
                       tax_table(PS107_merged))
#summarize
prevdf.tax.summary <- plyr::ddply(prevdf.tax, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})

#plot
prev_plot_phyl <- ggplot(prevdf.tax, aes(TotalAbundance, Prevalence / nsamples(PS107_merged),color=Phylum)) +
# Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")

#  Define prevalence threshold as 5% of total samples
prevalenceThreshold <- round(0.05 * nsamples(PS107_merged))
prevalenceThreshold

# Execute prevalence filter, using `prune_taxa()` function
PS107_merged.prev <-  prune_taxa((prevdf > prevalenceThreshold), PS107_merged)

#######################################
# Variance stabilized Transformation
#######################################
#without prevalance
PS107_merged.dds <- phyloseq_to_deseq2(PS107_merged, ~1)
varianceStabilizingTransformation(PS107_merged.dds, blind = TRUE, fitType = "parametric")
PS107_merged.dds <- estimateSizeFactors(PS107_merged.dds)
PS107_merged.dds <- estimateDispersions(PS107_merged.dds)
otu.vst <- getVarianceStabilizedData(PS107_merged.dds)

#make sure that the dimensions of the OTU table and the DEseq object are matching
dim(otu.vst)
dim(otu_table(PS107_merged))

PS107_merged.vst<-PS107_merged
otu_table(PS107_merged.vst)<- otu_table(otu.vst, taxa_are_rows = TRUE)

#after prevalance
PS107_merged.prev.dds <- phyloseq_to_deseq2(PS107_merged.prev, ~1)
varianceStabilizingTransformation(PS107_merged.prev.dds, blind = TRUE, fitType = "parametric")
PS107_merged.prev.dds <- estimateSizeFactors(PS107_merged.prev.dds)
PS107_merged.prev.dds <- estimateDispersions(PS107_merged.prev.dds)
otu.vst <- getVarianceStabilizedData(PS107_merged.prev.dds)

#make sure that the dimensions aof the OTU table and the DEseq object are matching
dim(otu.vst)
dim(otu_table(PS107_merged.prev))

PS107_merged.prev.vst<-PS107_merged.prev
otu_table(PS107_merged.prev.vst)<- otu_table(otu.vst, taxa_are_rows = TRUE)

#####################################
# export R-native serialized RDS file
#####################################
saveRDS(PS107_merged.prev.vst, "./Data/PS107_merged_prev_vst.rds")
saveRDS(PS107_merged.vst, "./Data/PS107_merged_vst.rds")
saveRDS(PS107_merged, "./Data/PS107_merged.rds")
saveRDS(PS107_merged.prev, "./Data/PS107_merged_prev.rds")

#####################################
#get session info and remove all objects and libraries
#####################################
sessionInfo()

rm(list = ls(all = TRUE))

#unload libraries
detach("package:phyloseq")
detach("package:vegan")
detach("package:tidyr")
detach("package:DESeq2")
detach("package:ggplot2")