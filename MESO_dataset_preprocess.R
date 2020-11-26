#load libraries
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")

#load functions
source('./scripts/color_palettes.R')
source('./scripts/extra_functions.R')

#####################################
#Parse for Phyloseq
#####################################
ASV_tab<- read.csv("./data/MESO_seqtab.txt", h=T, sep="\t")
#correct sample names in ASV table
names(ASV_tab)<- gsub("_F_filt.fastq.gz","",names(ASV_tab))

TAX<- as.matrix(read.csv("./data/MESO_taxonomy_table.txt", h=T,sep = "\t"))
#add unclassified taxonomy levels 
TAX<-add_uncl(TAX)

ENV <- read.csv("./data/MESO_samples_meta.csv", sep = "," , h = T, row.names = 1, fill = T, na.strings=c("","NA"))
#reorder samples
ENV <- ENV[names(ASV_tab),]


# Check order of samples
all.equal(rownames(ASV_tab), rownames(TAX))

#creating Phyloseq dataset
meso_ps <- phyloseq(otu_table(ASV_tab, taxa_are_rows=TRUE), 
                    sample_data(ENV), 
                    tax_table(TAX))

#add reference sequence and replace variants with ASVs
dna <- Biostrings::DNAStringSet(taxa_names(meso_ps))
names(dna) <- taxa_names(meso_ps)
meso_ps <- merge_phyloseq(meso_ps, dna)
taxa_names(meso_ps) <- paste0("ASV", seq(ntaxa(meso_ps)))

#remove unclassified on phylum level, chloroplast and Mitochondrial sequence variants
meso_ps0 <- subset_taxa(meso_ps, !Kingdom %in% c("Eukaryota") & !Phylum %in% c("NA_uncl") & !Order %in% c("Chloroplast") & !Family %in% c("Mitochondria") )


#####################################
#Fix categories 
#####################################
sample_data(meso_ps0)$Community<- factor(
  sample_data(meso_ps0)$Community, 
  levels = c("FL", "PA"))

sample_data(meso_ps0)$StationName<- factor(
  sample_data(meso_ps0)$StationName, 
  levels = c("T4","T1","T2","T5","T3"))

sample_data(meso_ps0)$Type <- factor(sample_data(meso_ps0)$Type, 
                                      levels = c("Surface-10","Chl.max-20-30","B.Chl.max-50",
                                                 "Epipelagic-100","Mesopelagic-200","Mesopelagic-400"))

#separation of surface and deep
sample_data(meso_ps0)$layers <- "down"
sample_data(meso_ps0)$layers[sample_data(meso_ps0)$Depth < 100] <- "up"

#####################################
#Plot total number of reads and OTUs per sample
#####################################
readsumsdf <- data.frame(nreads = sort(taxa_sums(meso_ps0), TRUE), sorted = 1:ntaxa(meso_ps0), 
                         type = "ASVs")
readsumsdf <- rbind(readsumsdf, data.frame(nreads = sort(sample_sums(meso_ps0), 
                                                         TRUE), sorted = 1:nsamples(meso_ps0), type = "Samples"))
title = "Total number of reads"
p <- ggplot(readsumsdf, aes(x = sorted, y = nreads)) + geom_bar(stat = "identity")
OTU_libs_overview <-  p + ggtitle(title) + scale_y_log10() + facet_wrap(~type, 1, scales = "free")

#####################################
# Compute prevalence of each feature
#####################################
prevdf <-  apply(X = otu_table(meso_ps0),
                 MARGIN = ifelse(taxa_are_rows(meso_ps0), yes = 1, no = 2),
                 FUN = function(x){sum(x > 0)})

# # Add taxonomy and total read counts to this data.frame
prevdf.tax  <-  data.frame(Prevalence = prevdf,
                           TotalAbundance = taxa_sums(meso_ps0),
                           tax_table(meso_ps0))
#summarize
prevdf.tax.summary <- plyr::ddply(prevdf.tax, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})

# 
# #plot
prev_plot_phyl <- ggplot(prevdf.tax, aes(TotalAbundance, Prevalence / nsamples(meso_ps0),color=Phylum)) +
  # # Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")

#  Define prevalence threshold as 5% of total samples
prevalenceThreshold <- round(0.05 * nsamples(meso_ps0))
prevalenceThreshold

# Execute prevalence filter, using `prune_taxa()` function
meso_ps0.prev <-  prune_taxa((prevdf > prevalenceThreshold), meso_ps0)

#####################################
# export R-native serialized RDS file
#####################################
saveRDS(meso_ps0, "./Data/MESO_ps0.rds")
saveRDS(meso_ps0.prev, "./Data/MESO_ps0_prev.rds")

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