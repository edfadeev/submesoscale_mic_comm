#load libraries
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")

#load functions
source('./scripts/color_palettes.R')
source('./scripts/extra_functions.R')

#####################################
#Parse for Phyloseq the PS107 HAUSGARTEN surface samples
#####################################
ASV_tab<- read.csv("./data/PS107_SRF_dada2/PS107_SRF_seqtab.txt", h=T, sep="\t")
#correct sample names in ASV table
names(ASV_tab)<- gsub("_F_filt.fastq.gz","",names(ASV_tab))

TAX<- as.matrix(read.csv("./data/PS107_SRF_dada2/PS107_SRF_taxonomy_table.txt", h=T,sep = "\t"))
#add unclassified taxonomy levels 
TAX<-add_uncl(TAX)

ENV <- read.csv("./data/PS107_SRF_dada2/PS107_SRF_samples_meta.csv", sep = "," , h = T, row.names = 1, fill = T, na.strings=c("","NA"))
#reorder samples
ENV <- ENV[names(ASV_tab),]
ENV$Type<- "Chl.max-20-30"
ENV$Community.PS107<- "FL"
ENV$Station.PS107<- "FL"

# Check order of samples
all.equal(rownames(ASV_tab), rownames(TAX))

#creating Phyloseq dataset
PS107_SRF_ps <- phyloseq(otu_table(ASV_tab, taxa_are_rows=TRUE), 
                    sample_data(ENV), 
                    tax_table(TAX))

#add reference sequence and replace variants with ASVs
dna <- Biostrings::DNAStringSet(taxa_names(PS107_SRF_ps))
names(dna) <- taxa_names(PS107_SRF_ps)
meso_ps <- merge_phyloseq(PS107_SRF_ps, dna)
taxa_names(PS107_SRF_ps) <- paste0("ASV", seq(ntaxa(PS107_SRF_ps)))

#remove unclassified on phylum level, chloroplast and Mitochondrial sequence variants
PS107_SRF_ps0 <- subset_taxa(PS107_SRF_ps, !Kingdom %in% c("Eukaryota") & !Phylum %in% c("NA_uncl") & !Order %in% c("Chloroplast") & !Family %in% c("Mitochondria") )


#####################################
#Merge datasets on genus level
#####################################
meso_ps0.prev <-  readRDS("./Data/MESO_ps0_prev.rds")
meso_ps0.up50m <- subset_samples(meso_ps0.prev, Type %in% c("Surface-10"))
sample_data(meso_ps0.up50m)$merging <- paste(sample_data(meso_ps0.up50m)$StationName,sample_data(meso_ps0.up50m)$Depth, sep ="-")
meso_ps0.up50m.merged<- merge_samples(meso_ps0.up50m, "merging",  fun=sum)
meso_ps0.up50m.merged.glom<- tax_glom(meso_ps0.up50m.merged, "Genus")


sample_names(PS107_SRF_ps0)<- sample_data(PS107_SRF_ps0)$StationName
PS107_SRF_ps0.glom<- tax_glom(PS107_SRF_ps0, "Genus")

merged_ps<- merge_phyloseq(meso_ps0.up50m.merged.glom,PS107_SRF_ps0.glom)
sample_data(merged_ps)$StationName<- sample_names(merged_ps)

#transform data
BAC_pruned.ra <- transform_sample_counts(merged_ps, function(x) x / sum(x))
BAC_pruned.ra.long <- psmelt(BAC_pruned.ra)

#calculate abundance for each Class and replace classes below 1% with "Other taxa"
BAC_pruned.ra.long.agg <- BAC_pruned.ra.long %>% select(StationName, Community, Type,Class,Abundance)%>%
  group_by(StationName,Community,Type,Class) %>%
  summarize(Abund.total= sum(Abundance)*100)  %>% 
  filter(Abund.total>0)

taxa_classes <- unique(BAC_pruned.ra.long.agg$Class)

BAC_pruned.ra.long.agg<- BAC_pruned.ra.long.agg %>% 
  mutate(Class = ifelse(Abund.total < 2, "Other taxa < 2%", Class))%>% 
  mutate(Class =factor(Class, levels = c(taxa_classes,"Other taxa < 2%")),
         StationName = factor(StationName, levels = c("EG4","EG1", "T4-10", "T4-30", "T3-10", "T3-30", "0E", "N5","N4", "T5-10","T5-25",
                                                      "T1-10","T1-20", "T2-10","T2-30","HG2","HG1",
                                                      "S3","SV4","SV2")))

#Plot 
barplots <- ggplot(BAC_pruned.ra.long.agg, aes(x = StationName, y = Abund.total, fill = Class)) + 
  #facet_grid(Type~Community, space= "fixed", scales = "free") +
  geom_col()+
  scale_fill_manual(values = class_col) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Sequence proportions (%) \n")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),#axis.text.x = element_blank(),
        text=element_text(size=14),legend.position = "bottom")

ggsave("./figures/PS107_meso.pdf", barplots, dpi = 300, 
       #width = 30, height = 30, 
       units = "cm")



#####################################
#NMDS plot of upper 50 m
#####################################
BAC_pruned.ra <- transform_sample_counts(merged_ps, function(x) x / sum(x))

#NMDS plot
merged_ps.gm.ord <- ordinate(BAC_pruned.ra, method = "NMDS", distance = "bray")
merged_ps.gm.df <- plot_ordination(merged_ps, merged_ps.gm.ord, axes = c(1,2,3),justDF = TRUE)

PS107.NMDS_up50m.p <- ggplot(data = merged_ps.gm.df, aes(x = NMDS1, y = NMDS2, label  = StationName))+
  geom_point(colour = "black", size = 5) +
  geom_point(size = 4) +
  geom_text(nudge_y= -0.02,size=3, colour = "black")+
  scale_colour_manual(values = c("in"="red1",
                                 "out"="blue1")) + 
  annotate(geom="text", x=-0.6, y=0.7, label= paste0("Stress = ", round(merged_ps.gm.ord$stress,2)),
           color="black", size = 5)+
  geom_hline(yintercept = 0.15, linetype = 2, alpha = 0.5)+
  coord_fixed()+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),#axis.text.x = element_blank(),
        text=element_text(size=14),legend.position = "bottom")

ggsave("./figures/NMDS_SP107_surf.pdf", PS107.NMDS_up50m.p, dpi = 300, 
       #width = 11.4, height = 23, 
       units = "cm")


#significance test
df <- as(sample_data(BAC_pruned.ra), "data.frame") %>% 
          mutate(Bloom = case_when(StationName %in% c("EG1","EG4","T1-10","T1-20","T2-10","T2-30","T5-10","T5-25") ~ "YES",
                                                      TRUE~"NO"),
                 Dataset = case_when(StationName %in% c("T1-10","T1-20","T2-10","T2-30","T5-10","T5-25","T3-10","T3-30","T4-10","T4-30") ~ "sub",
                 TRUE~"PS107"))


d <- phyloseq::distance(BAC_pruned.ra, "bray")
adonis_all <- adonis2(d ~ Bloom + Dataset, df)
adonis_all

 

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