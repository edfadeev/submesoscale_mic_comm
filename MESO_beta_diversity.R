#load libraries
library(ggplot2); packageVersion("ggplot2")
library(dplyr); packageVersion("dplyr")
library(phyloseq); packageVersion("phyloseq")
library(vegan); packageVersion("vegan")


library(cowplot); packageVersion("cowplot")
library(reshape2); packageVersion("reshape2")
library(ggsignif); packageVersion("ggsignif")
library(ggrepel); packageVersion("ggrepel")


#load functions
source('./scripts/color_palettes.R')
source('./scripts/extra_functions.R')

#####################################
#Load phyloseq object
####################################
meso_ps0.prev <-  readRDS("./Data/MESO_ps0_prev.rds")


#####################################
#Plot barplots of communities
#####################################
#transform data
BAC_pruned.ra <- transform_sample_counts(meso_ps0.prev, function(x) x / sum(x))
BAC_pruned.ra.long <- psmelt(BAC_pruned.ra)

#calculate abundance for each Class and replace classes below 1% with "Other taxa"
BAC_pruned.ra.long.agg <- BAC_pruned.ra.long %>% select(StationName,Community,Type,Class,Abundance)%>%
                        group_by(StationName,Community,Type,Class) %>%
                        summarize(Abund.total= sum(Abundance)*100)  %>% 
                        filter(Abund.total>0)

taxa_classes <- unique(BAC_pruned.ra.long.agg$Class)

BAC_pruned.ra.long.agg<- BAC_pruned.ra.long.agg %>% 
                        mutate(Class = ifelse(Abund.total < 2, "Other taxa < 2%", Class))%>% 
                              mutate(Class =factor(Class, levels = c(taxa_classes,"Other taxa < 2%")))

#Plot 
barplots <- ggplot(BAC_pruned.ra.long.agg, aes(x = StationName, y = Abund.total, fill = Class)) + 
  facet_grid(Type~Community, space= "fixed") +
  geom_col()+
  scale_fill_manual(values = tol21rainbow) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Sequence proportions (%) \n")+
  theme_classic()+
  theme(legend.position = "bottom")

ggsave("./figures/barplot-prev.pdf", barplots, dpi = 300, 
       width = 30, height = 30, 
       units = "cm")

#####################################
#NMDS plot for all samples
#####################################
#transform the ASV table using geometric mean
meso_ps0.gm_mean <- phyloseq_gm_mean_trans(meso_ps0.prev)

#NMDS plot
meso_ps0.gm_mean.ord <- ordinate(meso_ps0.gm_mean, method = "NMDS", distance = "euclidean")
meso_ps0.gm_mean.df <- plot_ordination(meso_ps0.gm_mean, meso_ps0.gm_mean.ord, axes = c(1,2,3),justDF = TRUE)

meso_ps0.NMDS.p <- ggplot(data = meso_ps0.gm_mean.df, aes(x = NMDS1, y = NMDS2, shape = Community, colour = Group))+
  geom_point(colour = "black", size = 5) +
  geom_point(size = 4) +
  geom_text(aes(x = NMDS1, y = NMDS2,label = paste(StationName, paste(Depth,"m",sep =""),sep="-")), 
            nudge_y= -8,size=3, colour = "black")+
  scale_colour_manual(values = c("in"="red1",
                                 "out"="blue1")) + 
  annotate(geom="text", x=-80, y=90, label= paste0("Stress = ", round(meso_ps0.gm_mean.ord$stress,2)),
           color="black", size = 5)+
  theme_classic()+
  theme(legend.position = "bottom")


ggsave("./figures/NMDS_prev.pdf", meso_ps0.NMDS.p, dpi = 300, 
       #width = 11.4, height = 23, 
       units = "cm")

#####################################
#NMDS plot of upper 50 m
#####################################
#subset the surface samples
meso_ps0.gm_mean_up50m <- subset_samples(meso_ps0.gm_mean, Type %in% c("Surface-10", 
                                                                    "Chl.max-20-30",
                                                                    "B.Chl.max-50"))
#NMDS plot
meso_ps0.gm_mean_up50m.ord <- ordinate(meso_ps0.gm_mean_up50m, method = "NMDS", distance = "euclidean")
meso_ps0.gm_mean_up50m.df <- plot_ordination(meso_ps0.gm_mean_up50m, meso_ps0.gm_mean_up50m.ord, axes = c(1,2,3),justDF = TRUE)

meso_ps0.NMDS_up50m.p <- ggplot(data = meso_ps0.gm_mean_up50m.df, aes(x = NMDS1, y = NMDS2, shape = Community, colour = Group))+
  geom_point(colour = "black", size = 5) +
  geom_point(size = 4) +
  geom_text(aes(x = NMDS1, y = NMDS2,label = paste(StationName, paste(Depth,"m",sep =""),sep="-")), 
            nudge_y= -2,size=2, colour = "black")+
  scale_colour_manual(values = c("in"="red1",
                                 "out"="blue1")) + 
  annotate(geom="text", x=-70, y=50, label= paste0("Stress = ", round(meso_ps0.gm_mean_up50m.ord$stress,2)),
           color="black", size = 5)+
  theme_classic()+
  theme(legend.position = "bottom")

ggsave("./figures/NMDS_prev_surf.pdf", meso_ps0.NMDS_up50m.p, dpi = 300, 
       #width = 11.4, height = 23, 
       units = "cm")

#####################################
#Statistics on comm. composition
#####################################
#significance test
df <- as(sample_data(meso_ps0.gm_mean), "data.frame")
d <- phyloseq::distance(meso_ps0.gm_mean, "euclidean")
adonis_all <- adonis2(d ~ Community + Group + Type, df)
adonis_all


#significance test on only surface samples
df.up50m <- as(sample_data(meso_ps0.gm_mean_up50m), "data.frame")
d.up50m <- phyloseq::distance(meso_ps0.gm_mean_up50m, "euclidean")
adonis_up50m <- adonis2(d.up50m ~ Community + Group + Type, df.up50m)
adonis_up50m


#significance test below 50 m
meso_ps0.gm_mean_below50m <- subset_samples(meso_ps0.gm_mean, Type %in% c("Epipelagic-100","Mesopelagic-200","Mesopelagic-400"))
df.below50m <- as(sample_data(meso_ps0.gm_mean_below50m), "data.frame")
d.below50m <- phyloseq::distance(meso_ps0.gm_mean_below50m, "euclidean")
adonis_below50m <- adonis2(d.below50m ~ Community + Group + Type, df.below50m)
adonis_below50m

#####################################
#Comm. composition enrichments
#####################################
res_all <- data.frame()

for (frac in c("FL","PA")){
  PS107_merged.prev.frac <- subset_samples(meso_ps0.prev, layers =="up")
  PS107_merged.prev.frac <- subset_samples(PS107_merged.prev.frac, Community ==frac)
  PS107_merged.prev.frac <- prune_taxa(taxa_sums(PS107_merged.prev.frac)>0,PS107_merged.prev.frac)
  
  #run DEseq
  BAC_sed.ddsMat <- phyloseq_to_deseq2(PS107_merged.prev.frac, ~Group)
  geoMeans = apply(counts(BAC_sed.ddsMat), 1, gm_mean)
  BAC_sed.ddsMat = estimateSizeFactors(BAC_sed.ddsMat, geoMeans = geoMeans)
  BAC_sed.DEseq = DESeq(BAC_sed.ddsMat, fitType="local")
  
  BAC_sed.DEseq.res <- results(BAC_sed.DEseq)
  
  BAC_sed.DEseq.res = BAC_sed.DEseq.res[order(BAC_sed.DEseq.res$padj, na.last=NA), ]
  alpha = 0.1
  sigtab = BAC_sed.DEseq.res[(BAC_sed.DEseq.res$padj < alpha), ]
  sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(meso_ps0.prev)[rownames(sigtab), ], "matrix"))
  
  sigtab$frac <- frac
  
  res_all <- rbind(res_all,sigtab, make.row.names = TRUE)
}

#aggregate the results
res_all.agg.in <- do.call(data.frame, aggregate(log2FoldChange~ Class+Order+Family + frac, res_all[res_all$log2FoldChange<0,], function(x) c(mean = mean(x), se = se(x), n = length(x))))
res_all.agg.out <- do.call(data.frame, aggregate(log2FoldChange~ Class+Order+Family + frac, res_all[res_all$log2FoldChange>0,], function(x) c(mean = mean(x), se = se(x), n = length(x))))

res_all.agg <- rbind(res_all.agg.in,res_all.agg.out)
Family <- res_all.agg$Family[res_all.agg$log2FoldChange.n>2]
res_all_for_plot <- res_all[res_all$Family %in% Family,]

#plot
meso_prev.enr.p<- ggplot(res_all_for_plot, aes(y=log2FoldChange , x=Family, colour = Class))+ 
  geom_point(aes(shape = frac),size = 4, colour = "black")+
  geom_point(aes(shape = frac),size = 3)+
  geom_boxplot(data= res_all_for_plot[res_all_for_plot$log2FoldChange>0,], aes(y=log2FoldChange , x=Family), colour = "gray", fill = NA, alpha = 0.5, size = 1,outlier.size = 0)+ 
  geom_boxplot(data= res_all_for_plot[res_all_for_plot$log2FoldChange<0,], aes(y=log2FoldChange , x=Family),colour = "gray", fill = NA, alpha = 0.5, size = 1,outlier.size = 0)+ 
  scale_colour_manual(values = phyla.col)+ 
  #scale_x_discrete("Class",expand = waiver())+
  geom_hline(aes(yintercept=0), linetype="dashed")+
  theme_classic()+
  theme(legend.position = "bottom", axis.text.x = element_text(angle =90))+
  coord_flip()+
  ylim(-26,26)+
  geom_vline(aes(xintercept=Inf)) + 
  geom_hline(aes(yintercept=Inf))

ggsave("./figures/enrichment.pdf", meso_prev.enr.p, dpi = 300, 
       #width = 11.4, height = 23, 
       units = "cm")

write.table(res_all.agg, "data/enriched_taxa.txt")

#####################################
#get session info and remove all objects and libraries
#####################################
sessionInfo()

rm(list = ls(all = TRUE))
