#load libraries
library("ggplot2"); packageVersion("ggplot2")
library("phyloseq"); packageVersion("phyloseq")
library("iNEXT"); packageVersion("iNEXT")
library("dplyr"); packageVersion("dplyr")
library("reshape2"); packageVersion("reshape2")
library("venn"); packageVersion("venn")

#load functions
source('./scripts/color_palettes.R')
source('./scripts/extra_functions.R')

#####################################
#Load phyloseq object
####################################
meso_ps0 <-  readRDS("./Data/meso_ps0.rds")

#####################################
#Plot rarefaction
####################################
iNEXT.out <- iNEXT(as.data.frame(otu_table(meso_ps0)), q=0, datatype="abundance")
rare <-fortify(iNEXT.out, type=1)

meta <- as(sample_data(meso_ps0), "data.frame")
meta$site <- rownames(meta)
rare$Community <- meta$Community[match(rare$site, meta$site)] 
rare$StationName <- meta$StationName[match(rare$site, meta$site)] 
rare$Type <- meta$Type[match(rare$site, meta$site)]
rare$label <- paste(rare$StationName,rare$Type, sep = "-")

rare.point <- rare[which(rare$method == "observed"),]
rare.line <- rare[which(rare$method != "observed"),]
rare.line$method <- factor (rare.line$method,
                            c("interpolated", "extrapolated"),
                            c("interpolation", "extrapolation"))

rare.p <- ggplot(rare, aes(x=x, y=y, colour = site))+
  geom_line(aes(linetype = method), lwd = 0.5, data= rare.line)+
  #geom_ribbon(aes(ymin=y.lwr, ymax= y.upr, colour = NULL), alpha = 0.2)+
  geom_point(aes(shape=Type, fill=Community), size =3, data= rare.point)+
  #geom_text(aes(label=label), size =2, data= rare.point, colour = "black", nudge_y = -100)+
  scale_colour_discrete(guide = FALSE)+
  labs(x = "Sample size", y = "Species richness")+
  facet_grid(~Community)+
  xlim(0,3e5)+
  #ylim(0,6000)+
  theme_classic(base_size = 12)+theme(legend.position="bottom")

ggsave("./figures/rarefaction.pdf", 
       plot = rare.p,
       units = "cm",
       width = 30, height = 30, 
       #scale = 1,
       dpi = 300)

#####################################
#Alpha diversity table
####################################
#dada2 workflow results
reads.tab <- read.csv2("./data/MESO_libs_summary.csv",header = TRUE, sep = "\t")
row.names(reads.tab) <-paste("X",row.names(reads.tab), sep ="")

#alpha div indeces
meso_ps0_alpha <- estimate_richness(meso_ps0, measures = c("Observed", "Chao1","Shannon", "InvSimpson"))
meso_ps0_alpha$Row.names<-rownames(meso_ps0_alpha)

#generate table
meso_ps0_summary_table <- merge(sample_data(meso_ps0),reads.tab,by =0) %>%
  merge(meso_ps0_alpha,by ="Row.names")%>%
  mutate_if(is.numeric, round, 2) %>%
  mutate(Seq.prop = round(tabled/input,2)) %>%
  select("Event_ID","StationName","Longitude","Latitude","Depth","Community", "input","tabled", "Seq.prop","Observed","Chao1","Shannon","InvSimpson")%>%
  rename("Event ID" = "Event_ID",
         "Fraction" = "Community",
         "Sampling depth [m]" = "Depth",
         "Raw seq." = "input",
         "Final seq." = "tabled",
         "Seq. proportions" = "Seq.prop",
         "Observed ASVs" = "Observed",
         "Chao1 richness est." = "Chao1",
         "Shannon Index" = "Shannon")

write.csv(meso_ps0_summary_table, "./Data/alpha_table.csv")

#####################################
#ASV distribution
#####################################
for (frac in c("FL","PA")){
  #mean number of ASV per sample
  meso_merged_sub <- subset_samples(meso_ps0, Community== frac)
  meso_merged_sub <- prune_taxa(taxa_sums(meso_merged_sub)>0, meso_merged_sub)
  
  vectorx <- vector()
  for (i in sample_names(meso_merged_sub)){
    meso_merged_station <- prune_samples(i, meso_merged_sub)
    meso_merged_station <- prune_taxa(taxa_sums(meso_merged_station)>0,meso_merged_station)
    vectorx <- c(vectorx,dim(otu_table(meso_merged_station))[1])
  }
  meso_bacteria <- subset_taxa(meso_merged_sub, Kingdom == "Bacteria")
  meso_bacteria <- prune_taxa(taxa_sums(meso_bacteria)>0,meso_bacteria)
  
  meso_archaea <- subset_taxa(meso_merged_sub, Kingdom == "Archaea")
  meso_archaea <- prune_taxa(taxa_sums(meso_archaea)>0,meso_archaea)
  
  assign(paste("df", frac, sep = "."), t(data.frame(total.reads= sum(sample_sums(meso_merged_sub)),
                                                    total.ASV= dim(otu_table(meso_merged_sub))[1],
                                                    total.BAC.ASV=dim(otu_table(meso_bacteria))[1],
                                                    total.ARC.ASV=dim(otu_table(meso_archaea))[1],
                                                    mean.ASV.per.sample=mean(vectorx),
                                                    sd.ASV.per.sample=se(vectorx))))
}

meso_tax.overview <- cbind(df.FL,df.PA)

#####################################
#Alpha diversity statistics
####################################
# Create new ps object with diversity estimates added to sample_data
meso_ps0_div <- merge_phyloseq(meso_ps0, sample_data(meso_ps0_alpha))
sampledata_DF <- data.frame(sample_data(meso_ps0_div))

# Run Shapiro normality test 
shapiro_test_Observed <- shapiro.test(sampledata_DF$Observed)
shapiro_test_Chao1 <- shapiro.test(sampledata_DF$Chao1)
shapiro_test_Shan <- shapiro.test(sampledata_DF$Shannon)

#ANOVA with Depth factor and Tukey’s HSD post-hoc test to determine which pairwise comparisons are different
aov.shannon <- aov(Shannon ~ Type, data = sampledata_DF)
TukeyHSD(aov.shannon)

aov.Chao1 <- aov(Chao1 ~ Type, data = sampledata_DF)
TukeyHSD(aov.Chao1)

aov.Observed <- aov(Observed  ~ Type, data = sampledata_DF)
TukeyHSD(aov.Observed)

#compare richness inside and outside
sampledata_DF_sub <- sampledata_DF[sampledata_DF$Depth< 60,]
t.test(Observed ~ Group, data = sampledata_DF_sub)

t.test(Chao1 ~ Group, data = sampledata_DF_sub)

#below 30 m
sampledata_DF_sub <- sampledata_DF[sampledata_DF$Depth> 60,]
t.test(Observed ~ Group, data = sampledata_DF_sub)

t.test(Chao1 ~ Group, data = sampledata_DF_sub)

#####################################
#Alpha diversity plots
####################################
sampledata_DF_long<- sampledata_DF %>% mutate(layer = case_when(Depth >60 ~ "Below 50m", Depth <60 ~ "Upper 50 m"),
                                              Group= case_when(Group == "in" ~ "Inside the filament",
                                                               Group == "out" ~ "Outside the filament")) %>% 
                  mutate(layer = factor(layer,level = c("Upper 50 m","Below 50m")),
                         Group = factor(Group,level = c("Inside the filament","Outside the filament"))) %>% 
                    select(Group, layer, StationName, Depth, Type, Community, Observed, Shannon, Chao1) %>% 
                      melt(measure.vars = c("Observed","Chao1","Shannon")) %>% 
                        mutate(variable = factor(variable, levels =c("Observed","Chao1","Shannon")))


rich.plot.type<- ggplot(sampledata_DF_long, aes (x = layer, y = value, group = interaction(layer,Group,Community), shape = Community, colour = Group))+
  geom_boxplot()+
  geom_jitter()+
  facet_wrap(variable~., scales = "free", ncol = 3)+
  scale_color_manual(values=c("Inside the filament"="red1","Outside the filament"="blue1"))+
  labs(x = "Water depth",y = "Index values")+
  #geom_signif(comparisons = list(c("in", "out")),
  #            map_signif_level=TRUE, test = "wilcox.test", color = "black")+
  theme_classic()+
  theme(legend.position = "bottom")

ggsave("./figures/alpha_div.pdf", 
       plot = rich.plot.type,
       units = "cm",
       width = 30, height = 15, 
       #scale = 1,
       dpi = 300)

#####################################
#overlaps of communities on ASV level
####################################
#subset upper water column
meso_ps0_up <- subset_samples(meso_ps0, layers == "up")
meso_ps0_up_merged <- merge_samples(meso_ps0_up, "Group")
meso_ps0_up_merged <- prune_taxa(taxa_sums(meso_ps0_up_merged)>0,meso_ps0_up_merged)
sample_names(meso_ps0_up_merged)<- c("up_in","up_out")

meso_ps0_down <- subset_samples(meso_ps0, layers == "down")
meso_ps0_down_merged <- merge_samples(meso_ps0_down, "Group")
meso_ps0_down_merged <- prune_taxa(taxa_sums(meso_ps0_down_merged)>0,meso_ps0_down_merged)
sample_names(meso_ps0_down_merged)<- c("down_in","down_out")

#merge phyloseq objects togather
meso_ps0_merged <- merge_phyloseq(meso_ps0_up_merged,meso_ps0_down_merged)

#down
meso_ps0_up_in <- prune_samples("up_in", meso_ps0_merged)
meso_ps0_up_in <- prune_taxa(taxa_sums(meso_ps0_up_in)>0,meso_ps0_up_in)

meso_ps0_up_out <- prune_samples("up_out", meso_ps0_merged)
meso_ps0_up_out <- prune_taxa(taxa_sums(meso_ps0_up_out)>0,meso_ps0_up_out)

meso_ps0_down_in <- prune_samples("down_in", meso_ps0_merged)
meso_ps0_down_in <- prune_taxa(taxa_sums(meso_ps0_down_in)>0,meso_ps0_down_in)

meso_ps0_down_out <- prune_samples("down_out", meso_ps0_merged)
meso_ps0_down_out <- prune_taxa(taxa_sums(meso_ps0_down_out)>0,meso_ps0_down_out)

z <- list()
z[["up_in"]] <- as.character(colnames(otu_table(meso_ps0_up_in, taxa_are_rows = FALSE)))
z[["up_out"]] <- as.character(colnames(otu_table(meso_ps0_up_out, taxa_are_rows = FALSE)))
z[["down_in"]] <- as.character(colnames(otu_table(meso_ps0_down_in, taxa_are_rows = FALSE)))
z[["down_out"]] <- as.character(colnames(otu_table(meso_ps0_down_out, taxa_are_rows = FALSE)))

#plot
venn(z, snames = names(z), ilab=TRUE, zcolor = "bw")

#generate overlap matrix and explore the results
ASVs_overlaps <- pres_abs_matrix(z)    
ASVs_overlaps$ASV <- rownames(ASVs_overlaps)

taxonomy <- as.data.frame(tax_table(meso_ps0))
taxonomy$ASV <- rownames(taxonomy)

ASVs_overlaps <- full_join(taxonomy,ASVs_overlaps, by = c("ASV"))

#explore results
ASVs_up_all <- ASVs_overlaps %>% 
  filter(up_in =="1",
         up_out =="1")%>%
  group_by(Phylum,Class,Order)%>%
  summarize(Total_ASVs=n())

ASVs_up_in <- ASVs_overlaps %>% 
  filter(up_in =="1",
         up_out =="0")%>%
  group_by(Phylum,Class,Order)%>%
  summarize(Total_ASVs=n())

ASVs_up_out <- ASVs_overlaps %>% 
  filter(up_in =="0",
         up_out =="1")%>%
  group_by(Phylum,Class,Order)%>%
  summarize(Total_ASVs=n())

#below 50 m
ASVs_down_all <- ASVs_overlaps %>% 
  filter(down_in =="1",
         down_out =="1")%>%
  group_by(Phylum,Class,Order)%>%
  summarize(Total_ASVs=n())

ASVs_down_in <- ASVs_overlaps %>% 
  filter(down_in =="1",
         down_out =="0")%>%
  group_by(Phylum,Class,Order)%>%
  summarize(Total_ASVs=n())


ASVs_down_out <- ASVs_overlaps %>% 
  filter(down_in =="0",
         down_out =="1")%>%
  group_by(Phylum,Class,Order)%>%
  summarize(Total_ASVs=n())


#####################################
#get session info and remove all objects and libraries
#####################################
sessionInfo()

rm(list = ls(all = TRUE))
