#set working directory in Windows
# wd <- dirname(rstudioapi::getActiveDocumentContext()$path)
# setwd(wd)

#load libraries
library("ggplot2"); packageVersion("ggplot2")
library("phyloseq"); packageVersion("phyloseq")
library("plyr"); packageVersion("plyr")
library("dplyr"); packageVersion("dplyr")
library("olsrr"); packageVersion("olsrr")
library("cowplot"); packageVersion("cowplot")
library("iNEXT"); packageVersion("iNEXT")

#set plots theme
theme_set(theme_classic())

#load colour palettes
source('./Scripts/color_palettes.R')

#####################################
#Load phyloseq object
####################################
PS107_merged <-  readRDS("./Data/PS107_merged_prev.rds")

#####################################
#Plot rarefaction
####################################
iNEXT.out <- iNEXT(as.data.frame(otu_table(PS107_merged.prev)), q=0, datatype="abundance")
rare <-fortify(iNEXT.out, type=1)

meta <- as(sample_data(PS107_merged.prev), "data.frame")
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
  #xlim(0,1e5)+ylim(0,6000)+
  theme_classic(base_size = 12)+theme(legend.position="bottom")

ggsave("./figures/rarefaction-prev.pdf", 
       plot = rare.p,
       units = "cm",
       width = 30, height = 30, 
       #scale = 1,
       dpi = 300)

#####################################
#Alpha diversity statistical tests
####################################
# Calculate richness
PS107_alpha.div <- estimate_richness(PS107_merged.prev, split = TRUE, measures = NULL)

#generate data set with all bacterial community characteristics
PS107_comm.char<- data.frame(Station = sample_data(PS107_merged.prev)$StationName,
                             Group = sample_data(PS107_merged.prev)$Group,
                             Depth = sample_data(PS107_merged.prev)$Depth,
                             Type = sample_data(PS107_merged.prev)$Type,
                             Community = sample_data(PS107_merged.prev)$Community,
                             Sequences= sample_sums(PS107_merged),
                             Observed = PS107_alpha.div$Observed,
                             Chao1 = PS107_alpha.div$Chao1,
                             Completness = round(100*PS107_alpha.div$Observed/PS107_alpha.div$Chao1, digits=2),
                             Shanonn = round(PS107_alpha.div$Shannon,digits=2),
                             Simpson = round(PS107_alpha.div$Simpson,digits=2),
                             Evenness = round(PS107_alpha.div$Shannon/log(PS107_alpha.div$Observed),digits=2))


write.csv(PS107_comm.char, "./Data/alpha_table_prev.csv")

#test siginifance of difference in alpha diversity in upper water column
PS107_comm.char.up <- PS107_comm.char[PS107_comm.char$Depth < 100, ]

wilcox.test(PS107_comm.char.up$Chao1[PS107_comm.char.up$Group =="in"],
            PS107_comm.char.up$Chao1[PS107_comm.char.up$Group =="out"])


Chao1_upper_50M.plot <- ggplot(PS107_comm.char.up, aes (x = Group, y = Chao1, fill = Group))+
  geom_boxplot(outlier.color = NULL, notch = FALSE)+
  scale_fill_manual(values = c("in"= "red","out" = "blue"))+
  geom_boxplot()+
  facet_grid(~Community)+
  theme_classic(base_size = 12)+
  theme(legend.position = "bottom")


#test siginifance of difference in alpha diversity below 100 m
PS107_comm.char.down <- PS107_comm.char[PS107_comm.char$Depth >= 100, ]

wilcox.test(PS107_comm.char.down$Chao1[PS107_comm.char.down$Group =="in"],
            PS107_comm.char.down$Chao1[PS107_comm.char.down$Group =="out"])

#Chao1 summary
PS107_comm.Chao1.agg <- do.call(data.frame, aggregate(Chao1~ Group+Type + Community, PS107_comm.char, function(x) c(mean = mean(x), se = se(x),median = median(x))))

  
#Shanonn summary
PS107_comm.Shanonn.agg <- do.call(data.frame, aggregate(Shanonn~ Group+Type + Community, PS107_comm.char, function(x) c(mean = mean(x), se = se(x),median = median(x))))

#Simpson summary
PS107_comm.Simpson.agg <- do.call(data.frame, aggregate(Simpson~ Group+Type + Community, PS107_comm.char, function(x) c(mean = mean(x), se = se(x),median = median(x))))


#####################################
#ASV distribution
#####################################
for (frac in c("FL","PA")){
  #mean number of ASV per sample
  PS107_merged_sub <- subset_samples(PS107_merged.prev, Community== frac)
  PS107_merged_sub <- prune_taxa(taxa_sums(PS107_merged_sub)>0, PS107_merged_sub)
  
  vectorx <- vector()
  for (i in sample_names(PS107_merged_sub)){
    PS107_merged_station <- prune_samples(i, PS107_merged_sub)
    PS107_merged_station <- prune_taxa(taxa_sums(PS107_merged_station)>0,PS107_merged_station)
    vectorx <- c(vectorx,dim(otu_table(PS107_merged_station))[1])
  }
  PS107_bacteria <- subset_taxa(PS107_merged_sub, Kingdom == "Bacteria")
  PS107_bacteria <- prune_taxa(taxa_sums(PS107_bacteria)>0,PS107_bacteria)
  
  PS107_archaea <- subset_taxa(PS107_merged_sub, Kingdom == "Archaea")
  PS107_archaea <- prune_taxa(taxa_sums(PS107_archaea)>0,PS107_archaea)
  
  assign(paste("df", frac, sep = "."), t(data.frame(total.reads= sum(sample_sums(PS107_merged_sub)),
                                                    total.ASV= dim(otu_table(PS107_merged_sub))[1],
                                                    total.BAC.ASV=dim(otu_table(PS107_bacteria))[1],
                                                    total.ARC.ASV=dim(otu_table(PS107_archaea))[1],
                                                    mean.ASV.per.sample=mean(vectorx),
                                                    sd.ASV.per.sample=se(vectorx))))
}

df.overview <- cbind(df.FL,df.PA)



#####################################
#get session info and remove all objects and libraries
#####################################
sessionInfo()

rm(list = ls(all = TRUE))
