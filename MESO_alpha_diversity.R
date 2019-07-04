#set working directory
wd <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd)

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
PS107_merged <-  readRDS("./Data/PS107_merged.rds")

#####################################
#Plot rarefaction
####################################
#PS107_merged <- subset_samples(PS107_merged,Type!="Sediment")
iNEXT.out <- iNEXT(as.data.frame(otu_table(PS107_merged)), q=0, datatype="abundance")
rare <-fortify(iNEXT.out, type=1)

meta <- as(sample_data(PS107_merged), "data.frame")
meta$site <- rownames(meta)
rare$Fraction <- meta$Fraction[match(rare$site, meta$site)] 
rare$StationName <- meta$StationName[match(rare$site, meta$site)] 
rare$Type <- meta$Type[match(rare$site, meta$site)]
rare$label <- paste(rare$StationName,rare$Type, sep = "-")

rare.point <- rare[which(rare$method == "observed"),]
rare.line <- rare[which(rare$method != "observed"),]
rare.line$method <- factor (rare.line$method,
                            c("interpolated", "extrapolated"),
                            c("interpolation", "extrapolation"))

#exclude sediment samples from the plot
#rare.line <- rare.line[rare.line$Type!="Sediment",]
#rare.point <- rare.point[rare.point$Type!="Sediment",]

rare.p <- ggplot(rare, aes(x=x, y=y, colour = site))+
  geom_line(aes(linetype = method), lwd = 0.5, data= rare.line)+
  #geom_ribbon(aes(ymin=y.lwr, ymax= y.upr, colour = NULL), alpha = 0.2)+
  geom_point(aes(shape=Type, fill=Fraction), size =3, data= rare.point)+
  #geom_text(aes(label=label), size =2, data= rare.point, colour = "black", nudge_y = -100)+
  scale_colour_discrete(guide = FALSE)+
  labs(x = "Sample size", y = "Species richness")+
  facet_grid(~Fraction)+
  #xlim(0,1e5)+ylim(0,6000)+
  theme_classic(base_size = 12)+theme(legend.position="bottom")

ggsave("./figures/Figure-rarefaction-sed.pdf", 
       plot = rare.p,
       units = "cm",
       #width = 17.8,
       scale = 1,
       #height = 17.4,
       dpi = 300)

#####################################
#Alpha diversity
####################################
#calculate alpha diversity indeces 100 times
# Initialize matrices to store richness and evenness estimates
#rare_size <- 5000
rare_size <- min(sample_sums(PS107_merged))
nsamp <-  nsamples(PS107_merged)
trials <- 100

#Observed number of OTUs
OTUs <- matrix(nrow = nsamp, ncol = trials)
row.names(OTUs) <- sample_names(PS107_merged)

richness <- matrix(nrow = nsamp, ncol = trials)
row.names(richness) <- sample_names(PS107_merged)

evenness <- matrix(nrow = nsamp, ncol = trials)
row.names(evenness) <- sample_names(PS107_merged)

#Shannon div. index
Shannon <- matrix(nrow = nsamp, ncol = trials)
row.names(Shannon) <- sample_names(PS107_merged)

#Simpson div.index
Simpson <- matrix(nrow = nsamp, ncol = trials)
row.names(Simpson) <- sample_names(PS107_merged)

#Rarefy the dataset by the smallest sample 100 times
set.seed(12345)

for (i in 1:100) {
PS107_merged.rare <- rarefy_even_depth(PS107_merged, sample.size = rare_size,
                                      rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)

# Calculate richness
PS107_alpha.div <- estimate_richness(PS107_merged.rare, split = TRUE, measures = NULL)

# OTU num.
OTUs[ ,i] <- PS107_alpha.div$Observed

#Chao1
richness[ ,i] <- PS107_alpha.div$Chao1

# Calculate Shannon
Shannon[ ,i] <- PS107_alpha.div$Shannon

# Calculate Simpson
Simpson[ ,i] <- PS107_alpha.div$Simpson

# Calculate evenness
evenness[ ,i] <- PS107_alpha.div$Shannon/log(PS107_alpha.div$Observed)
}

# Create a new dataframe to hold the means and standard deviations of observed OTUs
SampleID <- row.names(OTUs)
mean <- apply(OTUs, 1, mean)
sd <- apply(OTUs, 1, sd)
measure <- rep("OTUs", nsamp)
OTUs_stats <- data.frame(SampleID, mean, sd, measure)

# Create a new dataframe to hold the means and standard deviations of richness estimates
SampleID <- row.names(richness)
mean <- apply(richness, 1, mean)
sd <- apply(richness, 1, sd)
measure <- rep("Richness", nsamp)
rich_stats <- data.frame(SampleID, mean, sd, measure)

# Create a new dataframe to hold the means and standard deviations of evenness estimates
SampleID <- row.names(evenness)
mean <- apply(evenness, 1, mean)
sd <- apply(evenness, 1, sd)
measure <- rep("evenness J", nsamp)
even_stats <- data.frame(SampleID, mean, sd, measure)

# Create a new dataframe to hold the means and standard deviations of Shannon div. index estimates
SampleID <- row.names(Shannon)
mean <- apply(Shannon, 1, mean)
sd <- apply(Shannon, 1, sd)
measure <- rep("Shannon", nsamp)
Shannon_stats <- data.frame(SampleID, mean, sd, measure)

# Create a new dataframe to hold the means and standard deviations of Simpson div. index estimates
SampleID <- row.names(Simpson)
mean <- apply(Simpson, 1, mean)
sd <- apply(Simpson, 1, sd)
measure <- rep("Simpson", nsamp)
Simpson_stats <- data.frame(SampleID, mean, sd, measure)

#generate data set with all bacterial community characteristics
PS107_comm.char<- data.frame(Station = sample_data(PS107_merged)$StationName,
                           Depth = sample_data(PS107_merged)$Depth,
                           #Type = sample_data(PS107_merged)$Type,
                           Community = sample_data(PS107_merged)$Fraction,
                           #Region = sample_data(PS107_merged)$Region,
                           Sequences= sample_sums(PS107_merged),
                           Observed = paste(round(OTUs_stats$mean,digits=0),round(OTUs_stats$sd,digits=0),sep="\u00B1"),
                           Chao1 = paste(round(rich_stats$mean,digits=0),round(rich_stats$sd,digits=0),sep="\u00B1"),
                           Completness = round(100*OTUs_stats$mean/rich_stats$mean, digits=2),
                           Shanonn = paste(round(rich_stats$mean,digits=2),round(rich_stats$sd,digits=2),sep="\u00B1"),
                           Simpson = paste(round(Simpson_stats$mean,digits=2),round(Simpson_stats$sd,digits=2),sep="\u00B1"),
                           Evenness = paste(round(even_stats$mean,digits=2),round(even_stats$sd,digits=2),sep="\u00B1"))
                         
write.csv(PS107_comm.char, "./Data/alpha_table.csv")

#####################################
#Alpha diversity plots
####################################
#generate the same table with no SD
PS107_comm.char<- data.frame(Station = sample_data(PS107_merged)$StationName,
                            Longitude = sample_data(PS107_merged)$Longitude..degrees_east.,
                            Type = sample_data(PS107_merged)$Type,
                            Fraction = sample_data(PS107_merged)$Fraction,
                            Region = sample_data(PS107_merged)$Region,
                            Sequences= sample_sums(PS107_merged),
                            Observed = round(OTUs_stats$mean,digits=0),
                            Chao1 = round(rich_stats$mean,digits=1),
                            Completness = OTUs_stats$mean/rich_stats$mean,
                            Shanonn = round(rich_stats$mean,digits=2),
                            Simpson = round(Simpson_stats$mean,digits=2),
                            Evenness = round(even_stats$mean,digits=2))



#exclude sediment samples from the plots
PS107_comm.char.no.sed <- PS107_comm.char[PS107_comm.char$Type!="Sediment",]

#Chao1 plot
PS107_comm.Chao1.agg <- do.call(data.frame, aggregate(Chao1~ Type +Fraction, PS107_comm.char.no.sed, function(x) c(mean = mean(x), sd = sd(x),median = median(x))))

PS107_comm.Chao1.agg$Type<- factor(PS107_comm.Chao1.agg$Type, 
                                  levels = rev(c("SRF", "EPI", "MESO", "BATHY")))

PS107_comm.Chao1.agg$Type <- as.character(PS107_comm.Chao1.agg$Type)
#plot
PS107_Chao1 <-   ggplot(PS107_comm.char.no.sed, aes(y = Chao1, x = Type))+
  geom_point(data = PS107_comm.char.no.sed, aes(y = Chao1, x = Type, group = Fraction, shape = Fraction),colour="black",size = 5) +
  geom_point(data = PS107_comm.char.no.sed, aes(y = Chao1, x = Type, group = Fraction, colour = Region, shape = Fraction),  size =  4) +
  geom_point(data = PS107_comm.Chao1.agg, aes(y = Chao1.mean, x =Type, group = Fraction, shape = Fraction), size = 7) + 
  geom_line(data = PS107_comm.Chao1.agg, aes(y = Chao1.mean, x = Type, group = Fraction, linetype= Fraction), size = 1)+
  #geom_text(data = PS107_comm.char, aes(y = Chao1, x = Type, label = Station))+
  scale_x_discrete(limits = c("BATHY","MESO","EPI","SRF"))+
  coord_flip()+
  scale_color_manual(values = c("EGC" = "blue", "WSC"="red")) +
  xlab("Depth")+
  ylab("Chao1 Species Estimator")+
  theme_classic()+
  theme(legend.position = "none")


#Pielou's evenness (J)
PS107_comm.evenness.agg <- do.call(data.frame, aggregate(Evenness~ Type +Fraction, PS107_comm.char.no.sed, function(x) c(mean = mean(x), sd = sd(x), median = median(x))))

PS107_comm.evenness.agg$Type<- factor(
  PS107_comm.evenness.agg$Type, 
  levels = rev(c("SRF", "EPI", "MESO", "BATHY")))

#plot
PS107_evenness <-   ggplot()+
  geom_point(data = PS107_comm.char.no.sed, aes(y = Evenness, x = Type, group = Fraction, shape = Fraction),colour="black",size = 5)+
  geom_point(data = PS107_comm.char.no.sed, aes(y = Evenness, x = Type, group = Fraction, colour = Region, shape = Fraction), size = 4) +
  geom_point(data = PS107_comm.evenness.agg, aes(y = Evenness.mean, x = Type, group = Fraction, shape = Fraction), size = 7) +  
  geom_line(data = PS107_comm.evenness.agg, aes(y = Evenness.mean, x = Type, group = Fraction, linetype= Fraction), size = 1)+
  #geom_text(data = PS107_comm.char, aes(y = Evenness, x = Type, label = Station))+
  scale_x_discrete(limits = c("BATHY","MESO","EPI","SRF"))+
  coord_flip()+
  xlab("Depth")+
  scale_color_manual(values = c("EGC" = "blue", "WSC"="red")) +
  ylab("Pielou's evenness index (J)")+
  theme_classic()+
  theme(legend.position = "none")

#combined plot
plot_grid(PS107_Chao1, PS107_evenness, labels = c("A", "B"), ncol = 2, align = "hv")

ggsave("./figures/Figure-alpha.pdf", 
       plot = ggplot2::last_plot(),
       scale = 1,
       units = "cm",
       width = 17.8,
       #height = 17.4,
       dpi = 300)

#####################################
#Statistical tests
#####################################
PS107_comm.char_FL <- PS107_comm.char.no.sed[PS107_comm.char.no.sed$Fraction=="FL",]
PS107_comm.char_PA <- PS107_comm.char.no.sed[PS107_comm.char.no.sed$Fraction=="PA",]

#compare richness throughout the water column
#FL
wilcox.test(PS107_comm.char_FL[PS107_comm.char_FL$Type=="SRF",]$Chao1,
            PS107_comm.char_FL[PS107_comm.char_FL$Type=="EPI",]$Chao1)
wilcox.test(PS107_comm.char_FL[PS107_comm.char_FL$Type=="EPI",]$Chao1,
            PS107_comm.char_FL[PS107_comm.char_FL$Type=="MESO",]$Chao1)
wilcox.test(PS107_comm.char_FL[PS107_comm.char_FL$Type=="MESO",]$Chao1,
            PS107_comm.char_FL[PS107_comm.char_FL$Type=="BATHY",]$Chao1)

#PA
wilcox.test(PS107_comm.char_PA[PS107_comm.char_PA$Type=="SRF",]$Chao1,
            PS107_comm.char_PA[PS107_comm.char_PA$Type=="EPI",]$Chao1)
wilcox.test(PS107_comm.char_PA[PS107_comm.char_PA$Type=="EPI",]$Chao1,
            PS107_comm.char_PA[PS107_comm.char_PA$Type=="MESO",]$Chao1)
wilcox.test(PS107_comm.char_PA[PS107_comm.char_PA$Type=="MESO",]$Chao1,
            PS107_comm.char_PA[PS107_comm.char_PA$Type=="BATHY",]$Chao1)

#compare richness between fractions
wilcox.test(PS107_comm.char_FL[PS107_comm.char_FL$Type=="SRF",]$Chao1,
            PS107_comm.char_PA[PS107_comm.char_PA$Type=="SRF",]$Chao1)

wilcox.test(PS107_comm.char_FL[PS107_comm.char_FL$Type=="EPI",]$Chao1,
            PS107_comm.char_PA[PS107_comm.char_PA$Type=="EPI",]$Chao1)

wilcox.test(PS107_comm.char_FL[PS107_comm.char_FL$Type=="MESO",]$Chao1,
            PS107_comm.char_PA[PS107_comm.char_PA$Type=="MESO",]$Chao1)

wilcox.test(PS107_comm.char_FL[PS107_comm.char_FL$Type=="BATHY",]$Chao1,
            PS107_comm.char_PA[PS107_comm.char_PA$Type=="BATHY",]$Chao1)


#compare eveness throughout the water column
#FL
wilcox.test(PS107_comm.char_FL[PS107_comm.char_FL$Type=="SRF",]$Evenness,
            PS107_comm.char_FL[PS107_comm.char_FL$Type=="EPI",]$Evenness)
wilcox.test(PS107_comm.char_FL[PS107_comm.char_FL$Type=="EPI",]$Evenness,
            PS107_comm.char_FL[PS107_comm.char_FL$Type=="MESO",]$Evenness)
wilcox.test(PS107_comm.char_FL[PS107_comm.char_FL$Type=="MESO",]$Evenness,
            PS107_comm.char_FL[PS107_comm.char_FL$Type=="BATHY",]$Evenness)

#PA
wilcox.test(PS107_comm.char_PA[PS107_comm.char_PA$Type=="SRF",]$Evenness,
            PS107_comm.char_PA[PS107_comm.char_PA$Type=="EPI",]$Evenness)
wilcox.test(PS107_comm.char_PA[PS107_comm.char_PA$Type=="EPI",]$Evenness,
            PS107_comm.char_PA[PS107_comm.char_PA$Type=="MESO",]$Evenness)
wilcox.test(PS107_comm.char_PA[PS107_comm.char_PA$Type=="MESO",]$Evenness,
            PS107_comm.char_PA[PS107_comm.char_PA$Type=="BATHY",]$Evenness)

#####################################
#OTU distribution
#####################################
PS107_merged.no.sed <- subset_samples(PS107_merged, Type!="Sediment")
PS107_merged.no.sed <- prune_taxa(taxa_sums(PS107_merged.no.sed)>0,PS107_merged.no.sed)

for (frac in c("FL","PA")){
  #mean number of OTU per sample
  PS107_merged_sub <- subset_samples(PS107_merged.no.sed, Fraction== frac)
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
                                                    total.OTU= dim(otu_table(PS107_merged_sub))[1],
                                                    total.BAC.OTU=dim(otu_table(PS107_bacteria))[1],
                                                    total.ARC.OTU=dim(otu_table(PS107_archaea))[1],
                                                    mean.OTU.per.sample=mean(vectorx),
                                                    sd.OTU.per.sample=sd(vectorx)/sqrt(length(vectorx)))))
}

df.overview <- cbind(df.FL,df.PA)

#####################################
#get session info and remove all objects and libraries
#####################################
sessionInfo()

rm(list = ls(all = TRUE))
