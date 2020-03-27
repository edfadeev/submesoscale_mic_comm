# #set working directory
# wd <- dirname(rstudioapi::getActiveDocumentContext()$path)
# setwd(wd)

#load libraries
library("ggplot2"); packageVersion("ggplot2")
library("DESeq2"); packageVersion("DESeq2")
library("phyloseq"); packageVersion("phyloseq")
library("cowplot"); packageVersion("cowplot")
library("dplyr"); packageVersion("dplyr")
library("gage"); packageVersion("gage")

# calculate geometric means 
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

#set plots theme
theme_set(theme_classic())

#load colour palettes
source('./Scripts/color_palettes.R')

#####################################
#Load phyloseq object
####################################
PS107_merged.prev <-  readRDS("./Data/PS107_merged_prev.rds")


#create taxonomy db
BAC_tax <- as.data.frame(tax_table(PS107_merged.prev))
BAC_tax$OTU <- rownames(tax_table(PS107_merged.prev))

BAC_tax %>%
  mutate_if(is.factor, as.character) -> BAC_tax

#Extract the desired taxonomic level
Genera <- unique(BAC_tax$Order)
#generate list of OTU for each Taxa
OTU.gs <- list()
for (s in 1:length(Genera)){
  n <- Genera[s]
  Class_OTU <- subset(BAC_tax, Order == n)
  OTU.gs[[n]] <- Class_OTU$OTU
  
}

#################
frac <- c("FL","PA")
type <- c("Surface-10","Chl.max-20-30","B.Chl.max-50","Epipelagic-100","Mesopelagic-200","Mesopelagic-400")
layers <- c("up","down")

res_all <- data.frame()

for (i in 1:2){
PS107_merged.prev.comm <- subset_samples(PS107_merged.prev, Community ==frac[i])
for (m in 1:2){
PS107_merged.prev.frac <- subset_samples(PS107_merged.prev.comm, layers ==layers[m])
PS107_merged.prev.frac <- prune_taxa(taxa_sums(PS107_merged.prev.frac)>0,PS107_merged.prev.frac)

#run DEseq
BAC_sed.ddsMat <- phyloseq_to_deseq2(PS107_merged.prev.frac, ~Group)
varianceStabilizingTransformation(BAC_sed.ddsMat, blind = TRUE, fitType = "parametric")
BAC_sed.ddsMat <- estimateSizeFactors(BAC_sed.ddsMat)
BAC_sed.ddsMat <- estimateDispersions(BAC_sed.ddsMat)
BAC_sed.DEseq = DESeq(BAC_sed.ddsMat)
BAC_sed.DEseq.res <- results(BAC_sed.DEseq)

deseq2.fc <- BAC_sed.DEseq.res$log2FoldChange
names(deseq2.fc) <- rownames(BAC_sed.DEseq.res)
exp.fc=deseq2.fc

fc.Class.p <- gage(exp.fc, gsets = OTU.gs, same.dir=TRUE, ref = NULL, samp = NULL)

#plot
FL_enrch <- rbind(fc.Class.p$greater,fc.Class.p$less)
FL_enrch <-  data.frame(FL_enrch[FL_enrch[,"q.val"]<0.01 &
                                   !is.na(FL_enrch[,"q.val"]),])
FL_enrch$Order <- rownames(FL_enrch)
FL_enrch <- unique(merge(FL_enrch,BAC_tax[,c("Class","Order")]))


FL_enrch$frac <- frac[i]
FL_enrch$layers <- layers[m]

res_all <- rbind(res_all,FL_enrch, make.row.names = TRUE)
}
}

#plot
enr.FL.p <- ggplot(data=res_all[res_all$frac== "FL",], aes(y=stat.mean , x=Order, label = set.size, shape = frac))+ 
  geom_text(size = 4, aes(y=stat.mean , x=Order), nudge_y= 1.5, nudge_x= 0)+
  ylab("Mean log2foldchange")+ 
  geom_point(size = 5, aes(colour = Class))+
  scale_colour_manual(values = phyla.col)+ 
  ylim(-10,10)+
  scale_x_discrete("Class",expand = waiver())+
  geom_hline(aes(yintercept=0), linetype="dashed")+
  theme_classic()+
  theme(legend.position = "bottom", axis.text.x = element_text(angle =90))+
  coord_flip()+
  facet_grid(.~layers)

enr.PA.p <- ggplot(data=res_all[res_all$frac== "PA",], aes(y=stat.mean , x=Order, label = set.size, shape = frac))+ 
  geom_text(size = 4, aes(y=stat.mean , x=Order), nudge_y= 1.5, nudge_x= 0)+
  ylab("Mean log2foldchange")+ 
  geom_point(size = 5, aes(colour = Class))+
  scale_colour_manual(values = phyla.col)+ 
  ylim(-10,10)+
  scale_x_discrete("Class",expand = waiver())+
  geom_hline(aes(yintercept=0), linetype="dashed")+
  theme_classic()+
  theme(legend.position = "bottom", axis.text.x = element_text(angle =90))+
  coord_flip()+
  facet_grid(.~layers)

plot_grid(enr.FL.p,enr.PA.p, ncol = 1)


res_all$layers <- factor(res_all$layers, levels = c("up", "down"))
ggplot(data=res_all, aes(y=stat.mean , x=Order, label = set.size, shape= frac))+ 
  geom_text(size = 4, aes(y=stat.mean , x=Order), nudge_y= 1.5, nudge_x= 0)+
  geom_point(data=res_all[res_all$frac== "FL",],size = 7, shape = 16, colour = "black")+
  geom_point(data=res_all[res_all$frac== "PA",],size = 7, shape = 17, colour = "black")+
  geom_point(data=res_all[res_all$frac== "FL",],size = 5, shape = 16, aes(colour = Class))+
  geom_point(data=res_all[res_all$frac== "PA",],size = 5, shape = 17, aes(colour = Class))+
  ylab("Mean log2foldchange")+ 
  scale_colour_manual(values = phyla.col)+ 
  ylim(-10,10)+
  scale_x_discrete("Class",expand = waiver())+
  geom_hline(aes(yintercept=0), linetype="dashed")+
  theme_classic()+
  theme(legend.position = "bottom", axis.text.x = element_text(angle =90))+
  coord_flip()+
  facet_grid(.~layers)




df_shapes <- data.frame(shape = 0:24)
ggplot(df_shapes, aes(0, 0, shape = shape)) +
  geom_point(aes(shape = shape), size = 5, fill = 'red') +
  scale_shape_identity() +
  facet_wrap(~shape) +
  theme_void()
