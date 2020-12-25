#load library
library(RColorBrewer)


#####################################
#Color palettes for plots
#####################################
#large colours range
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

# The palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

tol21rainbow<- c("#771155", 
                 "#AA4488", 
                 "#CC99BB", 
                 "#114477", 
                 "#4477AA", 
                 "#77AADD", 
                 "#117777", 
                 "#44AAAA", 
                 "#77CCCC", 
                 "#117744", 
                 "#44AA77", 
                 "#88CCAA", 
                 "#777711", 
                 "#AAAA44", 
                 "#DDDD77", 
                 "#774411", 
                 "#AA7744", 
                 "#DDAA77", 
                 "#771122", 
                 "#AA4455", 
                 "#DD7788")


class_col <- c("Acidimicrobiia"="#DDDD77",
               "Alphaproteobacteria"= "#771155" ,
               "Bacteroidia"= "#77AADD", 
               "Gammaproteobacteria"="#117777", 
               "Nitrososphaeria" = "#E69F00",
               "Verrucomicrobiae" = "#AA7744",
               "Marinimicrobia_(SAR406_clade)_uncl"= "#34ABAA", 
               "Marinimicrobia (SAR406 clade)_uncl"= "#34ABAA", 
               "SAR324_clade(Marine_group_B)_uncl"= "#4477AA", 
               "Thermoplasmata" = "#0072B2",
               "Nitrospinia" = "#77CCCC", 
               "Dehalococcoidia"= "#AAAA44",
               "Planctomycetes"= "#777711", 
               "Bacteria_uncl" = "#AA4455",
               "OM190"="#117744",
               "NB1-j_uncl"= "#44AA77",
               "Phycisphaerae"= "#88CCAA", 
               "Cyanobacteriia" = "#771122",
               "Bdellovibrionia"="#AA4488",
               "Bacilli" = "#DDAA77",
               "Other taxa < 2%" = "gray50")
