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

phyla.col <- c("Alphaproteobacteria"= "#771155" ,
               "Acidimicrobiia"="#AA4488",
               "Arctic97B-4 marine group" = "#DDAA77",
               "Betaproteobacteria"="#CC99BB", 
               "BD2-11_terrestrial_group"= "#114477",
               "Cyanobacteria" = "#771122",
               "Clostridia" = "#AA4455",
               "Deltaproteobacteria"= "#4477AA", 
               "Bacteroidia"= "#77AADD", 
               "Gammaproteobacteria"="#117777", 
               "Epsilonproteobacteria" = "#FFFFFF",
               "Lentisphaeria" = "#DD7788",
               "Nitrososphaeria" = "#E69F00",
               "Marinimicrobia_(SAR406_clade)_uncl"= "#34ABAA",
               "Nitrospinia" = "#77CCCC",  
               "OM190"="#117744", 
               "Kiritimatiellae"= "#44AA77", 
               "Oligosphaeria" = "#F0E442",
               "Phycisphaerae"= "#88CCAA", 
               "Planctomycetes"= "#777711", 
               "SAR202 clade"= "#AAAA44", 
               "Dehalococcoidia"= "#AAAA44",
               "Sphingobacteriia" ="#DDDD77", 
               "SPOTSOCT00m83" = "#774411", 
               "Thermoplasmata" = "#0072B2", 
               "Verrucomicrobiae" = "#AA7744",
               "BD7-11" = "#CC79A7",
                "028H05-P-BN-P5"= "white",
                "WCHB1-41"= "red")




