This GitHub Repository contains code and data included in:
 
## [Submesoscale physicochemical dynamics directly shape bacterioplankton community structure in space and time](https://doi.org/10.1002/lno.11799)
Eduard Fadeev, Matthias Wietz, Wilken-Jon von Appen, Morten H. Iversen, Eva-Maria Nöthig, Anja Engel, Martin Graeve, Antje Boetius

### Abstract:
_Submesoscale eddies and fronts are important components of oceanic mixing and energy fluxes. These phenomena occur in the surface ocean for a period of several days, on scales between a few hundred meters and few tens of kilometers. Remote sensing and modeling suggest that eddies and fronts may influence marine ecosystem dynamics, but their limited temporal and spatial scales make them challenging for observation and in situ sampling. Here, the study of a submesoscale filament in summerly Arctic waters (depth 0–400 m) revealed enhanced mixing of Polar and Atlantic water masses, resulting in a ca. 4 km wide and ca. 50 km long filament with distinct physical and biogeochemical characteristics. Compared to the surrounding waters, the filament was characterized by a distinct phytoplankton bloom, associated with depleted inorganic nutrients, elevated chlorophyll a concentrations, as well as twofold higher phyto- and bacterioplankton cell abundances. High-throughput 16S rRNA gene sequencing of bacterioplankton communities revealed enrichment of typical phytoplankton bloom-associated taxonomic groups (e.g., Flavobacteriales) inside the filament. Furthermore, linked to the strong water subduction, the vertical export of organic matter to 400 m depth inside the filament was twofold higher compared to the surrounding waters. Altogether, our results show that physical submesoscale mixing can shape distinct biogeochemical conditions and microbial communities within a few kilometers of the ocean. Hence, the role of submesoscale features in polar waters for surface ocean biodiversity and biogeochemical processes need further investigation, especially with regard to the fate of sea ice in the warming Arctic Ocean._


### Content:
**data** - output files from the 16S dada2 workflow and metadata. \
**scripts** - supporting scripts for the statistical analysis. \
```MESO_dataset_preprocess.R``` - processing of the dada2 output into phyloseq. \
```MESO_alpha_diversity.R``` - Mic. communities alpha diversity calculations (included in supplementary material). \
```MESO_beta_diversity.R``` - Mic. communities beta diversity (Figure 3). \
```MESO_asv_overalps.R``` - Shared and unique ASVs between communities (Figure 4). \
```MESO_enrichment.R``` - taxonomic enrochment tests between mic. communities (Figure 5). \

### Author:
Eduard Fadeev([dr.eduard.fadeev@gmail.com](mailto:dr.eduard.fadeev@gmail.com)) 

### Software Versions:
R version 4.1.1 (2021-08-10)\
RStudio version: 1.4.1717
