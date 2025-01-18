## DNA sampling in the Gulf of Mexico on GOMECC-4

This repository has QIIME 2 compatible files, metadata, R Markdown files, and all other input files needed for analysis and figure generation for the following manuscript:

**Assessing the effects of warming and carbonate chemistry on marine microbes in the Gulf of Mexico through basin-scale DNA metabarcoding**<br/>
*Sean R. Anderson, Katherine Silliman, Leticia Barbero, Fabian A. Gomez, Beth A. Stauffer, Astrid Schnetzer, Christopher R. Kelble, and Luke R. Thompson, (2024)*

* [Preprint](https://www.biorxiv.org/content/10.1101/2024.07.30.605667v1)

### 1. Summary
Ocean acidification and warming threaten marine life, either directly or indirectly, which can lead to loss of biodiversity and negative impacts for marine ecosystems. The effects of carbonate chemistry parameters and temperature on diverse marine microbes remains unclear. Yet, it is important to consider microbes because they form the base of ocean food webs, mediate global carbon and nutrient cycles, and are potential indicators of environmental change. In large part, this disconnect is rooted in a lack of spatial sampling of microbes and environmental variables in many ocean basins, including the Gulf of Mexico (GOM). In this study, we collected DNA samples at the basin scale in the GOM as part of the fourth Gulf of Mexico Ecosystems and Carbon Cycle (GOMECC-4) cruise that sailed in the summer-fall of 2021. DNA samples were collected at 51 sites along 16 inshore-offshore transects and up to three depths per site that reflected the surface, deep chlorophyll maximum, and near bottom (481 total filters). DNA metabarcoding captured prokaryotes (16S V4-V5) and protists (18S V9) at previously unresolved spatial scales. Generalized linear models were used to reveal the effects of carbonate chemistry parameters, temperature, and other variables on group-specific relative abundance in the photic zone. Models supported prior physiological trends among certain microbes, like positive temperature effects on SAR11 and SAR86, as well as a negative response of *Prochlorococcus* to lower pH. New insights were observed for Syndiniales and Sagenista, ubiquitous protists that represent parasitic and herbivorous lifestyles. At the species level, picoeukaryotes like *Ostreococcus* sp. and *Emiliania huxleyi* were found to be indicator taxa of less buffered waters in the GOM at this time.

### 2. Bioinformatics
Code for 16S and 18S datasets are available in the `code` folder.
* 16S and 18S amplicons processed separately
* Remove primers with Cutadapt
* Process reads with [Tourmaline](https://github.com/aomlomics/tourmaline)
* Tourmaline uses a Snakemake workflow and wraps QIIME 2 and DADA2 to infer amplicon sequence variants (ASVs)
* Taxonomy assigned for 16S using SILVA (Version 138.1) and 18S with the PR2 database (Version 5.0.1)
* All Tourmaline output files and taxonomic reference databases for this project have been archived on [Zenodo](https://zenodo.org/records/13102580)

### 3. R data analysis and visualization
Code is available as markdown files for [16S](https://aomlomics.github.io/gomecc/GOMECC4_16S.html) and [18S](https://aomlomics.github.io/gomecc/GOMECC4_18S.html) datasets. Raw markdown files are also in the `code` folder.
* QIIME 2 files uploaded to R using [qiime2R](https://github.com/jbisanz/qiime2R)
* Hierarchical clustering of microbial composition to resolve spatial patterns
* Population dynamics observed in each cluster with stacked bar plots, PCoA plots, and diversity plots (richness and Shannon index)
* Generalized linear models (GLMs) constructed for major microbial groups (class level 18S; order level 16S) using relative abundances and non-collinear variables
* Models focused on the photic zone to mitigate collinearity
* Final GLMs applied to all GOMECC-4 sites to expand microbial distributions at the surface; plots made in [Ocean Data View](https://odv.awi.de/)
* ODV plots in the `odv-plots` folder.
* Indicator analysis performed to reveal ASVs indicative of less buffered surface waters (via TA:DIC ratio) using [indicspecies](https://emf-creaf.github.io/indicspecies/)
* All other input files necessary for analyses are in the `data-input` folder.

### 4. Links to associated data
* Raw sequence data for this project are available in NCBI SRA under Bioproject PRJNA887898
* Species observation data (counts) has been published on the Ocean Biodiversity Information System (OBIS) and the Global Biodiversity Information Facility [(GBIF)](https://www.gbif.org/dataset/9012def0-bd87-48a0-ac9e-e0e78dd37689)
* Biological data has also been submitted to the National Centers for Environmental Information [(NCEI)](https://www.ncei.noaa.gov/archive/accession/0250940/data/0-data/noaa-aoml-gomecc)
* Environmental measurements from the Niskin bottles and CTD profiles are also available at NCEI via Accession numbers [0294320](https://doi.org/10.25921/4twf-pp50) and [0284051](https://doi.org/10.25921/04h7-gv36)
* The GOMECC-4 [cruise report](https://doi.org/10.25923/rwx5-s713) is also available and details all sampling and processing on the cruise

**This work was funded in part through the NOAA Ocean Acidification Program (OAP) ROR #02bfn4816 under project numbers 21392 (Thompson) and 20708 (Barbero) and by awards NA16OAR4320199 and NA21OAR4320190 to the Northern Gulf Institute from NOAA’s Office of Oceanic and Atmospheric Research, U.S. Department of Commerce. This research was carried out in part under the auspices of the Cooperative Institute for Marine and Atmospheric Studies (CIMAS) and NOAA, cooperative agreement NA20OAR4320472.**
