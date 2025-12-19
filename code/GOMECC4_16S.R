## Code for processing 16S amplicon data from the GOMECC-4 cruise
# Updated 12/18/2025
# Sean R. Anderson

# Load in packages
library(tidyverse);library(vegan);library(qiime2R)
library(phyloseq);library(reshape2);library(corrplot)
library(fantaxtic);library(RColorBrewer);library(microbiome);library(MASS)
library(factoextra);library(microeco);library(data.table);library(patchwork)
library(file2meco);library(ggpubr);library(treemap);library(geosphere);library(compositions)
library(rcartocolor);library(indicspecies);library(performance);library(speedyseq)
library(dplyr);library(sjstats);library(lmtest);library(ranacapa);library(caret)
library(mgcv);library(gratia);library(scales);library(patchwork);library(marmap)

# Load 16S ASV count table
table <- read_qza(file="/Users/seananderson/Documents/NOAA NGI/GOMECC-4/data-input/16S/table-16S-merge.qza")
count_tab <- table$data %>% as.data.frame() # Convert to data frame 
#write.csv(count_tab, "Count16S_all.csv", row.names = T)

# Load taxonomy file
taxonomy <- read_qza(file="/Users/seananderson/Documents/NOAA NGI/GOMECC-4/data-input/16S/taxonomy-16S-merge.qza")
tax_tab_16S <- taxonomy$data %>% # Convert to data frame, tab separate and rename taxa levels, and remove row with confidence values
  as.data.frame() %>%
  separate(Taxon, sep = ";", c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")) %>% 
  column_to_rownames("Feature.ID") %>%
  dplyr::select(-Confidence)
tax_tab_16S$Kingdom <- gsub("^.{0,3}", "", tax_tab_16S$Kingdom) # Clean up the taxonomy names for each level
tax_tab_16S$Phylum <- gsub("^.{0,4}", "", tax_tab_16S$Phylum)
tax_tab_16S$Class <- gsub("^.{0,4}", "", tax_tab_16S$Class)
tax_tab_16S$Order <- gsub("^.{0,4}", "", tax_tab_16S$Order)
tax_tab_16S$Family <- gsub("^.{0,4}", "", tax_tab_16S$Family)
tax_tab_16S$Genus <- gsub("^.{0,4}", "", tax_tab_16S$Genus)
tax_tab_16S$Species<- gsub("^.{0,4}", "", tax_tab_16S$Species)
#write.csv(tax_tab_16S, "Taxonomy16S_all.csv", row.names = T)

# From the prior exported count and taxonomy .csv files, manually add a new column for 16S functional groups
# Re-upload these new .csv files for downstream processing. All files provided on GitHub. 
tax_new_16S = read.csv(file = "/Users/seananderson/Documents/NOAA NGI/GOMECC-4/data-input/16S/Taxonomy16S_all.csv.gz", header = T, row.names = NULL, check.names = F, fileEncoding = "UTF-8-BOM") # File compressed to reduce size
count_new_16S = read.csv(file = "/Users/seananderson/Documents/NOAA NGI/GOMECC-4/data-input/16S/Count16S_all.csv.gz", header = T, row.names = NULL, check.names = F, fileEncoding = "UTF-8-BOM") # File compressed to reduce size
count_new_16S = count_new_16S[ order(match(count_new_16S$ASV, tax_new_16S$ASV)), ] # Match order of ASVs for both files 
rownames(count_new_16S) <- count_new_16S$ASV # Rename row names
count_new_16S <- count_new_16S[ -c(1) ]
rownames(tax_new_16S) <- tax_new_16S$ASV # Rename row names to match count file
tax_new_16S <- tax_new_16S[ -c(1) ]

# Load metadata
sample_info_tab <- read.csv("/Users/seananderson/Documents/NOAA NGI/GOMECC-4/data-input/metadata_Aug2023.csv", header = T, row.names = NULL, check.names = F, fileEncoding = "UTF-8-BOM") 
row.names(sample_info_tab) <- sample_info_tab[,1]
sample_info_tab <- sample_info_tab[,-1]
sample_info_tab$dic = as.numeric(sample_info_tab$dic)
sample_info_tab$carbonate = as.numeric(sample_info_tab$carbonate)

# Create phyloseq object
ps1 <- phyloseq(tax_table(as.matrix(tax_new_16S)), otu_table(count_new_16S, taxa_are_rows = T), sample_data(sample_info_tab)) # Create phyloseq object
taxa_names(ps1) <- paste0("bASV", seq(ntaxa(ps1))) # Rename ASVs to be sequential

# Remove any unwanted groups
ps_new_16S = subset_taxa(ps1, Order !="Chloroplast" |is.na(Order)) # The following are groups we want to remove
ps_new_16S = subset_taxa(ps_new_16S, Family !="Mitochondria" |is.na(Family))
ps_new_16S = subset_taxa(ps_new_16S, Kingdom !="Eukaryota" |is.na(Kingdom))
ps_new_16S <- subset_taxa(ps_new_16S, Phylum !="Unassigned", Prune = T)

ps_new_16S = name_na_taxa(ps_new_16S) # Add an "unassigned" label to lowest annotation

# Remove controls for now
ps_sub_16S <- subset_samples(ps_new_16S, sample_type == "seawater") # Remove controls 
ps_sub_16S = prune_samples(sample_sums(ps_sub_16S) >=5000, ps_sub_16S) # Remove samples with very low sequence read numbers

# Remove singletons (ASVs observed once across the dataset)
ps_filt_16S = filter_taxa(ps_sub_16S, function (x) {sum(x) > 1}, prune = TRUE)

# Estimate number of reads
ps_min <- min(sample_sums(ps_filt_16S))
ps_mean <- mean(sample_sums(ps_filt_16S))
ps_max <- max(sample_sums(ps_filt_16S))

# Set color palettes used for certain figures
nb.cols <- 17
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)
group = carto_pal(12, "Safe")

# Plot rarefaction curves - Figure S1
rare_16S <- suppressWarnings(ggrare(ps_filt_16S, step = 100, plot = FALSE, parallel = FALSE, se = FALSE))
rare_16S$data$depth_category <- factor(rare_16S$data$depth_category, levels = c("Surface","DCM","Deep"))
rare_16S + theme(legend.position = "none") + theme_bw() + theme(legend.position = "right") + facet_wrap(~depth_category + distance,scales = "free_y")
#ggsave(filename = "16S_rare_v2.eps", plot = last_plot(), device = "eps", path = NULL, scale = 1, width = 8, height = 5, dpi = 150)

# Rarefy to even sampling depth
ps_rare_16S <- rarefy_even_depth(ps_filt_16S, sample.size = min(sample_sums(ps_filt_16S)), rngseed = 714, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)

# Convert to Aitchison distance and prepare for clustering
ps_clr_16S <- microbiome::transform(ps_rare_16S, "clr") # Centered log transform rarefied data
euc_16S = phyloseq::distance(ps_clr_16S, method = "euclidean") # Calculate Aitchison distance
euc.table <- as.matrix(dist(euc_16S)) # Convert to matrix

# Perform hierarchical clustering and observe clusters - Figure S6
spe.ward <- hclust(euc_16S, method = 'ward.D2') # Hierarchical clustering with Ward's method
fviz_nbclust(euc.table, factoextra::hcut, method = "silhouette") # Three clusters is optimal
#ggsave(filename = "16S_cluster_silhouette_v2.eps", plot = last_plot(), device = "eps", path = NULL, scale = 1, width = 5, height = 4, dpi = 150) 

# Split the samples into clusters
sub_grp <- cutree(spe.ward, k = 3) # Cut the data based on our clusters; cluster assignments already in metadata file
sub_grp = as.data.frame(sub_grp)
table(sub_grp)

# Remove three samples that were unexpectedly grouped
ps_subset_16S = subset_samples(ps_rare_16S, sample_names(ps_rare_16S) != "GOMECC4_CAPECORAL_Sta140_Deep_C" & sample_names(ps_rare_16S) != "GOMECC4_LA_Sta38_Deep_C"  & sample_names(ps_rare_16S) != "GOMECC4_FLSTRAITS_Sta123_Surface_B")

# Create new a distance matrix after removing three samples
ps_clr2_16S <- microbiome::transform(ps_subset_16S, "clr") 
euc_16S = phyloseq::distance(ps_clr2_16S, method = "euclidean") # Aitchison distance

# Run PERMANOVA to test for significance between treatments
metadata <- as(sample_data(ps_subset_16S), "data.frame") # Subset metadata
metadata$distance=as.factor(metadata$distance) # Convert to factor
metadata$region=as.factor(metadata$region) # Convert to factor
metadata$depth_category=as.factor(metadata$depth_category) # Convert to factor
adonis2(phyloseq::distance(ps_subset_16S, method ="euclidean")~distance, data = metadata, p.adjust.m = 'holm', perm = 9999) # Run with 9999 permutations
adonis2(phyloseq::distance(ps_subset_16S, method ="euclidean")~region, data = metadata, p.adjust.m = 'holm', perm = 9999)
adonis2(phyloseq::distance(ps_subset_16S, method ="euclidean")~depth_category, data = metadata, p.adjust.m = 'holm', perm = 9999)

# Plot PCoA and color samples by cluster - Figure 2A
ordu = ordinate(ps_subset_16S, "PCoA", distance = euc_16S) 
p = plot_ordination(ps_subset_16S, ordu, color = "cluster_16S")
p$data$cluster_16S <- factor(p$data$cluster_16S, levels = c("Cluster 1","Cluster 2","Cluster 3"))
p + theme_bw() + scale_fill_manual(values = c("#BB5566","#DDAA33","#004488")) +
  geom_point(aes(fill = cluster_16S), size = 5, shape = 21, colour = "black") + 
  theme(text = element_text(size = 14)) 
#ggsave(filename = "16S_pcoa_final.pdf", plot = last_plot(), device = "pdf", path = NULL, scale = 1, width = 7, height = 5, dpi = 150)

# Estimate diversity metrics
rich_16S <- estimate_richness(ps_subset_16S, measures = c("Observed", "Shannon")) # Estimate richness
clus = sample_data(ps_subset_16S)$cluster_16S # Subset out cluster
region = sample_data(ps_subset_16S)$region # Subset out region
depth = sample_data(ps_subset_16S)$depth_meters
rich_all <- data.frame(rich_16S, region, depth, clus) # Merge richness estimates with cluster and region
rich_all$depth = as.character(rich_all$depth)
df2 = reshape2::melt(rich_all) # Melt the data format for plotting
levels(df2$variable)[match("Observed", levels(df2$variable))] <- "# of 16S ASVs" # Change label
levels(df2$variable)[match("Shannon", levels(df2$variable))] <- "Shannon diversity index" # Change label

# Plot richness and diversity in each cluster - Figure 2B
p <- ggplot(df2, aes(x = factor(clus), y = value, fill = clus))
p$data$clus <- factor(p$data$clus, levels = c("Cluster 3","Cluster 2", "Cluster 1"))
p + geom_boxplot(alpha = 1, outlier.shape = NA, color = "black") + theme_bw() + 
  theme(text = element_text(size = 14)) + ylab("Diversity values") + theme(legend.position="right") + scale_fill_manual(values = c("#004488","#DDAA33","#BB5566")) +
  geom_point(aes(fill = clus), size = 5, shape = 21, alpha= 1,colour = "black", position = position_jitterdodge(jitter.width = 0.8)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black")) + theme(axis.title.x = element_blank()) +
  facet_wrap(~variable, scales = "free_x", nrow = 2) + coord_flip() +
  geom_pwc(aes(group = clus), tip.length = 0, method = "wilcox_test", label = "p.adj", p.adjust.method = "holm")
#ggsave(filename = "Diversity_16S_new.pdf", plot = last_plot(), device = "pdf", path = NULL, scale = 1, width = 6, height = 7, dpi = 150) 

# Plot diversity metrics with respect to transect - Figure S12
p <- ggplot(df2, aes(x = factor(region), y = value, fill = region))
p$data$region <- factor(p$data$region, levels = c("27N", "FLSTRAITS","CAPECORAL", "TAMPA", "PANAMACITY", "PENSACOLA", "LA", "GALVESTON", "PAISNP", "BROWNSVILLE", "TAMPICO", "VERACRUZ", "CAMPECHE", "MERIDA", "YUCATAN", "CATOCHE", "CANCUN")) # Set the transect order
p + geom_boxplot(alpha = 0.5, outlier.shape = NA, color = "black") + theme_bw() + geom_smooth(method = "loess", se = TRUE, color = "black", aes(group = 1)) +
  theme(text = element_text(size = 14)) + ylab("Diversity values") + theme(legend.position = "right") + scale_fill_manual(values = mycolors) +
  geom_point(aes(fill = region), size = 3, shape = 21, colour = "black", position=position_jitterdodge(), show.legend = FALSE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black")) + theme(axis.title.x = element_blank()) + facet_wrap(clus ~ variable, scales = "free_y")
#ggsave(filename = "Diversity_16S_region_v1.pdf", plot = last_plot(), device = "pdf", path = NULL, scale = 1, width = 20, height = 6, dpi = 150) 

# Plot diversity metrics with respect to absolute depth - Figure S10
df2$depth <- as.numeric(as.character(df2$depth)) # Convert to numeric for correct sorting
df2$depth <- factor(df2$depth, levels = rev(sort(unique(df2$depth)))) # Re-convert to factor with levels reversed (so smallest at top after flip)
p <- ggplot(df2, aes(x = depth, y = value, fill = clus)) + geom_point(aes(fill = clus), size = 3, shape = 21, colour = "black")  + 
scale_fill_manual(values = c("#BB5566","#DDAA33","#004488")) + geom_smooth(method = "loess", se = TRUE, color = "black", aes(group = 1)) +
scale_x_discrete(labels = label_wrap(0.1)) + coord_flip() + facet_wrap(~variable, scales = "free_x") + guides(fill = FALSE) + theme_bw()
p
#ggsave(filename = "16S_div_depth.pdf", plot = last_plot(), device = "pdf", path = NULL, scale = 1, width = 8, height = 6, dpi = 150)

# Plot relative abundance stacked bar charts - Figure 2C
tax_table(ps_subset_16S) <- tax_table(ps_subset_16S)[,2:8] # Subset out the functional category column and focus on taxonomy levels
barplot <- ps_subset_16S %>%
  tax_glom(taxrank = "Order", NArm = FALSE) %>% # Agglomerate to order level
  transform_sample_counts(function(OTU) 100* OTU/sum(OTU)) %>% # Transform to relative abundance
  psmelt() %>% # Melt data
  group_by(region,Order,cluster_16S) %>% # Group by cluster and transect
  summarise_at("Abundance", .funs = mean) # Summarize at the mean

focus <- c("SAR11_clade", "Synechococcales", "Marinimicrobia_(SAR406_clade)", "Nitrosopumilales", "SAR86_clade","Flavobacteriales","Actinomarinales","Marine_Group_II","Rhodospirillales","SAR324_clade(Marine_group_B)","Thiomicrospirales","Puniceispirillales") # Focus on top 12 groups
barplot$Order <- ifelse(barplot$Order %in% focus, barplot$Order, "Others") # Others category for plotting
barplot_16S = barplot # Rename
barplot_16S$Order<- as.character(barplot_16S$Order) # Convert to character

p <- ggplot(data = barplot_16S, aes(x = region, y=Abundance, fill = Order))
p$data$Order <- factor(p$data$Order, levels = c("Others","Puniceispirillales","Thiomicrospirales","SAR324_clade(Marine_group_B)","Rhodospirillales","Marine_Group_II","Actinomarinales","Flavobacteriales", "SAR86_clade","Nitrosopumilales","Marinimicrobia_(SAR406_clade)", "Synechococcales","SAR11_clade" )) # Set order of groups in the plot
p$data$region <- factor(p$data$region, levels = c("27N", "FLSTRAITS","CAPECORAL", "TAMPA", "PANAMACITY", "PENSACOLA", "LA", "GALVESTON", "PAISNP", "BROWNSVILLE", "TAMPICO", "VERACRUZ", "CAMPECHE", "MERIDA", "YUCATAN", "CATOCHE", "CANCUN")) # Set order of transects
p + geom_bar(aes(), stat = "identity", position = "fill", width = 0.9) +
  scale_y_continuous(expand = c(0, 0))+ geom_hline(yintercept=0) + theme_bw() + 
  scale_fill_manual(values = rev(c("#88CCEE", "#CC6677", "#44AA99" , "#882255" , "#DDCC77","#332288" , "#117733","#999933" , "#AA4499","#F5793A","#F7CDA4", "#A5CFCC","#757575"))) + 
  theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust = 1)) + 
  theme(legend.position = "right") + theme(text = element_text(size = 12)) + guides(fill = guide_legend(nrow = 14, ncol = 1)) + 
  facet_wrap(~cluster_16S,scales = "free_x") + labs(y = "Relative abundance (%)")
#ggsave(filename = "16S_stacked_class.pdf", plot = last_plot(), device = "pdf", path = NULL, scale = 1, width = 14, height = 5, dpi = 150) 

# Prepare for taxonomy tree maps - Figure S14
ps <- tax_glom(ps_subset_16S, "Genus", NArm = FALSE) # Agglomerate to genus level
x1 = speedyseq::psmelt(ps) # Melt the data
focus <- c("Proteobacteria", "Cyanobacteria", "Marinimicrobia_(SAR406_clade)", "Crenarchaeota", "Bacteroidota", "Actinobacteriota") # Focus on top phylum level groups
x1$Phylum <- ifelse(x1$Phylum %in% focus, x1$Phylum, "Others") # Others category for plotting
x1$Phylum <- factor(x1$Phylum, levels = c("Actinobacteriota", "Bacteroidota", "Crenarchaeota","Cyanobacteria","Marinimicrobia_(SAR406_clade)","Proteobacteria","Others")) # Set order for plotting
cluster1 <- x1[x1[["cluster_16S"]] == "Cluster 1", ] # Subset to Cluster 1
cluster2 <- x1[x1[["cluster_16S"]] == "Cluster 2", ] # Subset to Cluster 2
cluster3 <- x1[x1[["cluster_16S"]] == "Cluster 3", ] # Subset to Cluster 3

# Plot tree map for Cluster 1 (photic zone) 
#pdf("Cluster1_treemap_16S_v2.pdf", width = 12, height = 4)
treemap(dtf = cluster1,
        title = "Cluster 1 (2-99 m)", 
        algorithm = "pivotSize", border.lwds = c(2,0.5,0.1),
        border.col = c("black", "black", "black"),
        mapping = c(0,0,0),
        index = c("Order", "Genus"),
        vSize = "Abundance",
        vColor = "Phylum",
        palette = "Dark2",
        type="categorical",
        fontsize.labels = c(10,8),
        fontface.labels = c(2,1),
        fontcolor.labels = c("white","white"),
        bg.labels = 255, 
        position.legend = "bottom",
        align.labels = list(
          c("left", "top"), 
          c("center", "center")),
        inflate.labels = F,
        overlap.labels = 0.8,
        lowerbound.cex.labels = 0.2,
        force.print.labels = F) 
#dev.off()

# Plot tree map for Cluster 2 (offshore DCM)
#pdf("Cluster2_treemap_16S_v2.pdf", width = 12, height = 4)
treemap(dtf = cluster2,
        title = "Cluster 2 (2-124 m)", 
        algorithm = "pivotSize", border.lwds = c(2,0.5,0.1),
        border.col = c("black", "black", "black"),
        mapping = c(0,0,0),
        index = c("Order", "Genus"),
        vSize = "Abundance",
        vColor = "Phylum",
        palette = "Dark2",
        type="categorical",
        fontsize.labels = c(10,8),
        fontface.labels = c(2,1),
        fontcolor.labels = c("white","white"),
        bg.labels = 255, 
        position.legend = "bottom",
        align.labels=  list(
          c("left", "top"), 
          c("center", "center")),
        inflate.labels = F,
        overlap.labels = 0.8,
        lowerbound.cex.labels = 0.2,
        force.print.labels = F) 
#dev.off()

# Plot tree map for Cluster 3 (Aphotic zone)
#pdf("Cluster3_treemap_16S_v2.pdf", width = 12, height = 4)
treemap(dtf = cluster3,
        title = "Cluster 3 (135-3326 m)", 
        algorithm = "pivotSize", border.lwds = c(2,0.5,0.1),
        border.col = c("black", "black", "black"),
        mapping = c(0,0,0),
        index = c("Order", "Genus"),
        vSize = "Abundance",
        vColor = "Phylum",
        palette = "Dark2",
        type="categorical",
        fontsize.labels = c(10,8),
        fontface.labels = c(2,1),
        fontcolor.labels=  c("white","white"),
        bg.labels = 255, 
        position.legend = "bottom",
        align.labels = list(
          c("left", "top"), 
          c("center", "center")),
        inflate.labels = F,
        overlap.labels = 0.8,
        lowerbound.cex.labels = 0.2,
        force.print.labels = F) 
#dev.off()

# Perform centered log-ratio (CLR) transformation of the 16S ASV count data and save as separate data frame
otu <- as.data.frame(t(otu_table(ps_subset_16S)))
otu_clr <- clr(otu + 1) # CLR transform with pseudo-count 
otu_clr <- as.data.frame(otu_clr) # Convert to data frame

# Format data for generalized additive models (GAMs)
order_16S <- tax_glom(ps_subset_16S, taxrank = "Order", NArm = FALSE) # Aggregate at the order level
x1 = psmelt(order_16S) # Melt the data
otu_clr$Sample <- rownames(otu_clr) # Define sample row in CLR table
otu_clr_long <- reshape2::melt(otu_clr, id.vars = "Sample", variable.name = "OTU", value.name = "CLRAbundance") # Melt to format this table
x1 <- dplyr::left_join(x1, otu_clr_long, by = c("Sample", "OTU")) # Merge with melted phyloseq table; now CLR abundance is its own column

cluster1 <- x1[x1[["cluster_16S"]] == "Cluster 1", ] # Subset to Cluster 1
glm.1 = cluster1 %>% # Group the data based on class and preserve sample replication for the models
  dplyr::group_by(station, depth_category, Order, replicate, region) %>%
  summarise(temp = mean(temp),
            salinity= mean(salinity),
            oxygen = mean(oxygen),
            phosphate = mean(phosphate),
            nitrate = mean(nitrate),
            nitrite = mean(nitrite),
            silicate = mean(silicate),
            nh4 = mean(nh4),
            pH_corrected = mean(pH_corrected),
            total_alkalinity = mean(total_alkalinity),
            OMEGA_AR = mean(OMEGA_AR),
            dic = mean(dic),
            pCO2_corrected = mean(pCO2_corrected),
            carbonate = mean(carbonate),
            fluorescence = mean(fluorescence),
            Abundance = mean(Abundance),
            CLRAbundance = mean(CLRAbundance),
            .groups = 'drop')

# Repeat data format for genus level to include major cyanobacteria
genus_16S <- tax_glom(ps_subset_16S, taxrank = "Genus", NArm = FALSE) # Aggregate at the genus level
x2 = psmelt(genus_16S) # Melt the data
x2 <- dplyr::left_join(x2, otu_clr_long, by = c("Sample", "OTU")) # Merge with melted phyloseq table; now CLR abundance is its own column

cluster1 <- x2[x2[["cluster_16S"]] == "Cluster 1", ] # Subset to Cluster 1
glm.2 = cluster1 %>% # Group the data based on class and preserve sample replication for the models
  dplyr::group_by(station, depth_category, Genus, replicate, region) %>%
  summarise(temp = mean(temp),
            salinity = mean(salinity),
            oxygen = mean(oxygen),
            phosphate = mean(phosphate),
            nitrate = mean(nitrate),
            nitrite = mean(nitrite),
            silicate = mean(silicate),
            nh4 = mean(nh4),
            pH_corrected = mean(pH_corrected),
            total_alkalinity = mean(total_alkalinity),
            OMEGA_AR = mean(OMEGA_AR),
            dic = mean(dic),
            pCO2_corrected = mean(pCO2_corrected),
            carbonate = mean(carbonate),
            fluorescence = mean(fluorescence),
            Abundance = mean(Abundance),
            CLRAbundance = mean(CLRAbundance),
            .groups = 'drop')

# Select variables for initial models based on low collinearity
df1 = glm.1[,c(6:20)] # Subset to include environmental variables
df1 = na.omit(df1) # Omit any rows that have NA values
correlations = cor(df1, method = "spearman") # Perform Spearman correlations to assess collinearity
#write.csv(correlations, "Cluster1_corr_16S.csv",row.names = T) # Write .csv file - Table S3

df1_filt = glm.1[,c(3,6:8,9,10,13:14,17,21)] # Initial list of variables that were not collinear (Spearman < 0.7 or > -0.7)
df1_filt = na.omit(df1_filt) # Omit any rows that have NA values
model1 <- lm(Abundance ~., data = df1_filt) # Test model for VIFs
car::vif(model1) # Display VIFs from test model - VIFs should be < 10

df1_filt = glm.1[,c(3,5:8,9,10,13:14,17,21,22)] # Re-do to add back in the region column needed for outlier removal and CLR abundance column
df1_filt = na.omit(df1_filt) 

# Repeat this process for the genus level groups
df1_filt2 = glm.2[,c(3,6:8,9,10,13:14,17,21)]
df1_filt2 = na.omit(df1_filt2) # Omit any rows that have NA values
model1 <- lm(Abundance ~., data = df1_filt2)
car::vif(model1) # VIFs look good

df1_filt2 = glm.2[,c(3,5:8,9,10,13:14,17,21,22)] # Re-do to add back in the region column needed for outlier removal and CLR column
df1_filt2 = na.omit(df1_filt2) 

# Subset to the taxonomic groups of interest
df_sar <- df1_filt[df1_filt[["Order"]] == "SAR11_clade", ]
df_sar86 <- df1_filt[df1_filt[["Order"]] == "SAR86_clade", ]
df_flavo <- df1_filt[df1_filt[["Order"]] == "Flavobacteriales", ]
df_syne <- df1_filt2[df1_filt2[["Genus"]] == "Synechococcus_CC9902", ] # Synechococcus
df_pro <- df1_filt2[df1_filt2[["Genus"]] == "Prochlorococcus_MIT9313", ] # Prochlorococcus

## Workflow for group-specific 16S GAMs in the photic zone
# Remove outliers from the 27N line based on strict IQR
remove_region_outliers <- function(df, region_col = "region", target_region = "27N", vars, iqr_multiplier = 3) {
  df_clean <- df
  
  for (var in vars) {
    if (!var %in% colnames(df_clean)) next
    
# Subset 27N region
    region_data <- df_clean[df_clean[[region_col]] == target_region, ]
    
# Calculate IQR bounds 
    q1 <- quantile(region_data[[var]], 0.25, na.rm = TRUE)
    q3 <- quantile(region_data[[var]], 0.75, na.rm = TRUE)
    iqr <- q3 - q1
    lower <- q1 - iqr_multiplier * iqr
    upper <- q3 + iqr_multiplier * iqr
    
# Keep samples outside of outlier range
    keep_idx <- !(df_clean[[region_col]] == target_region & 
                    (df_clean[[var]] < lower | df_clean[[var]] > upper))
    df_clean <- df_clean[keep_idx, ]
  }
  
  return(df_clean)
}

# Variables to check as outliers in 27N line
predictor_vars <- c("temp", "salinity", "oxygen", "phosphate", "nh4", "nitrate", "pH_corrected", "dic")

# Clean the 27N extreme outliers only from all groups - 6 samples removed
df_sar_filt <- remove_region_outliers(df_sar, region_col = "region", target_region = "27N", vars = predictor_vars, iqr_multiplier = 3)
df_syne_filt <- remove_region_outliers(df_syne, region_col = "region", target_region = "27N", vars = predictor_vars, iqr_multiplier = 3)
df_pro_filt <- remove_region_outliers(df_pro, region_col = "region", target_region = "27N", vars = predictor_vars, iqr_multiplier = 3)
df_sar86_filt <- remove_region_outliers(df_sar86, region_col = "region", target_region = "27N", vars = predictor_vars, iqr_multiplier = 3)
df_flavo_filt <- remove_region_outliers(df_flavo, region_col = "region", target_region = "27N", vars = predictor_vars, iqr_multiplier = 3)

## Run through GAM workflow for each group - start with SAR11
# Convert to log abundance and log transform variables
df_sar_filt$LogAbundance <- log(df_sar_filt$Abundance + 1) # Log transform rarefied abundance
log_vars <- c("temp", "salinity", "phosphate", "nh4", "nitrate") # Select predictor variables to log transform

for (var in log_vars) {
  new_var <- paste0("log_", var)
  df_sar_filt[[new_var]] <- log(df_sar_filt[[var]] + 0.01)
} # Loop over the variables and transform

set.seed(1234) # Set for reproducibility

# Define predictors again with new log transformed variables
predictors <- c("log_temp", "log_salinity", "oxygen", "log_phosphate", "log_nh4", "log_nitrate", "pH_corrected", "dic")

# Stratified split of data into train and test sets
df_sar_filt$Strata <- interaction(df_sar_filt$region, cut(df_sar_filt$LogAbundance, breaks = 2)) # Cut data based on region and abundance
train_index <- suppressWarnings(createDataPartition(df_sar_filt$Strata, p = 0.8, list = FALSE)) # Define train set as 80%; avoid warning message, some samples no longer there after the outlier removal
train_data <- df_sar_filt[train_index, ] # Train set
test_data  <- df_sar_filt[-train_index, ] # Test set (20%)

# Create all non-empty predictor combinations
predictor_combos <- unlist(lapply(1:length(predictors), function(i) combn(predictors, i, simplify = FALSE)), recursive = FALSE)

# RMSE function to use with abundance data
rmse <- function(actual, predicted) sqrt(mean((actual - predicted)^2))

# Cross fold validation for Gaussian GAM (log abundance)
cv_gam_rmse <- function(vars, data, folds = 10) {
  set.seed(123)
  data <- data[sample(nrow(data)), ]  # Shuffle
  data$Fold <- sample(rep(1:folds, length.out = nrow(data)))
  rmses <- numeric(folds) # Use 10 folds
  
  for (i in 1:folds) {
    train_fold <- data[data$Fold != i, ] # Loops through each fold and partition into training and test
    test_fold  <- data[data$Fold == i, ]
    
# Create GAM formula based on chosen predictors
    formula_obj <- as.formula(paste("LogAbundance ~", paste(paste0("s(`", vars, "`, k = 10)"), collapse = " + ")))
    
# Fit GAM on training fold
    model <- try(gam(formula_obj, data = train_fold, family = gaussian(), select = TRUE, method = "REML"), silent = TRUE)
    if (inherits(model, "try-error")) {
      rmses[i] <- NA # Report NA for a fold if fitting fails
    } else {
      
# Predict on test fold and compute RMSE
      preds <- predict(model, newdata = test_fold, type = "response")
      rmses[i] <- sqrt(mean((test_fold$LogAbundance - preds)^2))
    }
  }
  
# Return the mean RMSE across all folds, ignoring failed folds (NA)
  return(mean(rmses, na.rm = TRUE))
}

# Run cross fold on training data
cv_results <- data.frame(Predictors = sapply(predictor_combos, paste, collapse = " + "), RMSE = sapply(predictor_combos, function(vars) cv_gam_rmse(vars, train_data)))

# Rank by RMSE
cv_results <- cv_results[order(cv_results$RMSE), ]

# Fit best GAM from above
best_vars <- unlist(strsplit(cv_results$Predictors[1], " + ", fixed = TRUE))
best_formula <- as.formula(paste("LogAbundance ~", paste(paste0("s(", best_vars, ", k = 10)"), collapse = " + ")))
final_gam_sar <- gam(best_formula, data = train_data, family = gaussian(), select = TRUE, method = "REML") 

# Extract p-values from smooth terms and adjust for multiple comparison to assess significance
gam_summary <- summary(final_gam_sar)
smooth_pvals <- gam_summary$s.table[, "p-value"]
adj_pvals <- p.adjust(smooth_pvals, method = "holm") # Holm correction

# Table of smooth results
smooth_results <- data.frame(
  Term = rownames(gam_summary$s.table),
  EDF = gam_summary$s.table[, "edf"],
  Stat = gam_summary$s.table[, if ("F" %in% colnames(gam_summary$s.table)) "F" else "Chi.sq"],
  RawP = smooth_pvals,
  AdjP = adj_pvals,
  Sig = adj_pvals < 0.05 # TRUE/FALSE significance
) 
print(smooth_results)

# Predict on train and test sets
train_data$Predicted <- predict(final_gam_sar, newdata = train_data, type = "response")
test_data$Predicted  <- predict(final_gam_sar, newdata = test_data, type = "response")

# Evaluate model performance with RMSE and spearman
train_rmse <- rmse(train_data$LogAbundance, train_data$Predicted)
test_rmse  <- rmse(test_data$LogAbundance, test_data$Predicted)
spearman_train <- cor(train_data$LogAbundance, train_data$Predicted, method = "spearman")
spearman_test  <- cor(test_data$LogAbundance, test_data$Predicted, method = "spearman")

# Estimate deviance explained for train and test sets
null_deviance_train <- sum((train_data$LogAbundance - mean(train_data$LogAbundance))^2)
resid_deviance_train <- sum((train_data$LogAbundance - train_data$Predicted)^2)
dev_expl_train <- 1 - resid_deviance_train / null_deviance_train

null_deviance_test <- sum((test_data$LogAbundance - mean(train_data$LogAbundance))^2)
resid_deviance_test <- sum((test_data$LogAbundance - test_data$Predicted)^2)
dev_expl_test <- 1 - resid_deviance_test / null_deviance_test

# Combine all test and training data for plotting
train_data$Set <- "Train"
test_data$Set <- "Test"
combined_data <- rbind(train_data, test_data)
combined_data$Fill <- ifelse(combined_data$Set == "Train", "#88CCEE", "gray70") # Set colors
combined_data$Shape <- ifelse(combined_data$Set == "Train", 21, 21) # Set plot parameters
combined_data$Alpha <- 1 # Set plot parameters
equation_text <- paste("Log(Abundance + 1) ~", paste(paste0("s(", best_vars, ")"), collapse = " + ")) # Set GAM equation text to display in plot

# Plot observed vs. predicted for train and test for Syndiniales - Figure 4A
ggplot(combined_data, aes(x = LogAbundance, y = Predicted)) +
  geom_point(aes(shape = factor(Shape), fill = Fill, alpha = Alpha), size = 5, stroke = 1, color = "black", show.legend = TRUE) +
  scale_shape_manual(values = c("21" = 21, "24" = 24)) + scale_fill_identity() + scale_alpha_identity() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") + annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.5,
  label = paste0("ρ (Train) = ", round(spearman_train, 2), "\nρ (Test) = ", round(spearman_test, 2)), size = 4.5) +
  annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -1.5, label = paste0("RMSE (Train) = ", round(train_rmse, 2), "\nRMSE (Test) = ", round(test_rmse, 2)), size = 4.5) + 
  labs(title = "Observed vs. Predicted Log Abundance (GAM)", subtitle = equation_text, x = "Observed log(Abundance + 1)", y = "Predicted log(Abundance + 1)", fill = "Data Set") +
  xlim(6.9, 8) + ylim(6.9, 8) + theme_bw()
#ggsave(filename = "Sar_gam_final.pdf", plot = last_plot(), device = "pdf", path = NULL, scale = 1, width = 6, height = 5, dpi = 150)

# Repeat GAM workflow for Prochlorococcus
df_pro_filt$LogAbundance <- log(df_pro_filt$Abundance + 1) # Log transform rarefied abundance
log_vars <- c("temp", "salinity", "phosphate", "nh4", "nitrate") # Select predictor variables to log transform

for (var in log_vars) {
  new_var <- paste0("log_", var)
  df_pro_filt[[new_var]] <- log(df_pro_filt[[var]] + 0.01)
} # Loop over the variables and transform

set.seed(1234) # Set for reproducibility

# Define predictors again with new log transformed variables
predictors <- c("log_temp", "log_salinity", "oxygen", "log_phosphate", "log_nh4", "log_nitrate", "pH_corrected", "dic")

# Stratified split of data into train and test sets
df_pro_filt$Strata <- interaction(df_pro_filt$region, cut(df_pro_filt$LogAbundance, breaks = 2)) # Cut data based on region and abundance
train_index <- suppressWarnings(createDataPartition(df_pro_filt$Strata, p = 0.8, list = FALSE)) # Define train set as 80%
train_data <- df_pro_filt[train_index, ] # Train set
test_data  <- df_pro_filt[-train_index, ] # Test set (20%)

# Create all non-empty predictor combinations
predictor_combos <- unlist(lapply(1:length(predictors), function(i) combn(predictors, i, simplify = FALSE)), recursive = FALSE)

# RMSE function to use with abundance data
rmse <- function(actual, predicted) sqrt(mean((actual - predicted)^2))

# Cross fold validation for Gaussian GAM (log abundance)
cv_gam_rmse <- function(vars, data, folds = 10) {
  set.seed(123)
  data <- data[sample(nrow(data)), ]  # Shuffle
  data$Fold <- sample(rep(1:folds, length.out = nrow(data)))
  rmses <- numeric(folds) # Use 10 folds
  
  for (i in 1:folds) {
    train_fold <- data[data$Fold != i, ] # Loops through each fold and partition into training and test
    test_fold  <- data[data$Fold == i, ]
    
# Create GAM formula based on chosen predictors
    formula_obj <- as.formula(paste("LogAbundance ~", paste(paste0("s(`", vars, "`, k = 10)"), collapse = " + ")))
    
# Fit GAM on training fold
    model <- try(gam(formula_obj, data = train_fold, family = gaussian(), select = TRUE, method = "REML"), silent = TRUE)
    if (inherits(model, "try-error")) {
      rmses[i] <- NA # Report NA for a fold if fitting fails
    } else {
      
# Predict on test fold and compute RMSE
      preds <- predict(model, newdata = test_fold, type = "response")
      rmses[i] <- sqrt(mean((test_fold$LogAbundance - preds)^2))
    }
  }
  
# Return the mean RMSE across all folds, ignoring failed folds (NA)
  return(mean(rmses, na.rm = TRUE))
}

# Run cross fold on training data
cv_results <- data.frame(Predictors = sapply(predictor_combos, paste, collapse = " + "), RMSE = sapply(predictor_combos, function(vars) cv_gam_rmse(vars, train_data)))

# Rank by RMSE
cv_results <- cv_results[order(cv_results$RMSE), ]

# Fit best GAM
best_vars <- unlist(strsplit(cv_results$Predictors[1], " + ", fixed = TRUE))
best_formula <- as.formula(paste("LogAbundance ~", paste(paste0("s(", best_vars, ", k = 10)"), collapse = " + ")))
final_gam_pro <- gam(best_formula, data = train_data, family = gaussian(), select = TRUE, method = "REML") 

# Extract p-values from smooth terms and adjust for multiple comparison to assess significance
gam_summary <- summary(final_gam_pro)
smooth_pvals <- gam_summary$s.table[, "p-value"]
adj_pvals <- p.adjust(smooth_pvals, method = "holm") # Holm correction

# Table of smooth results
smooth_results <- data.frame(
  Term = rownames(gam_summary$s.table),
  EDF = gam_summary$s.table[, "edf"],
  Stat = gam_summary$s.table[, if ("F" %in% colnames(gam_summary$s.table)) "F" else "Chi.sq"],
  RawP = smooth_pvals,
  AdjP = adj_pvals,
  Sig = adj_pvals < 0.05 # TRUE/FALSE significance
) 
print(smooth_results)

# Predict on train and test sets
train_data$Predicted <- predict(final_gam_pro, newdata = train_data, type = "response")
test_data$Predicted  <- predict(final_gam_pro, newdata = test_data, type = "response")

# Evaluate model performance with RMSE and spearman
train_rmse <- rmse(train_data$LogAbundance, train_data$Predicted)
test_rmse  <- rmse(test_data$LogAbundance, test_data$Predicted)
spearman_train <- cor(train_data$LogAbundance, train_data$Predicted, method = "spearman")
spearman_test  <- cor(test_data$LogAbundance, test_data$Predicted, method = "spearman")

# Estimate deviance explained for train and test sets
null_deviance_train <- sum((train_data$LogAbundance - mean(train_data$LogAbundance))^2)
resid_deviance_train <- sum((train_data$LogAbundance - train_data$Predicted)^2)
dev_expl_train <- 1 - resid_deviance_train / null_deviance_train

null_deviance_test <- sum((test_data$LogAbundance - mean(train_data$LogAbundance))^2)
resid_deviance_test <- sum((test_data$LogAbundance - test_data$Predicted)^2)
dev_expl_test <- 1 - resid_deviance_test / null_deviance_test

# Combine all test and training data for plotting
train_data$Set <- "Train"
test_data$Set <- "Test"
combined_data <- rbind(train_data, test_data)
combined_data$Fill <- ifelse(combined_data$Set == "Train", "#CC6677", "gray70")
combined_data$Shape <- ifelse(combined_data$Set == "Train", 21, 21)
combined_data$Alpha <- 1
equation_text <- paste("Log(Abundance + 1) ~", paste(paste0("s(", best_vars, ")"), collapse = " + "))

# Plot observed vs. predicted for train and test for Prochlorococcus - Figure 4B
ggplot(combined_data, aes(x = LogAbundance, y = Predicted)) +
  geom_point(aes(shape = factor(Shape), fill = Fill, alpha = Alpha), size = 5, stroke = 1, color = "black", show.legend = TRUE) +
  scale_shape_manual(values = c("21" = 21, "24" = 24)) + scale_fill_identity() + scale_alpha_identity() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") + annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.5,
  label = paste0("ρ (Train) = ", round(spearman_train, 2), "\nρ (Test) = ", round(spearman_test, 2)), size = 4.5) +
  annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -1.5, label = paste0("RMSE (Train) = ", round(train_rmse, 2), "\nRMSE (Test) = ", round(test_rmse, 2)), size = 4.5) + 
  labs(title = "Observed vs. Predicted Log Abundance (GAM)", subtitle = equation_text, x = "Observed log(Abundance + 1)", y = "Predicted log(Abundance + 1)", fill = "Data Set") +
  xlim(0, 7.9) + ylim(0, 7.9) + theme_bw()
#ggsave(filename = "Pro_gam_final.pdf", plot = last_plot(), device = "pdf", path = NULL, scale = 1, width = 6, height = 5, dpi = 150)

# Repeat GAM workflow for Synechococcus
df_syne_filt$LogAbundance <- log(df_syne_filt$Abundance + 1) # Log transform rarefied abundance
log_vars <- c("temp", "salinity", "phosphate", "nh4", "nitrate") # Select predictor variables to log transform

for (var in log_vars) {
  new_var <- paste0("log_", var)
  df_syne_filt[[new_var]] <- log(df_syne_filt[[var]] + 0.01)
} # Loop over the variables and transform

set.seed(1234) # Set for reproducibility

# Define predictors again with new log transformed variables
predictors <- c("log_temp", "log_salinity", "oxygen", "log_phosphate", "log_nh4", "log_nitrate", "pH_corrected", "dic")

# Stratified split of data into train and test sets
df_syne_filt$Strata <- interaction(df_syne_filt$region, cut(df_syne_filt$LogAbundance, breaks = 2)) # Cut data based on region and abundance
train_index <- suppressWarnings(createDataPartition(df_syne_filt$Strata, p = 0.8, list = FALSE)) # Define train set as 80%
train_data <- df_syne_filt[train_index, ] # Train set
test_data  <- df_syne_filt[-train_index, ] # Test set (20%)

# Create all non-empty predictor combinations
predictor_combos <- unlist(lapply(1:length(predictors), function(i) combn(predictors, i, simplify = FALSE)), recursive = FALSE)

# RMSE function to use with abundance data
rmse <- function(actual, predicted) sqrt(mean((actual - predicted)^2))

# Cross fold validation for Gaussian GAM (log abundance)
cv_gam_rmse <- function(vars, data, folds = 10) {
  set.seed(123)
  data <- data[sample(nrow(data)), ]  # Shuffle
  data$Fold <- sample(rep(1:folds, length.out = nrow(data)))
  rmses <- numeric(folds) # Use 10 folds
  
  for (i in 1:folds) {
    train_fold <- data[data$Fold != i, ] # Loops through each fold and partition into training and test
    test_fold  <- data[data$Fold == i, ]
    
# Create GAM formula based on chosen predictors
    formula_obj <- as.formula(paste("LogAbundance ~", paste(paste0("s(`", vars, "`, k = 10)"), collapse = " + ")))
    
# Fit GAM on training fold
    model <- try(gam(formula_obj, data = train_fold, family = gaussian(), select = TRUE, method = "REML"), silent = TRUE)
    if (inherits(model, "try-error")) {
      rmses[i] <- NA # Report NA for a fold if fitting fails
    } else {
      
# Predict on test fold and compute RMSE
      preds <- predict(model, newdata = test_fold, type = "response")
      rmses[i] <- sqrt(mean((test_fold$LogAbundance - preds)^2))
    }
  }
  
# Return the mean RMSE across all folds, ignoring failed folds (NA)
  return(mean(rmses, na.rm = TRUE))
}

# Run cross fold on training data
cv_results <- data.frame(Predictors = sapply(predictor_combos, paste, collapse = " + "), RMSE = sapply(predictor_combos, function(vars) cv_gam_rmse(vars, train_data)))

# Rank by RMSE
cv_results <- cv_results[order(cv_results$RMSE), ]

# Fit best GAM
best_vars <- unlist(strsplit(cv_results$Predictors[1], " + ", fixed = TRUE))
best_formula <- as.formula(paste("LogAbundance ~", paste(paste0("s(", best_vars, ", k = 10)"), collapse = " + ")))
final_gam_syne <- gam(best_formula, data = train_data, family = gaussian(), select = TRUE, method = "REML") 

# Extract p-values from smooth terms and adjust for multiple comparison to assess significance
gam_summary <- summary(final_gam_syne)
smooth_pvals <- gam_summary$s.table[, "p-value"]
adj_pvals <- p.adjust(smooth_pvals, method = "holm") # Holm correction

# Table of smooth results
smooth_results <- data.frame(
  Term = rownames(gam_summary$s.table),
  EDF = gam_summary$s.table[, "edf"],
  Stat = gam_summary$s.table[, if ("F" %in% colnames(gam_summary$s.table)) "F" else "Chi.sq"],
  RawP = smooth_pvals,
  AdjP = adj_pvals,
  Sig = adj_pvals < 0.05 # TRUE/FALSE significance
) 
print(smooth_results)

# Predict on train and test sets
train_data$Predicted <- predict(final_gam_syne, newdata = train_data, type = "response")
test_data$Predicted  <- predict(final_gam_syne, newdata = test_data, type = "response")

# Evaluate model performance with RMSE and spearman
train_rmse <- rmse(train_data$LogAbundance, train_data$Predicted)
test_rmse  <- rmse(test_data$LogAbundance, test_data$Predicted)
spearman_train <- cor(train_data$LogAbundance, train_data$Predicted, method = "spearman")
spearman_test  <- cor(test_data$LogAbundance, test_data$Predicted, method = "spearman")

# Estimate deviance explained for train and test sets
null_deviance_train <- sum((train_data$LogAbundance - mean(train_data$LogAbundance))^2)
resid_deviance_train <- sum((train_data$LogAbundance - train_data$Predicted)^2)
dev_expl_train <- 1 - resid_deviance_train / null_deviance_train

null_deviance_test <- sum((test_data$LogAbundance - mean(train_data$LogAbundance))^2)
resid_deviance_test <- sum((test_data$LogAbundance - test_data$Predicted)^2)
dev_expl_test <- 1 - resid_deviance_test / null_deviance_test

# Combine all test and training data for plotting
train_data$Set <- "Train"
test_data$Set <- "Test"
combined_data <- rbind(train_data, test_data)
combined_data$Fill <- ifelse(combined_data$Set == "Train", "#732633", "gray70")
combined_data$Shape <- ifelse(combined_data$Set == "Train", 21, 21)
combined_data$Alpha <- 1
equation_text <- paste("Log(Abundance + 1) ~", paste(paste0("s(", best_vars, ")"), collapse = " + "))

# Plot observed vs. predicted for train and test for Synechococcus - Figure 4C
ggplot(combined_data, aes(x = LogAbundance, y = Predicted)) +
  geom_point(aes(shape = factor(Shape), fill = Fill, alpha = Alpha), size = 5, stroke = 1, color = "black", show.legend = TRUE) +
  scale_shape_manual(values = c("21" = 21, "24" = 24)) + scale_fill_identity() + scale_alpha_identity() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") + annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.5,
  label = paste0("ρ (Train) = ", round(spearman_train, 2), "\nρ (Test) = ", round(spearman_test, 2)), size = 4.5) +
  annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -1.5, label = paste0("RMSE (Train) = ", round(train_rmse, 2), "\nRMSE (Test) = ", round(test_rmse, 2)), size = 4.5) + 
  labs(title = "Observed vs. Predicted Log Abundance (GAM)", subtitle = equation_text, x = "Observed log(Abundance + 1)", y = "Predicted log(Abundance + 1)", fill = "Data Set") +
  xlim(2, 8) + ylim(2, 8) + theme_bw()
#ggsave(filename = "Syne_gam_final.pdf", plot = last_plot(), device = "pdf", path = NULL, scale = 1, width = 6, height = 5, dpi = 150)

# Repeat GAM workflow for SAR86
df_sar86_filt$LogAbundance <- log(df_sar86_filt$Abundance + 1) # Log transform rarefied abundance
log_vars <- c("temp", "salinity", "phosphate", "nh4", "nitrate") # Select predictor variables to log transform

for (var in log_vars) {
  new_var <- paste0("log_", var)
  df_sar86_filt[[new_var]] <- log(df_sar86_filt[[var]] + 0.01)
} # Loop over the variables and transform

set.seed(1234) # Set for reproducibility

# Define predictors again with new log transformed variables
predictors <- c("log_temp", "log_salinity", "oxygen", "log_phosphate", "log_nh4", "log_nitrate", "pH_corrected", "dic")

# Stratified split of data into train and test sets
df_sar86_filt$Strata <- interaction(df_sar86_filt$region, cut(df_sar86_filt$LogAbundance, breaks = 2)) # Cut data based on region and abundance
train_index <- suppressWarnings(createDataPartition(df_sar86_filt$Strata, p = 0.8, list = FALSE)) # Define train set as 80%
train_data <- df_sar86_filt[train_index, ] # Train set
test_data <- df_sar86_filt[-train_index, ] # Test set (20%)

# Create all non-empty predictor combinations
predictor_combos <- unlist(lapply(1:length(predictors), function(i) combn(predictors, i, simplify = FALSE)), recursive = FALSE)

# RMSE function to use with abundance data
rmse <- function(actual, predicted) sqrt(mean((actual - predicted)^2))

# Cross fold validation for Gaussian GAM (log abundance)
cv_gam_rmse <- function(vars, data, folds = 10) {
  set.seed(123)
  data <- data[sample(nrow(data)), ]  # Shuffle
  data$Fold <- sample(rep(1:folds, length.out = nrow(data)))
  rmses <- numeric(folds) # Use 10 folds
  
  for (i in 1:folds) {
    train_fold <- data[data$Fold != i, ] # Loops through each fold and partition into training and test
    test_fold  <- data[data$Fold == i, ]
    
# Create GAM formula based on chosen predictors
    formula_obj <- as.formula(paste("LogAbundance ~", paste(paste0("s(`", vars, "`, k = 10)"), collapse = " + ")))
    
# Fit GAM on training fold
    model <- try(gam(formula_obj, data = train_fold, family = gaussian(), select = TRUE, method = "REML"), silent = TRUE)
    if (inherits(model, "try-error")) {
      rmses[i] <- NA # Report NA for a fold if fitting fails
    } else {
      
# Predict on test fold and compute RMSE
      preds <- predict(model, newdata = test_fold, type = "response")
      rmses[i] <- sqrt(mean((test_fold$LogAbundance - preds)^2))
    }
  }
  
# Return the mean RMSE across all folds, ignoring failed folds (NA)
  return(mean(rmses, na.rm = TRUE))
}

# Run cross fold on training data
cv_results <- data.frame(Predictors = sapply(predictor_combos, paste, collapse = " + "), RMSE = sapply(predictor_combos, function(vars) cv_gam_rmse(vars, train_data)))

# Rank by RMSE
cv_results <- cv_results[order(cv_results$RMSE), ]

# Fit best GAM
best_vars <- unlist(strsplit(cv_results$Predictors[1], " + ", fixed = TRUE))
best_formula <- as.formula(paste("LogAbundance ~", paste(paste0("s(", best_vars, ", k = 10)"), collapse = " + ")))
final_gam_sar86 <- gam(best_formula, data = train_data, family = gaussian(), select = TRUE, method = "REML") 

# Extract p-values from smooth terms and adjust for multiple comparison to assess significance
gam_summary <- summary(final_gam_sar86)
smooth_pvals <- gam_summary$s.table[, "p-value"]
adj_pvals <- p.adjust(smooth_pvals, method = "holm") # Holm correction

# Table of smooth results
smooth_results <- data.frame(
  Term = rownames(gam_summary$s.table),
  EDF = gam_summary$s.table[, "edf"],
  Stat = gam_summary$s.table[, if ("F" %in% colnames(gam_summary$s.table)) "F" else "Chi.sq"],
  RawP = smooth_pvals,
  AdjP = adj_pvals,
  Sig = adj_pvals < 0.05 # TRUE/FALSE significance
) 
print(smooth_results)

# Predict on train and test sets
train_data$Predicted <- predict(final_gam_sar86, newdata = train_data, type = "response")
test_data$Predicted  <- predict(final_gam_sar86, newdata = test_data, type = "response")

# Evaluate model performance with RMSE and spearman
train_rmse <- rmse(train_data$LogAbundance, train_data$Predicted)
test_rmse  <- rmse(test_data$LogAbundance, test_data$Predicted)
spearman_train <- cor(train_data$LogAbundance, train_data$Predicted, method = "spearman")
spearman_test  <- cor(test_data$LogAbundance, test_data$Predicted, method = "spearman")

# Estimate deviance explained for train and test sets
null_deviance_train <- sum((train_data$LogAbundance - mean(train_data$LogAbundance))^2)
resid_deviance_train <- sum((train_data$LogAbundance - train_data$Predicted)^2)
dev_expl_train <- 1 - resid_deviance_train / null_deviance_train

null_deviance_test <- sum((test_data$LogAbundance - mean(train_data$LogAbundance))^2)
resid_deviance_test <- sum((test_data$LogAbundance - test_data$Predicted)^2)
dev_expl_test <- 1 - resid_deviance_test / null_deviance_test

# Combine all test and training data for plotting
train_data$Set <- "Train"
test_data$Set <- "Test"
combined_data <- rbind(train_data, test_data)
combined_data$Fill <- ifelse(combined_data$Set == "Train", "#DDCC77", "gray70")
combined_data$Shape <- ifelse(combined_data$Set == "Train", 21, 21)
combined_data$Alpha <- 1
equation_text <- paste("Log(Abundance + 1) ~", paste(paste0("s(", best_vars, ")"), collapse = " + "))

# Plot observed vs. predicted for train and test for SAR86 - Figure 4D
ggplot(combined_data, aes(x = LogAbundance, y = Predicted)) +
  geom_point(aes(shape = factor(Shape), fill = Fill, alpha = Alpha), size = 5, stroke = 1, color = "black", show.legend = TRUE) +
  scale_shape_manual(values = c("21" = 21, "24" = 24)) + scale_fill_identity() + scale_alpha_identity() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") + annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.5,
  label = paste0("ρ (Train) = ", round(spearman_train, 2), "\nρ (Test) = ", round(spearman_test, 2)), size = 4.5) +
  annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -1.5, label = paste0("RMSE (Train) = ", round(train_rmse, 2), "\nRMSE (Test) = ", round(test_rmse, 2)), size = 4.5) + 
  labs(title = "Observed vs. Predicted Log Abundance (GAM)", subtitle = equation_text, x = "Observed log(Abundance + 1)", y = "Predicted log(Abundance + 1)", fill = "Data Set") +
  xlim(4.2, 6.7) + ylim(4.2, 6.7) + theme_bw()
#ggsave(filename = "SAR86_gam_final.pdf", plot = last_plot(), device = "pdf", path = NULL, scale = 1, width = 6, height = 5, dpi = 150)

# Repeat GAM workflow for Flavobacteriales
df_flavo_filt$LogAbundance <- log(df_flavo_filt$Abundance + 1) # Log transform rarefied abundance
log_vars <- c("temp", "salinity", "phosphate", "nh4", "nitrate") # Select predictor variables to log transform

for (var in log_vars) {
  new_var <- paste0("log_", var)
  df_flavo_filt[[new_var]] <- log(df_flavo_filt[[var]] + 0.01)
} # Loop over the variables and transform

set.seed(1234) # Set for reproducibility

# Define predictors again with new log transformed variables
predictors <- c("log_temp", "log_salinity", "oxygen", "log_phosphate", "log_nh4", "log_nitrate", "pH_corrected", "dic")

# Stratified split of data into train and test sets
df_flavo_filt$Strata <- interaction(df_flavo_filt$region, cut(df_flavo_filt$LogAbundance, breaks = 2)) # Cut data based on region and abundance
train_index <- suppressWarnings(createDataPartition(df_flavo_filt$Strata, p = 0.8, list = FALSE)) # Define train set as 80%
train_data <- df_flavo_filt[train_index, ] # Train set
test_data <- df_flavo_filt[-train_index, ] # Test set (20%)

# Create all non-empty predictor combinations
predictor_combos <- unlist(lapply(1:length(predictors), function(i) combn(predictors, i, simplify = FALSE)), recursive = FALSE)

# RMSE function to use with abundance data
rmse <- function(actual, predicted) sqrt(mean((actual - predicted)^2))

# Function for cross fold validation for Gaussian GAM (log abundance)
cv_gam_rmse <- function(vars, data, folds = 10) {
  set.seed(123)
  data <- data[sample(nrow(data)), ]  # Shuffle
  data$Fold <- sample(rep(1:folds, length.out = nrow(data)))
  rmses <- numeric(folds) # Use 10 folds
  
  for (i in 1:folds) {
    train_fold <- data[data$Fold != i, ] # Loops through each fold and partition into training and test
    test_fold  <- data[data$Fold == i, ]
    
# Create GAM formula based on chosen predictors
    formula_obj <- as.formula(paste("LogAbundance ~", paste(paste0("s(`", vars, "`, k = 10)"), collapse = " + ")))
    
# Fit GAM on training fold
    model <- try(gam(formula_obj, data = train_fold, family = gaussian(), select = TRUE, method = "REML"), silent = TRUE)
    if (inherits(model, "try-error")) {
      rmses[i] <- NA # Report NA for a fold if fitting fails
    } else {
      
# Predict on test fold and compute RMSE
      preds <- predict(model, newdata = test_fold, type = "response")
      rmses[i] <- sqrt(mean((test_fold$LogAbundance - preds)^2))
    }
  }
  
# Return the mean RMSE across all folds, ignoring failed folds (NA)
  return(mean(rmses, na.rm = TRUE))
}

# Run cross fold on training data
cv_results <- data.frame(Predictors = sapply(predictor_combos, paste, collapse = " + "), RMSE = sapply(predictor_combos, function(vars) cv_gam_rmse(vars, train_data)))

# Rank by RMSE
cv_results <- cv_results[order(cv_results$RMSE), ]

# Fit best GAM
best_vars <- unlist(strsplit(cv_results$Predictors[1], " + ", fixed = TRUE))
best_formula <- as.formula(paste("LogAbundance ~", paste(paste0("s(", best_vars, ", k = 10)"), collapse = " + ")))
final_gam_flavo <- gam(best_formula, data = train_data, family = gaussian(), select = TRUE, method = "REML") 

# Extract p-values from smooth terms and adjust for multiple comparison to assess significance
gam_summary <- summary(final_gam_flavo)
smooth_pvals <- gam_summary$s.table[, "p-value"]
adj_pvals <- p.adjust(smooth_pvals, method = "holm") # Holm correction

# Table of smooth results
smooth_results <- data.frame(
  Term = rownames(gam_summary$s.table),
  EDF = gam_summary$s.table[, "edf"],
  Stat = gam_summary$s.table[, if ("F" %in% colnames(gam_summary$s.table)) "F" else "Chi.sq"],
  RawP = smooth_pvals,
  AdjP = adj_pvals,
  Sig = adj_pvals < 0.05 # TRUE/FALSE significance
) 
print(smooth_results)

# Predict on train and test sets
train_data$Predicted <- predict(final_gam_flavo, newdata = train_data, type = "response")
test_data$Predicted  <- predict(final_gam_flavo, newdata = test_data, type = "response")

# Evaluate model performance with RMSE and spearman
train_rmse <- rmse(train_data$LogAbundance, train_data$Predicted)
test_rmse  <- rmse(test_data$LogAbundance, test_data$Predicted)
spearman_train <- cor(train_data$LogAbundance, train_data$Predicted, method = "spearman")
spearman_test  <- cor(test_data$LogAbundance, test_data$Predicted, method = "spearman")

# Estimate deviance explained for train and test sets
null_deviance_train <- sum((train_data$LogAbundance - mean(train_data$LogAbundance))^2)
resid_deviance_train <- sum((train_data$LogAbundance - train_data$Predicted)^2)
dev_expl_train <- 1 - resid_deviance_train / null_deviance_train

null_deviance_test <- sum((test_data$LogAbundance - mean(train_data$LogAbundance))^2)
resid_deviance_test <- sum((test_data$LogAbundance - test_data$Predicted)^2)
dev_expl_test <- 1 - resid_deviance_test / null_deviance_test

# Combine all test and training data for plotting
train_data$Set <- "Train"
test_data$Set <- "Test"
combined_data <- rbind(train_data, test_data)
combined_data$Fill <- ifelse(combined_data$Set == "Train", "#332288", "gray70")
combined_data$Shape <- ifelse(combined_data$Set == "Train", 21, 21)
combined_data$Alpha <- 1
equation_text <- paste("Log(Abundance + 1) ~", paste(paste0("s(", best_vars, ")"), collapse = " + "))

# Plot observed vs. predicted for train and test for Flavobacteriales - Figure 4E
ggplot(combined_data, aes(x = LogAbundance, y = Predicted)) +
  geom_point(aes(shape = factor(Shape), fill = Fill, alpha = Alpha), size = 5, stroke = 1, color = "black", show.legend = TRUE) +
  scale_shape_manual(values = c("21" = 21, "24" = 24)) + scale_fill_identity() + scale_alpha_identity() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") + annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.5,
  label = paste0("ρ (Train) = ", round(spearman_train, 2), "\nρ (Test) = ", round(spearman_test, 2)), size = 4.5) +
  annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -1.5, label = paste0("RMSE (Train) = ", round(train_rmse, 2), "\nRMSE (Test) = ", round(test_rmse, 2)), size = 4.5) + 
  labs(title = "Observed vs. Predicted Log Abundance (GAM)", subtitle = equation_text, x = "Observed log(Abundance + 1)", y = "Predicted log(Abundance + 1)", fill = "Data Set") +
  xlim(4, 7) + ylim(4, 7) + theme_bw()
#ggsave(filename = "Flavo_gam_final.pdf", plot = last_plot(), device = "pdf", path = NULL, scale = 1, width = 6, height = 5, dpi = 150)

## Loop through the above GAMs and plot residuals and QQ-plots for all groups in a single panel
final_gams <- list(SAR11 = final_gam_sar, Prochlorococcus = final_gam_pro, Synechococcus = final_gam_syne, SAR86 = final_gam_sar86, Flavobacteriales = final_gam_flavo) # Name models in a list 

# Function to calculate quantiles and plot residuals
make_diagnostic_plots <- function(model, group_name) {
  resid_vals <- residuals(model, type = "deviance")
  
# Calculate theoretical and sample quantiles
  qq_theory <- qqnorm(resid_vals, plot.it = FALSE)$x
  qq_sample <- qqnorm(resid_vals, plot.it = FALSE)$y
  
# Calculate slope and intercept from line through 1st & 3rd quartiles (like qqline)
  q_sample <- quantile(resid_vals, probs = c(0.25, 0.75), na.rm = TRUE)
  q_theory <- qnorm(c(0.25, 0.75))
  slope <- diff(q_sample) / diff(q_theory)
  intercept <- q_sample[1] - slope * q_theory[1]
  qq_data <- data.frame(theoretical = qq_theory, sample = qq_sample)
  
# Plot histogram of residuals
  hist_plot <- ggplot(data.frame(resid = resid_vals), aes(x = resid)) +
    geom_histogram(bins = 30, fill = "#004488", color = "black") + theme_minimal() +
    labs(title = paste(group_name, "- Residuals"), x = "Residuals", y = "Count")
  
# Plot quantiles
  qq_plot <- ggplot(qq_data, aes(x = theoretical, y = sample)) +
    geom_point(color = "#BB5566") + geom_abline(intercept = intercept, slope = slope, linetype = "dashed") +
    theme_minimal() + labs(title = paste(group_name, "- QQ Plot"), x = "Theoretical Quantiles", y = "Sample Quantiles")
  
  hist_plot + qq_plot # Combine histogram and q-q plots
}

# Generate a list of diagnostic plots
diagnostic_panels <- lapply(names(final_gams), function(name) {
  make_diagnostic_plots(final_gams[[name]], name)
})

# Combine all diagnostics into a single panel - Figure S3
wrap_plots(diagnostic_panels, ncol = 2)
#ggsave(filename = "GAM_residuals_16S.pdf", plot = last_plot(), device = "pdf", path = NULL, scale = 1, width = 10, height = 10, dpi = 150)

## Compare AIC values for all 16S groups with different GAM and GLM distributions
# Define data frames - make sure these have abundance, logAbundance, clrabundance, and log predictors
data_list <- list(SAR11 = df_sar_filt, Prochlorococcus = df_pro_filt, Synechococcus = df_syne_filt, SAR86 = df_sar86_filt, Flavobacteriales = df_flavo_filt)

# Define log-transformed predictors used in models
predictors <- c("log_temp", "log_salinity", "oxygen", "log_phosphate", "log_nh4", "log_nitrate", "pH_corrected", "dic")

# Initialize results table
aic_results <- data.frame(Group = character(), AIC_GLM_Poisson = numeric(), AIC_GLM_NB = numeric(), AIC_GAM_Gaussian_Log = numeric(), AIC_GAM_Gaussian_CLR = numeric(), AIC_GAM_NB = numeric(), stringsAsFactors = FALSE)

# Loop through each dataset and run different models
for (group_name in names(data_list)) {df <- data_list[[group_name]]
df <- df %>% drop_na(any_of(c("Abundance", "LogAbundance","CLRAbundance", predictors)))

# Model formulas
formula_raw <- as.formula(paste("Abundance ~", paste(predictors, collapse = " + ")))
formula_log <- as.formula(paste("LogAbundance ~", paste0("s(", predictors, ", bs='tp')", collapse = " + ")))
formula_clr <- as.formula(paste("CLRAbundance ~", paste0("s(", predictors, ", bs='tp')", collapse = " + ")))
formula_nb_gam <- as.formula(paste("Abundance ~", paste0("s(", predictors, ", bs='tp')", collapse = " + ")))

# Fit Models
glm_pois <- glm(formula_raw, data = df, family = poisson()) # GLM Poisson
glm_nb <- tryCatch(glm.nb(formula_raw, data = df), error = function(e) NULL) # GLM nb
gam_gauss <- gam(formula_log, data = df, family = gaussian()) # GAM guass
gam_gauss_clr <- gam(formula_clr, data = df, family = gaussian())
gam_nb <- gam(formula_nb_gam, data = df, family = nb()) # GAM nb

# Compute AICs for all models
aic_glm_pois <- AIC(glm_pois)
aic_glm_nb <- if (!is.null(glm_nb)) AIC(glm_nb) else NA
aic_gam_gauss <- AIC(gam_gauss)
aic_gam_gauss_clr <- AIC(gam_gauss_clr)
aic_gam_nb <- AIC(gam_nb)

# Combine results
aic_results <- rbind(aic_results, data.frame(Group = group_name, AIC_GLM_Poisson = round(aic_glm_pois, 2), AIC_GLM_NB = round(aic_glm_nb, 2), AIC_GAM_Gaussian_Log = round(aic_gam_gauss, 2), AIC_GAM_Gaussian_CLR = round(aic_gam_gauss_clr, 2), AIC_GAM_NB = round(aic_gam_nb, 2)))
}

# Export AIC results for all groups and models tested - Table S4
print(aic_results)
#write.csv(aic_results, "GAM_GLM_AIC_comparison_16S.csv", row.names = FALSE)

## Compare Gaussian GAMs with CLR vs. log transformed read counts for all 16S groups
# Function to build the plot of predicted vs. observed abundances - Figure S4
make_panel_plot <- function(df, obs, pred, r2, dev, title) {
  
# Equal axis limits to visualize 1:1 line
  rng <- range(c(df[[obs]], df[[pred]]), na.rm = TRUE)
  pad <- 0.05 * diff(rng)
  lims <- c(rng[1] - pad, rng[2] + pad)
  
ggplot(df, aes_string(x = obs, y = pred)) +
    geom_point(fill = "gray80", color = "black", shape = 21, size = 3, alpha = 0.8) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    coord_equal(xlim = lims, ylim = lims, expand = FALSE) + theme_bw(base_size = 12) +
    labs(title = title, subtitle = paste0("R² = ", sprintf("%.2f", r2),"  Dev = ", sprintf("%.1f", dev), "%"), x = "Observed", y = "Predicted")
}

# Specify groups
group_list <- list(SAR11 = df_sar_filt, SAR86 = df_sar86_filt, Flavo = df_flavo_filt, Pro = df_pro_filt, Syne = df_syne_filt)

# Specify predictor variables
predictors <- c("log_temp", "log_salinity", "oxygen", "log_phosphate", "log_nh4", "log_nitrate", "pH_corrected", "dic")

# Function to run paired log vs CLR GAMs for a given microbial group
run_gam_pair <- function(df, group_name, predictors) {
  
  df <- df %>% 
    dplyr::select(LogAbundance, CLRAbundance, all_of(predictors)) %>% 
    drop_na()
  
# Create formulas (each smooth has k = 10)
  smooth_terms <- paste0("s(", predictors, ", k = 10)")
  f_log <- as.formula(paste("LogAbundance ~", paste(smooth_terms, collapse = " + ")))
  f_clr <- as.formula(paste("CLRAbundance ~", paste(smooth_terms, collapse = " + ")))
  
# Fit the GAMs
  m_log <- gam(f_log, data = df, method = "REML", select = TRUE)
  m_clr <- gam(f_clr, data = df, method = "REML", select = TRUE)
  
# Predict
  df$Pred_Log <- predict(m_log, newdata = df)
  df$Pred_Clr <- predict(m_clr, newdata = df)
  
# Metrics
  r2_log  <- summary(m_log)$r.sq
  dev_log <- summary(m_log)$dev.expl * 100
  
  r2_clr  <- summary(m_clr)$r.sq
  dev_clr <- summary(m_clr)$dev.expl * 100
  
# Create log + CLR comparison panels
p1 <- make_panel_plot(df, "LogAbundance", "Pred_Log", r2_log, dev_log, paste0(group_name, ": Log GAM"))
  
p2 <- make_panel_plot(df, "CLRAbundance", "Pred_Clr", r2_clr, dev_clr, paste0(group_name, ": CLR GAM"))
p1 + p2
}

plot_list <- lapply(names(group_list), function(g) {
  run_gam_pair(group_list[[g]], g, predictors)
})

# Final figure panel for all groups
final_panel <- wrap_plots(plotlist = plot_list, ncol = 2)
final_panel
#ggsave(filename = "GAM_log_clr_16S.pdf", plot = last_plot(), device = "pdf", path = NULL, scale = 1, width = 10, height = 10, dpi = 150)

## Plot GAM partial effects for each 16S group - Figure 6
sar_pe <- draw(final_gam_sar, residuals = TRUE, resid_col = "black", ci_col = "#88CCEE", ci_alpha = 0.5, ncol = 4)
sar_pe & theme_bw()
#ggsave(filename = "Sar_partial_effects.pdf", plot = last_plot(), device = "pdf", path = NULL, scale = 1, width = 9, height = 5, dpi = 150)

pro_pe <- draw(final_gam_pro, residuals = TRUE, resid_col = "black", ci_col = "#CC6677", ci_alpha = 0.5, ncol = 4)
pro_pe & theme_bw()
#ggsave(filename = "Pro_partial_effects.pdf", plot = last_plot(), device = "pdf", path = NULL, scale = 1, width = 9, height = 5, dpi = 150)

syne_pe <- draw(final_gam_syne, residuals = TRUE, resid_col = "black", ci_col = "#732633", ci_alpha = 0.5, ncol = 4)
syne_pe & theme_bw()
#ggsave(filename = "Syne_partial_effects.pdf", plot = last_plot(), device = "pdf", path = NULL, scale = 1, width = 9, height = 5, dpi = 150)

sar86_pe <- draw(final_gam_sar86, residuals = TRUE, resid_col = "black", ci_col = "#DDCC77", ci_alpha = 0.5, ncol = 4)
sar86_pe & theme_bw()
#ggsave(filename = "Sar86_partial_effects.pdf", plot = last_plot(), device = "pdf", path = NULL, scale = 1, width = 9, height = 5, dpi = 150)

flavo_pe <- draw(final_gam_flavo, residuals = TRUE, resid_col = "black", ci_col = "#332288", ci_alpha = 0.5, ncol = 4)
flavo_pe & theme_bw()
#ggsave(filename = "Flavo_partial_effects.pdf", plot = last_plot(), device = "pdf", path = NULL, scale = 1, width = 9, height = 5, dpi = 150)

## Predict log abundance at all GOMECC-4 sites
# Models for 16S groups were used to predict log abundance at 135 GOMECC-4 sites at the surface (< 10 m)
# Data for 6 sites were not available at the surface. This resulted in 84 sites where DNA was not collected, allowing for an expanded view of microbial distribution in Gulf surface waters.
gomecc_all <- read.csv("/Users/seananderson/Documents/NOAA NGI/GOMECC-4/data-input/GOMECC-all-sites-final.csv", header=T, row.names = NULL, check.names=F,fileEncoding="UTF-8-BOM") # Load file of relevant metadata from all surface sites 
gomecc_all_pred <- gomecc_all # Rename

# Apply log(+ 0.01) transform to selected predictors
log_vars <- c("temp", "salinity", "phosphate", "nh4", "nitrate")
for (v in log_vars) {
  new_colname <- paste0("log_", v)
  gomecc_all_pred[[new_colname]] <- log(gomecc_all_pred[[v]] + 0.01)
}

gomecc_all_sar = gomecc_all_pred[,c(13,15,17,18,20:21)] # Subset to variables in final SAR11 model
new_val_sar = predict(final_gam_sar, newdata = gomecc_all_sar, type = "response") # Predict new values 
new_val_sar = as.data.frame(new_val_sar)
odv_sar = cbind(new_val_sar, gomecc_all_pred[4:6])

gomecc_all_sar86 = gomecc_all_pred[,c(9,13,15,17,19)] # Subset to variables in final SAR86 model
new_val_sar86 = predict(final_gam_sar86, newdata = gomecc_all_sar86, type = "response") # Predict new values
new_val_sar86 = as.data.frame(new_val_sar86)

gomecc_all_flavo = gomecc_all_pred[,c(9,13,17:18,20:21)] # Subset to variables in final Flavobacteriales model
new_val_flavo = predict(final_gam_flavo, newdata = gomecc_all_flavo, type = "response") # Predict new values
new_val_flavo = as.data.frame(new_val_flavo)

gomecc_all_pro = gomecc_all_pred[,c(9,15,17:19)] # Subset to variables in final Prochlorococcus model
new_val_pro = predict(final_gam_pro, newdata = gomecc_all_pro, type = "response") # Predict new values 
new_val_pro = as.data.frame(new_val_pro)

gomecc_all_syne = gomecc_all_pred[,c(13,17:19,21)] # Subset to variables in final Synechococcus model
new_val_syne = predict(final_gam_syne, newdata = gomecc_all_syne, type = "response") # Predict new values
new_val_syne = as.data.frame(new_val_syne)

# Combine new predicted log abundance values into one data frame. These file with predicted values will be uploaded and analyzed in ODV (ODV files provided on GitHub).
all_odv = cbind(odv_sar, new_val_sar86, new_val_flavo, new_val_pro, new_val_syne)
#write.csv(all_odv, "ODV_microbes_abund_16S.csv", ow.names = T) # Write .csv file of log abundance
#write.csv(gomecc_all_pred, "ODV_factors_16S.csv", row.names = T) # Write .csv file of log factors

## Species indicator analysis
# Explore microbes at the ASV level that are indicative of different acidic conditions in the GOM at this time. As with the models, focus on samples from the photic zone (Cluster 1). 
# We used the TA:DIC (total alkalinity:dissolved inorganic carbon) ratio as it is a good proxy of OA, indicating the buffering capacity of seawater 
# A lower TA:DIC = less buffered waters (more acidic), while a higher TA:DIC = more buffered waters

# Run the core indicator analysis program using 16S ASVs
ps_dsq <- subset_samples(ps_subset_16S, cluster_16S == "Cluster 1") # Subset to Cluster 1
ps_dsq <- subset_samples(ps_dsq, alk_dic_oa == "High" | alk_dic_oa =="Low") # Split samples based on high or low TA:DIC (see above)
ASV = otu_table(ps_dsq) # Subset ASVs
ASV =as.data.frame(ASV) # Convert to data frame
if(taxa_are_rows(ps_dsq)){ASV <- t(ASV)} # Flip table for analysis
ASV = as.data.frame(ASV) # Convert to df again

data <- as(sample_data(ps_dsq), "data.frame") # Subset metadata
data = data[!is.na(data$alk_dic_oa),] # Remove TA:DIC NA values 
bin = data$alk_dic_oa # Bin based on low vs. high categories
inv = multipatt(ASV, bin, func = "r.g", control = how(nperm = 999)) # Core function using 999 permutations

# Correct p-values for multiple comparisons
indisp.sign <- data.table::as.data.table(inv$sign, keep.rownames = TRUE) # Extract table of stats
indisp.sign[ ,p.value.bh:= p.adjust(p.value, method = "BH")] # Add adjusted p-value with Benjamini-Hochberg
indisp.sign2 = indisp.sign[p.value <= 0.05,] # Filter to include only p-values < 0.05; we will have adjusted values as well
#write.csv(indisp.sign2, "16S_indic_corr2.csv", row.names = T) # Write .csv file

# Results of indicator analysis (ASV, stat, and corrected p-value) combined in a single .csv file
# Upload this file into R (provided on GitHub)
# Match indicator values table for each ASV with the relative abundance of the same ASV in the photic zone via the ASV count table
# Data frames need to match

# Upload table with all indicator metadata for high and low TA:DIC categories
indic_OA <- read.csv("/Users/seananderson/Documents/NOAA NGI/GOMECC-4/data-input/16S/OA_indic_16S_new.csv", header = T, row.names = NULL, check.names = F,fileEncoding = "UTF-8-BOM") # Indicator ASVs file
colnames(indic_OA)[1]  <- "OTU" # Change column names
indic.high = indic_OA[indic_OA[["OA"]] == "High", ] # Subset to new data frame with only high TA:DIC
indic.low = indic_OA[indic_OA[["OA"]] == "Low", ] # Subset to new df with only low TA:DIC

ps4 <-  transform_sample_counts(ps_dsq, function(x) 100*x / sum(x) ) # Transform to relative abundance
x1 = speedyseq::psmelt(ps4) # Melt the data
OA.1 = x1 %>% # Group to the mean relative abundance at the ASV level
  dplyr::group_by(OTU, alk_dic_oa, Phylum, Genus, Species) %>%
  dplyr::summarise_at(.vars = c("Abundance"), .funs = mean)

OA.high = OA.1[OA.1[["alk_dic_oa"]] == "High", ] # Subset ASV relative abundance based on high TA:DIC
OA.low = OA.1[OA.1[["alk_dic_oa"]] == "Low", ] # Subset ASV based on low TA:DIC

df_h <- OA.high[OA.high$OTU %in% indic.high$OTU,] # Keep ASVs if they match indicator dataset
df_L <- OA.low[OA.low$OTU %in% indic.low$OTU,] # Keep ASVs if they match indicator dataset

#write.csv(df_h, "ASVs_16S_abund_high_new.csv",row.names=T) # Write .csv file
#write.csv(df_L, "ASVs_16S_abund_low_new.csv",row.names=T) # Write .csv file

## No prokaryote ASVs were significant indicators after correcting p-values - no further analysis or indicator plots

## Estimate spearman values of variables for Clusters 2-3 
# Spearman correlation matrices for all clusters and variables
cluster2 <- x1[x1[["cluster_16S"]] == "Cluster 2", ] # Subset to Cluster 2
clus.2 = cluster2 %>% # Group the data based on class and preserve sample replication for the models
  dplyr::group_by(station, depth_category, Order, replicate) %>%
  summarise(temp = mean(temp),
            salinity = mean(salinity),
            oxygen = mean(oxygen),
            phosphate = mean(phosphate),
            nitrate = mean(nitrate),
            nitrite = mean(nitrite),
            silicate = mean(silicate),
            nh4 = mean(nh4),
            pH_corrected = mean(pH_corrected),
            total_alkalinity = mean(total_alkalinity),
            OMEGA_AR = mean(OMEGA_AR),
            dic = mean(dic),
            pCO2_corrected = mean(pCO2_corrected),
            carbonate = mean(carbonate),
            fluorescence = mean(fluorescence),
            Abundance = mean(Abundance),
            .groups = 'drop')

cluster3 <- x1[x1[["cluster_16S"]] == "Cluster 3", ] # Subset to Cluster 3
clus.3 = cluster3 %>% # Group the data based on class and preserve sample replication for the models
  dplyr::group_by(station, depth_category, Order, replicate) %>%
  summarise(temp = mean(temp),
            salinity = mean(salinity),
            oxygen = mean(oxygen),
            phosphate = mean(phosphate),
            nitrate = mean(nitrate),
            nitrite = mean(nitrite),
            silicate = mean(silicate),
            nh4 = mean(nh4),
            pH_corrected = mean(pH_corrected),
            total_alkalinity = mean(total_alkalinity),
            OMEGA_AR = mean(OMEGA_AR),
            dic = mean(dic),
            pCO2_corrected = mean(pCO2_corrected),
            carbonate = mean(carbonate),
            fluorescence = mean(fluorescence),
            Abundance = mean(Abundance),
            .groups = 'drop')

# Spearman Cluster 2
df2 = clus.2[,c(5:19)] # Subset to include environmental variables
df2 = na.omit(df2) # Omit any rows that have NA values
correlations2 = cor(df2, method = "spearman") # Perform Spearman correlations to assess collinearity
#write.csv(correlations2, "Cluster2_corr_16S.csv", row.names = T) # Write .csv file - Table S3

# Spearman Cluster 3
df3 = clus.3[,c(5:19)] # Subset to include environmental variables
df3 = na.omit(df3) # Omit any rows that have NA values
correlations3 = cor(df3, method = "spearman") # Perform Spearman correlations to assess collinearity
#write.csv(correlations3, "Cluster3_corr_16S.csv", row.names = T) # Write .csv file - Table S3

## GOM map with points filled based on 16S cluster assignments at each depth - Figure S7
b = getNOAA.bathy(lon1 = -100, lon2 = -79, lat1 = 30, lat2 = 22, resolution = 3, keep = TRUE) # Import bathymetry data
sites_16S = read.csv(file="/Users/seananderson/Documents/NOAA NGI/GOMECC-4/data-input/18S/map_clus_16S.csv", header = T, row.names = NULL, check.names = F, fileEncoding = "UTF-8-BOM")
sites_16S$Latitude <- as.character(sites_16S$Latitude) # Format latitude
sites_16S$Longitude <- as.character(sites_16S$Longitude) # Format longitude
sites_16S$Depth <- as.character(sites_16S$Depth) # Format depth zones
sites_16S$Cluster <- as.factor(sites_16S$Cluster) # Format cluster assignments

# Partition based on depth and set map colors
sites_16S_surf <- sites_16S[sites_16S[["Depth"]] == "Surface", ] # Subset to surface
sites_16S_dcm <- sites_16S[sites_16S[["Depth"]] == "DCM", ] # Subset to DCM
sites_16S_bottom <- sites_16S[sites_16S[["Depth"]] == "Deep", ] # Subset to near bottom
blues <- c("lightsteelblue4", "lightsteelblue3", "lightsteelblue2", "lightsteelblue1")
greys <- c(grey(0.7), grey(0.9), grey(0.95))
ap <- c("#BB5566" ,"#DDAA33","#004488") # Color of clusters

# Build the map
plot(b, image = TRUE, land = TRUE, n=1, bpal = list(c(0, max(b), greys), c(min(b), 0, blues))) # Bathymetry colors
plot(b, lwd = 1, deep = 0, shallow = 0, step = 0, add = TRUE) # Highlight coastline with solid black line
plot(b, deep= -50, shallow= -50, lwd = 0.4, drawlabels = FALSE, add = TRUE) # Add -50m isobath
plot(b, deep= -200, shallow= -200, lwd = 0.4, drawlabels = FALSE, add = TRUE) # Add -200m isobath
plot(b, deep= -1000, shallow= -1000, lwd = 0.4, drawlabels = FALSE, add = TRUE) # Add -1000m isobath
plot(b, deep= -3000, shallow= -3000, lwd = 0.4, drawlabels = FALSE, add = TRUE) # Add -3000m isobath

coord <- par("usr") # Store coordinate system for plot
b2 <- getNOAA.bathy(-100, -79, 36, 18, res = 3.5) # Re-run import
#pdf("map_sites_surf_16S.pdf") # Name pdf file to store the surface map - store and run through each one at a time
#pdf("map_sites_dcm_16S.pdf") # Name pdf file to store DCM map
#pdf("map_sites_bottom_16S.pdf") # Name pdf file to store bottom map

plot(b2, image = TRUE, land = TRUE, n = 1, bpal = list(c(0, max(b2), greys), c(min(b2), 0, blues))) # Replot map colors
plot(b2, lwd = 1, deep = 0, shallow = 0, step = 0, add = TRUE) # Highlight coastline with solid black line
plot(b2, deep= -50, shallow= -50, lwd = 0.4, drawlabels = FALSE, add = TRUE) # Add -50m isobath
plot(b2, deep= -200, shallow= -200, lwd = 0.4, drawlabels = FALSE, add = TRUE) # Add -200m isobath
plot(b2, deep= -1000, shallow= -1000, lwd = 0.4, drawlabels = FALSE, add = TRUE) # Add -1000m isobath
plot(b2, deep= -3000, shallow= -3000, lwd = 0.4, drawlabels = FALSE, add = TRUE) # Add -3000m isobath

points(sites_16S_surf$Longitude, sites_16S_surf$Latitude, pch = 21, col = "black", cex = 3, asp = 1, bg = ap[sites_16S_surf$Cluster]) # Add sampling points
points(sites_16S_dcm$Longitude, sites_16S_dcm$Latitude, pch = 21, col = "black", cex = 3, asp = 1, bg = ap[sites_16S_dcm$Cluster]) # Add sampling points
points(sites_16S_bottom$Longitude, sites_16S_bottom$Latitude, pch = 21, col = "black", cex = 3, asp = 1, bg = ap[sites_16S_bottom$Cluster]) # Add sampling points
scaleBathy(b2, deg = 1, x = "bottomleft", inset = 10, y = NULL, angle = 90)
#dev.off() # Print pdf one at a time (see above)

# R Session information
sessionInfo()