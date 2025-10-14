## Code for processing 18S amplicon data from the GOMECC-4 cruise
# Updated 10/14/2025
# Sean R. Anderson

# Load in packages
library(tidyverse);library(vegan);library(qiime2R)
library(phyloseq);library(reshape2);library(corrplot)
library(fantaxtic);library(RColorBrewer);library(microbiome)
library(factoextra);library(microeco);library(data.table)
library(file2meco);library(ggpubr);library(treemap);library(geosphere)
library(rcartocolor);library(indicspecies);library(performance);library(speedyseq)
library(dplyr);library(sjstats);library(lmtest);library(ranacapa);library(caret)
library(mgcv);library(gratia);library(scales);library(patchwork);library(marmap)

# Load 18S count file
table <- read_qza(file = "/Users/seananderson/Documents/NOAA NGI/GOMECC-4/data-input/18S/table-18S.qza")
count_tab <- table$data %>% as.data.frame() # Convert to data frame 
#write.csv(count_tab, "Count18S.csv",row.names=T) # Write csv file

# Load 18S taxonomy file
taxonomy <- read_qza(file = "/Users/seananderson/Documents/NOAA NGI/GOMECC-4/data-input/18S/taxonomy-18S.qza")
tax_tab <- taxonomy$data %>% # Convert to data frame, tab separate and rename taxa levels, and remove row with confidence values
  as.data.frame() %>%
  separate(Taxon, sep = ";", c("Domain","Supergroup","Division","Subdivision", "Class","Order","Family", "Genus", "Species")) %>% 
  column_to_rownames("Feature.ID") %>%
  dplyr::select(-Confidence)
#write.csv(tax_tab, "Taxonomy18S.csv",row.names=T) # Write csv file

# From the prior exported count and taxonomy .csv files, manually add a new column for 18S functional groups
# Re-upload these new .csv files for downstream processing. All files provided on GitHub. 
tax_new = read.csv(file = "/Users/seananderson/Documents/NOAA NGI/GOMECC-4/data-input/18S/Taxonomy18S_new.csv.gz", header = T, row.names = NULL, check.names = F,fileEncoding = "UTF-8-BOM") # File is compressed to reduce size
count_new = read.csv(file = "/Users/seananderson/Documents/NOAA NGI/GOMECC-4/data-input/18S/Count18S_new.csv.gz", header = T, row.names = NULL, check.names = F,fileEncoding = "UTF-8-BOM") # File is compressed to reduce size
count_new = count_new[ order(match(count_new$ASV, tax_new$ASV)), ] # Match order of ASVs for both files 
rownames(count_new) <- count_new$ASV # Rename row names
count_new <- count_new[ -c(1) ]
rownames(tax_new) <- tax_new$ASV # Rename row names to match count file
tax_new <- tax_new[ -c(1) ]

# Load metadata
sample_info_tab <- read.csv("/Users/seananderson/Documents/NOAA NGI/GOMECC-4/data-input/metadata_Aug2023.csv", header=T, row.names = NULL, check.names = F, fileEncoding = "UTF-8-BOM") 
row.names(sample_info_tab) <- sample_info_tab[,1]
sample_info_tab <- sample_info_tab[,-1]
sample_info_tab$dic = as.numeric(sample_info_tab$dic)
sample_info_tab$carbonate = as.numeric(sample_info_tab$carbonate)

# Create phyloseq object
ps <- phyloseq(tax_table(as.matrix(tax_new)), otu_table(count_new, taxa_are_rows = T), sample_data(sample_info_tab)) # Create phyloseq object
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps))) # Rename ASVs to be sequential

# Remove any unwanted groups
ps_new = subset_taxa(ps, Class != "Craniata" |is.na(Class)) # The following are groups we want to remove
ps_new = subset_taxa(ps, Class != "Basidiomycota" |is.na(Class))
ps_new = subset_taxa(ps, Class != "Ascomycota" |is.na(Class))
ps_new = subset_taxa(ps_new, Division !="Streptophyta" |is.na(Division))
ps_new = subset_taxa(ps_new, Domain !="Bacteria" |is.na(Domain))
ps_new = subset_taxa(ps_new, Division !="Rhodophyta" |is.na(Division))
ps_new = subset_taxa(ps_new, Family !="Insecta" |is.na(Family))
ps_new = subset_taxa(ps_new, Family !="Archosauria" |is.na(Family))
ps_new <- subset_taxa(ps_new, Subdivision !="Unassigned", Prune = T)

ps_new = name_na_taxa(ps_new) # Add an "unassigned" label to lowest annotation

# Subset to include only protists
ps_prot = subset_taxa(ps_new, Subdivision !="Metazoa", Prune = T)
ps_meta = subset_taxa(ps_new, Subdivision =="Metazoa", Prune = T) # Subset out any remaining metazoans

# Remove controls for now
ps_sub <- subset_samples(ps_prot, sample_type == "seawater") # Remove controls  
ps_sub = prune_samples(sample_sums(ps_sub) >=3000, ps_sub) # Remove samples with very low sequence read numbers

# Remove singletons (ASVs observed once across the dataset)
ps_filt = filter_taxa(ps_sub, function (x) {sum(x) > 1}, prune = TRUE)

# Estimate number of reads
ps_min <- min(sample_sums(ps_filt))
ps_mean <- mean(sample_sums(ps_filt))
ps_max <- max(sample_sums(ps_filt))

# Set color palettes used for certain figures
nb.cols <- 17
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)
group = carto_pal(12, "Safe")

# Plot rarefaction curves - Figure S1
rare_18S <- suppressWarnings(ggrare(ps_filt, step = 100, plot = FALSE, parallel = FALSE, se = FALSE))
rare_18S$data$depth_category <- factor(rare_18S$data$depth_category, levels = c("Surface","DCM","Deep"))
rare_18S + theme(legend.position = "none") + theme_bw() + theme(legend.position = "right") + facet_wrap(~depth_category + distance,scales = "free_y")
#ggsave(filename = "18S_rare_v2.eps", plot = last_plot(), device = "eps", path = NULL, scale = 1, width = 8, height = 5, dpi = 150)

# Rarefy to even sampling depth
ps_rare <- rarefy_even_depth(ps_filt, sample.size = min(sample_sums(ps_filt)), rngseed = 714, replace = TRUE, trimOTUs = TRUE, verbose = FALSE)

# Convert to Aitchison distance and prepare for clustering
ps_clr <- microbiome::transform(ps_rare, "clr") # Center log transform rarefied data
euc = phyloseq::distance(ps_clr, method = "euclidean") # Calculate Aitchison distance
euc.table <- as.matrix(dist(euc)) # Convert to matrix

# Perform hierarchical clustering and observe clusters - Figure S4
spe.ward <- hclust(euc, method = "ward.D2") # Hierarchical clustering with Ward's method
fviz_nbclust(euc.table, factoextra::hcut, method = "silhouette") # Three clusters is optimal
#ggsave(filename = "18S_cluster_silhouette_v1.eps", plot = last_plot(), device = "eps", path = NULL, scale = 1, width = 5, height = 4, dpi = 150) 

# Split the samples into clusters
sub_grp <- cutree(spe.ward, k = 3) # Cut the data based on our clusters; cluster assignments already in metadata file
sub_grp = as.data.frame(sub_grp)
table(sub_grp)

# Remove three samples that were unexpectedly grouped
ps_subset = subset_samples(ps_rare, sample_names(ps_rare) != "GOMECC4_CAPECORAL_Sta140_Deep_C" & sample_names(ps_rare) != "GOMECC4_LA_Sta38_Deep_C"  & sample_names(ps_rare) != "GOMECC4_FLSTRAITS_Sta123_Surface_B")

# Plot environmental factors across the three clusters - Figure S6
dataset <- phyloseq2meco(ps_subset) # Create microeco object
dataset$sample_table$cluster_18S <- as.factor(dataset$sample_table$cluster_18S) # Convert cluster to factor
t1 <- trans_env$new(dataset = dataset, env_cols = c(2,30:31,36:40,42:44,47:50,75)) # Subset to environmental factors of interest
t1$cal_diff(method = "wilcox", group = "cluster_18S") # Calculate statistical difference between clusters for each factor
head(t1$res_diff)
valid_measures <- intersect(unique(t1$res_diff$Measure), colnames(t1$data_env)) # Place all the plots into a list
tmp <- list() # Loop through each factor and plot
for(i in valid_measures){
  tmp[[i]] <- t1$plot_diff(measure = i, add_sig_text_size = 5, xtext_size = 12, add = "jitter", color_values = c("#BB5566","#DDAA33","#004488")) +
  theme(plot.margin = unit(c(0.1, 0, 0, 1), "cm")) + theme_bw() + theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black")) 
}
p1 = gridExtra::grid.arrange(grobs = tmp, ncol = 3) # Arrange panel
#ggsave(filename = "18S_factors_clusters.pdf", plot = p1, device = "pdf", path = NULL, scale = 1, width =6, height = 10, dpi = 150)

# Create new a distance matrix after removing three samples
ps_clr2 <- microbiome::transform(ps_subset, "clr") 
euc2 = phyloseq::distance(ps_clr2, method = "euclidean") # Aitchison distances

# Run PERMANOVA to test for significance between treatments
metadata <- as(sample_data(ps_subset), "data.frame") # Subset metadata
metadata$distance = as.factor(metadata$distance) # Convert to factor
metadata$region = as.factor(metadata$region) # Convert to factor
metadata$depth_category = as.factor(metadata$depth_category) # Convert to factor
adonis2(phyloseq::distance(ps_subset, method = "euclidean")~distance, data = metadata, p.adjust.m = 'holm', perm = 9999) # Run with 9999 permutations
adonis2(phyloseq::distance(ps_subset, method = "euclidean")~region, data = metadata,  p.adjust.m = 'holm', perm = 9999) # Repeat with transect
adonis2(phyloseq::distance(ps_subset, method = "euclidean")~depth_category, data = metadata,  p.adjust.m = 'holm', perm = 9999) # Repeat with depth

# Plot PCoA and color samples by cluster - Figure 3A
ordu = ordinate(ps_subset, "PCoA", distance = euc2) # Ordination
p = plot_ordination(ps_subset, ordu, color = "cluster_18S")
p$data$depth_meters <- as.factor(p$data$depth_meters) # Convert sampling depth to factor
p + theme_bw() + scale_fill_manual(values = c("#BB5566","#DDAA33","#004488")) +
  geom_point(aes(fill = cluster_18S), size = 5, shape = 21, colour = "black")  + 
  theme(text = element_text(size = 14)) 
#ggsave(filename = "18S_pcoa_final.pdf", plot = last_plot(), device = "pdf", path = NULL, scale = 1, width = 7, height = 5, dpi = 150)

# Estimate diversity metrics
rich_18S <- estimate_richness(ps_subset, measures = c("Observed", "Shannon")) # Estimate richness
clus = sample_data(ps_subset)$cluster_18S # Subset out cluster
region = sample_data(ps_subset)$region # Subset out region
depth = sample_data(ps_subset)$depth_meters # Subset out absolute depth
rich_all <- data.frame(rich_18S, region, depth, clus) # Merge richness estimates with cluster, region, and depth
rich_all$depth = as.character(rich_all$depth)
df2 = melt(rich_all) # Melt the data format for plotting
levels(df2$variable)[match("Observed", levels(df2$variable))] <- "# of 18S ASVs" # Change label
levels(df2$variable)[match("Shannon", levels(df2$variable))] <- "Shannon diversity index" # Change label

# Plot richness and diversity in each cluster - Figure 3B
p <- ggplot(df2, aes(x = factor(clus), y = value, fill = clus))
p$data$clus <- factor(p$data$clus, levels = c("Cluster 3","Cluster 2", "Cluster 1"))
p + geom_boxplot(alpha = 1, outlier.shape = NA, color = "black") + theme_bw() + 
  theme(text = element_text(size = 14)) + ylab("Diversity values") + theme(legend.position = "right") + scale_fill_manual(values = c("#004488","#DDAA33","#BB5566"))+
  geom_point(aes(fill = clus), size = 5, shape = 21, alpha = 1,colour = "black", position = position_jitterdodge(jitter.width = 0.8)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black")) + theme(axis.title.x = element_blank()) +
  facet_wrap(~variable, scales = "free_x",nrow = 2) + coord_flip() +
  geom_pwc(aes(group = clus), tip.length = 0, method = "wilcox_test", label = "p.adj", p.adjust.method = "holm")
#ggsave(filename = "Diversity_18S_new.pdf", plot = last_plot(), device = "pdf", path = NULL, scale = 1, width = 6, height = 7, dpi = 150) 

# Plot diversity metrics with respect to transect - Figure S11
p <- ggplot(df2, aes(x = factor(region), y = value, fill = region))
p$data$region <- factor(p$data$region, levels = c("27N", "FLSTRAITS","CAPECORAL", "TAMPA", "PANAMACITY", "PENSACOLA", "LA", "GALVESTON", "PAISNP", "BROWNSVILLE", "TAMPICO", "VERACRUZ", "CAMPECHE", "MERIDA", "YUCATAN", "CATOCHE", "CANCUN")) # Set the transect order
p + geom_boxplot(alpha = 0.5, outlier.shape = NA, color = "black") + theme_bw() +  geom_smooth(method = "loess", se = TRUE, color = "black", aes(group = 1)) + 
  theme(text = element_text(size = 14)) + ylab("Diversity values") + theme(legend.position = "right") + 
  scale_fill_manual(values = mycolors)+ geom_point(aes(fill = region), size = 3, shape = 21, colour = "black", position = position_jitterdodge(), show.legend = FALSE) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black")) + theme(axis.title.x = element_blank()) + facet_wrap(clus ~ variable, scales = "free_y")
#ggsave(filename = "Diversity_18S_region_v1.pdf", plot = last_plot(), device = "pdf", path = NULL, scale = 1, width = 20, height =6, dpi = 150) 

# Plot diversity metrics with respect to absolute depth - Figure S9
df2$depth <- as.numeric(as.character(df2$depth)) # Convert to numeric for correct sorting
df2$depth <- factor(df2$depth, levels = rev(sort(unique(df2$depth)))) # Re-convert to factor with levels reversed (so smallest at top after flip)
p <- ggplot(df2, aes(x = depth, y = value, fill = clus)) + geom_point(aes(fill = clus),size = 3, shape = 21,colour = "black")  + 
scale_fill_manual(values = c("#BB5566","#DDAA33","#004488")) + geom_smooth(method = "loess", se = TRUE, color = "black", aes(group = 1)) + scale_x_discrete(labels = label_wrap(0.1)) +
coord_flip() + facet_wrap(~variable, scales = "free_x") + guides(fill = FALSE) + theme_bw()
p
#ggsave(filename = "18S_div_depth.pdf", plot = last_plot(), device = "pdf", path = NULL, scale = 1, width = 6, height = 8, dpi = 150)

# Plot relative abundance stacked bar charts - Figure 3C
tax_table(ps_subset) <- tax_table(ps_subset)[,2:10] # Subset out the functional category column and focus on taxonomy levels assigned via PR2
barplot <- ps_subset %>%
  tax_glom(taxrank = "Class", NArm=FALSE) %>% # Agglomerate to class level
  transform_sample_counts(function(OTU) 100* OTU/sum(OTU)) %>% # Transform to relative abundance
  psmelt() %>% # Melt data
  group_by(region, Class, cluster_18S) %>% # Group by cluster and transect
  summarise_at("Abundance", .funs = mean) # Summarize at the mean

focus <- c("Syndiniales", "Dinophyceae", "Polycystinea", "Diplonemea", "Prymnesiophyceae","Mamiellophyceae","Sagenista","Bacillariophyceae","Opalozoa","Mediophyceae","Acantharea","RAD-B") # Focus on top 12 groups
barplot$Class <- ifelse(barplot$Class %in% focus, barplot$Class, "Others") # Others category
barplot_18S = barplot # Rename data frame
barplot_18S$Class<- as.character(barplot_18S$Class) # Convert to character

p <- ggplot(data = barplot_18S, aes(x = region, y = Abundance, fill = Class))
p$data$Class <- factor(p$data$Class, levels = c("Others","Bacillariophyceae","Opalozoa","Mediophyceae","Acantharea","Sagenista","Prymnesiophyceae","RAD-B", "Mamiellophyceae","Diplonemea","Polycystinea", "Dinophyceae","Syndiniales" )) # Set order of groups in the plot
p$data$region <- factor(p$data$region, levels = c("27N", "FLSTRAITS","CAPECORAL", "TAMPA", "PANAMACITY", "PENSACOLA", "LA", "GALVESTON", "PAISNP", "BROWNSVILLE", "TAMPICO", "VERACRUZ", "CAMPECHE", "MERIDA", "YUCATAN", "CATOCHE", "CANCUN")) # Set order of transects
p + geom_bar(aes(), stat="identity", position = "fill", width = 0.9)+
  scale_y_continuous(expand = c(0, 0))+  geom_hline(yintercept = 0) + theme_bw() + 
  scale_fill_manual(values = rev(c("#88CCEE", "#CC6677", "#44AA99", "#882255", "#999933", "#117733", "#DDCC77", "#332288", "#AA4499", "#F5793A", "#F7CDA4", "#A5CFCC", "#757575"))) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + theme(legend.position ="right") + theme(text = element_text(size = 12)) +
  guides(fill = guide_legend(nrow = 14, ncol = 1)) + facet_wrap(~cluster_18S, scales = "free_x") + labs(y = "Relative abundance (%)")
#ggsave(filename = "18S_stacked_class.pdf", plot = last_plot(), device = "pdf", path = NULL, scale = 1, width = 14, height = 5, dpi = 150)

# Prepare for taxonomy tree maps - Figure S13
ps <- tax_glom(ps_rare, "Genus", NArm = FALSE) # Agglomerate to genus level
x1 = speedyseq::psmelt(ps) # Melt the data
focus <- c("Alveolata", "Discoba", "Chlorophyta", "Rhizaria", "Stramenopiles", "Haptophyta") # Focus on top division level groups
x1$Division <- ifelse(x1$Division %in% focus, x1$Division, "Others") # Others category for plotting
x1$Division <- factor(x1$Division, levels = c("Alveolata", "Chlorophyta", "Discoba","Haptophyta","Rhizaria","Stramenopiles","Others")) # Set order for plotting
cluster1 <- x1[x1[["cluster_18S"]] == "Cluster 1", ] # Subset to Cluster 1
cluster2 <- x1[x1[["cluster_18S"]] == "Cluster 2", ] # Subset to Cluster 2
cluster3 <- x1[x1[["cluster_18S"]] == "Cluster 3", ] # Subset to Cluster 3

# Plot tree map for Cluster 1 (photic zone) 
#pdf("Cluster1_treemap_18S_v2.pdf", width = 12, height = 4)
treemap(dtf = cluster1,
        title = "Cluster 1 (2-52 m)", 
        algorithm = "pivotSize", border.lwds = c(2,0.5,0.1),
        border.col = c("black", "black", "black"),
        mapping = c(0,0,0),
        index = c("Class", "Genus"),
        vSize = "Abundance",
        vColor = "Division",
        palette = "Dark2",
        type = "categorical",
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
#pdf("Cluster2_treemap_18S_v2.pdf", width = 12, height = 4)
treemap(dtf = cluster2,
        title = "Cluster 2 (2-124 m)", 
        algorithm = "pivotSize", border.lwds = c(2,0.5,0.1),
        border.col = c("black", "black", "black"),
        mapping = c(0,0,0),
        index = c("Class", "Genus"),
        vSize = "Abundance",
        vColor = "Division",
        palette = "Dark2",
        type = "categorical",
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

# Plot tree map for Cluster 3 (Aphotic zone)
#pdf("Cluster3_treemap_18S_v2.pdf", width = 12, height = 4)
treemap(dtf = cluster3,
        title = "Cluster 3 (135-3326 m)", 
        algorithm = "pivotSize", border.lwds = c(2,0.5,0.1),
        border.col = c("black", "black", "black"),
        mapping = c(0,0,0),
        index = c("Class", "Genus"),
        vSize = "Abundance",
        vColor = "Division",
        palette = "Dark2",
        type = "categorical",
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

# Format data for generalized additive models (GAMs)
class_18S <- tax_glom(ps_subset, taxrank = "Class", NArm = FALSE) # Aggregate to the class level
x1 = psmelt(class_18S) # Melt the data
cluster1 <- x1[x1[["cluster_18S"]] == "Cluster 1", ] # Subset to Cluster 1
glm.1 = cluster1 %>% # Group the data based on class and preserve sample replication for the models
  dplyr::group_by(station, depth_category, Class, replicate, region) %>%
  dplyr::summarise(temp = mean(temp),
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

# Select variables for initial models based on low collinearity
df1 = glm.1[,c(6:20)] # Subset to include environmental variables
df1 = na.omit(df1) # Omit any rows that have NA values
correlations = cor(df1, method = "spearman") # Perform Spearman correlations to assess collinearity
#write.csv(correlations, "Cluster1_corr_18S.csv",row.names=T) # Write .csv file - Table S3

df1_filt = glm.1[,c(3,6:8,9,10,13:14,17,21)] # Initial list of variables that were not collinear (Spearman < 0.7 or > -0.7)
df1_filt = na.omit(df1_filt) # Omit any rows that have NA values
model1 <- lm(Abundance ~., data = df1_filt) # Test model for VIFs
car::vif(model1) # Display VIFs from test model - VIFs should be < 10

df1_filt = glm.1[,c(3,5:8,9,10,13:14,17,21)] # Re-do to add back in the region column needed for outlier removal
df1_filt = na.omit(df1_filt) 

# Subset to the top 4 most relatively abundant protist groups in photic zone
df_syn <- df1_filt[df1_filt[["Class"]] == "Syndiniales", ] # Subset to Syndiniales
df_sag <- df1_filt[df1_filt[["Class"]] == "Sagenista", ] # Subset to Sagenista
df_dino<- df1_filt[df1_filt[["Class"]] == "Dinophyceae", ] # Subset to Dinophyceae
df_prym <- df1_filt[df1_filt[["Class"]] == "Prymnesiophyceae", ] # Subset to Prymnesiophyceae

## Workflow for group-specific GAMs in the photic zone

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

# Clean the 27N extreme outliers only from all groups with above function - 2 samples are removed
df_syn_filt <- remove_region_outliers(df_syn, region_col = "region", target_region = "27N", vars = predictor_vars, iqr_multiplier = 3)
df_dino_filt <- remove_region_outliers(df_dino, region_col = "region", target_region = "27N", vars = predictor_vars, iqr_multiplier = 3)
df_sag_filt <- remove_region_outliers(df_sag, region_col = "region", target_region = "27N", vars = predictor_vars, iqr_multiplier = 3)
df_prym_filt <- remove_region_outliers(df_prym, region_col = "region", target_region = "27N", vars = predictor_vars, iqr_multiplier = 3)

## Run through GAM workflow for each group - start with Syndiniales 
# Convert to log abundance and log transform variables that were slightly right skewed
df_syn_filt$LogAbundance <- log(df_syn_filt$Abundance + 1) # Log transform rarefied abundance
log_vars <- c("temp", "salinity", "phosphate", "nh4", "nitrate") # Select predictor variables to log transform

for (var in log_vars) {
  new_var <- paste0("log_", var)
df_syn_filt[[new_var]] <- log(df_syn_filt[[var]] + 0.01)
} # Loop over the variables and transform

set.seed(1234) # Set for reproducibility

# Define predictors again with new log transformed variables
predictors <- c("log_temp", "log_salinity", "oxygen", "log_phosphate", "log_nh4", "log_nitrate", "pH_corrected", "dic")

# Stratified split of data into train and test sets that takes into account spatial dynamics and a gradient in log abundance
df_syn_filt$Strata <- interaction(df_syn_filt$region, cut(df_syn_filt$LogAbundance, breaks = 2)) # Cut data based on region and abundance
train_index <- createDataPartition(df_syn_filt$Strata, p = 0.8, list = FALSE) # Define train set as 80%
train_data <- df_syn_filt[train_index, ] # Train set
test_data  <- df_syn_filt[-train_index, ] # Test set (20%)

# Create all non-empty predictor combinations
predictor_combos <- unlist(lapply(1:length(predictors), function(i) combn(predictors, i, simplify = FALSE)), recursive = FALSE)

# RMSE function to use with abundance data
rmse <- function(actual, predicted) sqrt(mean((actual - predicted)^2))

# Cross fold validation function for Gaussian GAM (log abundance)
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
model <- try(gam(formula_obj, data = train_fold, family = gaussian(), select = TRUE, method = "REML"),silent = TRUE)
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
final_gam_syn <- gam(best_formula, data = train_data, family = gaussian(), select = TRUE, method = "REML") 

# Extract p-values from smooth terms and adjust for multiple comparison to assess significance
gam_summary <- summary(final_gam_syn)
smooth_pvals <- gam_summary$s.table[, "p-value"]
adj_pvals <- p.adjust(smooth_pvals, method = "holm") # Holm correction

# Table of smooth results
smooth_results <- data.frame(
  Term = rownames(gam_summary$s.table),
  EDF = gam_summary$s.table[, "edf"],
  Stat = gam_summary$s.table[, if ("F" %in% colnames(gam_summary$s.table)) "F" else "Chi.sq"],
  RawP = smooth_pvals,
  AdjP = adj_pvals
) 
print(smooth_results)

# Predict on train and test sets
train_data$Predicted <- predict(final_gam_syn, newdata = train_data, type = "response")
test_data$Predicted  <- predict(final_gam_syn, newdata = test_data, type = "response")

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

# Plot observed vs. predicted for train and test for Syndiniales - Figure 5A
ggplot(combined_data, aes(x = LogAbundance, y = Predicted)) +
  geom_point(aes(shape = factor(Shape), fill = Fill, alpha = Alpha), size = 5, stroke = 1, color = "black", show.legend = TRUE) +
  scale_shape_manual(values = c("21" = 21, "24" = 24)) + scale_fill_identity() + scale_alpha_identity() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") + annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.5,
  label = paste0("ρ (Train) = ", round(spearman_train, 2), "\nρ (Test) = ", round(spearman_test, 2)), size = 4.5) +
  annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -1.5, label = paste0("RMSE (Train) = ", round(train_rmse, 2), "\nRMSE (Test) = ", round(test_rmse, 2)), size = 4.5) + 
  labs(title = "Observed vs. Predicted Log Abundance (GAM)", subtitle = equation_text, x = "Observed log(Abundance + 1)", y = "Predicted log(Abundance + 1)", fill = "Data Set") +
  xlim(4.4, 7.5) + ylim(4.4, 7.5) + theme_bw()
#ggsave(filename = "Syn_gam_final.pdf", plot = last_plot(), device = "pdf", path = NULL, scale = 1, width = 6, height = 5, dpi = 150)

# Repeat GAM workflow for Dinophyceae
df_dino_filt$LogAbundance <- log(df_dino_filt$Abundance + 1) # Log transform rarefied abundance
log_vars <- c("temp", "salinity", "phosphate", "nh4", "nitrate") # Select predictor variables to log transform

for (var in log_vars) {
  new_var <- paste0("log_", var)
  df_dino_filt[[new_var]] <- log(df_dino_filt[[var]] + 0.01)
} # Loop over the variables and transform

set.seed(1234) # Set for reproducibility

# Define predictors again with new log transformed variables
predictors <- c("log_temp", "log_salinity", "oxygen", "log_phosphate", "log_nh4", "log_nitrate", "pH_corrected", "dic")

# Stratified split of data into train and test sets
df_dino_filt$Strata <- interaction(df_dino_filt$region, cut(df_dino_filt$LogAbundance, breaks = 2)) # Cut data based on region and abundance
train_index <- createDataPartition(df_dino_filt$Strata, p = 0.8, list = FALSE) # Define train set as 80%
train_data <- df_dino_filt[train_index, ] # Train set
test_data  <- df_dino_filt[-train_index, ] # Test set (20%)

# Create all non-empty predictor combinations
predictor_combos <- unlist(lapply(1:length(predictors), function(i) combn(predictors, i, simplify = FALSE)), recursive = FALSE)

# RMSE function to use with abundant data
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
    model <- try(gam(formula_obj, data = train_fold, family = gaussian(), select = TRUE, method = "REML"),silent = TRUE)
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
final_gam_dino <- gam(best_formula, data = train_data, family = gaussian(), select = TRUE, method = "REML") 

# Extract p-values from smooth terms and adjust for multiple comparison to assess significance
gam_summary <- summary(final_gam_dino)
smooth_pvals <- gam_summary$s.table[, "p-value"]
adj_pvals <- p.adjust(smooth_pvals, method = "holm") # Holm correction

# Table of smooth results
smooth_results <- data.frame(
  Term = rownames(gam_summary$s.table),
  EDF = gam_summary$s.table[, "edf"],
  Stat = gam_summary$s.table[, if ("F" %in% colnames(gam_summary$s.table)) "F" else "Chi.sq"],
  RawP = smooth_pvals,
  AdjP = adj_pvals
) 
print(smooth_results)

# Predict on train and test sets
train_data$Predicted <- predict(final_gam_dino, newdata = train_data, type = "response")
test_data$Predicted  <- predict(final_gam_dino, newdata = test_data, type = "response")

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

# Plot observed vs. predicted for train and test for Dinophyceae - Figure 5B
ggplot(combined_data, aes(x = LogAbundance, y = Predicted)) +
  geom_point(aes(shape = factor(Shape), fill = Fill, alpha = Alpha), size = 5, stroke = 1, color = "black", show.legend = TRUE) +
  scale_shape_manual(values = c("21" = 21, "24" = 24)) + scale_fill_identity() + scale_alpha_identity() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") + annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.5,
  label = paste0("ρ (Train) = ", round(spearman_train, 2), "\nρ (Test) = ", round(spearman_test, 2)), size = 4.5) +
  annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -1.5, label = paste0("RMSE (Train) = ", round(train_rmse, 2), "\nRMSE (Test) = ", round(test_rmse, 2)), size = 4.5) + 
  labs(title = "Observed vs. Predicted Log Abundance (GAM)", subtitle = equation_text, x = "Observed log(Abundance + 1)", y = "Predicted log(Abundance + 1)", fill = "Data Set") +
  xlim(6, 7.5) + ylim(6, 7.5)+ theme_bw()
#ggsave(filename = "Dino_gam_final.pdf", plot = last_plot(), device = "pdf", path = NULL, scale = 1, width = 6, height = 5, dpi = 150)

# Repeat GAM workflow for Prymnesiophyceae
df_prym_filt$LogAbundance <- log(df_prym_filt$Abundance + 1) # Log transform rarefied abundance
log_vars <- c("temp", "salinity", "phosphate", "nh4", "nitrate") # Select predictor variables to log transform

for (var in log_vars) {
  new_var <- paste0("log_", var)
  df_prym_filt[[new_var]] <- log(df_prym_filt[[var]] + 0.01)
} # Loop over the variables and transform

set.seed(1234) # Set for reproducibility

# Define predictors again with new log transformed variables
predictors <- c("log_temp", "log_salinity", "oxygen", "log_phosphate", "log_nh4", "log_nitrate", "pH_corrected", "dic")

# Stratified split of data into train and test sets
df_prym_filt$Strata <- interaction(df_prym_filt$region, cut(df_prym_filt$LogAbundance, breaks = 2)) # Cut data based on region and abundance
train_index <- createDataPartition(df_prym_filt$Strata, p = 0.8, list = FALSE) # Define train set as 80%
train_data <- df_prym_filt[train_index, ] # Train set
test_data  <- df_prym_filt[-train_index, ] # Test set (20%)

# Create all non-empty predictor combinations
predictor_combos <- unlist(lapply(1:length(predictors), function(i) combn(predictors, i, simplify = FALSE)), recursive = FALSE)

# RMSE function to use with abundant data
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
    model <- try(gam(formula_obj, data = train_fold, family = gaussian(), select = TRUE, method = "REML"),silent = TRUE)
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
final_gam_prym <- gam(best_formula, data = train_data, family = gaussian(), select = TRUE, method = "REML") 

# Extract p-values from smooth terms and adjust for multiple comparison to assess significance
gam_summary <- summary(final_gam_prym)
smooth_pvals <- gam_summary$s.table[, "p-value"]
adj_pvals <- p.adjust(smooth_pvals, method = "holm") # Holm correction

# Table of smooth results
smooth_results <- data.frame(
  Term = rownames(gam_summary$s.table),
  EDF = gam_summary$s.table[, "edf"],
  Stat = gam_summary$s.table[, if ("F" %in% colnames(gam_summary$s.table)) "F" else "Chi.sq"],
  RawP = smooth_pvals,
  AdjP = adj_pvals
) 
print(smooth_results)

# Predict on train and test sets
train_data$Predicted <- predict(final_gam_prym, newdata = train_data, type = "response")
test_data$Predicted  <- predict(final_gam_prym, newdata = test_data, type = "response")

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

# Plot observed vs. predicted for train and test for Prymnesiophyceae - Figure 5C
ggplot(combined_data, aes(x = LogAbundance, y = Predicted)) +
  geom_point(aes(shape = factor(Shape), fill = Fill, alpha = Alpha), size = 5, stroke = 1, color = "black", show.legend = TRUE) +
  scale_shape_manual(values = c("21" = 21, "24" = 24)) + scale_fill_identity() + scale_alpha_identity() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") + annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.5,
  label = paste0("ρ (Train) = ", round(spearman_train, 2), "\nρ (Test) = ", round(spearman_test, 2)), size = 4.5) +
  annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -1.5, label = paste0("RMSE (Train) = ", round(train_rmse, 2), "\nRMSE (Test) = ", round(test_rmse, 2)), size = 4.5) + 
  labs(title = "Observed vs. Predicted Log Abundance (GAM)", subtitle = equation_text, x = "Observed log(Abundance + 1)", y = "Predicted log(Abundance + 1)", fill = "Data Set") +
  xlim(0, 6) + ylim(0, 6)+ theme_bw()
#ggsave(filename = "Prym_gam_final.pdf", plot = last_plot(), device = "pdf", path = NULL, scale = 1, width = 6, height = 5, dpi = 150)

# Repeat GAM workflow for Sagenista
df_sag_filt$LogAbundance <- log(df_sag_filt$Abundance + 1) # Log transform rarefied abundance
log_vars <- c("temp", "salinity", "phosphate", "nh4", "nitrate") # Select predictor variables to log transform

for (var in log_vars) {
  new_var <- paste0("log_", var)
  df_sag_filt[[new_var]] <- log(df_sag_filt[[var]] + 0.01)
} # Loop over the variables and transform

set.seed(1234) # Set for reproducibility

# Define predictors again with new log transformed variables
predictors <- c("log_temp", "log_salinity", "oxygen", "log_phosphate", "log_nh4", "log_nitrate", "pH_corrected", "dic")

# Stratified split of data into train and test sets
df_sag_filt$Strata <- interaction(df_sag_filt$region, cut(df_sag_filt$LogAbundance, breaks = 2)) # Cut data based on region and abundance
train_index <- createDataPartition(df_sag_filt$Strata, p = 0.8, list = FALSE) # Define train set as 80%
train_data <- df_sag_filt[train_index, ] # Train set
test_data  <- df_sag_filt[-train_index, ] # Test set (20%)

# Create all non-empty predictor combinations
predictor_combos <- unlist(lapply(1:length(predictors), function(i) combn(predictors, i, simplify = FALSE)), recursive = FALSE)

# RMSE function to use with abundant data
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
    model <- try(gam(formula_obj, data = train_fold, family = gaussian(), select = TRUE, method = "REML"),silent = TRUE)
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
final_gam_sag <- gam(best_formula, data = train_data, family = gaussian(), select = TRUE, method = "REML") 

# Extract p-values from smooth terms and adjust for multiple comparison to assess significance
gam_summary <- summary(final_gam_sag)
smooth_pvals <- gam_summary$s.table[, "p-value"]
adj_pvals <- p.adjust(smooth_pvals, method = "holm") # Holm correction

# Table of smooth results
smooth_results <- data.frame(
  Term = rownames(gam_summary$s.table),
  EDF = gam_summary$s.table[, "edf"],
  Stat = gam_summary$s.table[, if ("F" %in% colnames(gam_summary$s.table)) "F" else "Chi.sq"],
  RawP = smooth_pvals,
  AdjP = adj_pvals
) 
print(smooth_results)

# Predict on train and test sets
train_data$Predicted <- predict(final_gam_sag, newdata = train_data, type = "response")
test_data$Predicted  <- predict(final_gam_sag, newdata = test_data, type = "response")

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

# Plot observed vs. predicted for train and test for Sagenista - Figure 5D
ggplot(combined_data, aes(x = LogAbundance, y = Predicted)) +
  geom_point(aes(shape = factor(Shape), fill = Fill, alpha = Alpha), size = 5, stroke = 1, color = "black", show.legend = TRUE) +
  scale_shape_manual(values = c("21" = 21, "24" = 24)) + scale_fill_identity() + scale_alpha_identity() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") + annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.5,
  label = paste0("ρ (Train) = ", round(spearman_train, 2), "\nρ (Test) = ", round(spearman_test, 2)), size = 4.5) +
  annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -1.5, label = paste0("RMSE (Train) = ", round(train_rmse, 2), "\nRMSE (Test) = ", round(test_rmse, 2)), size = 4.5) + 
  labs(title = "Observed vs. Predicted Log Abundance (GAM)", subtitle = equation_text, x = "Observed log(Abundance + 1)", y = "Predicted log(Abundance + 1)", fill = "Data Set") +
  xlim(2, 5.5) + ylim(2, 5.5)+ theme_bw()
#ggsave(filename = "Sag_gam_final.pdf", plot = last_plot(), device = "pdf", path = NULL, scale = 1, width = 6, height = 5, dpi = 150)

## Loop through above GAMs and plot residuals and QQ-plots for all groups in a single panel
final_gams <- list(Syndiniales = final_gam_syn, Dinophyceae = final_gam_dino, Prymnesiophyceae = final_gam_prym, Sagenista = final_gam_sag) # Name models in a list 

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

# Combine diagnostics into a single panel - Figure S3
wrap_plots(diagnostic_panels, ncol = 2)
#ggsave(filename = "GAM_residuals_18S.pdf", plot = last_plot(), device = "pdf", path = NULL, scale = 1, width = 10, height = 10, dpi = 150)

## Compare AIC values for all 18S groups with different GAM and GLM distributions
# Define data frames - make sure these have abundance, logAbundance, and log predictors
data_list <- list(Dinophyceae = df_dino_filt, Prymnesiophyceae = df_prym_filt, Syndiniales = df_syn_filt, Sagenista = df_sag_filt)

# Define log-transformed predictors used in models
predictors <- c("log_temp", "log_salinity", "oxygen", "log_phosphate", "log_nh4", "log_nitrate", "pH_corrected", "dic")

# Initialize results table
aic_results <- data.frame(Group = character(), AIC_GLM_Poisson = numeric(), AIC_GLM_NB = numeric(), AIC_GAM_Gaussian_Log = numeric(), AIC_GAM_NB = numeric(), stringsAsFactors = FALSE)

# Loop through each dataset and run different models
for (group_name in names(data_list)) {df <- data_list[[group_name]]
  df <- df %>% drop_na(any_of(c("Abundance", "LogAbundance", predictors)))

# Model formulas
formula_raw <- as.formula(paste("Abundance ~", paste(predictors, collapse = " + ")))
formula_log <- as.formula(paste("LogAbundance ~", paste0("s(", predictors, ", bs='tp')", collapse = " + ")))
formula_nb_gam <- as.formula(paste("Abundance ~", paste0("s(", predictors, ", bs='tp')", collapse = " + ")))

# Fit Models
glm_pois <- glm(formula_raw, data = df, family = poisson()) # GLM Poisson
glm_nb <- tryCatch(glm.nb(formula_raw, data = df), error = function(e) NULL) # GLM nb
gam_gauss <- gam(formula_log, data = df, family = gaussian()) # GAM guass
gam_nb <- gam(formula_nb_gam, data = df, family = nb()) # GAM nb

# Compute AICs for all models
aic_glm_pois <- AIC(glm_pois)
aic_glm_nb <- if (!is.null(glm_nb)) AIC(glm_nb) else NA
aic_gam_gauss <- AIC(gam_gauss)
aic_gam_nb <- AIC(gam_nb)

# Combine results
aic_results <- rbind(aic_results, data.frame(Group = group_name, AIC_GLM_Poisson = round(aic_glm_pois, 2), AIC_GLM_NB = round(aic_glm_nb, 2), AIC_GAM_Gaussian_Log = round(aic_gam_gauss, 2), AIC_GAM_NB = round(aic_gam_nb, 2)))
}

# Export AIC results for all groups and models tested - Table S4
print(aic_results)
#write.csv(aic_results, "GAM_GLM_AIC_comparison.csv", row.names = FALSE)

## Plot GAM partial effects for each 18S group - Figure 7
syn_pe <- draw(final_gam_syn, residuals = TRUE, resid_col = "black", ci_col = "#88CCEE", ci_alpha = 0.5, ncol = 4)
syn_pe & theme_bw()
#ggsave(filename = "Syn_partial_effects.pdf", plot = last_plot(), device = "pdf", path = NULL, scale = 1, width = 9, height = 5, dpi = 150)

dino_pe <- draw(final_gam_dino, residuals = TRUE, resid_col = "black", ci_col = "#CC6677", ci_alpha = 0.5, ncol = 4)
dino_pe & theme_bw()
#ggsave(filename = "Dino_partial_effects.pdf", plot = last_plot(), device = "pdf", path = NULL, scale = 1, width = 9, height = 5, dpi = 150)

prym_pe <- draw(final_gam_prym, residuals = TRUE, resid_col = "black", ci_col = "#DDCC77", ci_alpha = 0.5, ncol = 4)
prym_pe & theme_bw()
#ggsave(filename = "Prym_partial_effects.pdf", plot = last_plot(), device = "pdf", path = NULL, scale = 1, width = 9, height = 5, dpi = 150)

sag_pe <- draw(final_gam_sag, residuals = TRUE, resid_col = "black", ci_col = "#332288", ci_alpha = 0.5, ncol = 4)
sag_pe & theme_bw()
#ggsave(filename = "Sag_partial_effects.pdf", plot = last_plot(), device = "pdf", path = NULL, scale = 1, width = 8, height = 2.5, dpi = 150)

## Predict log abundance at all GOMECC-4 sites
# Models for 18S groups were used to predict log abundance at 135 GOMECC-4 sites at the surface (< 10 m)
# Data for 6 sites were not available at the surface. This resulted in 84 sites where DNA was not collected, allowing for an expanded view of microbial distribution in Gulf surface waters.
gomecc_all <- read.csv("/Users/seananderson/Documents/NOAA NGI/GOMECC-4/data-input/GOMECC-all-sites-final.csv", header=T, row.names = NULL, check.names=F,fileEncoding="UTF-8-BOM") # Load file of relevant metadata from all surface sites 
gomecc_all_trans <- gomecc_all # Rename

# Apply log(+ 0.01) transform to selected skewed predictors
log_vars <- c("temp", "salinity", "phosphate", "nh4", "nitrate")
for (v in log_vars) {
  new_colname <- paste0("log_", v)
  gomecc_all_trans[[new_colname]] <- log(gomecc_all_trans[[v]] + 0.01)
}

gomecc_all_syn = gomecc_all_trans[,c(9,15,18:21)] # Subset to variables in final Syndiniales model
new_val_syn = predict(final_gam_syn, newdata = gomecc_all_syn,type = "response") # Predict new values based on GAM
new_val_syn = as.data.frame(new_val_syn)
odv_syn = cbind(new_val_syn,gomecc_all_trans[4:6])

gomecc_all_dino = gomecc_all_trans[,c(9,13,15,17:18,20:21)] # Subset to variables in final Dinophyceae model
new_val_dino = predict(final_gam_dino, newdata=gomecc_all_dino,type="response") # Predict new values
new_val_dino = as.data.frame(new_val_dino)

gomecc_all_prym = gomecc_all_trans[,c(9,13,15,17:19)] # Subset to variables in final Prymnesiophyceae model
new_val_prym = predict(final_gam_prym, newdata = gomecc_all_prym,type = "response") # Predict new values
new_val_prym = as.data.frame(new_val_prym)

gomecc_all_sag = gomecc_all_trans[,c(17,20:21)] # Subset to variables in final Sagenista model
new_val_sag = predict(final_gam_sag, newdata=gomecc_all_sag,type="response") # Predict new values
new_val_sag = as.data.frame(new_val_sag)

# Combine new predicted log abundance values into one data frame. This file with predicted values will be uploaded and analyzed in ODV (ODV files provided on GitHub).
all_odv = cbind(odv_syn, new_val_sag, new_val_dino, new_val_prym)
#write.csv(all_odv, "ODV_microbes_abund.csv",row.names = T) # Write .csv file of log abundance
#write.csv(gomecc_all_trans, "ODV_factors.csv",row.names = T) # Write .csv file of log factors

## Species indicator analysis
# Explore microbes at the ASV level that are indicative of different acidic conditions in the GOM at this time. As with the models, focus on samples from the photic zone (Cluster 1). 
# We used the TA:DIC (total alkalinity:dissolved inorganic carbon) ratio as it is a good proxy of OA, indicating the buffering capacity of seawater 
# A lower TA:DIC = less buffered waters (more acidic), while a higher TA:DIC = more buffered waters

# Plot values of TA:DIC in the photic zone to manually define cutoff for low vs. high TA:DIC - Figure 9A
ps_dsq <- subset_samples(ps_subset, cluster_18S == "Cluster 1") # Subset to Cluster 1
data <- as(sample_data(ps_dsq), "data.frame") # Convert to data frame 
data = data[!is.na(data$alk_dic_oa),] # Remove NA values
histogram <- ggplot(data, aes(x = alk_dic_ratio)) # Density histogram of TA:DIC
histogram + geom_histogram(aes(y = ..density..), 
colour = "black", fill = "white") + geom_density(alpha = 0.2, fill = "#FF6666") + 
geom_vline(aes(xintercept = 1.16), linetype = "dashed", size = 0.6) + theme_bw() # Set cutoff line at 1.16 based on plot
#ggsave(filename = "18S_dic_alk_hist_final.pdf", plot = last_plot(), device = "pdf", path = NULL, scale = 1, width = 5, height = 4, dpi = 150)

# Plot TA:DIC ratio vs. corrected in situ pH values and color samples by transect in the photic zone - Figure 9B
p <- ggscatter(data, x = "alk_dic_ratio", y = "pH_corrected", cor.method = "pearson", cor.coef =T, point = "FALSE", add.params = list(color = "black", fill = "darkgray"), col = "region",add = "reg.line", conf.int = TRUE, ylab = "pH", xlab = "TA:DIC",title = "TA:DIC ratio GOM") + 
geom_point(aes(fill = region), size = 5, shape = 21, colour = "black", position="jitter") + scale_fill_manual(values = mycolors) + theme(legend.position = "right") + theme(text = element_text(size = 18)) + theme_bw()+ geom_vline(aes(xintercept = 1.16), linetype = "dashed", size = 0.6)
p$data$region <- factor(p$data$region, levels = c("27N", "FLSTRAITS","CAPECORAL", "TAMPA", "PANAMACITY", "PENSACOLA", "LA", "GALVESTON", "PAISNP", "BROWNSVILLE", "TAMPICO", "VERACRUZ", "CAMPECHE", "MERIDA", "YUCATAN", "CATOCHE", "CANCUN")) # Set order of transects
p
#ggsave(filename = "18S_dic_alk_pH_final.pdf", plot = last_plot(), device = "pdf", path = NULL, scale = 1, width = 6, height = 4, dpi = 150)

# Run the core indicator analysis program using 18S ASVs
ps_dsq <- subset_samples(ps_subset, cluster_18S == "Cluster 1") # Subset to Cluster 1
ps_dsq <- subset_samples(ps_dsq, alk_dic_oa == "High" | alk_dic_oa == "Low") # Split samples based on high or low TA:DIC (see above)
ASV = otu_table(ps_dsq) # Subset ASVs
ASV = as.data.frame(ASV) # Convert to data frame
if(taxa_are_rows(ps_dsq)){ASV <- t(ASV)} # Flip table for analysis
ASV = as.data.frame(ASV) # Convert to df again

data <- as(sample_data(ps_dsq), "data.frame") # Subset metadata
data = data[!is.na(data$alk_dic_oa),] # Remove TA:DIC NA values 
bin = data$alk_dic_oa # Bin based on low vs. high categories
inv = multipatt(ASV, bin, func = "r.g", control = how(nperm = 999)) # Core function using 999 permutations

# Correct p-values for multiple comparisons
indisp.sign <- data.table::as.data.table(inv$sign, keep.rownames = TRUE) # Extract table of stats
indisp.sign[ ,p.value.bh:=p.adjust(p.value, method = "BH")] # Add adjusted p-value with Benjamini-Hochberg (more strict)
indisp.sign2 = indisp.sign[p.value <= 0.05,] # Filter to include only p-values < 0.05; we will have adjusted values as well
#write.csv(indisp.sign2, "18S_indic_corr2.csv",row.names = T) # Write .csv file

# Results of indicator analysis (ASV, stat, and corrected p-value) combined in a single .csv file
# Upload this file into R (provided on GitHub)
# Match indicator values table for each ASV with the relative abundance of the same ASV in the photic zone via the ASV count table
# Data frames need to match

# Upload table with all indicator metadata for high and low TA:DIC categories
indic_OA <- read.csv("/Users/seananderson/Documents/NOAA NGI/GOMECC-4/data-input/18S/OA_indic_new.csv", header = T, row.names = NULL, check.names = F,fileEncoding = "UTF-8-BOM") # Indicator ASVs file
colnames(indic_OA)[1]  <- "OTU" # Change column names
indic.high = indic_OA[indic_OA[["OA"]] == "High", ] # Subset to new data frame with only high TA:DIC
indic.low = indic_OA[indic_OA[["OA"]] == "Low", ] # Subset to new df with only low TA:DIC

ps4 <-  transform_sample_counts(ps_dsq, function(x) 100*x / sum(x) ) # Transform to relative abundance
x1 = speedyseq::psmelt(ps4) # Melt the data
OA.1 = x1 %>% # Group to the mean relative abundance at the ASV level
  dplyr::group_by(OTU,alk_dic_oa, Division, Genus, Species) %>%
  dplyr::summarise_at(.vars = c("Abundance"), .funs = mean)

OA.high = OA.1[OA.1[["alk_dic_oa"]] == "High", ] # Subset ASV relative abundance based on high TA:DIC
OA.low = OA.1[OA.1[["alk_dic_oa"]] == "Low", ] # Subset ASV based on low TA:DIC

df_h <- OA.high[OA.high$OTU %in% indic.high$OTU,] # Keep ASVs if they match indicator dataset
df_L <- OA.low[OA.low$OTU %in% indic.low$OTU,] # Keep ASVs if they match indicator dataset

#write.csv(df_h, "ASVs_18S_abund_high_new.csv",row.names = T) # Write .csv file
#write.csv(df_L, "ASVs_18S_abund_low_new.csv",row.names = T) # Write .csv file

# The previous two files were used in combination with the indicator file to create a master .csv file that has ASV IDs, indicator values, and relative abundance for each specific ASV in the photic zone
# This file to be loaded into R (provided on GitHub)

# Upload final indicator file and plot results - Figure 9C-D
indic_OA2 <- read.csv("/Users/seananderson/Documents/NOAA NGI/GOMECC-4/data-input/18S/18S_indic_corrected.csv", header = T, row.names = NULL, check.names = F,fileEncoding = "UTF-8-BOM") # Upload file
p <- ggplot(indic_OA2, aes(x = stat, y = Abundance)) + geom_point(aes(fill = Division), shape = 21,size = 5) + theme_bw() + theme(text = element_text(size = 14)) +
scale_y_continuous(limits = c(0, NA)) + facet_wrap(~OA, scales = "free") + scale_fill_brewer(palette = "Dark2") 
p$data$Division <- factor(p$data$Division, levels = c("Alveolata", "Chlorophyta", "Discoba", "Haptophyta", "Rhizaria", "Stramenopiles", "Others")) # Set order
p
#ggsave(filename = "18S_indicator_final.pdf", plot = last_plot(), device = "pdf", path = NULL, scale = 1, width = 12, height = 5, dpi = 150)

## Estimate spearman values of variables for Clusters 2-3 and plot correlation matrices for all clusters
# Spearman values for Clusters 2-3
cluster2 <- x1[x1[["cluster_18S"]] == "Cluster 2", ] # Subset to Cluster 2
clus.2 = cluster2 %>% # Group the data based on class and preserve sample replication for the models
  dplyr::group_by(station, depth_category, Class, replicate) %>%
  dplyr::summarise(temp = mean(temp),
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

cluster3 <- x1[x1[["cluster_18S"]] == "Cluster 3", ] # Subset to Cluster 3
clus.3 = cluster3 %>% # Group the data based on class and preserve sample replication for the models
  dplyr::group_by(station, depth_category, Class, replicate) %>%
  dplyr::summarise(temp = mean(temp),
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
#write.csv(correlations2, "Cluster2_corr_18S.csv",row.names=T) # Write .csv file - Table S3

# Spearman Cluster 3
df3 = clus.3[,c(5:19)] # Subset to include environmental variables
df3 = na.omit(df3) # Omit any rows that have NA values
correlations3 = cor(df3, method = "spearman") # Perform Spearman correlations to assess collinearity
#write.csv(correlations3, "Cluster3_corr_18S.csv",row.names=T) # Write .csv file - Table S3

# Correlation matrices - start with Cluster 1
df1 = glm.1[,c(6:20)] # Subset to include environmental variables
df1 = na.omit(df1) # Omit any rows that have NA values
factors = df1[!duplicated(df1$carbonate), ] # Remove duplicates
correlations = cor(df1, method = "spearman") # Perform Spearman correlations to assess collinearity

res2 <- Hmisc::rcorr(as.matrix(factors), (as.matrix(factors)), type = "spearman") 
coeff = res2$r # Subset the coefficients
coeff = coeff[, c(1:15)]
coeff= coeff[-c(1:15), ] 

pvalue = res2$P # Subset the p-values
pvalue[is.na(pvalue)] <- 0
pvalue = pvalue[, c(1:15)]
pvalue= pvalue[-c(1:15), ]

# Plot correlation matrix for Cluster 1
corrplot(coeff, p.mat = pvalue, outline = TRUE, type = "full", insig = "blank", sig.level = 0.05, pch.cex = .9, tl.col = "black", tl.srt = 90, method = "color", addgrid.col = "black", order = "hclust") 

# Repeat for Cluster 2 and 3
df1 = clus.2[,c(5:19)] # Subset to include environmental variables
df1 = na.omit(df1) # Omit any rows that have NA values
factors = df1[!duplicated(df1$carbonate), ] # Remove duplicates
correlations = cor(df1, method = "spearman") # Perform Spearman correlations to assess collinearity

res2 <- Hmisc::rcorr(as.matrix(factors), (as.matrix(factors)), type = "spearman") 
coeff = res2$r # Subset the coefficients
coeff = coeff[, c(1:15)]
coeff= coeff[-c(1:15), ] 

pvalue = res2$P # Subset the p-values
pvalue[is.na(pvalue)] <- 0
pvalue = pvalue[, c(1:15)]
pvalue= pvalue[-c(1:15), ]

# Plot correlation matrix for Cluster 2
corrplot(coeff, p.mat = pvalue, outline = TRUE, type = "full", insig = "blank", sig.level = 0.05, pch.cex = .9, tl.col = "black", tl.srt = 90, method = "color", addgrid.col = "black", order = "hclust") 

# Spearman Cluster 3
df1 = clus.3[,c(5:19)] # Subset to include environmental variables
df1 = na.omit(df1) # Omit any rows that have NA values
factors = df1[!duplicated(df1$carbonate), ] # Remove duplicates
correlations = cor(df1, method = "spearman") # Perform Spearman correlations to assess collinearity

res2 <- Hmisc::rcorr(as.matrix(factors), (as.matrix(factors)), type = "spearman") 
coeff = res2$r # Subset the coefficients
coeff = coeff[, c(1:15)]
coeff= coeff[-c(1:15), ] 

pvalue = res2$P # Subset the p-values
pvalue[is.na(pvalue)] <- 0
pvalue = pvalue[, c(1:15)]
pvalue= pvalue[-c(1:15), ]

# Plot correlation matrix for Cluster 3
corrplot(coeff, p.mat = pvalue, outline = TRUE, type = "full", insig = "blank", sig.level = 0.05, pch.cex = .9, tl.col = "black", tl.srt = 90, method = "color", addgrid.col = "black", order = "hclust") 

## GOM map of DNA samples collected - Figure 1A
b = getNOAA.bathy(lon1 = -100, lon2 = -79, lat1 = 30, lat2 = 22, resolution = 3, keep=TRUE) # Import bathymetry data
sites = read.csv(file="/Users/seananderson/Documents/NOAA NGI/GOMECC-4/data-input/18S/map_coord_full.csv", header=T, row.names = NULL, check.names=F,fileEncoding="UTF-8-BOM") # Import site information
sites$Latitude <- as.character(sites$Latitude) # Format latitude
sites$Longitude <- as.character(sites$Longitude) # Format longitude
sites$Distance <- as.factor(sites$Distance) # Format distance to shore
sites$Station <- as.character(sites$Station) # Format station

# Set map colors
blues <- c("lightsteelblue4", "lightsteelblue3","lightsteelblue2", "lightsteelblue1") # Set colors on map
greys <- c(grey(0.7), grey(0.9), grey(0.95)) # Set greys
ap <- c("#9ebcda","#8856a7") # Define color for sample points

# Build the initial map
plot(b, image = TRUE, land = TRUE, n=1, bpal = list(c(0, max(b), greys), c(min(b), 0, blues))) # Bathymetry colors
plot(b, lwd = 1, deep = 0, shallow = 0, step = 0, add = FALSE) # Highlight coastline with solid black line
plot(b, deep=-50, shallow=-50, lwd = 0.4, drawlabels = FALSE, add = TRUE) # Add -50m isobath
plot(b, deep=-200, shallow=-200, lwd = 0.4, drawlabels = FALSE, add = TRUE) # Add -200m isobath
plot(b, deep=-1000, shallow=-1000, lwd = 0.4, drawlabels = FALSE, add = TRUE) # Add -1000m isobath
plot(b, deep=-3000, shallow=-3000, lwd = 0.4, drawlabels = FALSE, add = TRUE) # Add -3000m isobath

coord <- par("usr") # Store coordinate system for plot
b2 <- getNOAA.bathy(-100, -79, 36, 18, res = 3) # Re-run import
#pdf("map_sites_all_final.pdf") # Name pdf file to store the map
plot(b2, image = TRUE, land = TRUE, n = 1, bpal = list(c(0, max(b2), greys), c(min(b2), 0, blues))) # Replot map colors
plot(b2, lwd = 1, deep = 0, shallow = 0, step = 0, add = TRUE) # Highlight coastline with solid black line
plot(b2, deep = -50, shallow = -50, lwd = 0.4, drawlabels = FALSE, add = TRUE) # Add -50m isobath
plot(b2, deep = -200, shallow = -200, lwd = 0.4, drawlabels = FALSE, add = TRUE) # Add -200m isobath
plot(b2, deep = -1000, shallow = -1000, lwd = 0.4, drawlabels = FALSE, add = TRUE) # Add -1000m isobath
plot(b2, deep = -3000, shallow = -3000, lwd = 0.4, drawlabels = FALSE, add = TRUE) # Add -3000m isobath

# Adds dotted line between sampling sites that are along the same transect (16 transects)
lines(sites$Longitude[sites$Station=="27N"], sites$Latitude[sites$Station=="27N"],col="black",lwd=1.5,lty=2)
lines(sites$Longitude[sites$Station=="TAMPA"], sites$Latitude[sites$Station=="TAMPA"],col="black",lwd=1.5,lty=2)
lines(sites$Longitude[sites$Station=="PANAMACITY"], sites$Latitude[sites$Station=="PANAMACITY"],col="black",lwd=1.5,lty=2)
lines(sites$Longitude[sites$Station=="PENSACOLA"], sites$Latitude[sites$Station=="PENSACOLA"],col="black",lwd=1.5,lty=2)
lines(sites$Longitude[sites$Station=="LA"], sites$Latitude[sites$Station=="LA"],col="black",lwd=1.5,lty=2)
lines(sites$Longitude[sites$Station=="GALVESTON"], sites$Latitude[sites$Station=="GALVESTON"],col="black",lwd=1.5,lty=2)
lines(sites$Longitude[sites$Station=="BROWNSVILLE"], sites$Latitude[sites$Station=="BROWNSVILLE"],col="black",lwd=1.5,lty=2)
lines(sites$Longitude[sites$Station=="TAMPICO"], sites$Latitude[sites$Station=="TAMPICO"],col="black",lwd=1.5,lty=2)
lines(sites$Longitude[sites$Station=="VERACRUZ"], sites$Latitude[sites$Station=="VERACRUZ"],col="black",lwd=1.5,lty=2)
lines(sites$Longitude[sites$Station=="CAMPECHE"], sites$Latitude[sites$Station=="CAMPECHE"],col="black",lwd=1.5,lty=2)
lines(sites$Longitude[sites$Station=="MERIDA"], sites$Latitude[sites$Station=="MERIDA"],col="black",lwd=1.5,lty=2)
lines(sites$Longitude[sites$Station=="YUCATAN"], sites$Latitude[sites$Station=="YUCATAN"],col="black",lwd=1.5,lty=2)
lines(sites$Longitude[sites$Station=="CATOCHE"], sites$Latitude[sites$Station=="CATOCHE"],col="black",lwd=1.5,lty=2)
lines(sites$Longitude[sites$Station=="CANCUN"], sites$Latitude[sites$Station=="CANCUN"],col="black",lwd=1.5,lty=2)
lines(sites$Longitude[sites$Station=="FLSTRAITS"], sites$Latitude[sites$Station=="FLSTRAITS"],col="black",lwd=1.5,lty=2)
lines(sites$Longitude[sites$Station=="CAPECORAL"], sites$Latitude[sites$Station=="CAPECORAL"],col="black",lwd=1.5,lty=2)

# Add sampling points to the map
points(sites$Longitude, sites$Latitude, pch = 21, col = "black", cex = 3, asp = 1, bg = ap[sites$Distance]) # Color specified above
scaleBathy(b2, deg = 1, x = "bottomleft", inset = 10, y = NULL, angle = 90) # Add distance scale bar in km
#dev.off() # Print off the above pdf

# Same map with points filled based on 18S cluster assignments at each depth - Figure 1B-D
b = getNOAA.bathy(lon1 = -100, lon2 = -79, lat1 = 30, lat2 = 22, resolution = 3, keep = TRUE) # Import bathymetry data
sites = read.csv(file="/Users/seananderson/Documents/NOAA NGI/GOMECC-4/data-input/18S/map_clus_18S.csv", header = T, row.names = NULL, check.names = F, fileEncoding="UTF-8-BOM") # Import sites
sites$Latitude <- as.character(sites$Latitude) # Format latitude
sites$Longitude <- as.character(sites$Longitude) # Format longitude
sites$Depth <- as.character(sites$Depth) # Format depth zones
sites$Cluster <- as.factor(sites$Cluster) # Format cluster assignments

# Partition based on depth and set map colors
sites_surf <- sites[sites[["Depth"]] == "Surface", ] # Subset to surface
sites_dcm <- sites[sites[["Depth"]] == "DCM", ] # Subset to DCM
sites_bottom <- sites[sites[["Depth"]] == "Deep", ] # Subset to near bottom
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
#pdf("map_sites_surf.pdf") # Name pdf file to store the surface map - store and run through each one at a time
#pdf("map_sites_dcm.pdf") # Name pdf file to store DCM map
#pdf("map_sites_bottom.pdf") # Name pdf file to store bottom map

plot(b2, image = TRUE, land = TRUE, n=1, bpal = list(c(0, max(b2), greys), c(min(b2), 0, blues))) # Replot map colors
plot(b2, lwd = 1, deep = 0, shallow = 0, step = 0, add = TRUE) # Highlight coastline with solid black line
plot(b2, deep= -50, shallow= -50, lwd = 0.4, drawlabels = FALSE, add = TRUE) # Add -50m isobath
plot(b2, deep= -200, shallow= -200, lwd = 0.4, drawlabels = FALSE, add = TRUE) # Add -200m isobath
plot(b2, deep= -1000, shallow= -1000, lwd = 0.4, drawlabels = FALSE, add = TRUE) # Add -1000m isobath
plot(b2, deep= -3000, shallow= -3000, lwd = 0.4, drawlabels = FALSE, add = TRUE) # Add -3000m isobath

# Adds filled sampling points based on cluster assignments - run one at a time and repeat
points(sites_surf$Longitude, sites_surf$Latitude, pch = 21, col = "black",cex = 3,asp = 1, bg = ap[sites_surf$Cluster]) # Add sampling points
points(sites_dcm$Longitude, sites_dcm$Latitude, pch = 21, col = "black",cex = 3,asp = 1, bg = ap[sites_dcm$Cluster]) # Add sampling points
points(sites_bottom$Longitude, sites_bottom$Latitude, pch = 21, col = "black",cex = 3,asp = 1, bg = ap[sites_bottom$Cluster]) # Add sampling points
scaleBathy(b2, deg = 1, x = "bottomleft", inset = 10, y = NULL, angle = 90) # Add scale bar
#dev.off() # Print pdf one at a time (see above)

# R Session information
sessionInfo()