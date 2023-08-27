
# Install necessary packages if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2")

install.packages("dplyr")
install.packages("tibble")

# Load the dada2 package
library(dada2); packageVersion("dada2")

# Download required files and organize them into a folder
# Set the working directory
# Click "Session" -> "Set Working Directory" -> "Choose Directory"
setwd("")
path <- getwd()
list.files(path)

# Read in sequencing files

# Get forward and reverse fastq filenames
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

# Extract sample names from filenames
sample.names <- sapply(strsplit(basename(fnFs), "_S"), `[`, 1)

# Generate and save quality profile plots for forward reads
png("Fwd_Read_Quality_Profiles.png")
plotQualityProfile(fnFs[1:2])
dev.off()

# Generate and save quality profile plots for reverse reads
png("Rev_Read_Quality_Profiles.png")
plotQualityProfile(fnRs[1:2])
dev.off()

# Filter and Trim

# Assign filenames and paths for filtered fastq.gz files
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# Filter and trim reads
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft=c(19,20), rm.phix=TRUE)
head(out)

# Save filter results to a text file
write.table(out, "Read_Filter_In_Out.txt", sep = "\t", quote=F)

# Generate and save quality profile plots for filtered forward reads
png("FiltFwd_Read_Quality_Profiles.png")
plotQualityProfile(filtFs[1:2])
dev.off()

# Generate and save quality profile plots for filtered reverse reads
png("FiltRev_Read_Quality_Profiles.png")
plotQualityProfile(filtRs[1:2])
dev.off()


# Generate an error model of our data
# Learns the specific error-signature of our dataset for later steps
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

### Dereplication
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

### Inferring ASVs - Sample Inference Algorithm
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

### Merge Paired Reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE, justConcatenate=TRUE)

### Construct an amplicon sequence variant table (ASV) table
seqtab <- makeSequenceTable(mergers)

### Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

### Track reads through the pipeline
getN <- function(x) sum(getUniques(x))

track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim), round(rowSums(seqtab.nochim)/out[,1]*100,1))
colnames(track) <- c("input", "filtered", "dadaF", "dadaR", "merged", "nonchim", "final_perc_reads_retained")
rownames(track) <- sample.names

write.table(track, "read-count-tracking.tsv", quote=FALSE, sep="\t", col.names=NA)

### Assign taxonomy with latest Silva database, v.138.1
taxa <- assignTaxonomy(seqtab.nochim, "silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)

write.table(taxa, "Taxonomic_Table.txt", sep = "\t", quote=F)



# CREATE ABUNDANCE TABLE #

### Make matrix table
# Transpose the sequence abundance table
t_seqtab <- as.data.frame(t(seqtab.nochim))
t_seqtab_ab <- t_seqtab
# Calculate the sum of abundances for each ASV
t_seqtab_ab$abund_sum <- rowSums(t_seqtab_ab)

# Merge taxonomy information and abundance data
taxabun <- merge(taxa, t_seqtab_ab, by = 'row.names', type = "full", match = "all")

# Sort rows by most abundant taxa
taxabun <- taxabun[with(taxabun, order(abund_sum, decreasing=TRUE)),]

# Remove the row sums column
taxabun = taxabun[,!(names(taxabun) %in% "abund_sum")]

# Write the absolute abundance table to a CSV file
write.csv(taxabun, "Taxa Absolute Abundance Table.csv", row.names = FALSE, quote = F)

## Import merged metaphlan table
Taxa_Absolute_Abundance_Table <- read_csv("Taxa Absolute Abundance Table.csv")

# Store the imported table in otumat
otumat <- Taxa_Absolute_Abundance_Table

# Create taxa table
taxmat <- data.frame(otumat[2:7])
taxmat <- replace(taxmat, is.na(taxmat), "Unclassified")

# Make column 1 row names and remove columns K01-K48
otumat <- data.frame(otumat[8:52])

# Convert row names of both tables into OTU1...OTU[rownumber]
rownames(otumat) <- paste0("OTU", 1:nrow(otumat))
rownames(taxmat) <- paste0("OTU", 1:nrow(taxmat))
otumat <- as.matrix(otumat)
taxmat <- as.matrix(taxmat)

# Create otu_table and tax_table for phyloseq object
OTU <- otu_table(otumat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)

# Construct sample data.frame as a phyloseq object
material <- list(rep("Non-Plastic", 17), rep("Plastic", 17), rep("Water_Column", 9), rep("PCRNeg", 1), rep("ExtrNeg", 1))
samdf <- data.frame(label = unlist(material))
rownames(samdf) <- colnames(OTU)

# Create the phyloseq object with the updated sample data
ps <- phyloseq(OTU, TAX, sample_data(samdf))

# Filter out the "Chloroplast" order from the dataset
ps_filtered <- subset_taxa(ps, Order != "Chloroplast")
families_to_remove <- c("Mitochondria")

# Remove "PCRNeg" and "ExtrNeg" samples
ps_filtered <- subset_samples(ps_filtered, !(sample_data(ps_filtered)$label %in% c("PCRNeg", "ExtrNeg")))

# Identify the taxa present in PCRNeg and ExtrNeg samples
taxa_to_remove <- taxa_names(ps_filtered)

# Get the taxa names associated with PCRNeg and ExtrNeg samples
taxa_to_remove <- intersect(taxa_to_remove, taxa_names(subset_taxa(ps, sample_data(ps)$label %in% c("PCRNeg", "ExtrNeg"))))

# Remove the identified taxa from the rest of the dataset
ps_filtered_pruned <- prune_taxa(setdiff(taxa_names(ps_filtered), taxa_to_remove), ps_filtered)

# Remove specified samples from the phyloseq object
samples_to_remove <- c("K47", "K48", "K16", "K15", "K17", "K08")
ps <- subset_samples(ps_filtered_pruned, !sample_names(ps_filtered))
                                                       
## Calculate Shared and Unique OTUs ##

# Get sample names for each category
plastic_samples <- sample_names(ps)[sample_data(ps)$label == "Plastic"]
non_plastic_samples <- sample_names(ps)[sample_data(ps)$label == "Non-Plastic"]
water_column_samples <- sample_names(ps)[sample_data(ps)$label == "Water_Column"]

# Subset the OTU table directly based on sample labels
plastic_otu <- otu_table(ps, taxa_are_rows = TRUE)[, plastic_samples]
non_plastic_otu <- otu_table(ps, taxa_are_rows = TRUE)[, non_plastic_samples]
water_column_otu <- otu_table(ps, taxa_are_rows = TRUE)[, water_column_samples]

# Calculate the number of shared and unique OTUs for each comparison
shared_plastic_non_plastic <- sum(taxa_sums(plastic_otu) > 0 & taxa_sums(non_plastic_otu) > 0)
shared_plastic_water_column <- sum(taxa_sums(plastic_otu) > 0 & taxa_sums(water_column_otu) > 0)
shared_non_plastic_water_column <- sum(taxa_sums(non_plastic_otu) > 0 & taxa_sums(water_column_otu) > 0)

# Calculate the number of unique OTUs for each sample type
unique_plastic_otus <- sum(taxa_sums(plastic_otu) > 0) - shared_plastic_non_plastic - shared_plastic_water_column
unique_non_plastic_otus <- sum(taxa_sums(non_plastic_otu) > 0) - shared_plastic_non_plastic - shared_non_plastic_water_column
unique_water_column_otus <- sum(taxa_sums(water_column_otu) > 0) - shared_plastic_water_column - shared_non_plastic_water_column

# Calculate the number of OTUs shared between all three sample types
shared_all_three <- sum(taxa_sums(plastic_otu) > 0 & taxa_sums(non_plastic_otu) > 0 & taxa_sums(water_column_otu) > 0)

# Create a data frame to store the counts, including the new row
otu_counts <- data.frame(
  Category = c("Shared: Plastic & Non-Plastic", "Shared: Plastic & Water_Column", "Shared: Non-Plastic & Water_Column",
               "Unique to Plastic", "Unique to Non-Plastic", "Unique to Water_Column", "Shared: All Three"),
  Count = c(shared_plastic_non_plastic, shared_plastic_water_column, shared_non_plastic_water_column,
            unique_plastic_otus, unique_non_plastic_otus, unique_water_column_otus, shared_all_three))

# Print the data table
print(otu_counts)



#### Statistical test on significance in amount of sample specific families

# Calculate the counts of shared families (based on your previous code)
# Create a contingency table for the comparison between exclusive OTUs in plastic vs. non-plastic
contingency_table_exclusive_plastic_nonplastic <- matrix(c(unique_plastic_otus - shared_plastic_non_plastic, unique_non_plastic_otus - shared_plastic_non_plastic), ncol = 2, byrow = TRUE)

# Run a Chi-squared test
chi_result_exclusive_plastic_nonplastic <- chisq.test(contingency_table_exclusive_plastic_nonplastic)

# Print the results of the Chi-squared test
print("Comparison between Exclusive to Plastic vs. Exclusive to Non-Plastic:")
print(chi_result_exclusive_plastic_nonplastic)

# Create a contingency table for the comparison between exclusive OTUs in plastic vs. water
contingency_table_exclusive_plastic_water <- matrix(c(unique_plastic_otus - shared_plastic_water_column, unique_water_column_otus - shared_plastic_water_column), ncol = 2, byrow = TRUE)

# Run a Chi-squared test
chi_result_exclusive_plastic_water <- chisq.test(contingency_table_exclusive_plastic_water)

# Print the results of the Chi-squared test
print("Comparison between Exclusive to Plastic vs. Exclusive to Water_Column:")
print(chi_result_exclusive_plastic_water)

# Create a contingency table for the comparison between exclusive OTUs in non-plastic vs. water
contingency_table_exclusive_nonplastic_water <- matrix(c(unique_non_plastic_otus - shared_non_plastic_water_column, unique_water_column_otus - shared_non_plastic_water_column), ncol = 2, byrow = TRUE)

# Run a Chi-squared test
chi_result_exclusive_nonplastic_water <- chisq.test(contingency_table_exclusive_nonplastic_water)

# Print the results of the Chi-squared test
print("Comparison between Exclusive to Non-Plastic vs. Exclusive to Water_Column:")
print(chi_result_exclusive_nonplastic_water)

#### Statistical test on significance in amount of families shared in substrates vs substrate + water

# Perform a proportion test comparing Plastic vs. Non-Plastic shared proportions
prop_test_result_plastic_nonplastic <- prop.test(
  x = c(shared_plastic_non_plastic, shared_plastic_water_column),
  n = c(shared_plastic_non_plastic + unique_plastic_otus, shared_plastic_water_column + unique_plastic_otus),
  alternative = "greater",
  correct = FALSE
)

# Perform a proportion test comparing Non-Plastic vs. Water_Column shared proportions
prop_test_result_non_plastic_water <- prop.test(
  x = c(shared_non_plastic_water_column, shared_plastic_water_column),
  n = c(shared_non_plastic_water_column + unique_non_plastic_otus, shared_plastic_water_column + unique_plastic_otus),
  alternative = "greater",
  correct = FALSE)

# Print the results of the proportion tests for the comparisons
print("Comparison between Plastic vs. Non-Plastic and Plastic vs. Water_Column:")
print(prop_test_result_plastic_nonplastic)

print("Comparison between Non-Plastic vs. Water_Column and Non-Plastic vs. Plastic:")
print(prop_test_result_non_plastic_water)





##Create a bar plot of family abundances across samples
# This plot visualizes the distribution of family abundances within different sample types

# Transform the sample counts to percentages
ps.rel = transform_sample_counts(ps, function(x) x/sum(x)*100)

# Agglomerate taxa at the Family level
glom <- tax_glom(ps.rel, taxrank = 'Family', NArm = FALSE)
ps.melt <- psmelt(glom)
ps.melt$Family <- as.character(ps.melt$Family)  # Convert Family to character

# Calculate median abundances for each Family within each sample type
ps.melt <- ps.melt %>%
  group_by(label, Family) %>%
  mutate(median=median(Abundance))

# Keep only Families with median abundance > 1
keep <- unique(ps.melt$Family[ps.melt$median > 1])
ps.melt$Family[!(ps.melt$Family %in% keep)] <- "< 1%"

# Summarize abundances by Sample, label, and Family
ps.melt_sum <- ps.melt %>%
  group_by(Sample, label, Family) %>%
  summarise(Abundance=sum(Abundance))

# Create a grouped bar plot using ggplot
# Create a grouped bar plot to visualize sample abundance by family

# Set the data source and aesthetics for x (Sample), y (Abundance), and fill (Family)
ggplot(ps.melt_sum, aes(x = Sample, y = Abundance, fill = Family)) + 
  
  # Add bars with heights based on the 'Abundance' variable
  geom_bar(stat = "identity") + 
  
  # Set the color scale using viridis colors with specific range
  scale_fill_viridis_d(option = "plasma", begin = 0, end = 0.8) +
  
  # Customize axis labels by removing them
  labs(x = "", y = "%") +
  
  # Create facets (subplots) for each 'label' with free x-axis scales
  facet_wrap(~label, scales = "free_x", nrow = 1) +
  
  # Use the classic theme for the plot
  theme_classic() + 
  
  # Customize specific theme elements
  theme(
    strip.background = element_blank(),  # Remove background of facet labels
    axis.text.x.bottom = element_text(angle = -90),  # Rotate x-axis text
    legend.text = element_text(size = 15),  # Customize legend text size
    axis.text = element_text(size = 13),  # Customize axis text size
    strip.text = element_text(size = 13)  # Customize facet label text size
  ) +
  
  # Customize the appearance of the legend
  guides(fill = guide_legend(nrow = 23))  # Set the number of rows in the legend


# Calculate alpha diversity indices (Shannon and Simpson) at the Family level
# This code calculates alpha diversity for each sample type and performs statistical tests

# Subset the phyloseq object to include only non-zero taxa
Family_physeq <- prune_taxa(taxa_sums(Family_physeq) > 0, Family_physeq)

# Create a data frame to store alpha diversity data
alpha_diversity_df <- data.frame(Subset_Names = character(),
                                 Shannon_Index = numeric(),
                                 Simpson_Index = numeric(),
                                 stringsAsFactors = FALSE)

# Loop through each sample type and calculate alpha diversity
for (sample_type in sample_types) {
  # Subset the phyloseq object to the current sample type
  sub_physeq <- subset_samples(Family_physeq, label == sample_type)
  
  # Calculate alpha diversity (Shannon and Simpson indices) for the current sample type
  alpha_diversity <- estimate_richness(sub_physeq, measures = c("Shannon", "Simpson"))
  
  # Add the alpha diversity data to the data frame
  alpha_diversity_df <- rbind(alpha_diversity_df,
                              data.frame(Subset_Names = sample_type,
                                         Shannon_Index = alpha_diversity$Shannon,
                                         Simpson_Index = alpha_diversity$Simpson))}

# Print the combined alpha diversity data
print(alpha_diversity_df)

# Perform Kruskal-Wallis tests and pairwise Wilcoxon rank-sum tests
# This code performs statistical tests to compare alpha diversity among sample types

# Perform Kruskal-Wallis test on Shannon Index
kruskal.test(Shannon_Index ~ Subset_Names, data = alpha_diversity_df)

# Perform Kruskal-Wallis test on Simpson Index
kruskal.test(Simpson_Index ~ Subset_Names, data = alpha_diversity_df)

# Perform pairwise Wilcoxon rank-sum tests for Shannon Index
posthoc_shannon <- pairwise.wilcox.test(alpha_diversity_df$Shannon_Index, alpha_diversity_df$Subset_Names, p.adjust.method = "holm")

# Print the post hoc test results for Shannon Index
print(posthoc_shannon)

# Perform pairwise Wilcoxon rank-sum tests for Simpson Index
posthoc_simpson <- pairwise.wilcox.test(alpha_diversity_df$Simpson_Index, alpha_diversity_df$Subset_Names, p.adjust.method = "holm")

# Print the post hoc test results for Simpson Index
print(posthoc_simpson)

###
# Beta diversity NMDS plot and supporting statistics
# This code calculates and visualizes the beta diversity using NMDS and associated statistics

# Transform the OTU table using centered log-ratio transformation
otu_z_Family <- cmultRepl(as.matrix(Family_physeq@otu_table), method = "GBM", output = 'p-counts', z.warning = 0.99)
clr <- function(x) sweep(log(x), 1, rowMeans(log(x)), "-")
otu_tx_Family <- data.frame(t(clr(t(otu_z_Family))))
otu_m_Family <- as.matrix(t(otu_tx_Family))

# Perform NMDS analysis on the transformed OTU table
nmds_Family <- metaMDS(otu_m_Family, distance = "euclidean")
nmds_scores_Family <- as.data.frame(nmds_Family$points)

# Remove specific samples from the metadata data frame
samples_to_remove <- c("K47", "K48", "K16", "K15", "K17", "K08")
samdf <- subset(samdf, !(rownames(samdf) %in% samples_to_remove))
sample_names_samdf <- rownames(samdf)
nmds_scores_Family$Samples <- samdf$label
nmds_scores_Family$Label <- sample_names_samdf
colnames(nmds_scores_Family) <- c("NMDS1", "NMDS2", "Label", "Samples")

# Factorize the Label variable for color mapping
nmds_scores_Family$Label <- factor(nmds_scores_Family$Label)

# Define a color-blind friendly palette for sample types
color_palette <- c("Plastic" = "tan", "Non-Plastic" = "forestgreen", "Water_Column" = "blue")

# Create the NMDS plot using ggplot and ggforce
# Create a scatter plot with NMDS scores and ellipses for the "Family" data

# Set the data source and initial plot title
plotted_e_Family <- ggplot(data = nmds_scores_Family) +
  ggtitle("Family") +
  
  # Adjust the position of the plot title
  theme(plot.title = element_text(hjust = 0.5)) +
  
  # Add ellipses based on NMDS scores and fill/color according to the 'Label' variable
  geom_mark_ellipse(aes(x = NMDS1, y = NMDS2, fill = Label, color = Label), expand = unit(0.5, "mm")) +
  
  # Add points based on NMDS scores and adjust shape/color based on the 'Label' variable
  geom_point(aes(x = NMDS1, y = NMDS2, shape = Label, color = Label)) +
  
  # Add text labels for points with 'Label' as "Plastic"
  geom_text(data = subset(nmds_scores_Family, Label == "Plastic"), aes(x = NMDS1, y = NMDS2, label = Samples),
            nudge_x = 0.02, nudge_y = -1, size = 4, color = "black") +
  
  # Set the manual fill colors using the defined 'color_palette'
  scale_fill_manual(values = color_palette) +
  
  # Set the manual color scale using the defined 'color_palette'
  scale_color_manual(values = color_palette) +
  
  # Set the manual shape scale based on the 'Label' variable values
  scale_shape_manual(values = c("Plastic" = 16, "Non-Plastic" = 17, "Water_Column" = 18)) +
  
  # Use a minimal theme for the plot
  theme_minimal() +
  
  # Customize theme elements like font sizes for legend and axis labels
  theme(
    legend.text = element_text(size = 11),
    axis.text = element_text(size = 11),
    strip.text = element_text(size = 11)
  ) +
  
  # Customize the legend appearance with a specified number of rows
  guides(fill = guide_legend(nrow = 3))


# Print the NMDS plot
print(plotted_e_Family)
# Perform ANOSIM and PERMANOVA on the entire dataset
ano <- anosim(otu_m_Family, samdf$label, distance = "euclidean", permutations = 99)

# Function to perform ANOSIM and PERMANOVA on subsets
perform_anosim_permanova <- function(data, labels) {
  result_anosim <- anosim(data, labels, distance = "euclidean", permutations = 99)
  result_permanova <- adonis2(formula = data ~ labels, permutations = 99, method = "euclidean")
  return(list(anosim = result_anosim, permanova = result_permanova))}

# Perform ANOSIM and PERMANOVA on different subsets of the data
# Replace subsetX_labels with actual subset labels
subsetX_labels <- c("Label1", "Label2", ...)  
indices_subsetX <- which(samdf$label %in% subsetX_labels)
result_subsetX <- perform_anosim_permanova(otu_m_Family[indices_subsetX, ], samdf$label[indices_subsetX])

# Create an empty data frame to store the p-values
p_values_df <- data.frame(
  Subset = character(),
  Test = character(),
  P_Value = numeric(),
  stringsAsFactors = FALSE)

# List of subset label-result pairs
subset_label_results <- list(
  subset1_labels = result_subset1,
  subset2_labels = result_subset2,
  subset3_labels = result_subset3,
  subset4_labels = result_subset4)

# Iterate through the list of subset label-result pairs
for (subset_label in names(subset_label_results)) {
  result_subset <- subset_label_results[[subset_label]]
  
  # Extract p-values from ANOSIM and PERMANOVA results
  p_value_anosim <- result_subset$anosim$signif
  p_value_permanova <- result_subset$permanova$`Pr(>F)`[1] # Assuming the first row is the overall p-value
  
  # Append the p-values to the data frame
  p_values_df <- rbind(
    p_values_df,
    data.frame(
      Subset = rep(subset_label, 2),
      Test = c("ANOSIM", "PERMANOVA"),
      P_Value = c(p_value_anosim, p_value_permanova),
      stringsAsFactors = FALSE))
}

# Export the p-values data frame as a CSV file
write.csv(p_values_df, "p_values_table_final.csv", row.names = FALSE)





## Begin heat map
# Load the DESeq2 package
library(DESeq2)
packageVersion("DESeq2")

# Factorize the label variable for DESeq2 analysis
sample_data(ps)$label <- as.factor(sample_data(ps)$label)

# Combine Plastic, Non-Plastic, and Water_Column samples for comparison
ps.taxa.sub <- subset_samples(ps, label %in% c("Plastic", "Non-Plastic", "Water_Column"))

# Filter out sparse features with > 90% zeros
ps.taxa.pse.sub <- prune_taxa(rowSums(otu_table(ps.taxa.sub) == 0) < ncol(otu_table(ps.taxa.sub)) * 0.9, ps.taxa.sub)
ps_ds = phyloseq_to_deseq2(ps.taxa.pse.sub, ~ label)

# Use alternative size factor estimator for genes containing a zero in every sample
ds <- estimateSizeFactors(ps_ds, type = "poscounts")
ds = DESeq(ds, test = "Wald", fitType = "parametric")
alpha = 0.05 
res = results(ds, alpha = alpha)
res = res[order(res$padj, na.last = NA), ]
taxa_sig = rownames(res[1:50, ]) # Select bottom 50 taxa with lowest p.adj values

# Create a new phyloseq object with only the significant taxa
ps.taxa.rel <- transform_sample_counts(ps, function(x) x / sum(x) * 100)
ps.taxa.rel.sig <- prune_taxa(taxa_sig, ps.taxa.rel)
ps.taxa.rel.sig <- prune_samples(colnames(otu_table(ps.taxa.pse.sub)), ps.taxa.rel.sig)

# Create a matrix of the significant taxa's relative abundances
matrix <- as.matrix(data.frame(otu_table(ps.taxa.rel.sig)))
rownames(matrix) <- as.character(tax_table(ps.taxa.rel.sig)[, "Family"])
metadata_sub <- data.frame(sample_data(ps.taxa.rel.sig))

# Define annotation colors for columns and rows
annotation_col <- data.frame(
  `Sample Type` = as.factor(metadata_sub$label),
  check.names = FALSE)

rownames(annotation_col) <- rownames(metadata_sub)

annotation_row <- data.frame(
  Order = as.factor(tax_table(ps.taxa.rel.sig)[, "Order"]))

# Make row names unique in the annotation_row data frame
annotation_row$Order <- factor(annotation_row$Order, levels = unique(annotation_row$Order))
rownames(annotation_row) <- make.unique(as.character(annotation_row$Order))

# Define annotation colors using named vectors
num_levels <- length(levels(annotation_row$Order))
Order_col <- sample(colors(10), num_levels)
names(Order_col) <- levels(annotation_row$Order)
ann_colors = list(
  "Sample Type" = c(Plastic = "tan", `Non-Plastic` = "forestgreen", "Water_Column" = "blue"),
  Order = Order_col)

# Load the 'viridis' package for color scale
library(viridis)
packageVersion("viridis")

# Define a custom color scale using the 'viridis' palette
custom_colors <- viridis(length(unique(samdf$label)), option = "viridis", direction = -1)

# Load the 'ComplexHeatmap' package for creating heatmaps
library(ComplexHeatmap)
packageVersion("Complexheatmap")

# Load the 'pheatmap' package for more heatmap features
library(pheatmap)
packageVersion("pheatmap")

# Create the heatmap using ComplexHeatmap
ComplexHeatmap::pheatmap(matrix, 
                         annotation_col = annotation_col, 
                         annotation_row = annotation_row, 
                         annotation_colors = ann_colors,
                         scale = "row",
                         color = custom_colors,
                         heatmap_legend_param = list(title = "Standard Deviations From Mean"))



# Load required libraries
library(dplyr)
packageVersion("dplyr")
library(tidyr)
packageVersion("tidyr")
library(stringr)
packageVersion("stringr")

# Group and calculate median for each Family within each label
ps.melt <- ps.melt %>%
  group_by(label, Family) %>%
  mutate(median = median(Abundance))

# Identify Families with median abundance > 1
keep <- unique(ps.melt$Family[ps.melt$median > 1])
ps.melt$Family[!(ps.melt$Family %in% keep)] <- "< 1%"

# Summarize Family abundance at the Sample, label, and Family level
ps.m.melt_sum <- ps.melt %>%
  group_by(Sample, label, Family) %>%
  summarise(Abundance = sum(Abundance))

# Get the top 6 Families in abundance for each Sample
Family_All_top6 <- ps.m.melt_sum %>%
  group_by(Sample) %>%
  arrange(desc(Abundance)) %>%
  slice_head(n = 6) %>%
  ungroup()

# Filter the top Families for plastic samples
Family_All_Data_subsetP_top6 <- Family_All_top6 %>%
  filter(label == "Plastic")

# Remove the "< 1%" Family from the filtered data
Family_All_Data_subsetP_top5 <- Family_All_Data_subsetP_top6 %>%
  filter(Family != "< 1%")

# Save the top 5 Families data for plastic samples
write.csv(Family_All_Data_subsetP_top5, "Family_All_Data_subsetP_top5", row.names = FALSE)
# Save the top 6 Families data for all samples
write.csv(Family_All_top6, "Family_All_top6", row.names = FALSE)

# Create a bubble plot with viridis color palette
library(ggplot2)
packageVersion("ggplot2")
library(viridis)
packageVersion("viridis")

# Create a scatter plot using ggplot to visualize abundance data for specific families in plastic samples

# Set the data and aesthetic mappings
ggplot(Family_All_Data_subsetP_top5, aes(x = Family, y = Sample, size = Abundance, color = Abundance)) +
  
  # Add points for each data point
  geom_point() +
  
  # Set the size of the points within a specified range
  scale_size(range = c(5, 20)) +
  
  # Set color scale for points based on the 'Abundance' variable
  scale_color_viridis(option = "cividis") +
  
  # Customize plot titles and axis labels
  labs(title = "Prevalence and Abundance of Non-specific Families in Plastics",
       x = "Family", y = "Sample", size = "Abundance", color = "Abundance") +
  
  # Apply a minimal theme to the plot
  theme_minimal() +
  
  # Customize theme elements
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis text
    panel.background = element_rect(fill = "antiquewhite3"),  # Set panel background color
    text = element_text(size = 25),  # Customize text size
    panel.grid = element_blank()  # Remove grid lines
  )





break

# Read in the Taxa Absolute Abundance Table
Taxa_Absolute_Abundance_Table <- read_csv("Taxa Absolute Abundance Table.csv")

# Create the OTU matrix and taxa table
otumat <- Taxa_Absolute_Abundance_Table
taxmat <- data.frame(otumat[2:7])
taxmat <- replace(taxmat, is.na(taxmat), "Unclassified")
otumat <- data.frame(otumat[8:52])
rownames(otumat) <- paste0("OTU", 1:nrow(otumat))
rownames(taxmat) <- paste0("OTU", 1:nrow(taxmat))
otumat <- as.matrix(otumat)
taxmat <- as.matrix(taxmat)

# Create the OTU and TAX objects for the phyloseq object
OTU <- otu_table(otumat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)

# Create the sample data frame as a phyloseq object
material <- list(rep("Non-Plastic", 17), rep("Plastic", 17), rep("Water_Column", 9), rep("PCRNeg", 1), rep("ExtrNeg", 1))
samdf <- data.frame(label = unlist(material))
rownames(samdf) <- colnames(OTU)
ps <- phyloseq(OTU, TAX, sample_data(samdf))

# Filter out the "Chloroplast" order and remove specified samples
ps_filtered <- subset_taxa(ps, Order != "Chloroplast")
families_to_remove <- c("Mitochondria")
ps_filtered <- subset_samples(ps_filtered, !(sample_data(ps_filtered)$label %in% c("PCRNeg", "ExtrNeg", "Water_Column", "Non-Plastic")))

# Identify and remove specific taxa associated with control samples
taxa_to_remove <- taxa_names(ps_filtered)
taxa_to_remove <- intersect(taxa_to_remove, taxa_names(subset_taxa(ps, sample_data(ps)$label %in% c("PCRNeg", "ExtrNeg", "Water_Column", "Non-Plastic"))))
ps_filtered_pruned <- prune_taxa(setdiff(taxa_names(ps_filtered), taxa_to_remove), ps_filtered)

# Remove specified samples from the phyloseq object
samples_to_remove <- c("K47", "K48", "K16", "K15", "K17", "K08")
ps <- subset_samples(ps_filtered_pruned, !sample_names(ps_filtered_pruned) %in% samples_to_remove)

# Transform sample counts to relative percentages
ps.rel = transform_sample_counts(ps, function(x) x/sum(x)*100)

# Agglomerate taxa at the 'Family' level
glom <- tax_glom(ps.rel, taxrank = 'Family', NArm = FALSE)

# Melt the data to long format for plotting
ps.melt <- psmelt(glom)

# Convert 'Family' column to character for adjustment
ps.melt$Family <- as.character(ps.melt$Family)

# Group the melted data by 'label' and 'Family', and calculate median Abundance
ps.melt <- ps.melt %>%
  group_by(label, Family) %>%
  mutate(median = median(Abundance))

# Select families with median Abundance > 1
keep <- unique(ps.melt$Family[ps.melt$median > 1])
ps.melt$Family[!(ps.melt$Family %in% keep)] <- "< 1%"

# Group melted data by 'Sample', 'label', and 'Family', and calculate summed Abundance
ps.melt_sum <- ps.melt %>%
  group_by(Sample, label, Family) %>%
  summarise(Abundance = sum(Abundance))

# Create a grouped bar plot using ggplot
# Create a bar plot using ggplot to visualize abundance data across samples

# Set the data and aesthetic mappings
ggplot(ps.melt_sum, aes(x = Sample, y = Abundance, fill = Family)) + 
  
  # Add bars with heights corresponding to the 'Abundance' values
  geom_bar(stat = "identity") + 
  
  # Set color scale for fill based on the 'Family' variable
  scale_fill_viridis_d(option = "plasma", begin = 0, end = 0.8) +
  
  # Customize axis labels
  labs(x = "", y = "%") +
  
  # Facet the plot by the 'label' variable, arranging facets in a single row
  facet_wrap(~label, scales = "free_x", nrow = 1) +
  
  # Apply a classic theme to the plot
  theme_classic() + 
  
  # Customize theme elements
  theme(
    strip.background = element_blank(),  # Remove strip background
    axis.text.x.bottom = element_text(angle = -90),  # Rotate x-axis text
    legend.text = element_text(size = 15),  # Customize legend text size
    axis.text = element_text(size = 13),  # Customize axis text size
    strip.text = element_text(size = 13)  # Customize strip text size
  ) +
  
  # Customize the legend appearance
  guides(fill = guide_legend(nrow = 24))  # Set number of rows in the legend




library(dplyr)  # Load the dplyr package for data manipulation
library(tidyr)  # Load the tidyr package for data reshaping
library(stringr)  # Load the stringr package for string manipulation

# Group the melted data by 'label' and 'Family', and calculate median Abundance
ps.melt <- ps.melt %>%
  group_by(label, Family) %>%
  mutate(median = median(Abundance))

# Identify unique families with median Abundance > 1
keep <- unique(ps.melt$Family[ps.melt$median > 1])

# Set families with median Abundance <= 1 to "< 1%"
ps.melt$Family[!(ps.melt$Family %in% keep)] <- "< 1%"

# Group melted data by 'Sample', 'label', and 'Family', and calculate summed Abundance
ps.m.melt_sum <- ps.melt %>%
  group_by(Sample, label, Family) %>%
  summarise(Abundance = sum(Abundance))

# Group samples by 'Sample', arrange in descending Abundance, and keep top 6
Family_Plastics_top6 <- ps.m.melt_sum %>%
  group_by(Sample) %>%
  arrange(desc(Abundance)) %>%
  slice_head(n = 6) %>%
  ungroup()

# Filter 'Family_All_top6' to include only 'Plastic' labeled samples
Family_All_Data_subsetP_top6 <- Family_All_top6 %>%
  filter(label == "Plastic")

# Filter 'Family_All_Data_subsetP_top6' to exclude families with "< 1%" Abundance
Family_PlasticSpecific_top5 <- Family_All_Data_subsetP_top6 %>%
  filter(Family != "< 1%")

# Write the filtered data to CSV files
write.csv(Family_PlasticSpecific_top5, "Family_PlasticSpecific_top5", row.names = FALSE)




# Create a bubble plot using ggplot with specified aesthetics
ggplot(Family_PlasticSpecific_top5, aes(x = Family, y = Sample, size = Abundance, color = Abundance)) +
  
  # Add points to the plot
  geom_point() +
  
  # Set the range of point sizes
  scale_size(range = c(5, 20)) +
  
  # Set the color scale using the 'viridis' palette
  scale_color_viridis(option = "cividis") +
  
  # Set plot labels and titles
  labs(title = "Prevalence and Abundance of Plastic-specific Families",
       x = "Family", y = "Sample", size = "Abundance", color = "Abundance") +
  
  # Set the minimal theme
  theme_minimal() +
  
  # Customize theme settings
  theme(
    axis.text.x = element_text(angle = 20, hjust = 1),  # Rotate x-axis text
    panel.background = element_rect(fill = "antiquewhite3"),  # Set panel background color
    text = element_text(size = 25),  # Set text size
    panel.grid = element_blank())  # Remove grid lines



knitr::write_bib(c(.packages(), "bookdown"), "packages.bib")
# Generate a BibTeX bibliography file containing package citations

# Load the knitr package to use the write_bib function
library(knitr)

# Get a vector of names of all currently loaded packages
loaded_packages <- .packages()

# Add "bookdown" to the vector of loaded packages
loaded_packages <- c(loaded_packages, "bookdown")

# Specify the output BibTeX file name
output_file <- "packages.bib"

# Use the write_bib function to create the BibTeX file
write_bib(packages = loaded_packages, file = output_file)


