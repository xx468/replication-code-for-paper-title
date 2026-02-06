# replication-code-for-paper-title
Code and data to reproduce the findings of Oxygen Availability Shapes Eukaryotic Plankton Assembly and Coalescences in River-Lake Ecotones

Corresponding Paper: Oxygen Availability Shapes Eukaryotic Plankton Assembly and Coalescences in River-Lake Ecotones
Authors:Fuchao Zheng a, b, Wendian Zhu a, Jinghan Zhang d, Shenglai Yin e, Qinghui You f, Qiwu Hu a, Chaoyang Fang a, Yong Ge a, *, Kunlin Yang b, Guofang Xu c, *
Journal/Status: Ecological Informatics

Brief Description: This repository contains the complete code and data used to generate all analyses, models, figures, and statistical results in the paper " Oxygen Availability Shapes Eukaryotic Plankton Assembly and Coalescences in River-Lake Ecotones". This study primarily investigates [Describe the core research content in 1-2 sentences, e.g., the impact of climate change on the distribution of a certain species].


├── README.md                   # This file
├── LICENSE                     # Open-source license (MIT License)
├── .gitignore                  # Files ignored by Git
│
├── data/                       # Data directory
│   ├── processed/              # **Analysis-ready data** (directly used by scripts)
│   │   ├── DO_data.csv         # Main dataset for analysis
│   │   ├── zoo_diversity indice.txt    # Diverse data
│   │   └── phy_diversity.txt           # Diverse data
│   └── raw/                    # **Raw data** (see instructions below for access)
│       ├── phy_OTU.txt         # OTU data
│       ├── zoo_OTU.txt         # OTU data 
│       ├── phy_taxa.csv        # tax data
│       └── zoo_taxa.csv        # tax data


# The code used in the article

# Single Factor Bar Plot Code
# Read local data
dt <- read.csv("DO_data.csv", header = TRUE)

# View first 12 rows of data
head(dt, 12)

# Install ggbeeswarm package
# install.packages("ggbeeswarm")

# Load ggbeeswarm and ggplot2 packages
library(ggbeeswarm)
library(ggplot2)

# Define custom plot theme for fine adjustments
top.mar <- 0.2
right.mar <- 0.2
bottom.mar <- 0.2
left.mar <- 0.2

# Custom theme 1
mytheme1 <- theme(
  panel.background = element_blank(),
  axis.ticks.length = unit(1.6, "mm"),
  plot.margin = unit(x = c(top.mar, right.mar, bottom.mar, left.mar), units = "inches")
)

# Custom theme 2
# Hide vertical axis and define font style, axis thickness, color, and tick length
mytheme2 <- theme_bw() +
  theme(
    text = element_text(family = "serif", colour = "black", size = 12),
    axis.line = element_line(size = 0.8, colour = "black"),
    axis.ticks = element_line(size = 0.8, colour = "gray10"),
    axis.ticks.length = unit(1.5, units = "mm"),
    plot.margin = unit(x = c(top.mar, right.mar, bottom.mar, left.mar), units = "inches")
  )

# No need to calculate mean and standard deviation in advance, directly create bar plot
dt$Group <- factor(dt$Group, levels = c("DO6", "DO8", "DO10"))

p1 <- ggplot(data = dt, mapping = aes(x = Group, y = DO, color = Group, fill = Group)) +
  scale_y_continuous(breaks = seq(0, 13, 4), limits = c(0, 13)) +
  stat_summary(fun = mean, geom = "bar", alpha = 0.0, width = 0.54, color = "black", show.legend = FALSE)

# Add error bars
# mult is multiplier, here we shrink error bars to half of their original size
p2 <- p1 + 
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1/2),
               geom = "errorbar", width = 0.2, color = "black", show.legend = FALSE)  # Error bars forced to black

# Add scatter points on bar plot
p3 <- p2 +
  geom_jitter(shape = 16, alpha = 0.8, size = 2, width = 0.15, show.legend = FALSE)

# Apply custom theme 1
p3 +
  scale_color_manual(values = c("DO6" = "#c1272d", "DO8" = "#0000a7", "DO10" = "#008176")) +
  scale_y_continuous(expand = expansion(add = c(0.1, 0.3))) +
  theme_bw() +
  xlab("") + 
  ylab("Dissolved Oxygen (mg/L)") +
  theme(text = element_text(family = "serif", colour = "black", size = 15))


# Single Factor Beeswarm Plot Code
rm(list = ls())
setwd("C:/Users/Think/Desktop/Data Analysis/Alpha Diversity/Update")

dt <- read.csv('Shannon.csv', row.names = 1, check.names = FALSE)

# Beeswarm plot
library(ggplot2)
library(ggbeeswarm)

dt$time <- factor(dt$time, levels = c("DO6", "DO8", "DO10"))

# Modified beeswarm plot code
p3 <- ggplot(dt, aes(x = site, y = w, color = time)) +
  geom_boxplot(position = position_dodge(0.9), width = 0.6, outlier.shape = NA) +
  
  geom_beeswarm(dodge.width = 0.9, size = 2, alpha = 1) +
  
  # Optionally add mean points and error bars
  # stat_summary(fun = mean, geom = "point", shape = 18, size = 4,
  #              position = position_dodge(0.9), color = "black") +
  # stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3,
  #              position = position_dodge(0.9), color = "black") +
  
  # Optionally add significance letters
  # geom_text(aes(label = Tukey, y = w + se), family = "serif", size = 0,
  #           vjust = -0.5, position = position_dodge(0.9)) +
  
  labs(x = "", y = "Shannon Index", color = "black") +
  ggpubr::theme_pubr(legend = c("right")) +  # Legend position
  scale_color_manual(values = c("#c1272d", "#0000a7", "#008176")) +  # Using color instead of fill
  theme(text = element_text(family = "serif", colour = "black", size = 15))
p3


### Linear Regression Code ###
library(ggplot2)
library(ggpmisc)
library(MASS)

data <- read.csv("Shannon_phy.csv")
str(data)

# Linear regression key values
# Check normality
hist(data$Shannon, probability = TRUE)
shapiro.test(data$Shannon)

# Linear regression
lm_regression <- lm(Shannon ~ DO, data = data)
summary(lm_regression)

# Quadratic regression
data$DO2 <- data$DO^2
quadraticModel <- lm(Shannon ~ DO + DO2, data = data)
summary(quadraticModel)

# Plot
data$treatment <- factor(data$treatment, levels = c("DO6", "DO8", "DO10"))

plotsample <- ggplot(data, aes(x = DO, y = Shannon)) +
  geom_point(size = 5, aes(color = factor(treatment)), alpha = 0.95) +
  scale_colour_manual(values = c("DO6" = "#c1272d",
                               "DO8" = "#0000a7",
                               "DO10" = "#008176")) +
  theme_bw() +
  geom_smooth(formula = y ~ x,
              color = 'black', fill = 'gray40',
              method = "lm", se = TRUE, fullrange = FALSE) +
  labs(x = "DO", y = "Shannon", fill = "") +
  theme(text = element_text(family = "serif")) +
  scale_x_continuous(limits = c(6, 13)) +
  scale_y_continuous(limits = c(0, 400)) +
  annotate(geom = 'text', x = 12, y = 8, label = "")
Plotsample

### PCoA Code ###
library(vegan)
library(ggplot2)
library(ggrepel)

otu <- read.delim('phy_OTU.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
otu <- data.frame(t(otu))

bray_dis <- vegdist(otu, method = 'bray')
bray <- as.matrix(bray_dis)

pcoa <- cmdscale(bray_dis, k = (nrow(otu) - 1), eig = TRUE)
pcoa_eig <- pcoa$eig
pcoa_exp <- pcoa_eig/sum(pcoa_eig)`

site <- scores(pcoa)
pcoa1 <- paste('PCoA1 (', round(100 * pcoa_exp[1], 2), '%)')
pcoa2 <- paste('PCoA2 (', round(100 * pcoa_exp[2], 2), '%)')

site <- data.frame(pcoa$point)[1:2]
site$name <- rownames(site)
site$group <- c(rep('DO6', 18), rep('DO8', 18), rep('DO10', 18))
site$group <- factor(site$group, levels = c("DO6", "DO8", "DO10"))

library(ggpubr)
p <- ggplot(data = site, aes(X1, X2)) +
  geom_point(aes(color = group), size = 2) +
  geom_text_repel(aes(X1, X2, label = '')) +
  stat_ellipse(aes(fill = group), geom = 'polygon', level = 0.95, alpha = 0.1, show.legend = FALSE) +
  scale_color_manual(values = c("#c1272d", "#0000a7", "#008176")) +
  scale_fill_manual(values = c("#c1272d", "#0000a7", "#008176")) +
  theme(panel.background = element_rect(color = 'black', fill = 'transparent'),
        plot.title = element_text(hjust = 0.5), legend.position = 'right') +
  theme(text = element_text(size = 12, family = "serif")) +
  geom_vline(xintercept = 0, color = 'gray', size = 0.5) +
  geom_hline(yintercept = 0, color = 'gray', size = 0.5) +
  labs(x = pcoa1, y = pcoa2, title = '') +
  theme(panel.grid = element_blank())
P


### ANOSIM Analysis Code ###
library(vegan)
library(ggplot2)
library(ggsignif)

data <- read.delim("zoo_OTU.txt", row.names = 1)
data <- data.frame(t(data))
data$group <- c(rep('DO6', 18), rep('DO8', 18), rep('DO10', 18))
data$group <- factor(data$group, levels = c("DO6", "DO8", "DO10"))

anosim_result <- anosim(data[, 1:124], data$group, distance = "bray", permutations = 999)
summary(anosim_result)

anosim <- data.frame(vec = anosim_result$class.vec, rank = anosim_result$dis.rank)
p_val <- anosim_result$signif
p_val <- p.adjust(p_val, method = "BH")
p_label <- paste("p-value:", p_val, "")
r_val <- anosim_result$statistic

compaired = list(c("DO6", "DO8"), c("DO8", "DO10"), c("DO6", "DO10"))

p1 <- ggplot(anosim) + 
  geom_boxplot(aes(x = vec, y = rank, fill = vec)) +
  scale_fill_manual(values = c("#c1272d", "#0000a7", "#008176")) +
  geom_signif(aes(x = vec, y = rank), comparisons = compaired, family = "serif",
              step_increase = 0.1, y_position = c(800, 850, 900),
              map_signif_level = TRUE) +
  theme_bw() +
  labs(y = "Bray-Curtis ANOSIM") +
  scale_y_continuous(limits = c(0, 1000)) +
  annotate("text", label = p_label, x = 1.2, y = 750, family = "serif") +
  annotate("text", label = sprintf("italic(R)==%.4f", r_val), x = 1.2, y = 850, size = 4, parse = TRUE, family = "serif") +
  theme(panel.background = element_blank(), legend.background = element_blank(), legend.position = "none") +
  theme(legend.title = element_blank(), axis.title.x = element_blank(),
        strip.background = element_rect(fill = "black", size = 0.1),
        strip.text = element_text(colour = "white")) +
  theme(text = element_text(family = "serif"))
p1

### Sankey Diagram Code ###
library(reshape2)
library(dplyr)
library(ggplot2)
library(ggalluvial)
library(ggsci)

# Import data
otu_abundance <- read.table(file = "phy_OTU.txt", sep = "\t", header = TRUE, check.names = FALSE)
sample_group <- read.table(file = "DO_grouping.txt", sep = "\t", header = TRUE, check.names = FALSE)
taxonomy <- read.table(file = "Phy_tax.txt", sep = "\t", header = TRUE, check.names = FALSE)

# Convert OTU abundance data to long format and merge data
otu_long <- melt(otu_abundance, id.vars = "OTU_id", variable.name = "SampleID", value.name = "Abundance") %>%
  merge(taxonomy, by = "OTU_id") %>%
  merge(sample_group, by.x = "SampleID", by.y = "samples")

# Convert DO groups to numeric values
otu_long$DO_value <- as.numeric(gsub("DO", "", otu_long$Group))

# Function to calculate Spearman correlation and filter significant taxa
get_significant_taxa <- function(data, taxa_level) {
  results <- data.frame()
  for(taxa in unique(data[[taxa_level]])) {
    subset_data <- data[data[[taxa_level]] == taxa, ]
    if(nrow(subset_data) > 1) {
      cor_test <- cor.test(subset_data$Abundance, subset_data$DO_value, 
                           method = "spearman", exact = FALSE)
      if(!is.na(cor_test$p.value) && cor_test$p.value < 0.05) {
        results <- rbind(results, data.frame(
          Taxa = taxa,
          Level = taxa_level,
          Correlation = cor_test$estimate,
          Pvalue = cor_test$p.value
        ))
      }
    }
  }
  return(results)
}

# Calculate significant correlated taxa at Phylum and Genus levels
phylum_cor <- get_significant_taxa(otu_long, "Phylum") %>% 
  mutate(Correlation_direction = ifelse(Correlation > 0, "Positive", "Negative"))

genus_cor <- get_significant_taxa(otu_long, "Genus") %>% 
  mutate(Correlation_direction = ifelse(Correlation > 0, "Positive", "Negative"))

# Merge results and remove duplicates
taxa_colors <- rbind(
  phylum_cor %>% select(Taxa, Correlation_direction, Level),
  genus_cor %>% select(Taxa, Correlation_direction, Level)
) %>% distinct()

# Get names of significant Phyla and Genera
significant_phyla <- unique(phylum_cor$Taxa)
significant_genera <- unique(genus_cor$Taxa)

# Filter otu_long data, keeping only significant Phyla AND Genera
otu_long_filtered <- otu_long %>%
  filter(Phylum %in% significant_phyla & Genus %in% significant_genera)

otu_long_filtered$Group <- factor(otu_long_filtered$Group, levels = c("DO6", "DO8", "DO10"))

p1 <- ggplot(data = otu_long_filtered,
             aes(axis1 = Group, axis2 = Phylum, axis3 = Genus, y = NULL)) +
  scale_x_discrete(limits = c("Group", "Phylum", "Genus"), expand = c(.1, .2)) +
  scale_y_continuous(expand = c(0, 0)) +
  geom_alluvium(aes(fill = Phylum)) +
  geom_stratum(width = 0.35) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 5) +
  scale_fill_npg() +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 15, color = "black"),
        legend.key.size = unit(1.5, 'lines'),
        legend.text = element_text(size = 20)) +
  labs(title = "", x = NULL, y = NULL)
p1


### Niche Breadth Code ###
library(vegan)

env <- read.table('env1.txt', header = TRUE)
env.stand <- decostand(env, 'standardize', MARGIN = 2)
env.dm <- vegdist(env.stand, method = 'euclidean')
env.KM <- cascadeKM(env.dm, inf.gr = 2, sup.gr = 20, iter = 100, criterion = 'ssi')
plot(env.KM, sortg = TRUE, gridcol = 'darkgrey')
write.csv(env.KM$partition, 'env_KM_partition.csv')

RS <- read.csv('phy_grouping.csv', row.names = 1)
rowSumRS <- rowSums(RS)
a <- RS / rowSumRS
a2 <- a * a
rowSuma2 <- rowSums(a2)
niche <- 1 / rowSuma2
niche.A <- (niche - 1) / (8 - 1)
write.csv(niche.A, 'species_niche_A.csv')

otu <- read.table('DO6_phy.txt', row.names = 1)
otu[otu > 0] <- 1
niche <- read.csv('species_niche_A.csv', row.names = 1)
niche <- as.matrix(niche)
niche.matrix <- otu * niche
Bcom <- colSums(niche.matrix) / colSums(otu)
write.csv(Bcom, 'Bcom.csv')



### Niche Overlap Code ###
library(spaa)

otu <- read.table('DO_phy_zoo_combined.txt', row.names = 1)
niche_overlap <- niche.overlap(otu, method = 'levins')
niche_overlap <- as.matrix(niche_overlap)
write.table(niche_overlap, 'niche_overlap.txt', sep = '\t', col.names = NA, quote = FALSE)

# Bootstrap confidence intervals
set.seed(123)
niche_overlap_boot <- niche.overlap.boot(otu, method = 'levins', times = 1000, quant = c(0.025, 0.975))

# Heatmap
library(ggplot2)
library(corrplot)
M = cor(niche_overlap)
corrplot.mixed(M, order = 'AOE')
corrplot(M, tl.col = "black")


### Community Assembly Code ###
# Null model
library(NST)
library(vegan)

# Read files
index <- read.delim('zoo_transposed.txt', header = TRUE, row.names = 1)
groups <- read.delim('DO_grouping.txt', header = TRUE, row.names = 1)

# Calculate NST
set.seed(123)
tnst <- tNST(comm = index, group = groups, dist.method = 'bray', null.model = 'PF', 
             rand = 1000, nworker = 1)

# View NST for within-group sample pairs
nst_groups <- tnst$index.pair.grp
nst_groups
pmax(0, pmin(1, nst_groups))

# Output main statistical results
write.table(nst_groups, 'nst_group.txt', sep = '\t', row.names = FALSE, quote = FALSE)


# Neutral Community Model
library(Hmisc)
library(minpack.lm)
library(stats4)

# Load species abundance data
# spp: Species or taxon abundance table, rows represent taxa, columns represent samples
spp <- read.csv('DO10_zoo.txt', header = TRUE, stringsAsFactors = FALSE, row.names = 1, sep = "\t")
spp <- t(spp)

## Fit the neutral model by Sloan et al. (2006) to a species abundance table
## and return several fit statistics. Alternatively, return the predicted 
## occurrence frequencies of each taxon based on their abundance in the metacommunity.
## Non-linear least squares (NLS) is used to fit the model parameters.

N <- mean(apply(spp, 1, sum))  # Calculate mean total abundance per sample
p.m <- apply(spp, 2, mean)     # Calculate mean relative abundance for each species
p.m <- p.m[p.m != 0]           # Remove zero-abundance species
p <- p.m / N                   # Calculate normalized relative abundance
spp.bi <- 1 * (spp > 0)        # Convert to presence-absence matrix
freq <- apply(spp.bi, 2, mean) # Calculate occurrence frequency for each species
freq <- freq[freq != 0]        # Remove species with zero frequency

# Merge abundance and frequency data
C <- merge(p, freq, by = 0)
C <- C[order(C[, 2]),]         # Sort by abundance
C <- as.data.frame(C)
C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),]  # Remove rows with zeros
p <- C.0[, 2]                  # Extract relative abundance
freq <- C.0[, 3]               # Extract occurrence frequency
names(p) <- C.0[, 1]           # Assign names
names(freq) <- C.0[, 1]

d <- 1 / N  # Calculate immigration parameter

# Fit neutral model using non-linear least squares
m.fit <- nlsLM(freq ~ pbeta(d, N * m * p, N * m * (1 - p), lower.tail = FALSE), 
               start = list(m = 0.1))
m.fit  # Display model fit results

# Calculate confidence intervals for m
m.ci <- confint(m.fit, 'm', level = 0.95)

# Predict occurrence frequencies
freq.pred <- pbeta(d, N * coef(m.fit) * p, N * coef(m.fit) * (1 - p), lower.tail = FALSE)

# Calculate confidence intervals using Wilson method
pred.ci <- binconf(freq.pred * nrow(spp), nrow(spp), alpha = 0.05, 
                   method = "wilson", return.df = TRUE)

# Calculate R-squared
R2 <- 1 - (sum((freq - freq.pred) ^ 2)) / (sum((freq - mean(freq)) ^ 2))
R2  # Display R-squared value

# Export three statistical tables: mean relative abundance (p.csv), 
# observed occurrence frequency (freq.csv), and predicted occurrence frequency (freq.pred.csv)
# write.csv(p, file = "p.csv")
# write.csv(freq, file = "freq.csv")
# write.csv(freq.pred, file = "freq.pred.csv")

# p: mean relative abundance
# freq: observed occurrence frequency  
# freq.pred: predicted occurrence frequency (neutral model fit)

# Create visualization
bacnlsALL <- data.frame(p, freq, freq.pred, pred.ci[, 2:3])

# Color code points based on their position relative to neutral prediction
inter.col <- rep('#0000a7', nrow(bacnlsALL))
inter.col[bacnlsALL$freq <= bacnlsALL$Lower] <- '#008176'  # Below neutral prediction
inter.col[bacnlsALL$freq >= bacnlsALL$Upper] <- '#c1272d'  # Above neutral prediction

library(grid)
grid.newpage()
pushViewport(viewport(h = 0.6, w = 0.6))
pushViewport(dataViewport(xData = range(log10(bacnlsALL$p)), 
                         yData = c(0, 1.02), extension = c(0.02, 0)))
grid.rect()
grid.points(log10(bacnlsALL$p), bacnlsALL$freq, pch = 20, 
           gp = gpar(col = inter.col, cex = 0.7))
grid.yaxis()
grid.xaxis()

# Add prediction lines
grid.lines(log10(bacnlsALL$p), bacnlsALL$freq.pred, 
          gp = gpar(col = 'black', lwd = 2), default = 'native')
grid.lines(log10(bacnlsALL$p), bacnlsALL$Lower, 
          gp = gpar(col = 'black', lwd = 2, lty = 2), default = 'native') 
grid.lines(log10(bacnlsALL$p), bacnlsALL$Upper, 
          gp = gpar(col = 'black', lwd = 2, lty = 2), default = 'native')  

# Add axis labels
grid.text(y = unit(0, 'npc') - unit(2.5, 'lines'), 
         label = 'Mean Relative Abundance (log10)', gp = gpar(fontface = 2)) 
grid.text(x = unit(0, 'npc') - unit(3, 'lines'), 
         label = 'Frequency of Occurrence', gp = gpar(fontface = 2), rot = 90) 

# Function to add model statistics
draw.text <- function(just, i, j) {
  grid.text(paste0("R\u00B2 = ", round(R2, 3), "\n", "Nm = ", round(coef(m.fit) * N)), 
           x = x[j], y = y[i], just = just)
}

x <- unit(1:4 / 5, "npc")
y <- unit(1:4 / 5, "npc")
draw.text(c("centre", "bottom"), 1, 4)  # Position of statistics label


### Co-occurrence Network and Keystone Species Code ###

# Read OTU data
otu_table <- read.delim('DO10_phy.txt', row.names = 1, check.names = FALSE)

# Calculate relative abundance
rel_abundance <- apply(otu_table, 2, function(x) x/sum(x))
rel_abundance

# Export relative abundance table
# write.csv(rel_abundance, 'non_GC_rel_abundance.csv')

# Filter OTUs based on relative abundance
mean_rel_abundance <- rowMeans(rel_abundance)
low_rel_abundance_otu <- rownames(otu_table)[mean_rel_abundance < 0.0001]  # OTUs with mean relative abundance < 0.01%
otu_table_filtered <- otu_table[!(rownames(otu_table) %in% low_rel_abundance_otu), ]  # Remove low abundance OTUs

# Filter OTUs based on occurrence frequency
freq <- apply(otu_table_filtered, 1, function(x) sum(x > 0)/length(x))
keep <- freq >= 1/5  # Adjust frequency threshold as needed
otu_table_filt <- otu_table_filtered[keep, ]  # Keep only OTUs with frequency above threshold

# Export filtered OTU table for network analysis
write.csv(otu_table_filt, 'non_GC_otu_table_filt.csv')
otu <- otu_table_filt

# Load required libraries
library(WGCNA)
library(psych)
library(reshape2)
library(igraph)

# Calculate Spearman correlation coefficients and p-values between OTUs
cor <- corAndPvalue(t(otu), y = NULL, use = "pairwise.complete.obs", 
                   alternative = 'two.sided', method = 'spearman')
r <- cor$cor  # Correlation coefficients
p <- cor$p    # P-values

# Adjust p-values using Benjamini-Hochberg method
p <- p.adjust(p, method = 'BH')

# Filter correlations: set to 0 if p > 0.001 or |r| < 0.60
r[p > 0.001 | abs(r) < 0.60] <- 0

# Write correlation matrix to CSV
write.csv(data.frame(r, check.names = FALSE), 'corr_matrix.csv')

# Create weighted undirected graph from correlation matrix
g <- graph_from_adjacency_matrix(r, mode = "undirected", weighted = TRUE, diag = FALSE)
g

# Remove isolated nodes (degree = 0)
g <- delete_vertices(g, names(degree(g)[degree(g) == 0]))

# Set edge attributes
E(g)$correlation <- E(g)$weight  # Correlation values
E(g)$weight <- abs(E(g)$weight)  # Absolute values for weight

# Read taxonomic information
tax <- read.delim('phy_zoo_grouping.txt', row.names = 1, check.names = FALSE)

# Add taxonomic information to vertices
tax <- tax[as.character(V(g)$name), ]
V(g)$Kingdom <- tax$Kingdom
V(g)$Phylum <- tax$Phylum
V(g)$Class <- tax$Class
V(g)$Order <- tax$Order
V(g)$Family <- tax$Family
V(g)$Genus <- tax$Genus
V(g)$Species <- tax$Species

# Create node list
node_list <- data.frame(
  label = names(V(g)),
  kingdom = V(g)$Kingdom,
  phylum = V(g)$Phylum,
  class = V(g)$Class,
  order = V(g)$Order,
  family = V(g)$Family,
  genus = V(g)$Genus,
  species = V(g)$Species
)
head(node_list)
write.csv(node_list, 'network_node_list.csv')

# Create edge list
edge <- data.frame(as_edgelist(g))
edge_list <- data.frame(
  source = edge[[1]],
  target = edge[[2]],
  weight = E(g)$weight,
  correlation = E(g)$correlation
)
head(edge_list)
write.csv(edge_list, 'network_edge_list.csv')

# Export graph for visualization in Gephi
write_graph(g, 'network.graphml', format = 'graphml')

###### Calculate Common Network Topological Parameters #####
nodes_num <- length(V(g))                    # Number of nodes
edges_num <- length(E(g))                    # Number of edges
positive_cor_num <- sum(E(g)$correlation > 0)  # Number of positive correlations
negative_cor_num <- sum(E(g)$correlation < 0)  # Number of negative correlations
average_degree <- mean(degree(g))            # Average degree
average_path_length <- mean_distance(g, directed = FALSE)  # Average path length
network_diameter <- diameter(g, directed = FALSE)          # Network diameter
network_density <- edge_density(g)                        # Network density
clustering_coefficient <- transitivity(g)                 # Clustering coefficient

# Compile network parameters
network_parameter <- data.frame(
  nodes_num, 
  edges_num, 
  positive_cor_num, 
  negative_cor_num, 
  average_degree,
  average_path_length,
  network_diameter, 
  network_density,
  clustering_coefficient
)
network_parameter
write.csv(network_parameter, 'network_parameter.csv')

# Create adjacency matrix
otu1 <- otu
otu1[otu1 > 0] <- 1
write.csv(otu1, 'adjacency_matrix.csv')

### Calculate Zi-Pi Values and Plot ###
library(igraph)
library(ggplot2)

# Read graph and detect communities
igraph <- read_graph("network.graphml", format = "graphml")
wtc <- cluster_louvain(igraph)
modularity(wtc)
V(igraph)$module <- membership(wtc)

# Within-module degree z-score function
within_module_deg_z_score <- function(g, A = NULL, weighted = FALSE) {
  stopifnot(is_igraph(g))
  if (is.null(A)) {
    if (isTRUE(weighted)) {
      A <- as_adj(g, sparse = FALSE, names = TRUE, attr = 'weight')
    } else {
      A <- as_adj(g, sparse = FALSE, names = TRUE)
    }
  }
  memb <- vertex_attr(g, "module")
  N <- max(memb)
  nS <- tabulate(memb)
  z <- Ki <- rep.int(0, dim(A)[1L])
  Ksi <- sigKsi <- rep.int(0, N)
  names(z) <- names(Ki) <- rownames(A)
  for (S in seq_len(N)) {
    x <- rowSums(A[memb == S, memb == S, drop = FALSE])
    Ki[memb == S] <- x
    Ksi[S] <- sum(x) / nS[S]
    sigKsi[S] <- sqrt(sum((x - Ksi[S])^2) / (nS[S] - 1))
  }
  z <- (Ki - Ksi[memb]) / sigKsi[memb]
  z[is.infinite(z)] <- 0
  df <- data.frame(Ki, z, row.names = names(Ki))
  return(df)
}

# Participation coefficient function
part_coeff <- function(g, A = NULL, weighted = FALSE) {
  stopifnot(is_igraph(g))
  if (is.null(A)) {
    if (isTRUE(weighted)) {
      A <- as_adj(g, sparse = FALSE, attr = 'weight')
    } else {
      A <- as_adj(g, sparse = FALSE)
    }
  }
  memb <- vertex_attr(g, "module")
  Ki <- colSums(A)
  N <- max(memb)
  Kis <- t(rowsum(A, memb))
  pi <- 1 - ((1 / Ki^2) * rowSums(Kis^2))
  names(pi) <- rownames(A)
  return(pi)
}

# Calculate Zi and Pi values
zi <- within_module_deg_z_score(igraph)
pi <- part_coeff(igraph)
zi_pi <- data.frame(zi, pi)
zi_pi <- na.omit(zi_pi)

# Classify nodes based on Zi-Pi thresholds
zi_pi[which(zi_pi$z < 2.5 & zi_pi$pi < 0.62), 'type'] <- 'Peripheral'
zi_pi[which(zi_pi$z < 2.5 & zi_pi$pi > 0.62), 'type'] <- 'Connectors'
zi_pi[which(zi_pi$z > 2.5 & zi_pi$pi < 0.62), 'type'] <- 'Module hubs'
zi_pi[which(zi_pi$z > 2.5 & zi_pi$pi > 0.62), 'type'] <- 'Network hubs'

write.csv(zi_pi, file = "DO10_zi_pi.csv")

## Zi-Pi Plot
ggplot(zi_pi, aes(pi, z)) +
  geom_point(aes(color = type), alpha = 0.5, size = 4) +
  scale_color_manual(values = c("#15B59D", "#7B3396", "#F8C58E", '#eecc16'),
                     limits = c('Peripheral', 'Connectors', 'Module hubs', 'Network hubs')) +
  theme(panel.grid = element_blank(), 
        axis.line = element_line(colour = 'black'),
        panel.background = element_blank(), 
        legend.key = element_blank()) +
  labs(x = 'Participation Coefficient (pi)', y = 'Within-module Degree (zi)', color = '') +
  geom_vline(xintercept = 0.62) +
  geom_hline(yintercept = 2.5)
