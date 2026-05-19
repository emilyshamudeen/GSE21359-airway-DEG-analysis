# Load libraries 

library(GEOquery)
library(dplyr)
library(stringr)
library(limma)
library(pheatmap)
library(ggplot2)
library(ggrepel)

# Import the data 

gse <- getGEO("GSE21359")[[1]]

pData(gse)   ## sample information
fData(gse)   ## gene annotation
exprs(gse)   ## expression matrix

# Check normalization and scale of expression values

# Check what processing has already been applied
pData(gse)$data_processing[1]

# Summarise expression values — RMA data should be log2 (range 0–16)
summary(exprs(gse))

# Values look fine — data is already RMA normalized and log2 transformed
# If values exceeded 16, log2 transform here: exprs(gse) <- log2(exprs(gse))

# Boxplot to confirm distributions are similar across samples (normalized)
boxplot(exprs(gse), outline = FALSE, main = "Expression distribution per sample")

# Inspect the clinical variables

sampleInfo <- pData(gse)
head(sampleInfo)

# Titles follow the pattern: "small airways, <type> <id>"
# where type is one of: "non-smoker", "smoker", "COPD"
table(sampleInfo$title)

# Build group variable
# ** non-smoker must be checked before smoker because "non-smoker"
# contains the substring "smoker", wrong order misclassifies 53 samples

sampleInfo$group <- ""

for(i in 1:nrow(sampleInfo)){
  
  if(str_detect(sampleInfo$title[i], regex("non-smoker", ignore_case = TRUE))){
    sampleInfo$group[i] <- "non_smoker"
  }
  
  if(str_detect(sampleInfo$title[i], regex("COPD", ignore_case = TRUE))){
    sampleInfo$group[i] <- "copd_smoker"
  }
  
  if(str_detect(sampleInfo$title[i], regex("smoker", ignore_case = TRUE)) &&
     !str_detect(sampleInfo$title[i], regex("non-smoker", ignore_case = TRUE)) &&
     !str_detect(sampleInfo$title[i], regex("COPD", ignore_case = TRUE))){
    sampleInfo$group[i] <- "healthy_smoker"
  }
}

table(sampleInfo$group)

# Extract sex and age as covariables
sampleInfo$sex <- ""

for(i in 1:nrow(sampleInfo)){
  
  if(str_detect(sampleInfo$`characteristics_ch1.1`[i], regex("^M$|male",   ignore_case = TRUE))){
    sampleInfo$sex[i] <- "male"
  }
  
  if(str_detect(sampleInfo$`characteristics_ch1.1`[i], regex("^F$|female", ignore_case = TRUE))){
    sampleInfo$sex[i] <- "female"
  }
}

sampleInfo$age <- as.numeric(sampleInfo$`age:ch1`)

# Keep only the columns we need
sampleInfo <- sampleInfo %>%
  select(title, geo_accession, group, sex, age)

head(sampleInfo)

# Align expression matrix with metadata 

# Row order of sampleInfo must match column order of exprs(gse)
# A mismatch here silently corrupts all downstream results
sampleInfo           <- sampleInfo[match(colnames(exprs(gse)), sampleInfo$geo_accession), ]
rownames(sampleInfo) <- sampleInfo$geo_accession

# Must return TRUE
all(colnames(exprs(gse)) == sampleInfo$geo_accession)


# Sample clustering — correlation heatmap 

corMatrix <- cor(exprs(gse), use = "c")
pheatmap(corMatrix)

# Print row names and column names to confirm they match
rownames(sampleInfo)
colnames(corMatrix)

# Add group and sex annotation sidebar
pheatmap(corMatrix, annotation_col = sampleInfo %>% select(group, sex))


# Principal Components Analysis

# Transpose the expression matrix so samples are rows
pca <- prcomp(t(exprs(gse)))

# Join PCs to sample information and plot
cbind(sampleInfo, pca$x) %>%
  ggplot(aes(x = PC1, y = PC2, col = group, label = geo_accession)) +
  geom_point() +
  geom_text_repel()

#  Differential expression analysis

design <- model.matrix(~ 0 + sampleInfo$group)
design

# Rename columns to match group names
colnames(design) <- c("copd_smoker", "healthy_smoker", "non_smoker")
design

# Filter lowly expressed probes
# Keep probes expressed above the median in more than 3 samples
cutoff       <- median(exprs(gse))
is_expressed <- exprs(gse) > cutoff
keep         <- rowSums(is_expressed) > 3

table(keep)

gse <- gse[keep, ]

# Calculate array weights to downweight any outlier arrays
aw <- arrayWeights(exprs(gse), design)
aw

# Fit the linear model
fit <- lmFit(exprs(gse), design, weights = aw)
head(fit$coefficients)

# Define contrasts — three biologically meaningful pairwise comparisons:
#   1. smoker vs non-smoker     — effect of smoking alone
#   2. COPD vs non-smoker       — smoking + disease combined
#   3. COPD vs healthy smoker   — disease effect within smokers
contrasts <- makeContrasts(
  healthy_smoker - non_smoker,
  copd_smoker    - non_smoker,
  copd_smoker    - healthy_smoker,
  levels = design
)

fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2)

topTable(fit2)
topTable(fit2, coef = 1)
topTable(fit2, coef = 2)
topTable(fit2, coef = 3)

# How many genes are differentially expressed overall?
summary(decideTests(fit2))


# Gene annotation

anno <- fData(gse)
head(anno)

anno       <- select(anno, ID, GB_ACC)
fit2$genes <- anno

topTable(fit2)


# Volcano plots

full_results1 <- topTable(fit2, coef = 1, number = Inf)
full_results2 <- topTable(fit2, coef = 2, number = Inf)
full_results3 <- topTable(fit2, coef = 3, number = Inf)

p_cutoff <- 0.05
fc_cutoff <- 1

# Contrast 1: healthy smoker vs non-smoker
full_results1 %>%
  mutate(Significant = adj.P.Val < p_cutoff & abs(logFC) > fc_cutoff) %>%
  ggplot(aes(x = logFC, y = -log10(adj.P.Val), col = Significant)) +
  geom_point() +
  labs(title = "Healthy Smoker vs Non-Smoker",
       x = "log2 Fold Change",
       y = "-log10(adjusted p-value)")

# Contrast 2: COPD vs non-smoker
full_results2 %>%
  mutate(Significant = adj.P.Val < p_cutoff & abs(logFC) > fc_cutoff) %>%
  ggplot(aes(x = logFC, y = -log10(adj.P.Val), col = Significant)) +
  geom_point() +
  labs(title = "COPD vs Non-Smoker",
       x = "log2 Fold Change",
       y = "-log10(adjusted p-value)")

# Contrast 3: COPD vs healthy smoker
full_results3 %>%
  mutate(Significant = adj.P.Val < p_cutoff & abs(logFC) > fc_cutoff) %>%
  ggplot(aes(x = logFC, y = -log10(adj.P.Val), col = Significant)) +
  geom_point() +
  labs(title = "COPD vs Healthy Smoker",
       x = "log2 Fold Change",
       y = "-log10(adjusted p-value)")


# Gene of interest lookup 

# CXCL14 (NM_004887) — primary gene studied in GSE21359
# NKX3-1 (NM_006167) — lung epithelial transcription factor

filter(full_results1, GB_ACC == "NM_004887")
filter(full_results2, GB_ACC == "NM_004887")
filter(full_results3, GB_ACC == "NM_004887")

filter(full_results1, GB_ACC == "NM_006167")
filter(full_results2, GB_ACC == "NM_006167")
filter(full_results3, GB_ACC == "NM_006167")
