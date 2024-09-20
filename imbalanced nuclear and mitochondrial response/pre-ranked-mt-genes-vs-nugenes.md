pre-ranked-mt-genes-vs-nugenes
================

``` r
# Import the mitochondrial genes
genelist <- read.table("/Users/yuanyuan/Downloads/github/hek293_rna_seq/mt-genes-vs-nu-genes/mtgenes.txt", header=T, quote = "", stringsAsFactors=F)
genelist <- genelist$KEGG_OXIDATIVE_PHOSPHORYLATION
MTgenes <- genelist[grep("MT-",genelist)]
nuMTgenes <- genelist[-grep("MT-",genelist)]
```

## plot

``` r
#the list of significance was imported from GLM analysis in edgeR, with gene expression to heteroplasmy levels

merged_results_imported <- read.csv("/Users/yuanyuan/Downloads/github/hek293_rna_seq/intermediate_results/edgeR_YM2_RM_correct_heteroplasmy_34_25.csv", header=TRUE,row.names=1)

merged_results_imported$signed_log_pvalue <- -log10(merged_results_imported$PValue) * sign(merged_results_imported$logFC) #using raw p value

#ranked
merged_results_imported <- merged_results_imported[order(merged_results_imported$signed_log_pvalue, decreasing = TRUE), ]

ranked_genes_values <- merged_results_imported$signed_log_pvalue
names(ranked_genes_values) <- merged_results_imported$hgnc_symbol
```

## test plot

``` r
# Load necessary libraries
library(ggplot2)

# Import the mitochondrial genes
genelist <- read.table("/Users/yuanyuan/Downloads/github/hek293_rna_seq/mt-genes-vs-nu-genes/mtgenes.txt", header=T, quote = "", stringsAsFactors=F)
genelist <- genelist$KEGG_OXIDATIVE_PHOSPHORYLATION
MTgenes <- genelist[grep("MT-",genelist)]
nuMTgenes <- genelist[-grep("MT-",genelist)]

# Import the ranked genes with signed p values
merged_results_imported <- read.csv("/Users/yuanyuan/Downloads/github/hek293_rna_seq/intermediate_results/edgeR_YM2_RM_correct_heteroplasmy_34_25.csv", header=TRUE, row.names=1)

# Calculate signed log p value
merged_results_imported$signed_log_pvalue <- -log10(merged_results_imported$PValue) * sign(merged_results_imported$logFC) # using raw p value

# Rank the genes
merged_results_imported <- merged_results_imported[order(merged_results_imported$signed_log_pvalue, decreasing = TRUE), ]

# Create a named vector for ranked genes values
ranked_genes_values <- merged_results_imported$signed_log_pvalue
names(ranked_genes_values) <- merged_results_imported$hgnc_symbol

# Prepare data for plotting
df <- data.frame(
  gene = names(ranked_genes_values),
  signed_log_pvalue = ranked_genes_values
)

# Label the genes as mtDNA or nuDNA
df$gene_type <- ifelse(df$gene %in% MTgenes, "mtDNA", ifelse(df$gene %in% nuMTgenes, "nuDNA", "Other"))

# Plot the data
ggplot(df, aes(x = gene, y = signed_log_pvalue, color = gene_type)) +
  geom_point() +
  scale_color_manual(values = c("mtDNA" = "blue", "nuDNA" = "red", "Other" = "grey")) +
  labs(title = "Signed Log P Values of Genes",
       x = "Genes",
       y = "Signed Log P Value",
       color = "Gene Type") +
  theme(axis.text.x = element_blank(), # Hides the x-axis text for readability
        axis.ticks.x = element_blank()) # Hides the x-axis ticks for readability
```

![](pre-ranked-mt-genes-vs-nugenes_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

``` r
# Load necessary libraries
library(ggplot2)
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(ggrepel)  # for labeling

# Import the mitochondrial genes
genelist <- read.table("/Users/yuanyuan/Downloads/github/hek293_rna_seq/mt-genes-vs-nu-genes/mtgenes.txt", header=TRUE, quote = "", stringsAsFactors=FALSE)
genelist <- genelist$KEGG_OXIDATIVE_PHOSPHORYLATION
MTgenes <- genelist[grep("MT-", genelist)]
nuMTgenes <- genelist[-grep("MT-", genelist)]

# Import the merged results
merged_results_imported <- read.csv("/Users/yuanyuan/Downloads/github/hek293_rna_seq/intermediate_results/edgeR_YM2_RM_correct_heteroplasmy_34_25.csv", header=TRUE, row.names=1)

# Calculate signed log p-value
merged_results_imported$signed_log_pvalue <- -log10(merged_results_imported$FDR) * sign(merged_results_imported$logFC) # using raw p value

# Rank genes by signed log p-value in decreasing order
merged_results_imported <- merged_results_imported[order(merged_results_imported$signed_log_pvalue, decreasing = TRUE), ]

# Filter genes to include only those in genelist
filtered_genes <- merged_results_imported %>% filter(hgnc_symbol %in% genelist)

# Add a new column to differentiate between MT genes and nu genes
filtered_genes$gene_type <- ifelse(filtered_genes$hgnc_symbol %in% MTgenes, "mtDNA", "nuDNA")

# Plot the data
ggplot(filtered_genes, aes(x = reorder(hgnc_symbol, -signed_log_pvalue), y = signed_log_pvalue, color = gene_type)) +
  geom_point() +
  geom_text_repel(aes(label = ifelse(gene_type == "mtDNA", hgnc_symbol, "")), color = "red", size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_color_manual(values = c("mtDNA" = "red", "nuDNA" = "blue")) +
  labs(x = "Rank sorted genes", y = "Signed log FDR", color = "Gene Type") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),   # Remove minor grid lines
    axis.line = element_line(color = "black")  # Retain axis lines
  ) +
  ggtitle("Signed log p-value of Pre-ranked Genes")
```

    ## Warning: ggrepel: 8 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

![](pre-ranked-mt-genes-vs-nugenes_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
# Load necessary libraries
library(ggplot2)
library(dplyr)
library(ggrepel)  # for labeling

# Import the mitochondrial genes
genelist <- read.table("/Users/yuanyuan/Downloads/github/hek293_rna_seq/mt-genes-vs-nu-genes/mtgenes.txt", header=TRUE, quote = "", stringsAsFactors=FALSE)
genelist <- genelist$KEGG_OXIDATIVE_PHOSPHORYLATION
MTgenes <- genelist[grep("MT-", genelist)]
nuMTgenes <- genelist[-grep("MT-", genelist)]

# Import the merged results
merged_results_imported <- read.csv("/Users/yuanyuan/Downloads/github/hek293_rna_seq/intermediate_results/edgeR_YM2_RM_correct_heteroplasmy_34_25.csv", header=TRUE, row.names=1)

# Calculate signed log p-value
merged_results_imported$signed_log_pvalue <- -log10(merged_results_imported$FDR) * sign(merged_results_imported$logFC) # using raw p value

# Rank genes by signed log p-value in decreasing order
merged_results_imported <- merged_results_imported[order(merged_results_imported$signed_log_pvalue, decreasing = TRUE), ]

# Filter genes to include only those in genelist
filtered_genes <- merged_results_imported %>% filter(hgnc_symbol %in% genelist)

# Add a new column to differentiate between MT genes and nu genes
filtered_genes$gene_type <- ifelse(filtered_genes$hgnc_symbol %in% MTgenes, "mtDNA", "nuDNA")

# Plot the data
plot <-  ggplot(filtered_genes, aes(x = reorder(hgnc_symbol, -signed_log_pvalue), y = signed_log_pvalue, color = gene_type)) +
  geom_point() +
  geom_text_repel(aes(label = ifelse(gene_type == "mtDNA", hgnc_symbol, "")), color = "red", size = 3, max.overlaps = Inf) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_color_manual(values = c("mtDNA" = "red", "nuDNA" = "blue")) +
  labs(x = "Rank sorted genes", y = "Signed log FDR", color = "Gene Type") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),   # Remove minor grid lines
    axis.line = element_line(color = "black")  # Retain axis lines
  ) +
  ggtitle("Signed log p-value of Pre-ranked Genes")

print(plot)
```

![](pre-ranked-mt-genes-vs-nugenes_files/figure-gfm/unnamed-chunk-2-2.png)<!-- -->

``` r
# ggsave("/Users/yuanyuan/Downloads/github/hek293_rna_seq/plots and tables to present/mtDNA_nuDNA_signed_logFDR.pdf", plot = plot, width = 8, height = 5)
```

## Question:

Does the length of the MT-DNA genes affect the gene expression level
post mtND5 mutation?

``` r
# Load necessary libraries
library(ggplot2)
library(dplyr)
library(ggpubr)  # for stat_cor function
library(ggrepel)  # for geom_text_repel

# Import the mitochondrial genes
genelist <- read.table("/Users/yuanyuan/Downloads/github/hek293_rna_seq/mt-genes-vs-nu-genes/mtgenes.txt", header=TRUE, quote = "", stringsAsFactors=FALSE)
genelist <- genelist$KEGG_OXIDATIVE_PHOSPHORYLATION
MTgenes <- genelist[grep("MT-", genelist)]

# Import the merged results
merged_results_imported <- read.csv("/Users/yuanyuan/Downloads/github/hek293_rna_seq/intermediate_results/edgeR_YM2_RM_correct_heteroplasmy_34_25.csv", header=TRUE, row.names=1)

# Filter to include only MT genes
mt_genes_data <- merged_results_imported %>% filter(hgnc_symbol %in% MTgenes)

# Create the mt_gene_lengths dataframe with lengths (in base pairs)
mt_gene_lengths <- data.frame(
  gene = c("MT-ND1", "MT-ND2", "MT-CO1", "MT-CO2", "MT-ATP6", "MT-ATP8", "MT-CO3", "MT-ND3", "MT-ND4L", "MT-ND4", "MT-ND5", "MT-ND6", "MT-CYB"),
  length = c(956, 1042, 1541, 684, 681, 207, 784, 351, 297, 1378, 1812, 525, 1140)
)

# Merge the lengths with the expression data
mt_genes_data <- mt_genes_data %>% 
  inner_join(mt_gene_lengths, by = c("hgnc_symbol" = "gene"))

# Plot the data
ggplot(mt_genes_data, aes(x = length, y = logFC, label = hgnc_symbol)) +
  geom_point(color = "red") +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  stat_cor(method = "pearson", label.x = 1200, label.y = max(mt_genes_data$logFC) + 0.5) +  # Adjusted position for p-value and r-value
  labs(x = "Gene Length (bp)", y = "Log Fold Change", title = "Correlation between mtDNA Gene Length and Expression Level") +
  theme_classic() +
  geom_text_repel()  # Use geom_text_repel to avoid overlapping annotations
```

    ## `geom_smooth()` using formula = 'y ~ x'

    ## Warning: The following aesthetics were dropped during statistical transformation: label
    ## ℹ This can happen when ggplot fails to infer the correct grouping structure in
    ##   the data.
    ## ℹ Did you forget to specify a `group` aesthetic or to convert a numerical
    ##   variable into a factor?

![](pre-ranked-mt-genes-vs-nugenes_files/figure-gfm/gene%20length%20correlation%20to%20gene%20expression-1.png)<!-- -->

\##Could the result above be artifact? Include gene length when
conducting GLM analysis in RNA seq data

``` r
##this part can be ignored by reading in the exported table directly


# Load necessary library
library(biomaRt)
library(dplyr)

##select genes
rawCounts<-read.table('/Users/yuanyuan/Downloads/github/hek293_rna_seq/5022D_rawCounts.txt',header=T,row.names=1)
ensembl_ids <- row.names(rawCounts)

# Specify the Ensembl database
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Retrieve canonical transcript information
canonical_transcripts <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "transcript_length"),
                               filters = "ensembl_gene_id",
                               values = ensembl_ids,
                               mart = ensembl)

# Print the first few rows to check the results
head(canonical_transcripts)
```

    ##   ensembl_gene_id hgnc_symbol transcript_length
    ## 1 ENSG00000012817       KDM5D              5501
    ## 2 ENSG00000012817       KDM5D              6825
    ## 3 ENSG00000012817       KDM5D              5579
    ## 4 ENSG00000012817       KDM5D              5134
    ## 5 ENSG00000012817       KDM5D              5235
    ## 6 ENSG00000012817       KDM5D              4701

``` r
# Retrieve other gene information
biomart_data <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "description"), 
                      filters = "ensembl_gene_id", 
                      values = ensembl_ids, 
                      mart = ensembl)

# Merge the canonical transcript lengths with the other gene information
biomart_data <- merge(biomart_data, canonical_transcripts, by = c("ensembl_gene_id", "hgnc_symbol"))



# Print the first few rows to check the results
#head(biomart_data)

# ##Export the data to avoid biomart querying in the future
# # Specify the file path and name
# output_file <- "/Users/yuanyuan/Downloads/github/hek293_rna_seq/intermediate_results/biomart_data_with_transcript_length.csv"
# 
# # Export the biomart_data dataframe to a CSV file
# write.csv(biomart_data, file = output_file, row.names = FALSE)
```

Next do the fitting with transcript length

``` r
# Load necessary libraries
library(edgeR)
```

    ## Loading required package: limma

``` r
library(dplyr)

# Load the raw counts
rawCounts<-read.table('/Users/yuanyuan/Downloads/github/hek293_rna_seq/5022D_rawCounts.txt',header=T,row.names=1)
# Remove YW2 column
rawCounts <- rawCounts[,-2]

# Prepare the DGEList object
dge <- DGEList(counts = rawCounts)

# Filter lowly expressed genes
keep <- filterByExpr(dge)
```

    ## Warning in filterByExpr.DGEList(dge): All samples appear to belong to the same
    ## group.

``` r
dge <- dge[keep, ]

# Normalize the data
dge <- calcNormFactors(dge)

# Define heteroplasmy levels (should match the number of columns in rawCounts)
heteroplasmy_levels_corr <- c(0, 1, 1, 0, 0, 0.6, 0.6)

# Prepare the design matrix for heteroplasmy levels
design <- model.matrix(~ heteroplasmy_levels_corr)

# Estimate dispersion
dge <- estimateDisp(dge, design)

# Fit the GLM
fit <- glmFit(dge, design)

# Perform likelihood ratio test
lrt <- glmLRT(fit)

# Extract results
results <- topTags(lrt, n = Inf)$table

# Print top results
head(results)
```

    ##                      logFC   logCPM       LR       PValue          FDR
    ## ENSG00000241553 -0.7738171 4.919317 41.86266 9.791489e-11 7.760861e-07
    ## ENSG00000074842 -0.6471856 6.690900 41.83892 9.911067e-11 7.760861e-07
    ## ENSG00000126768 -0.5514007 6.180053 39.58295 3.144169e-10 1.641361e-06
    ## ENSG00000160799 -0.6066968 6.117878 35.89204 2.085592e-09 8.165613e-06
    ## ENSG00000148291 -0.6896875 4.627463 34.31380 4.690435e-09 1.469138e-05
    ## ENSG00000142507 -0.6680732 6.116228 33.20531 8.292364e-09 1.834993e-05

``` r
# Specify the file path for biomart data
input_file <- "/Users/yuanyuan/Downloads/github/hek293_rna_seq/intermediate_results/biomart_data_with_transcript_length.csv"

# Read the CSV file into a dataframe
biomart_data_imported <- read.csv(input_file, header = TRUE, stringsAsFactors = FALSE)

# Prioritize by keeping the row with a non-NA hgnc_symbol (if any)
biomart_data_imported <- biomart_data_imported[order(biomart_data_imported$ensembl_gene_id, !is.na(biomart_data_imported$hgnc_symbol)), ]

# Remove duplicated Ensembl gene IDs, keeping the first occurrence
biomart_data_unique <- biomart_data_imported[!duplicated(biomart_data_imported$ensembl_gene_id), ]

# Ensure row names in rawCounts match the Ensembl IDs in biomart_data_unique
gene_lengths <- biomart_data_unique %>%
  filter(ensembl_gene_id %in% rownames(rawCounts)) %>%
  dplyr::select(ensembl_gene_id, transcript_length)

# Merge results with gene lengths
results$ensembl_gene_id <- rownames(results)
results_merged <- merge(results, gene_lengths, by = "ensembl_gene_id")



# 
# 
# # Plot logFC vs. transcript_length with x-axis limit set to 0-5000
# library(ggplot2)
# ggplot(results_merged, aes(x = transcript_length, y = logFC)) +
#   geom_point(color = "red") +
#   geom_smooth(method = "lm", se = TRUE, color = "blue") +
#   stat_cor(method = "pearson", label.x = 3500, label.y = max(results_merged$logFC) * 0.9) +
#   labs(x = "Transcript Length (bp)", y = "Log Fold Change", title = "Correlation between Gene Length and Log Fold Change") +
#   theme_minimal() +
#   xlim(0, 5000)


# Fit the linear model
fit <- lm(logFC ~ transcript_length, data = results_merged)

# Extract the coefficients
coefficients <- coef(fit)
intercept <- round(coefficients[1], 3)
slope <- round(coefficients[2], 3)

# Create the formula string
fit_formula <- paste("y =", intercept, "+", slope, "* x")


# Load necessary libraries
library(ggplot2)
library(ggpubr)

# Plot with the fitted line and formula
ggplot(results_merged, aes(x = transcript_length, y = logFC)) +
  geom_point(color = "red") +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  stat_cor(method = "pearson", label.x = 3500, label.y = max(results_merged$logFC) * 0.9) +
  annotate("text", x = 3500, y = max(results_merged$logFC) * 0.8, label = fit_formula, color = "black") +
  labs(x = "Transcript Length (bp)", y = "Log Fold Change", title = "Correlation between Gene Length and Log Fold Change") +
  theme_minimal() +
  xlim(0, 5000)
```

    ## `geom_smooth()` using formula = 'y ~ x'

    ## Warning: Removed 1387 rows containing non-finite values (`stat_smooth()`).

    ## Warning: Removed 1387 rows containing non-finite values (`stat_cor()`).

    ## Warning: Removed 1387 rows containing missing values (`geom_point()`).

![](pre-ranked-mt-genes-vs-nugenes_files/figure-gfm/GLM%20fitting-1.png)<!-- -->

``` r
####the fitted model

# Fit the linear model
fit <- lm(logFC ~ transcript_length, data = results_merged)

# Extract the coefficients without rounding
coefficients <- coef(fit)
intercept <- coefficients[1]
slope <- coefficients[2]

# Print the coefficients
print(coefficients)
```

    ##       (Intercept) transcript_length 
    ##     -2.195039e-02      1.801518e-05

``` r
# Get the summary of the linear model
summary(fit)
```

    ## 
    ## Call:
    ## lm(formula = logFC ~ transcript_length, data = results_merged)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -3.3547 -0.1657 -0.0004  0.1683  1.5060 
    ## 
    ## Coefficients:
    ##                     Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)       -2.195e-02  2.856e-03  -7.687  1.6e-14 ***
    ## transcript_length  1.802e-05  8.371e-07  21.522  < 2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.2823 on 15587 degrees of freedom
    ## Multiple R-squared:  0.02886,    Adjusted R-squared:  0.0288 
    ## F-statistic: 463.2 on 1 and 15587 DF,  p-value: < 2.2e-16

Conclusion: There is no systematic correlation for all genes. The mtDNA
gene logFC to transcript length is true finding. Next question: Does
this correlation applies to other nuclear expressed oxidative
phosphorylation genes?

``` r
# Load required libraries
library(dplyr)
library(readr)
library(ggplot2)
library(ggrepel)

# Specify the file path for biomart data
input_file <- "/Users/yuanyuan/Downloads/github/hek293_rna_seq/intermediate_results/biomart_data_with_transcript_length.csv"

# Read the CSV file into a dataframe
biomart_data_imported <- read_csv(input_file, show_col_types = FALSE)

# Prioritize by keeping the row with a non-NA hgnc_symbol (if any)
biomart_data_imported <- biomart_data_imported[order(biomart_data_imported$ensembl_gene_id, !is.na(biomart_data_imported$hgnc_symbol)), ]

# Remove duplicated Ensembl gene IDs, keeping the first occurrence
biomart_data_unique <- biomart_data_imported[!duplicated(biomart_data_imported$ensembl_gene_id), ]

# Import the mitochondrial genes
genelist <- read_table("/Users/yuanyuan/Downloads/github/hek293_rna_seq/mt-genes-vs-nu-genes/mtgenes.txt", col_types = cols())
genelist <- genelist$KEGG_OXIDATIVE_PHOSPHORYLATION
MTgenes <- genelist[grep("MT-", genelist)]
nuMTgenes <- genelist[-grep("MT-", genelist)]

# Import the merged results
merged_results_imported <- read_csv("/Users/yuanyuan/Downloads/github/hek293_rna_seq/intermediate_results/edgeR_YM2_RM_correct_heteroplasmy_34_25.csv", col_types = cols(), col_names = TRUE)
```

    ## New names:
    ## • `` -> `...1`

``` r
# Rename the first column to 'ensembl_gene_id'
colnames(merged_results_imported)[1] <- "ensembl_gene_id"

# Calculate signed log p-value
merged_results_imported <- merged_results_imported %>%
  mutate(signed_log_pvalue = -log10(FDR) * sign(logFC))

# Rank genes by signed log p-value in decreasing order
merged_results_imported <- merged_results_imported %>%
  arrange(desc(signed_log_pvalue))

# Filter genes to include only those in genelist
filtered_genes <- merged_results_imported %>%
  filter(hgnc_symbol %in% genelist)

# Add a new column to differentiate between MT genes and nu genes
filtered_genes <- filtered_genes %>%
  mutate(gene_type = ifelse(hgnc_symbol %in% MTgenes, "mtDNA", "nuDNA"))

# Merge filtered genes with biomart data to add transcript length
final_genes <- filtered_genes %>%
  left_join(biomart_data_unique %>% dplyr::select(ensembl_gene_id, transcript_length), by = "ensembl_gene_id")

# Convert transcript length to kiloBase (kb)
final_genes <- final_genes %>%
  mutate(transcript_length_kb = transcript_length / 1000)

# Calculate the Fold Change (10^logFC)
final_genes <- final_genes %>%
  mutate(FoldChange = 10^logFC)

# Ensure data is ready for plotting
final_genes <- final_genes %>%
  mutate(color = ifelse(gene_type == "mtDNA", "#ed3479", "black"))

# Separate the data for mtDNA and nuDNA genes
mtDNA_genes <- final_genes %>% filter(gene_type == "mtDNA")
nuDNA_genes <- final_genes %>% filter(gene_type == "nuDNA")

# Run regression models for mtDNA and nuDNA genes
model_mtDNA <- lm(FoldChange ~ transcript_length_kb, data = mtDNA_genes)
summary_model_mtDNA <- summary(model_mtDNA)
r_squared_mtDNA <- summary_model_mtDNA$r.squared
p_value_mtDNA <- summary_model_mtDNA$coefficients[2, 4]

model_nuDNA <- lm(FoldChange ~ transcript_length_kb, data = nuDNA_genes)
summary_model_nuDNA <- summary(model_nuDNA)
r_squared_nuDNA <- summary_model_nuDNA$r.squared
p_value_nuDNA <- summary_model_nuDNA$coefficients[2, 4]

# Format p-values for display
formatted_p_value_mtDNA <- format(p_value_mtDNA, scientific = TRUE)
formatted_p_value_nuDNA <- format(p_value_nuDNA, scientific = TRUE)

# Plot
plot <- ggplot(final_genes, aes(x = transcript_length_kb, y = FoldChange, color = color)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  scale_color_identity() +
  ylim(0, 3) +  # Set y-axis limits
  geom_text_repel(data = subset(final_genes, gene_type == "mtDNA"), aes(label = hgnc_symbol), color = "#ed3479") +
  annotate("text", x = -Inf, y = Inf, label = paste0("mtDNA R^2: ", round(r_squared_mtDNA, 3), "\nP: ", formatted_p_value_mtDNA), hjust = -0.6, vjust = 1.1, size = 5, color = "#ed3479") +
  annotate("text", x = Inf, y = Inf, label = paste0("nuDNA R^2: ", round(r_squared_nuDNA, 3), "\nP: ", formatted_p_value_nuDNA), hjust = 1.1, vjust = 1.1, size = 5, color = "black") +
  labs(title = "Fold Change vs. Transcript Length for mtDNA and nuDNA Genes", x = "Transcript Length (kb)", y = "Fold Change") +
  theme_classic()

print(plot)
```

    ## `geom_smooth()` using formula = 'y ~ x'

    ## Warning: Removed 2 rows containing non-finite values (`stat_smooth()`).

    ## Warning: Removed 2 rows containing missing values (`geom_point()`).

![](pre-ranked-mt-genes-vs-nugenes_files/figure-gfm/r%20correlation%20in%20nuDNA%20mtDNA%20genes-1.png)<!-- -->

``` r
# # Save the plot as an SVG file
# ggsave("/Users/yuanyuan/Downloads/github/hek293_rna_seq/plots and tables to present/mtDNA_nuDNA_transcript_length_to_FoldChange.pdf", plot = plot, width = 8, height = 5)
```

Here is the equation for the fold change in the expression of *MT-ND5*
relative to the control:

$$
\text{Fold change (MT-ND5/Control)} = k \times \text{transcript length (kp)}
$$

\##the fitted model parameters for nuDNA and mtDNA

``` r
##work with previous trunk
print("13 mtDNA encoded genes:")
```

    ## [1] "13 mtDNA encoded genes:"

``` r
summary(model_mtDNA) ### 13 mtDNA encoded genes ##Fold Change ## Kilobases in x axis
```

    ## 
    ## Call:
    ## lm(formula = FoldChange ~ transcript_length_kb, data = mtDNA_genes)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.17546 -0.07708 -0.04870  0.03451  0.47553 
    ## 
    ## Coefficients:
    ##                      Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)          -0.06152    0.10877  -0.566    0.583    
    ## transcript_length_kb  0.93350    0.10904   8.561 3.41e-06 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.1872 on 11 degrees of freedom
    ## Multiple R-squared:  0.8695, Adjusted R-squared:  0.8576 
    ## F-statistic: 73.29 on 1 and 11 DF,  p-value: 3.41e-06

``` r
# print("total nucelar oxidative phosphotylated genes:")
# summary(model_nuDNA) ##total nucelar oxidative phosphotylated genes

# print("total genes sequenced:")
# summary(fit) ##total genes sequenced
# 
```

## 3243A\>G heteroplasmy

``` r
##source: https://www.pnas.org/doi/full/10.1073/pnas.1414028111
# Load necessary library
library(data.table)

file <- "/Users/yuanyuan/Downloads/github/hek293_rna_seq/mtDNA gene expression and transcript length/GSE56158_RPKM_cybrids_normalized_readgroups.txt"

data <- read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)


data <- t(data)

# Convert back to dataframe if necessary
data <- as.data.frame(data)

# Assign the first row as the column names
colnames(data) <- data[1, ]
data <- data[-1, ]

# Import the mitochondrial genes
genelist <- read_table("/Users/yuanyuan/Downloads/github/hek293_rna_seq/mt-genes-vs-nu-genes/mtgenes.txt", col_types = cols())
genelist <- genelist$KEGG_OXIDATIVE_PHOSPHORYLATION
MTgenes <- genelist[grep("MT-", genelist)]
nuMTgenes <- genelist[-grep("MT-", genelist)]

# Convert the data to numeric (since it is currently in character format)
row_names <- rownames(data)
# Convert data to numeric while preserving the structure
data_numeric <- data
data_numeric[] <- lapply(data_numeric, as.numeric)

# Restore the row names
rownames(data_numeric) <- row_names
#data <- data.frame(lapply(data, as.numeric))

data_filtered <- data.frame(Gene = rownames(data_numeric), data_numeric, stringsAsFactors = FALSE)

# Filter for MT genes
MT_data <- data_filtered[data_filtered$Gene %in% MTgenes, ]
MT_data$Gene_Type <- "MT"

# Filter for nuMT genes
nuMT_data <- data_filtered[data_filtered$Gene %in% nuMTgenes, ]
nuMT_data$Gene_Type <- "nuMT"

# Combine the MT and nuMT datasets
filtered_data <- rbind(MT_data, nuMT_data)



# the control columns are the first 3 columns after transposing
control <- data[, 1:3]
hetero_20 <- data[,4:6]
hetero_30 <- data[,7:9]
hetero_50 <- data[,10:12]
hetero_60 <- data[,13:15]
hetero_90 <- data[,16:18]
hetero_100 <- data[,19:21]



control_mean <- rowMeans(control)
hetero_20_mean <- rowMeans(hetero_20)
hetero_30_mean <- rowMeans(hetero_30)
hetero_50_mean <- rowMeans(hetero_50)
hetero_60_mean <- rowMeans(hetero_60)
hetero_90_mean <- rowMeans(hetero_90)
hetero_100_mean <- rowMeans(hetero_100)
```
