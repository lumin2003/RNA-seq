countdata <- read.table("AV.txt", header = TRUE, sep = "\t", row.names = 1)
colnames(countdata) <- factor(c("A1","A2", "A3","AV1","AV2", "AV3"))

Treat <- factor(c("A","A","A","AV","AV","AV"))

batch <- factor(substring(colnames(countdata), 
                       nchar(colnames(countdata)), 
                       nchar(colnames(countdata))))

# Create DGEList object
y <- DGEList(counts = countdata, group = Treat)

# Filter out lowly expressed genes
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes = FALSE]

# Normalize the data
y <- calcNormFactors(y)

# Plot MDS to visualize sample relationships
plotMDS(y, col = rep(1:2, each = 3))

# Design matrix including Time and interaction between Time and Treatment
design <- model.matrix(~ Treat+batch,data =countdata)

# Estimate dispersion
y <- estimateDisp(y, design, robust = TRUE)

# Plot Biological Coefficient of Variation (BCV)
plotBCV(y)

# Fit the model using Quasi-likelihood (QL) F-tests
fit <- glmQLFit(y, design, robust = TRUE)

# Perform QL F-test for all coefficients
qlf_DE <- glmQLFTest(fit)
topTags(qlf_DE)

####forDEGs
results <- glmQLFTest(fit, coef=2)
                  
resutls4 <- results$table

EnhancedVolcano(
  resutls4,
  lab = rownames(resutls4),
  x = 'logFC',
  y = 'PValue',
  xlim = c(-8, 8),
  ylim = c(0, -log10(min(resutls4$PValue, na.rm = TRUE))),
  title = 'A_vs_AV_sorting_volcano',
  selectLab = c("Cd68",  "Il21r", "Csf2rb", "Cd80", "Cd86", 
                "Stat1", "Stat2"),
  pCutoff = 0.05,
  FCcutoff = 1.0,
  col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
  gridlines.major = FALSE, gridlines.minor = FALSE) 

fit$counts
normalized_counts <- as.data.frame(fit$counts)
lognorm <- log2(normalized_counts_2+1)

heatmap_AV <- lognorm[c("Cd68","Il21r", "Cd80", "Cd86", 
                     "Stat1", "Stat2"),]
library(stringr)
heatmap_AV <- heatmap_AV %>% 
  pivot_longer(-1,names_to = "Sample",values_to = "Exp")%>%
  group_by(gene_name)%>%
  mutate(Exp=(Exp-mean(Exp))/sd(Exp))%>%
  ungroup()
library("tidyHeatmap")
tidyHeatmap::heatmap(.data = heatmap_in,
                     .row = gene_name,
                     .column=Sample,
                     .value = Exp,
                     cluster_columns = FALSE,
                     show_row_names=TRUE,
                     row_names_gp = gpar(fontsize = 10),
                     column_names_gp = gpar(fontsize = 12),
                     palette_value = c("blue", "white", "red"))

