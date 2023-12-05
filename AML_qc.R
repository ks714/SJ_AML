head(pca.df)

pca.df1 = filter(pca.df, pca.df$PC1<=50)
pca.df2 = filter(pca.df, pca.df$PC1>50) 

# load package
library(edgeR)
library(limma)

pca.df1$Group = "Normal"
pca.df2$Group = "Outliers"

counts_with_ID = cbind(all_4_counts_protein, all_4_meta)

normal_counts = filter(counts_with_ID, counts_with_ID$Sample %in% pca.df1$Sample)
#mean(normal_counts$ENSG00000155511.18)
outlier_counts = filter(counts_with_ID, counts_with_ID$Sample %in% pca.df2$Sample)

normal_counts$Group = "Normal"
outlier_counts$Group = "Outliers"
write.table(outlier_counts[,20005:20009], file = "outlier_counts_meta", sep = "\t", col.names = TRUE, row.names = TRUE, quote = F)
write.table(normal_counts[,20005:20009], file = "normal_counts_meta", sep = "\t", col.names = TRUE, row.names = TRUE, quote = F)


qc = rbind(normal_counts,outlier_counts)

qc_df = rbind(pca.df1, pca.df2)

# define group
groupLabels <- qc$Group
TS <- factor(groupLabels, levels = c("Normal", "Outliers"))

#batch correction: batch can be added as a covariate if needed
#this code does not do batch correction
design <- model.matrix(~ 0 + TS )
colnames(design) <- levels(TS)
#counts <- read.table(countFile, sep = "\t", header = TRUE, row.names = 1, check.names = F)
dge <- DGEList(counts = t(qc[,1:20004]), group = groupLabels)
write.table(qc[,1:20004], file = "count_matrix", sep = "\t", col.names = TRUE, row.names = TRUE, quote = F)


# filter out low expressed genes
count_cutoff=10
cpm_cutoff <- as.vector(cpm(count_cutoff,mean(dge$samples$lib.size) ) )
keep <- rowSums(cpm(dge) > cpm_cutoff) >= min(as.numeric(table(groupLabels)))
dge <- dge[keep, keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)

v <- voom(dge, design, plot = F)
fit <- lmFit(v, design)
cont.matrix <- makeContrasts(NormalvsOutliers = (Normal - Outliers), levels = design)
fitcon <- contrasts.fit(fit, cont.matrix)
fitcon <- eBayes(fitcon)

results1 <- topTable(fitcon, n = Inf, sort.by="P", coef="NormalvsOutliers")
results1$geneSymbol = geneInfo$geneSymbol[match(rownames(results1), geneInfo$geneID)]

outFile1 <- paste0("Normal_vs_Outliers1", "_diff.txt")
write.table(results1, file = outFile1, sep = "\t", col.names = TRUE, row.names = TRUE, quote = F)


## Graphing
volcano_plot <- with(results1, plot(logFC, -log10(adj.P.Val), pch=20, main="Volcano Plot", xlim=c(-5,5)))
N <- 5  # Adjust N based on your preference
top_genes <- head(results1, N)
points(results1$logFC[results1$adj.P.Val < 0.05], -log10(results1$adj.P.Val[results1$adj.P.Val < 0.05]), col="red", pch=20)
text(
  x = top_genes$logFC,
  y = -log10(top_genes$adj.P.Val),
  labels = top_genes$geneSymbol,
  pos = 4,
  col = "red",
  offset = 0.5
)

# Install and load the plotly package
# install.packages("plotly")
library(plotly)

# Create a scatter plot using plotly
plot_ly(
  x = results1$logFC,
  y = -log10(results1$adj.P.Val),
  type = "scatter",
  mode = "markers",
  text = results1$geneSymbol,  # Text to be displayed on hover
  marker = list(
    size = 10,
    color = ifelse(results1$adj.P.Val < 0.05, "red", "black")
  )
) %>%
  layout(
    title = "Interactive Volcano Plot",
    xaxis = list(title = "logFC"),
    yaxis = list(title = "-log10(Adjusted P-Value)")
  )

## Batch Correction
all_4_meta$LibraryPrep <- gsub("Stranded-Reverse", "Stranded", all_4_meta$LibraryPrep)

normal_meta = filter(all_4_meta, all_4_meta$Sample %in% normal_counts$Sample)

dds <- DESeqDataSetFromMatrix(t(all_4_counts %>% data.matrix()),colData = all_4_meta, design = ~ LibraryPrep + Source)
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep,]

vsd <- vst(dds, blind=FALSE)

mat <- assay(vsd)

design <- model.matrix(~ LibraryPrep + Source, data=all_4_meta)
bb <- design[,grep("LibraryPrep|Source",colnames(design))]
#mm <- design[,grep("LibraryPrep",colnames(design),invert = TRUE)]

mat <- limma::removeBatchEffect(x = mat, covariates = bb)
assay(vsd) <- mat

mat2 = mat
mat2[mat2<0]= 0

## PCA post-batch correction
top=3000
vsd_top_var = getTopMAD(mat2,topgenes = top)
pca_data = prcomp(t(log10(vsd_top_var+0.01)),retx = T,scale. = T,center = T)

dist_matrix <- dist(t(vsd_top_var))
mds_result <- cmdscale(dist_matrix, eig = TRUE)

# Extract eigenvalues
eigenvalues <- mds_result$eig

proportion_variance_2d <- sum(eigenvalues[1:2]) / sum(eigenvalues)

print(proportion_variance_2d)

plot(mds_result$points, pch = 19, col = "blue", main = "MDS Plot", xlab = "Dimension 1", ylab = "Dimension 2")


percentVar <- pca_data$sdev^2 / sum(pca_data$sdev^2)

pca.df = cbind(pca_data$x[,1:15], all_4_meta) %>% as_tibble()

p1 = ggplot(pca.df, aes(x = PC1,y=PC2, color = Source))+ 
  geom_point(size=1.5) + 
  #geom_vline(xintercept = 50, linetype = "dashed", color = "red") +
  xlab(paste0("PC1: ",round(percentVar[1] * 100,digits = 2),"% variance")) +
  ylab(paste0("PC2: ",round(percentVar[2] * 100,digits = 2),"% variance")) +
  #ggrepel::geom_text_repel(force=4,size=3,show.legend = F) +
  theme_bw() + 
  ggsci::scale_color_d3() +
  ggtitle("Count-based PCA",subtitle = "Top 1000 Protein Coding Genes - Normal Samples") +
  theme(plot.title = element_text(hjust = 0.5,size = 10),
        plot.subtitle = element_text(hjust = 0.5,size = 8))

p1

p2 = ggplot(pca.df, aes(x = PC1,y=PC2, color = LibraryPrep))+ 
  geom_point(size=1.5) + 
  geom_vline(xintercept = 50, linetype = "dashed", color = "red") +
  xlab(paste0("PC1: ",round(percentVar[1] * 100,digits = 2),"% variance")) +
  ylab(paste0("PC2: ",round(percentVar[2] * 100,digits = 2),"% variance")) +
  #ggrepel::geom_text_repel(force=4,size=3,show.legend = F) +
  theme_bw() + 
  ggsci::scale_color_d3() +
  ggtitle("Count-based PCA",subtitle = "20004 Protein Coding Genes") +
  theme(plot.title = element_text(hjust = 0.5,size = 10),
        plot.subtitle = element_text(hjust = 0.5,size = 8))

p2

gene_selector <- nearZeroVar(t(mat2), allowParallel = TRUE)

cor_matrix <- cor(t(mat2))
high_correlation_features <- findCorrelation(cor_matrix, cutoff = 0.8)
filtered_data <- gene_expression_data[, -high_correlation_features]


