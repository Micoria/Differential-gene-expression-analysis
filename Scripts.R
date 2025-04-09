if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("limma", "clusterProfiler", "org.Hs.eg.db", "hgu133plus2.db", "enrichplot", "pheatmap", "ggplot2", "VennDiagram"))
library(limma)
library(clusterProfiler)
library(org.Hs.eg.db)
library(hgu133plus2.db)
library(enrichplot)
library(pheatmap)
library(ggplot2)
library(VennDiagram)
if (!requireNamespace("AnnotationDbi", quietly = TRUE))
  BiocManager::install("AnnotationDbi")
library(AnnotationDbi)
if (!requireNamespace("limma", quietly = TRUE))
  BiocManager::install("limma")
library(limma)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
if (!requireNamespace("msigdbr", quietly = TRUE)) {
  install.packages("msigdbr")
}
library(msigdbr)
library(dplyr)


setwd("/Users/micoria/Documents/work/MedSci/0407")

lines <- readLines("GSE39582_series_matrix.csv")
lines <- readLines("GSE39582_series_matrix.csv")
expr_start <- grep("^ID_REF,", lines)  
expr_end <- grep("^!series_matrix_table_end", lines)  
expr_data <- read.table("GSE39582_series_matrix.csv",
                        skip = expr_start - 1,  
                        nrows = expr_end - expr_start - 1,
                        header = TRUE,
                        sep = ",",
                        stringsAsFactors = FALSE,
                        check.names = FALSE)  
head(expr_data)
str(expr_data)
grep("Sample_source_name_ch1", lines, value = TRUE)

#extract information of samples
lines <- readLines("GSE39582_series_matrix.csv")
sample_ids_line <- grep("^!Sample_geo_accession", lines, value = TRUE)
sample_ids <- unlist(strsplit(sample_ids_line, ","))[-1]  
group_line <- grep("^!Sample_source_name_ch1", lines, value = TRUE)
group_info <- unlist(strsplit(group_line, ","))[-1]  
sample_info <- data.frame(sample = sample_ids, group = group_info, stringsAsFactors = FALSE)
head(sample_info)
sample_info$group <- ifelse(grepl("non tumoral", sample_info$group),
                            "Normal",
                            "Tumor")
head(sample_info)
table(sample_info$group)

#####transfer probe ID into gene symbol
probe_ids <- expr_data$ID_REF
BiocManager::install("hgu133plus2.db")
library(hgu133plus2.db)

gene_symbols <- mapIds(hgu133plus2.db,
                       keys = probe_ids,
                       column = "SYMBOL",
                       keytype = "PROBEID",
                       multiVals = "first")
expr_data$ID_REF <- gene_symbols
expr_data <- expr_data[!is.na(expr_data$ID_REF), ]
head(expr_data)

######Limma
library(limma)
str(expr_data)
expr_data <- expr_data[!duplicated(expr_data$ID_REF), ]
rownames(expr_data) <- expr_data$ID_REF
expr_matrix <- expr_data[, -1]  
expr_matrix <- expr_matrix[, sample_info$sample]
group_list <- factor(sample_info$group, levels = c("Normal", "Tumor"))
design <- model.matrix(~ group_list)

fit <- lmFit(expr_matrix, design)
fit <- eBayes(fit)
deg <- topTable(fit, coef = 2, number = Inf, adjust.method = "BH")
head(deg)
write.csv(deg, "DEG_results.csv")

######go/kegg
# set threshold：logFC > 1 or < -1, adj.P.Val < 0.05
deg_up <- deg[deg$logFC > 1 & deg$adj.P.Val < 0.05, ]
count(deg_up)
count(deg_down)
deg_down <- deg[deg$logFC < -1 & deg$adj.P.Val < 0.05, ]
deg_all <- rbind(deg_up, deg_down)
cat("Upregulate_gene：", nrow(deg_up), "\n")
cat("Downregulate_gene：", nrow(deg_down), "\n")
write.csv(deg_up, "DEG_up.csv")
write.csv(deg_down, "DEG_down.csv")
write.csv(deg_all, "DEG_all.csv")

deg_genes <- rownames(deg_all)
entrez_ids <- mapIds(org.Hs.eg.db,
                     keys = deg_genes,
                     column = "ENTREZID",
                     keytype = "SYMBOL",
                     multiVals = "first")
entrez_ids <- na.omit(entrez_ids)
cat("the number of transfering ENTREZ ID sucessfully：", length(entrez_ids), "\n")

# GO richment
library(clusterProfiler)

ego <- enrichGO(gene = entrez_ids,
                OrgDb = org.Hs.eg.db,
                ont = "ALL",  
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.2,
                readable = TRUE)
head(ego)
write.csv(as.data.frame(ego), "GO_enrichment_results.csv")
library(ggplot2)
dotplot(ego, showCategory = 20, title = "GO Enrichment Dotplot") +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
ggsave("GO_dotplot.png", width = 10, height = 8)
barplot(ego, showCategory = 20, title = "GO Enrichment Barplot") +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
ggsave("GO_barplot.png", width = 10, height = 8)

#kegg
kegg <- enrichKEGG(gene = entrez_ids,
                   organism = 'hsa',
                   pvalueCutoff = 0.05)

head(kegg)
write.csv(as.data.frame(kegg), "KEGG_enrichment_results.csv")
dotplot(kegg, showCategory = 20, title = "KEGG Enrichment Dotplot") +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
ggsave("KEGG_dotplot.png", width = 10, height = 8)
barplot(kegg, showCategory = 20, title = "KEGG Enrichment Barplot") +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
ggsave("KEGG_barplot.png", width = 10, height = 8)

#########Gene intersection
#download circadian genes
library(msigdbr)
# check all of collection and subcollection
#regulator gene
collections_info <- msigdbr_collections()
print(collections_info)
library(dplyr)
msigdb_data <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG_LEGACY")
circadian_genes <- msigdb_data %>%
  filter(grepl("CIRCADIAN", gs_name, ignore.case = TRUE)) %>%
  pull(gene_symbol) %>%
  unique()
write.table(circadian_genes, file = "circadian_genes.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
#download CLDN gene
cldn_family_genes <- paste0("CLDN", 1:23)
cldn_genes_in_data <- cldn_family_genes[cldn_family_genes %in% rownames(expr_matrix)]
write.table(cldn_genes_in_data, file = "cldn_family_genes.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

deg_gene_symbols <- rownames(deg_all)

circadian_deg <- intersect(deg_gene_symbols, circadian_genes)
cldn_deg <- intersect(deg_gene_symbols, cldn_genes_in_data)

print(circadian_deg)
print(cldn_deg)
write.table(circadian_deg, file = "circadian_deg_genes.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(cldn_deg, file = "cldn_deg_genes.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

######venn figure
# make sure that you have these two variables：
# 1. deg_genes ：差异表达基因 SYMBOL
# 2. circadian_genes ：节律基因 SYMBOL



if (!requireNamespace("VennDiagram", quietly = TRUE)) {
  install.packages("VennDiagram")
}
library(VennDiagram)

gene_list <- list(
  DEG = deg_gene_symbols,
  Circadian = circadian_genes
)


library(ggVennDiagram)
library(ggplot2)

venn_plot <- ggVennDiagram(gene_list, label_alpha = 0, label_geom = "label") +
  scale_fill_gradient(low = "white", high = "skyblue") +
  theme_void() +
  theme(
    text = element_text(size = 12, face = "bold"),  
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    plot.margin = margin(2, 2, 2, 2, "cm") 
  ) +
  ggtitle("Overlap Between DEG and Circadian Genes")

ggsave("Venn_DEG_Circadian.png", plot = venn_plot, width = 8, height = 8, dpi = 300)


dev.off()
intersect_genes <- intersect(deg_genes, circadian_genes)
cat("交集基因数量：", length(intersect_genes), "\n")
write.csv(intersect_genes, "Intersection_DEG_Circadian.csv")







