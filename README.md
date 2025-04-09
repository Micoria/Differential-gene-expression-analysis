# Differential-gene-expression-analysis

CLDN Family and Circadian Genes in Colorectal Cancer: Differential Expression and Functional Enrichment Analysis
📖 Project Overview
This project aims to explore the expression patterns and potential biological significance of the CLDN (Claudin) family genes and circadian rhythm-related genes in colorectal cancer (CRC).
Using the GSE39582 dataset from the GEO database, we performed differential expression analysis between CRC samples and normal controls, followed by functional enrichment analysis (GO and KEGG), and gene intersection analysis to reveal candidate genes potentially linking circadian regulation and the CLDN family in CRC.

🎯 Research Background
Colorectal cancer is one of the most common malignant tumors worldwide. Increasing evidence indicates that circadian rhythm disruptions may play a role in cancer development and progression. Meanwhile, the CLDN gene family, as key components of tight junctions, is crucial for maintaining epithelial barrier integrity, and their dysregulation is associated with cancer metastasis and progression.

This project aims to integrate these two important gene sets:

Circadian rhythm regulation genes from the MSigDB database.

CLDN family genes (CLDN1 to CLDN23).

We systematically analyze their expression profiles in CRC and explore their functional involvement through enrichment analysis.

📂 Data Source
GSE39582 Dataset

Downloaded from the GEO database.

Platform: Affymetrix Human Genome U133 Plus 2.0 Array (GPL570).

Samples: 566 CRC patients, 19 healthy controls.

Circadian rhythm genes

Extracted from MSigDB v7.5.1, KEGG Legacy collection.

CLDN family genes

Curated manually: CLDN1 ~ CLDN23.

🧩 Methods
Data Processing
Raw data were preprocessed: background correction, normalization, and probe ID conversion using hgu133plus2.db.

Removed probes without gene symbols and duplicates.

Annotated gene symbols and prepared expression matrix.

Differential Expression Analysis
Tool: limma R package.

Model: Linear model fitted, with empirical Bayes moderation.

Thresholds:

|log2(FC)| > 1

adj. P-value < 0.05

Result:

Upregulated genes: 934

Downregulated genes: 867

Functional Enrichment Analysis
Tool: clusterProfiler

Databases: GO (BP, CC, MF) and KEGG.

Thresholds:

P-value < 0.05

Q-value < 0.2

Visualization:

Dotplots and Barplots

GO Results Example:

KEGG Results Example:

Gene Set Intersection Analysis
Compared DEGs with circadian genes and CLDN family genes.

Tool: intersect() function and ggVennDiagram.

Intersection Results:

Circadian-related DEG: BHLHE41

CLDN-related DEGs: CLDN1, CLDN2, CLDN7, CLDN8

Venn Diagram Example: (You can upload your Venn diagram image and link it here.)

📊 Results Summary
Differential Expression Analysis

1801 DEGs identified (934 upregulated, 867 downregulated).

Robust dataset with high-quality annotations.

Functional Enrichment

GO analysis highlighted processes such as:

Mitotic nuclear division

Chromosome segregation

Regulation of cell cycle

KEGG pathways enriched:

Cell cycle

Focal adhesion

Chemical carcinogenesis - DNA adducts

Intersection Genes

CLDN family and circadian gene sets showed meaningful overlaps with DEGs.

CLDN1, CLDN2, CLDN7, CLDN8 were identified as DEGs.

BHLHE41 emerged as a core circadian gene of interest, showing differential expression in CRC.

💡 Interpretation
The results suggest that BHLHE41, a critical circadian rhythm gene, might have potential regulatory roles in colorectal cancer development.
Moreover, the identification of multiple CLDN genes as DEGs reinforces the significance of tight junction dysregulation in CRC progression. These findings provide a valuable foundation for further studies investigating the co-regulatory networks between circadian genes and CLDN family genes in CRC.

🔧 Tools & Environment
R version 4.2.3

R packages:

limma

clusterProfiler

org.Hs.eg.db

hgu133plus2.db

enrichplot

pheatmap

ggplot2

VennDiagram

ggVennDiagram

msigdbr

📌 Next Steps (For Future Work)
Conduct co-expression analysis between BHLHE41 and CLDN genes.

Explore survival analysis using CRC patient data.

Investigate protein-protein interaction (PPI) networks.

Validate findings through external datasets or experimental verification.

🙏 Acknowledgements
This study utilizes publicly available data from the GEO database and MSigDB, and makes extensive use of open-source R packages. We thank the contributors of these valuable resources.

