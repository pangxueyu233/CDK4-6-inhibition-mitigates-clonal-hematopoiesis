# **CDK4/6 inhibition mitigates chemotherapy-induced expansion of  TP53-mutant clonal hematopoiesis**

This page recorded the codes and data used and mentioned in [*xxx*](XXX). And you could downloaded this paper by clicking [here](pdf/XXX)

**ABSTRACT**

Therapy-related myeloid neoplasms (tMN) represent a fatal consequence of exposure to cytotoxic therapy administered in the treatment of cancer. Individuals with pre-existing *TP53* clonal hematopoiesis (CH) are at high risk for development of tMN but currently avoidance of therapy is the only strategy to reduce tMN risk. Here, we show in four randomized clinical trials of the CDK4/6 inhibitor trilaciclib given in conjunction with a variety of chemotherapeutic regimens and across diverse cancer patient populations that trilaciclib mitigates expansion of chemotherapy-related CH clones with mutations in DNA damage response genes (including *TP53*, *PPM1D*, and *CHEK2* mutations). This finding was also observed in a syngeneic murine model of *TP53* mutant CH demonstrating that trilaciclib blocks platinum-induced *TP53* competitive repopulation through promoting hematopoietic stem and progenitor quiescence and decreasing *TP53* mutant stemness advantage. This represents the first demonstration of a pharmacologic strategy to block chemotherapy-induced expansion of pre-leukemic *TP53*-mutant clones.

**1. Codes of analyzing and visualization**

- **[Figure making](Figuremaking.md): Visualization and Analysis of Hematopoietic Cell Populations in scRNA-seq Data**
  - This R script processes and visualizes single-cell RNA sequencing (scRNA-seq) data to investigate hematopoietic cell populations. It includes cell annotation refinement, dimensional reduction plots (openTSNE, UMAP), gene expression analysis, stemness and cell cycle profiling, and differential gene expression analysis. The script leverages various bioinformatics packages such as Seurat, monocle, tricycle, and clusterProfiler to analyze hematopoietic stem and progenitor cells (HSPCs), myeloid, and lymphoid lineages under different conditions, including chemotherapy and CDK inhibition. The figures generated illustrate cellular distributions, gene expression dynamics, and functional implications of stemness and differentiation.

# **2. Raw data download**

- **Description**: This section includes all the raw FASTQ files from our study. These files are crucial for in-depth data analysis and understanding the sequencing results from single-cell RNA.

- **Download**: You can access and download these files from the [GEO database](https://chat.openai.com/c/link-to-download).

Below is a detailed annotation of the file structure and contents:

```shell
[4.0K]  .
├── [4.0K]  scRNA
│   ├── [ 44G]  ACC10_RNA_S1_L001_R1_001.fastq.gz
3 directories, 42 files
```

Contents

- Each `_RNA.tar.gz` file includes the filtered_feature_bc_matrix output and loupe file from the Cell Ranger count model.

# **3. Processed Data Download**

## 3.1. CellRanger  Output

- **Description**: This section includes the output files from Cell Ranger and Space Ranger, essential for the initial data processing and analysis of single-cell RNA data.
- **Download**: These files are available for access and download from the [Zenodo1XXX](https://zenodo.org/records/).

Below is a detailed annotation of the file structure and contents:

```shell
tree -lh
[4.0K]  .
├── [219M]  All.CytoTRACE.DS.pool2.rds # Contains CytoTRACE scores for Pool2 dataset, estimating differentiation potential across cells.
├── [2.3G]  All_filter_miQC.pool1.merge.Sv4.rds #Merged and quality-filtered single-cell RNA-seq dataset for Pool1, processed using miQC for doublet and low-quality cell removal.
├── [2.5G]  All_filter_miQC.pool2.merge.Sv4.rds # Merged and quality-filtered single-cell RNA-seq dataset for Pool2, similar to Pool1 but a separate dataset.
├── [ 12K]  CD45.1.and.CD45.2.cells.propo.v2.rds #Contains proportions of CD45.1 and CD45.2 cells across different experimental conditions.
├── [1.7K]  CD45.2.cells.propo.v2.tricycleGroup.v2.rds #CD45.2 cell proportions categorized by tricycle-inferred cell cycle phases.
├── [107K]  Fig_New_T00_obj.markers.rds #Differential expression markers identified in the T00 object, used for figure generation.
├── [ 47K]  MPP.GO_v2.rds #Gene ontology enrichment results for multipotent progenitor (MPP) cells.
├── [155K]  Mus_musculus_TF.txt #List of transcription factors in Mus musculus (mouse), used for regulatory analysis.
├── [ 534]  mouse_cell_cycle_genes.rds #Predefined list of mouse cell cycle genes for cell cycle scoring in Seurat.
├── [ 23M]  msigdb.v7.1.symbols.gmt #MSigDB v7.1 gene set collections, used for gene set enrichment analysis.
├── [2.2G]  only_Myeolid.pool2.merge.Sv4.rds #Merged and filtered single-cell RNA-seq dataset containing only myeloid lineage cells from Pool2.
├── [3.2K]  only_T.cells.propo.v2.rds #Proportion data for T cell subtypes across different experimental conditions.
├── [ 80M]  only_T_cells.pool1.rds #Single-cell RNA-seq dataset containing only T cells from Pool1.
└── [ 11M]  pool2.v1_cell_anno.DEGs.rds #Differentially expressed genes (DEGs) annotated in Pool2 dataset.

0 directories, 14 files
```

# **Citation**

Our paper has been published in [*XXX Journal*](https://). For further reference and details, you can access the publication at the provided link.

The raw data supporting the findings of this study can be downloaded from the following repositories:

- **GEO Database**: Access our dataset by visiting [GSEXXX](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSEXXX). This link will take you directly to the dataset's page.
- **Zenodo**: Additional data files are available on Zenodo. Download them at [Zenodo1XXX](https://zenodo.org/records/).