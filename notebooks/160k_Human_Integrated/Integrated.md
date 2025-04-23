# Analysis of Single-cell Gene Expression Data <span style="font-size:20px">(integration) v2.0.0</span>

## Bioinformatics Core Facility, University of Manchester

1. [Prepare workspace](#1.-Prepare-workspace)
2. [Import data](#2.-Import-data)
3. [Prepare for data integration](#3.-Prepare-for-data-integration)
4. [Data integration](#4.-Data-integration)
5. [Clustering of cells](#5.-Clustering-of-cells)
6. [Expression of manually selected genes](#6.-Expression-of-manually-selected-genes)
7. [Define per-cluster cell types](#7.-Define-per-cluster-cell-types)
8. [Marker gene detection (per-cluster)](#8.-Marker-gene-detection-(per-cluster))
9. [DE analysis between conditions](#9.-DE-analysis-between-conditions)
10. [DE analysis between conditions in each cluster/cell type](#10.-DE-analysis-between-conditions-in-each-cluster/cell-type)
11. [Functional analysis using `enrichR`](#11.-Functional-analysis-using-enrichR)

# Project summary

*Add experimental details here. For example:*

<div class="alert alert-info">
    <b>Kidney and Lung DTCs Matched PBMC MultiPro® Human Discovery Panel</b> <a href ="https://www.10xgenomics.com/datasets/160k_DTC_Matched_PBMC_MultiPro_Human_Discovery_Panel">Link</a>
    <br /><br />
    DTCs and matched PBMCs were obtained from Discovery Life Sciences, including Stage I kidney cancer DTCs, Stage I kidney cancer matched PBMCs, Stage III kidney cancer DTCs, Stage III kidney cancer matched PBMCs, Stage III lung cancer DTCs, and Stage III lung cancer matched PBMCs. Control PBMCs were obtained from AllCells.
    <br /><br />
    Samples were prepared in accordance with the Cell Preparation Guide Handbook (CG00053) for Tips & Best Practices on handling and counting cells. Consult Cell Thawing Protocols for Single Cell Assays (CG000447) for guidance on thawing dissociated tumor cells.
    <br /><br />
    Antibody staining was performed as described in Intracellular Protein Labeling Protocol section of the Demonstrated Protocol Cell Surface & Intracellular Protein Labeling for GEM-X Flex Gene Expression (CG000781, Rev A) and Flex assay was performed as described in the GEM-X Flex Gene Expression Reagent Kits for Multiplexed Samples with Feature Barcode technology for Protein using Barcode Oligo Capture User Guide (CG000789). Libraries were sequenced on an Illumina NovaSeq 6000 with approximately 40,000 paired end reads/cell for Gene Expression and 20,000 paired end reads/cell for Feature Barcode libraries.
    <br /><br />
    Dual indexing libraries were sequenced with the following scheme:
    <ul>
        <li>28 cycles Read 1</li>
        <li>10 cycles i7</li>
        <li>10 cycles i5</li>
        <li>90 cycles Read 2</li>
    </ul>
</div>

## Sample summary

<a href="#Show-run-info">See below</a>.

# 1. Prepare workspace

**Tested on R version 4.4 and Bioconductor version 3.20**

- Set sample names and paths to the directories where the HDF5-based object were saved
- Load R libraries
- Set up colour palettes

## Define sample names and HDF5 paths


```R
# Set sample name
sample_names <- c("Control1",
                  "KidneyCancer",
                  "LungCancer")

# Set the paths to the directories where the HDF5-based object were saved, in the same order as sample names
h5_files <- c("../160k_Human_Control1_PBMC/Control1_h5_sce", 
              "../160k_Human_Kidney_Cancer_PBMC/KidneyCancer_h5_sce",
              "../160k_Human_Lung_Cancer_PBMC/LungCancer_h5_sce")

# Show info
h5_df <- data.frame(Sample = sample_names, Type = "Processed SingleCellExperiment", HDF5 = h5_files)
h5_df$Sample <- factor(h5_df$Sample, levels = sample_names)
h5_df
```


<table class="dataframe">
<caption>A data.frame: 3 × 3</caption>
<thead>
	<tr><th scope=col>Sample</th><th scope=col>Type</th><th scope=col>HDF5</th></tr>
	<tr><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><td>Control1    </td><td>Processed SingleCellExperiment</td><td>../160k_Human_Control1_PBMC/Control1_h5_sce         </td></tr>
	<tr><td>KidneyCancer</td><td>Processed SingleCellExperiment</td><td>../160k_Human_Kidney_Cancer_PBMC/KidneyCancer_h5_sce</td></tr>
	<tr><td>LungCancer  </td><td>Processed SingleCellExperiment</td><td>../160k_Human_Lung_Cancer_PBMC/LungCancer_h5_sce    </td></tr>
</tbody>
</table>



## Load libraries

The required R packages are listed below. This workflow has been tested on **R version 4.4 (Bioconductor version 3.20)** and latest versions of the packages supported in this R environment (see <a href=#Session-Info>Session Info</a>).  

Commented line below are packages that are required but we are not loading and attaching them.

<div class="alert alert-info">
    The <strong>scRUtils</strong> R package is current available only on <a href="https://github.com/ycl6/scRUtils" target="_blank">GitHub</a>. It can be installed using <code>devtools::install_github("ycl6/scRUtils")</code>.
</div>


```R
suppressPackageStartupMessages({
    # R-4.4.2
    library(batchelor)     # Single-cell batch correction methods
    library(BiocNeighbors) # AnnoyParam
    library(BiocParallel)  # MulticoreParam
    library(bluster)       # clusterRows, NNGraphParam, mergeCommunities
    library(cowplot)       # plotColData, plotRowData, plot_grid
    library(DESeq2)
    library(edgeR)
    library(enrichR)
    library(ggforce)       # gather_set_data, geom_parallel_sets*
    library(ggplot2)
    library(pheatmap)
    library(scales)
    library(scater)
    library(scran)
    library(viridis)       # scale_color_viridis

#    Access the exact function with "::" without load and attach package
#    library(BiocSingular)  # RandomParam & IrlbaParam
#    library(dplyr)         # select
#    library(ggplotify)     # as.ggplot
#    library(gtools)        # mixedsort
#    library(HDF5Array)     # loadHDF5SummarizedExperiment, saveHDF5SummarizedExperiment
#    library(igraph)        # cut_at
#    library(limma)         # vennDiagram, plotMDS
#    library(pals)          # cols25
#    library(RColorBrewer)  # brewer.pal
#    library(stringi)       # stri_join_list
#    library(stringr)       # str_to_title

    library(scRUtils)       # with customised functions
    library(tidyverse)
})
```

## Set default options


```R
# Set output window width
options(width = 110) # default 80; getOption("width")

# Set figure size; width = 10, height = 7, res = 120
fig()
```

Choose a `BiocParallel` Interface:
- **SerialParam**: Supported on all platforms. Evaluate BiocParallel-enabled code with parallel evaluation disabled.
- **MulticoreParam**: Supported on Unix and Mac. Evaluate BiocParallel-enabled code using multiple cores on a single computer.
- **SnowParam**: Supported on all platforms. Evaluate BiocParallel-enabled code across several distinct instances, on one or several computers.
- **BatchtoolsParam**: Applicable to clusters with formal schedulers. Evaluate BiocParallel-enabled code by submitting to a cluster scheduler like SGE.

<div class="alert alert-info">
    HDF5-based arrays only allow serial reading within a process. Therefore, some of the processing steps in this workflow are using <code>dgCMatrix</code> class matrices to improve performance, while keeping the assays in the <code>SingleCellExperiment</code> object in <code>DelayedMatrix</code> class.
</div>


```R
# Set number of threads for parallelisation
nthreads <- 8

# Choose BiocParallel Interface
bpp <- MulticoreParam(nthreads)
```


```R
# Initialise c30() and c40() palettes
c30 <- c30()
c40 <- c40()

fig(width = 14, height = 7)
par(mfrow = c(1,2), mar=c(0,0,1,0))
pie(rep(1,30), col = c30, radius = 1.5, main = "c30")
pie(rep(1,40), col = c40, radius = 1.5, main = "c40")
reset.fig()
```


    
![png](Integrated_files/Integrated_10_0.png)
    



```R
# Set colours for samples
c_sample_col <- choosePalette(h5_df$Sample, c30[11:13])
c_sample_col

# Set colours for cell cycle phases
c_phase_col <- setNames(c30[c(20, 3, 9)], c("G1", "S", "G2M"))
c_phase_col

# Set colours for heatmaps
c_heatmap_col1 <- plasma(256, direction = -1) # logcounts
breaks <- seq(-3, 3, by = 0.05) # 121 breaks
c_heatmap_col2 <- colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdYlBu")))(length(breaks)) # scaled/z-score
```


<style>
.dl-inline {width: auto; margin:0; padding: 0}
.dl-inline>dt, .dl-inline>dd {float: none; width: auto; display: inline-block}
.dl-inline>dt::after {content: ":\0020"; padding-right: .5ex}
.dl-inline>dt:not(:first-of-type) {padding-left: .5ex}
</style><dl class=dl-inline><dt>Control1</dt><dd>'#8fbc8f'</dd><dt>KidneyCancer</dt><dd>'#483d8b'</dd><dt>LungCancer</dt><dd>'#cd853f'</dd></dl>




<style>
.dl-inline {width: auto; margin:0; padding: 0}
.dl-inline>dt, .dl-inline>dd {float: none; width: auto; display: inline-block}
.dl-inline>dt::after {content: ":\0020"; padding-right: .5ex}
.dl-inline>dt:not(:first-of-type) {padding-left: .5ex}
</style><dl class=dl-inline><dt>G1</dt><dd>'#dc143c'</dd><dt>S</dt><dd>'#0000ff'</dd><dt>G2M</dt><dd>'#00ff00'</dd></dl>



# 2. Import data

## Load HDF5 objects

Here, we use the `loadHDF5SummarizedExperiment` function from the `HDF5Array` R package to import HDF5 files.


```R
all.sce <- list()
barcodes <- data.frame(Sample = NULL, Barcodes = NULL)
for(i in 1:length(h5_df$Sample)) {
    sample_name <- as.character(h5_df$Sample[i])
    message(paste("Read HDF5 files:", sample_name))
    sce <- HDF5Array::loadHDF5SummarizedExperiment(h5_df$HDF5[i])
    levels(sce$Sample) <- sample_name

    # Add Sample ID to Barcode as prefix to make cell barcodes unique
    colnames(sce) <- paste0(sce$Sample, "_", sce$Barcode)
    sce$Barcode <- colnames(sce)

    print(sce)
    flush.console()

    all.sce[[sample_name]] <- sce
}
```

    Read HDF5 files: Control1
    


    class: SingleCellExperiment 
    dim: 18129 18460 
    metadata(8): Samples cyclone ... DoubletDensity runInfo
    assays(2): counts logcounts
    rownames(18129): SAMD11 NOC2L ... MT-ND6 MT-CYB
    rowData names(12): ID Symbol ... pct_dropout is_hvg
    colnames(18460): Control1_AAACAAGCAACAAGTTACTTTAGG-1 Control1_AAACAAGCAACTAGTGACTTTAGG-1 ...
      Control1_TTTGTGAGTGGCGTAGACTTTAGG-1 Control1_TTTGTGAGTTAATTCGACTTTAGG-1
    colData names(23): Sample Barcode ... label log10Sum
    reducedDimNames(3): PCA TSNE UMAP
    mainExpName: Gene Expression
    altExpNames(1): Antibody Capture


    Read HDF5 files: KidneyCancer
    


    class: SingleCellExperiment 
    dim: 18129 17168 
    metadata(8): Samples cyclone ... DoubletDensity runInfo
    assays(2): counts logcounts
    rownames(18129): SAMD11 NOC2L ... MT-ND6 MT-CYB
    rowData names(12): ID Symbol ... pct_dropout is_hvg
    colnames(17168): KidneyCancer_AAACAAGCAAATACCGATGTTGAC-1 KidneyCancer_AAACAAGCAACAGATTATGTTGAC-1
      ... KidneyCancer_TTTGTGAGTGTCCTTCATGTTGAC-1 KidneyCancer_TTTGTGAGTTGGATGAATGTTGAC-1
    colData names(23): Sample Barcode ... label log10Sum
    reducedDimNames(3): PCA TSNE UMAP
    mainExpName: Gene Expression
    altExpNames(1): Antibody Capture


    Read HDF5 files: LungCancer
    


    class: SingleCellExperiment 
    dim: 18129 14037 
    metadata(8): Samples cyclone ... DoubletDensity runInfo
    assays(2): counts logcounts
    rownames(18129): SAMD11 NOC2L ... MT-ND6 MT-CYB
    rowData names(12): ID Symbol ... pct_dropout is_hvg
    colnames(14037): LungCancer_AAACAAGCAAGGCCTGAGCTGTGA-1 LungCancer_AAACAAGCACCTTTGGAGCTGTGA-1 ...
      LungCancer_TTTGTGAGTTGAGTCTAGCTGTGA-1 LungCancer_TTTGTGAGTTTACGACAGCTGTGA-1
    colData names(23): Sample Barcode ... label log10Sum
    reducedDimNames(3): PCA TSNE UMAP
    mainExpName: Gene Expression
    altExpNames(1): Antibody Capture


## Show run info


```R
rI <- metadata(all.sce[[1]])[["runInfo"]]
df1 <- data.frame(matrix(ncol = length(rI[[1]]), nrow = 0))
df2 <- data.frame(matrix(ncol = length(rI[[2]]), nrow = 0))

for(i in 1:length(all.sce))
{
    rI <- metadata(all.sce[[i]])[["runInfo"]]
    df1 <- rbind(df1, as.data.frame(t(matrix(rI[[1]]))))
    df2 <- rbind(df2, as.data.frame(t(matrix(rI[[2]]))))
}
colnames(df1) <- names(rI[[1]])
colnames(df2) <- names(rI[[2]])
df2 <- cbind(df1[, 1, drop = FALSE], df2) # Add Sample ID

# Add runInfo
runInfo <- list("Sample" = df1, "Data" = df2)

if (suppressPackageStartupMessages(requireNamespace("kableExtra")) && requireNamespace("IRdisplay")) {
    suppressPackageStartupMessages(library(kableExtra))
    df1 %>% kable("html") %>% kable_styling(font_size = 12, position = "left", full_width = FALSE) %>% 
    column_spec(1:ncol(df1), extra_css = "vertical-align: middle !important;") %>%
    add_header_above(c("Single-cell Sample Summary" = ncol(df1)), bold = TRUE, background = "blue", 
                     color = "white") %>% as.character() %>% IRdisplay::display_html()
    
    df2 %>% kable("html") %>% kable_styling(font_size = 12, position = "left", full_width = FALSE) %>% 
    column_spec(1:ncol(df2), extra_css = 'vertical-align: middle !important;') %>%
    add_header_above(c("Data Processing Summary" = ncol(df2)), bold = TRUE, background = "blue", 
                     color = "white") %>% as.character() %>% IRdisplay::display_html()    
} else {
    message("Single-cell Sample Summary")
    data.frame(df1)
    message("Data Processing Summary")
    data.frame(df2)
}
```


<table class="table" style="font-size: 12px; width: auto !important; ">
 <thead>
<tr><th style="border-bottom:hidden;padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; font-weight: bold; color: white !important;padding-right: 4px; padding-left: 4px; background-color: blue !important;" colspan="7"><div style="border-bottom: 1px solid #ddd; padding-bottom: 5px; ">Single-cell Sample Summary</div></th></tr>
  <tr>
   <th style="text-align:left;"> Sample ID </th>
   <th style="text-align:left;"> Platform </th>
   <th style="text-align:left;"> Chemistry </th>
   <th style="text-align:left;"> Reference </th>
   <th style="text-align:left;"> Transcriptome </th>
   <th style="text-align:left;"> Probe Set Name </th>
   <th style="text-align:left;"> Pipeline Version </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;vertical-align: middle !important;"> Control1 </td>
   <td style="text-align:left;vertical-align: middle !important;"> 10x Genomics </td>
   <td style="text-align:left;vertical-align: middle !important;"> Flex Gene Expression </td>
   <td style="text-align:left;vertical-align: middle !important;"> GRCh38-2024-A </td>
   <td style="text-align:left;vertical-align: middle !important;"> GRCh38-2024-A </td>
   <td style="text-align:left;vertical-align: middle !important;"> Chromium Human Transcriptome Probe Set v1.1.0 </td>
   <td style="text-align:left;vertical-align: middle !important;"> cellranger-9.0.0 </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align: middle !important;"> KidneyCancer </td>
   <td style="text-align:left;vertical-align: middle !important;"> 10x Genomics </td>
   <td style="text-align:left;vertical-align: middle !important;"> Flex Gene Expression </td>
   <td style="text-align:left;vertical-align: middle !important;"> GRCh38-2024-A </td>
   <td style="text-align:left;vertical-align: middle !important;"> GRCh38-2024-A </td>
   <td style="text-align:left;vertical-align: middle !important;"> Chromium Human Transcriptome Probe Set v1.1.0 </td>
   <td style="text-align:left;vertical-align: middle !important;"> cellranger-9.0.0 </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align: middle !important;"> LungCancer </td>
   <td style="text-align:left;vertical-align: middle !important;"> 10x Genomics </td>
   <td style="text-align:left;vertical-align: middle !important;"> Flex Gene Expression </td>
   <td style="text-align:left;vertical-align: middle !important;"> GRCh38-2024-A </td>
   <td style="text-align:left;vertical-align: middle !important;"> GRCh38-2024-A </td>
   <td style="text-align:left;vertical-align: middle !important;"> Chromium Human Transcriptome Probe Set v1.1.0 </td>
   <td style="text-align:left;vertical-align: middle !important;"> cellranger-9.0.0 </td>
  </tr>
</tbody>
</table>



<table class="table" style="font-size: 12px; width: auto !important; ">
 <thead>
<tr><th style="border-bottom:hidden;padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; font-weight: bold; color: white !important;padding-right: 4px; padding-left: 4px; background-color: blue !important;" colspan="9"><div style="border-bottom: 1px solid #ddd; padding-bottom: 5px; ">Data Processing Summary</div></th></tr>
  <tr>
   <th style="text-align:left;"> Sample ID </th>
   <th style="text-align:left;"> Lab/Facility </th>
   <th style="text-align:left;"> Analyst </th>
   <th style="text-align:left;"> Run ID </th>
   <th style="text-align:left;"> Run Name </th>
   <th style="text-align:left;"> Sequencer </th>
   <th style="text-align:left;"> User </th>
   <th style="text-align:left;"> PI </th>
   <th style="text-align:left;"> Organism </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;vertical-align: middle !important;"> Control1 </td>
   <td style="text-align:left;vertical-align: middle !important;"> UoM BCF </td>
   <td style="text-align:left;vertical-align: middle !important;"> I-Hsuan Lin </td>
   <td style="text-align:left;vertical-align: middle !important;"> SOME_RUN_ID </td>
   <td style="text-align:left;vertical-align: middle !important;"> SOME_RUN_NAME </td>
   <td style="text-align:left;vertical-align: middle !important;"> NovaSeq 6000 </td>
   <td style="text-align:left;vertical-align: middle !important;"> SOME_USER </td>
   <td style="text-align:left;vertical-align: middle !important;"> SOME_PI </td>
   <td style="text-align:left;vertical-align: middle !important;"> Human </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align: middle !important;"> KidneyCancer </td>
   <td style="text-align:left;vertical-align: middle !important;"> UoM BCF </td>
   <td style="text-align:left;vertical-align: middle !important;"> I-Hsuan Lin </td>
   <td style="text-align:left;vertical-align: middle !important;"> SOME_RUN_ID </td>
   <td style="text-align:left;vertical-align: middle !important;"> SOME_RUN_NAME </td>
   <td style="text-align:left;vertical-align: middle !important;"> NovaSeq 6000 </td>
   <td style="text-align:left;vertical-align: middle !important;"> SOME_USER </td>
   <td style="text-align:left;vertical-align: middle !important;"> SOME_PI </td>
   <td style="text-align:left;vertical-align: middle !important;"> Human </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align: middle !important;"> LungCancer </td>
   <td style="text-align:left;vertical-align: middle !important;"> UoM BCF </td>
   <td style="text-align:left;vertical-align: middle !important;"> I-Hsuan Lin </td>
   <td style="text-align:left;vertical-align: middle !important;"> SOME_RUN_ID </td>
   <td style="text-align:left;vertical-align: middle !important;"> SOME_RUN_NAME </td>
   <td style="text-align:left;vertical-align: middle !important;"> NovaSeq 6000 </td>
   <td style="text-align:left;vertical-align: middle !important;"> SOME_USER </td>
   <td style="text-align:left;vertical-align: middle !important;"> SOME_PI </td>
   <td style="text-align:left;vertical-align: middle !important;"> Human </td>
  </tr>
</tbody>
</table>


## Visualising with t-SNE and UMAP


```R
p <- list()
for(i in 1:length(all.sce)) {
    p[[i]] <- plotProjections(all.sce[[i]], "Sample", dimnames = c("TSNE", "UMAP"), point_alpha = 0.01, 
                              feat_color = c_sample_col, text_by = "label", text_size = 7, 
                              legend_pos = "none", show_subtitle = FALSE, theme_size = 14, 
                              titles = sprintf("%s (%s)", names(all.sce[i]), c("TSNE","UMAP")))
}

fig(width = 16, height = 8)
plot_grid(plotlist = p, ncol = 2, align = "h")
reset.fig()
```


    
![png](Integrated_files/Integrated_17_0.png)
    


# 3. Prepare for data integration

## Find common genes

Reduce the genes in each dataset to those common genes only


```R
universe1 <- Reduce(intersect, lapply(all.sce, rownames))
print(paste("Common genes:", length(universe1)))

all.common <- list()
for(i in 1:length(all.sce)) {
    all.common[[names(all.sce)[i]]] <- all.sce[[i]][match(universe1, rownames(all.sce[[i]])),]
}
all.common
```

    [1] "Common genes: 18129"



    $Control1
    class: SingleCellExperiment 
    dim: 18129 18460 
    metadata(8): Samples cyclone ... DoubletDensity runInfo
    assays(2): counts logcounts
    rownames(18129): SAMD11 NOC2L ... MT-ND6 MT-CYB
    rowData names(12): ID Symbol ... pct_dropout is_hvg
    colnames(18460): Control1_AAACAAGCAACAAGTTACTTTAGG-1 Control1_AAACAAGCAACTAGTGACTTTAGG-1 ...
      Control1_TTTGTGAGTGGCGTAGACTTTAGG-1 Control1_TTTGTGAGTTAATTCGACTTTAGG-1
    colData names(23): Sample Barcode ... label log10Sum
    reducedDimNames(3): PCA TSNE UMAP
    mainExpName: Gene Expression
    altExpNames(1): Antibody Capture
    
    $KidneyCancer
    class: SingleCellExperiment 
    dim: 18129 17168 
    metadata(8): Samples cyclone ... DoubletDensity runInfo
    assays(2): counts logcounts
    rownames(18129): SAMD11 NOC2L ... MT-ND6 MT-CYB
    rowData names(12): ID Symbol ... pct_dropout is_hvg
    colnames(17168): KidneyCancer_AAACAAGCAAATACCGATGTTGAC-1 KidneyCancer_AAACAAGCAACAGATTATGTTGAC-1
      ... KidneyCancer_TTTGTGAGTGTCCTTCATGTTGAC-1 KidneyCancer_TTTGTGAGTTGGATGAATGTTGAC-1
    colData names(23): Sample Barcode ... label log10Sum
    reducedDimNames(3): PCA TSNE UMAP
    mainExpName: Gene Expression
    altExpNames(1): Antibody Capture
    
    $LungCancer
    class: SingleCellExperiment 
    dim: 18129 14037 
    metadata(8): Samples cyclone ... DoubletDensity runInfo
    assays(2): counts logcounts
    rownames(18129): SAMD11 NOC2L ... MT-ND6 MT-CYB
    rowData names(12): ID Symbol ... pct_dropout is_hvg
    colnames(14037): LungCancer_AAACAAGCAAGGCCTGAGCTGTGA-1 LungCancer_AAACAAGCACCTTTGGAGCTGTGA-1 ...
      LungCancer_TTTGTGAGTTGAGTCTAGCTGTGA-1 LungCancer_TTTGTGAGTTTACGACAGCTGTGA-1
    colData names(23): Sample Barcode ... label log10Sum
    reducedDimNames(3): PCA TSNE UMAP
    mainExpName: Gene Expression
    altExpNames(1): Antibody Capture



## Organise data structure

Make all `sce` objects have identical `rowData` and `colData` column names and order.

<div class="alert alert-warning">
    <strong>Warning!</strong> You may add more columns to show additional features in the integrated object if the columns are present in <u>all the samples</u>.
</div>


```R
# Select rowData columns to keep
rowData_colnames <- c("ID","Symbol","Type","SEQNAME","is_mito")

# Select colData columns to keep
# Add 'coarse_cell_type' and 'fine_cell_type' if using Cell Ranger's Cell Annotation prediction
colData_colnames <- c("Sample","Barcode","sum","detected","subsets_Mt_percent","log10Sum","sizeFactor","CellCycle",
                      #"coarse_cell_type","fine_cell_type",
                      "CellType","walktrap","louvain","leiden","label")

# Make all sce have identical rowData and colData column names and order
for(i in 1:length(all.common)) {
    rowData(all.common[[i]]) <- rowData(all.common[[i]])[, rowData_colnames]
    colData(all.common[[i]]) <- colData(all.common[[i]])[, colData_colnames]

    # Clear old reducedDims
    reducedDims(all.common[[i]]) <- list()
}
```


```R
all.common
```


    $Control1
    class: SingleCellExperiment 
    dim: 18129 18460 
    metadata(8): Samples cyclone ... DoubletDensity runInfo
    assays(2): counts logcounts
    rownames(18129): SAMD11 NOC2L ... MT-ND6 MT-CYB
    rowData names(5): ID Symbol Type SEQNAME is_mito
    colnames(18460): Control1_AAACAAGCAACAAGTTACTTTAGG-1 Control1_AAACAAGCAACTAGTGACTTTAGG-1 ...
      Control1_TTTGTGAGTGGCGTAGACTTTAGG-1 Control1_TTTGTGAGTTAATTCGACTTTAGG-1
    colData names(13): Sample Barcode ... leiden label
    reducedDimNames(0):
    mainExpName: Gene Expression
    altExpNames(1): Antibody Capture
    
    $KidneyCancer
    class: SingleCellExperiment 
    dim: 18129 17168 
    metadata(8): Samples cyclone ... DoubletDensity runInfo
    assays(2): counts logcounts
    rownames(18129): SAMD11 NOC2L ... MT-ND6 MT-CYB
    rowData names(5): ID Symbol Type SEQNAME is_mito
    colnames(17168): KidneyCancer_AAACAAGCAAATACCGATGTTGAC-1 KidneyCancer_AAACAAGCAACAGATTATGTTGAC-1
      ... KidneyCancer_TTTGTGAGTGTCCTTCATGTTGAC-1 KidneyCancer_TTTGTGAGTTGGATGAATGTTGAC-1
    colData names(13): Sample Barcode ... leiden label
    reducedDimNames(0):
    mainExpName: Gene Expression
    altExpNames(1): Antibody Capture
    
    $LungCancer
    class: SingleCellExperiment 
    dim: 18129 14037 
    metadata(8): Samples cyclone ... DoubletDensity runInfo
    assays(2): counts logcounts
    rownames(18129): SAMD11 NOC2L ... MT-ND6 MT-CYB
    rowData names(5): ID Symbol Type SEQNAME is_mito
    colnames(14037): LungCancer_AAACAAGCAAGGCCTGAGCTGTGA-1 LungCancer_AAACAAGCACCTTTGGAGCTGTGA-1 ...
      LungCancer_TTTGTGAGTTGAGTCTAGCTGTGA-1 LungCancer_TTTGTGAGTTTACGACAGCTGTGA-1
    colData names(13): Sample Barcode ... leiden label
    reducedDimNames(0):
    mainExpName: Gene Expression
    altExpNames(1): Antibody Capture



## Recomputes log-normalized expression values

The `multiBatchNorm` function recomputes log-normalized expression values (using `logNormCounts`) after adjusting the size factors for systematic differences in coverage between `SingleCellExperiment` objects.


```R
set.seed(12345)
all.common.normed <- do.call(multiBatchNorm, all.common)
all.common.normed
```


    $Control1
    class: SingleCellExperiment 
    dim: 18129 18460 
    metadata(8): Samples cyclone ... DoubletDensity runInfo
    assays(2): counts logcounts
    rownames(18129): SAMD11 NOC2L ... MT-ND6 MT-CYB
    rowData names(5): ID Symbol Type SEQNAME is_mito
    colnames(18460): Control1_AAACAAGCAACAAGTTACTTTAGG-1 Control1_AAACAAGCAACTAGTGACTTTAGG-1 ...
      Control1_TTTGTGAGTGGCGTAGACTTTAGG-1 Control1_TTTGTGAGTTAATTCGACTTTAGG-1
    colData names(13): Sample Barcode ... leiden label
    reducedDimNames(0):
    mainExpName: Gene Expression
    altExpNames(1): Antibody Capture
    
    $KidneyCancer
    class: SingleCellExperiment 
    dim: 18129 17168 
    metadata(8): Samples cyclone ... DoubletDensity runInfo
    assays(2): counts logcounts
    rownames(18129): SAMD11 NOC2L ... MT-ND6 MT-CYB
    rowData names(5): ID Symbol Type SEQNAME is_mito
    colnames(17168): KidneyCancer_AAACAAGCAAATACCGATGTTGAC-1 KidneyCancer_AAACAAGCAACAGATTATGTTGAC-1
      ... KidneyCancer_TTTGTGAGTGTCCTTCATGTTGAC-1 KidneyCancer_TTTGTGAGTTGGATGAATGTTGAC-1
    colData names(13): Sample Barcode ... leiden label
    reducedDimNames(0):
    mainExpName: Gene Expression
    altExpNames(1): Antibody Capture
    
    $LungCancer
    class: SingleCellExperiment 
    dim: 18129 14037 
    metadata(8): Samples cyclone ... DoubletDensity runInfo
    assays(2): counts logcounts
    rownames(18129): SAMD11 NOC2L ... MT-ND6 MT-CYB
    rowData names(5): ID Symbol Type SEQNAME is_mito
    colnames(14037): LungCancer_AAACAAGCAAGGCCTGAGCTGTGA-1 LungCancer_AAACAAGCACCTTTGGAGCTGTGA-1 ...
      LungCancer_TTTGTGAGTTGAGTCTAGCTGTGA-1 LungCancer_TTTGTGAGTTTACGACAGCTGTGA-1
    colData names(13): Sample Barcode ... leiden label
    reducedDimNames(0):
    mainExpName: Gene Expression
    altExpNames(1): Antibody Capture




```R
p <- list()
for(i in 1:length(all.common)) {
    p[[i]] <- ggplot(data.frame(Original = all.common[[i]]$sizeFactor, Normed = all.common.normed[[i]]$sizeFactor),
                     aes(x = Original, y = Normed)) + geom_point(size = 0.5) + 
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = 2, linewidth = 1) +
    theme_bw(base_size = 16) + ggtitle(paste(names(all.common[i]), "sf"))
}

fig(width = 16, height = 4)
plot_grid(plotlist = p, ncol = 4, align = "h")
reset.fig()
```


    
![png](Integrated_files/Integrated_25_0.png)
    


## Identifying highly variable genes (HVGs)

Identifying HVG based on common genes using re-normalised `SingleCellExperiment` objects.

### Quantifying per-gene variation

Choose to model the variance of the log-expression profiles for each gene (`modelGeneVar`) or model the per-gene count variance with Poisson noise (`modelGeneVarByPoisson`).


```R
# Choose variance modelling method to use
hvg_model <- "modelGeneVarByPoisson" # Or modelGeneVar

message(paste0("Using '", hvg_model, "' method"))
if(hvg_model == "modelGeneVar") {
    set.seed(12345)
    all.var <- lapply(all.common.normed, modelGeneVar, assay.type = "logcounts", BPPARAM = bpp)
} else {
    set.seed(12345)
    all.var <- lapply(all.common.normed, modelGeneVarByPoisson, assay.type = "counts", BPPARAM = bpp)
}

# Print DataFrame
for(i in 1:length(h5_df$Sample)) {
    print(paste0(h5_df$Sample[i], ":"))
    all.var[[i]] %>% as.data.frame %>% arrange(FDR, desc(bio)) %>% dplyr::select(1:6) %>% DataFrame %>% print
}
```

    Using 'modelGeneVarByPoisson' method
    


    [1] "Control1:"
    DataFrame with 18129 rows and 6 columns
                 mean     total      tech       bio   p.value       FDR
            <numeric> <numeric> <numeric> <numeric> <numeric> <numeric>
    LYZ       1.45634   4.59072  0.632304   3.95841         0         0
    CD74      3.76971   4.04530  0.178336   3.86697         0         0
    GNLY      1.20764   4.16534  0.640817   3.52452         0         0
    NKG7      1.49643   3.89281  0.627994   3.26482         0         0
    S100A9    1.12832   3.82526  0.634779   3.19048         0         0
    ...           ...       ...       ...       ...       ...       ...
    TGIF2LY         0         0         0         0       NaN       NaN
    PCDH11Y         0         0         0         0       NaN       NaN
    AMELY           0         0         0         0       NaN       NaN
    TSPY1           0         0         0         0       NaN       NaN
    DAZ2            0         0         0         0       NaN       NaN
    [1] "KidneyCancer:"
    DataFrame with 18129 rows and 6 columns
                 mean     total      tech       bio   p.value       FDR
            <numeric> <numeric> <numeric> <numeric> <numeric> <numeric>
    GNLY     1.138577   3.42077  0.486239   2.93453         0         0
    FOS      2.051771   2.70202  0.380513   2.32151         0         0
    CCL5     1.212722   2.68365  0.486587   2.19706         0         0
    CD74     1.797954   2.55010  0.424426   2.12568         0         0
    NKG7     0.993822   2.47213  0.475249   1.99688         0         0
    ...           ...       ...       ...       ...       ...       ...
    PCDH11Y         0         0         0         0       NaN       NaN
    AMELY           0         0         0         0       NaN       NaN
    TBL1Y           0         0         0         0       NaN       NaN
    TSPY1           0         0         0         0       NaN       NaN
    DAZ2            0         0         0         0       NaN       NaN
    [1] "LungCancer:"
    DataFrame with 18129 rows and 6 columns
                 mean     total      tech       bio   p.value       FDR
            <numeric> <numeric> <numeric> <numeric> <numeric> <numeric>
    GNLY     1.720334   5.00735  0.452945   4.55441         0         0
    CD74     2.968303   4.47523  0.231526   4.24371         0         0
    CCL5     2.221786   4.47919  0.362607   4.11658         0         0
    IGHM     0.940615   4.10908  0.482037   3.62705         0         0
    NKG7     1.662264   3.66641  0.461710   3.20470         0         0
    ...           ...       ...       ...       ...       ...       ...
    PCDH11Y         0         0         0         0       NaN       NaN
    AMELY           0         0         0         0       NaN       NaN
    TBL1Y           0         0         0         0       NaN       NaN
    TSPY1           0         0         0         0       NaN       NaN
    DAZ2            0         0         0         0       NaN       NaN


### Selecting HVGs

We perform feature selection by averaging the variance components across all batches with the `combineVar` function. We compute the average as it is responsive to batch-specific HVGs while still preserving the within-batch ranking of genes. In contrast, approaches based on taking the **intersection** or **union** of HVGs across batches become increasingly conservative or liberal, respectively, with an increasing number of batches.


```R
# Select the top N genes with the highest biological components
#hvg_genes <- Reduce(intersect, lapply(all.var, getTopHVGs, n = 5000, var.field = "bio", var.threshold = 0))

# Select the top N% genes with the highest biological components
#hvg_genes <- Reduce(intersect, lapply(all.var, getTopHVGs, prop = 0.20, var.field = "bio", var.threshold = 0))
```

When integrating datasets of variable composition, it is generally safer to err on the side of **including more genes** than are used in a single dataset analysis, to ensure that markers are retained for any dataset-specific subpopulations that might be present. 


```R
combined.var <- do.call(combineVar, all.var)

# Select the top N genes with the highest biological components
hvg_genes <- getTopHVGs(combined.var, n = 3500, var.field = "bio", var.threshold = 0)

# Select genes with positive average biological components
#hvg_genes <- rownames(combined.var)[combined.var$bio > 0]

# Select the top N% genes with the highest biological components
#hvg_genes <- getTopHVGs(combined.var, prop = 0.20, var.field = "bio", var.threshold = 0)
```


```R
message(sprintf("Top %d HVGs selected (%s):", length(hvg_genes), hvg_model))
head(hvg_genes, 20)
```

    Top 3500 HVGs selected (modelGeneVarByPoisson):
    



<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'GNLY'</li><li>'CD74'</li><li>'NKG7'</li><li>'FOS'</li><li>'CCL5'</li><li>'IGHM'</li><li>'LYZ'</li><li>'TRBC2'</li><li>'S100A9'</li><li>'DUSP1'</li><li>'H1-10'</li><li>'IL32'</li><li>'EFHD2'</li><li>'VIM'</li><li>'JUN'</li><li>'IL7R'</li><li>'SRGN'</li><li>'TCF7'</li><li>'IFI30'</li><li>'SPI1'</li></ol>



<div class="alert alert-info">
    <strong>Tip!</strong> In the accompanied Shiny App, top N HVGs are used as example genes in multi-gene expression visualisation.
</div>


```R
# Add to runInfo
runInfo <- c(runInfo, list(
    "HVG" = list(
        "Method" = hvg_model, 
        "Estimates" = combined.var, 
        "Genes" = hvg_genes)
    )
)
```

# 4. Data integration

Large single-cell RNA sequencing (scRNA-seq) projects usually need to generate data across multiple batches due to logistical constraints. However, the processing of different batches is often subject to uncontrollable differences, and results in systematic differences in the observed expression in cells from different batches. We also often observe a strong donor effect in integrated datasets. This might be due to differences in cell type composition between donors, but the more likely explanation is that of a technical difference in sample preparation processing or uninteresting genotypic differences.

Computational correction of batch effects is critical for eliminating batch-to-batch variation, allowing data across multiple batches to be combined for common downstream analysis. See [OSCA reference](https://bioconductor.org/books/3.20/OSCA.multisample/integrating-datasets.html)

## Diagnosing batch effects

Before performing any correction, we check that there actually is a batch effect across these datasets by checking that they cluster separately. Here, we use the `cbind` function to combine `SingleCellExperiment` objects in the list without applying any correction, and we informally verify that cells from different batches are separated using a t-SNE plot.

It is also possible to combine `SingleCellExperiment` objects with the `correctExperiments` function, and without applying any correction using the `NoCorrectParam` flag. However, data stored in `metadata` and `altExps` is not retained.

<div class="alert alert-warning">
    <strong>Warning!</strong>
    <br />
    If it takes longer than 1 minute for the merging to complete, it means there is a problem with <code>cbind()</code> concerning the columns in either <code>rowData()</code> and/or <code>colData()</code>. Interrupt the kernel and check for consistency of the columns of the objects stored in <code>all.common.normed</code>. Alternatively, use <code>correctExperiments()</code> to merge data and it'll throw out warning messages if there's a problem with the data structure.
</div>

### Create a single `SingleCellExperiment` object


```R
#set.seed(12345)
#combined <- correctExperiments(all.common.normed, PARAM = NoCorrectParam())
# Use `cbind` to combine SingleCellExperiment to retain "metadata" and altExps"
combined <- do.call(cbind, all.common.normed)

# Sort meta data
combined$Sample <- factor(combined$Sample, levels = gtools::mixedsort(levels(combined$Sample)))
combined$CellType <- factor(combined$CellType, levels = gtools::mixedsort(levels(combined$CellType)))

# Remove "merged" matrix in the assays slot, as it is identical to the "logcounts" matrix
assay(combined, "merged") <- NULL
combined
```


    class: SingleCellExperiment 
    dim: 18129 49665 
    metadata(24): Samples cyclone ... DoubletDensity runInfo
    assays(2): counts logcounts
    rownames(18129): SAMD11 NOC2L ... MT-ND6 MT-CYB
    rowData names(5): ID Symbol Type SEQNAME is_mito
    colnames(49665): Control1_AAACAAGCAACAAGTTACTTTAGG-1 Control1_AAACAAGCAACTAGTGACTTTAGG-1 ...
      LungCancer_TTTGTGAGTTGAGTCTAGCTGTGA-1 LungCancer_TTTGTGAGTTTACGACAGCTGTGA-1
    colData names(13): Sample Barcode ... leiden label
    reducedDimNames(0):
    mainExpName: Gene Expression
    altExpNames(1): Antibody Capture


Fix duplicate `metadata` names after combining datasets


```R
# If metadata in single-sample data is identical
if(length(metadata(combined)) %% length(names(all.common.normed)) == 0) {
    n_meta <- length(metadata(combined)) / length(names(all.common.normed))
    for(i in 1:length(names(all.common.normed))) {
        name <- names(all.common.normed)[i]
#        print(name)
        to <- i * n_meta
        from <- to - n_meta + 1
        names(metadata(combined))[from:to] <- paste0(name, "_", names(metadata(combined))[from:to])
    }
} else {
    stop("Not fixed automatically!")
}
```


```R
# Set colours for cell types
c_celltype_col <- choosePalette(combined$CellType, c40)
#c_celltype_col
```

### Add additional variables that describes the samples  to `colData`

<div class="alert alert-warning">
    <strong>Warning!</strong> Need to add additional information about your samples here.
</div>


```R
combined$condition <- combined$Sample
levels(combined$condition) <- c("Control","Cancer","Cancer","Cancer")

# Set colours for condition
c_cond_col <- setNames(c("cyan", "magenta"), c("Control", "Cancer"))
c_cond_col
```


<style>
.dl-inline {width: auto; margin:0; padding: 0}
.dl-inline>dt, .dl-inline>dd {float: none; width: auto; display: inline-block}
.dl-inline>dt::after {content: ":\0020"; padding-right: .5ex}
.dl-inline>dt:not(:first-of-type) {padding-left: .5ex}
</style><dl class=dl-inline><dt>Control</dt><dd>'cyan'</dd><dt>Cancer</dt><dd>'magenta'</dd></dl>



Create a `dgCMatrix` logcounts matrix for faster computation in some steps.


```R
normed_l <- list()
for(i in names(all.common.normed)) {
    normed_l[[i]] <-  as(logcounts(all.common.normed[[i]], withDimnames = FALSE), "dgCMatrix")
    rownames(normed_l[[i]]) <- rownames(all.common.normed[[i]])
}
```

### Run `multiBatchPCA` on the *uncorrected* `logcounts` matrix

The `multiBatchPCA` function perform a principal components analysis across multiple gene expression matrices (i.e. batches) to project all cells to a common low-dimensional space. Default to `d = 50` dimensions.

**weights**: A numeric vector, logical scalar or list specifying the weighting scheme to use. See `?multiBatchPCA` for more details.


```R
set.seed(12345)

# Use Randomized SVD algorithm, faster for file-backed matrices
#pca <- multiBatchPCA(normed_l, subset.row = hvg_genes, get.variance = TRUE, 
#                     weights = rep(1, length(normed_l)), 
#                     BSPARAM = BiocSingular::RandomParam(), BPPARAM = bpp)

# Use IRLBA algorithm, default behavior is more accurate
pca <- multiBatchPCA(normed_l, subset.row = hvg_genes, get.variance = TRUE, 
                     #weights = rep(1, length(normed_l)), 
                     BSPARAM = BiocSingular::IrlbaParam(), BPPARAM = bpp)
```


```R
# Re-construct pca object and attributes
rotation <- metadata(pca)$rotation
colnames(rotation) <- paste0("PC", seq_len(ncol(rotation)))
varExplained <- metadata(pca)$var.explained
percentVar <- varExplained/metadata(pca)$var.total*100

pca <- unlist(pca, use.names = FALSE)
colnames(pca) <- paste0("PC", seq_len(ncol(pca)))
attr(pca, "varExplained") <- varExplained
attr(pca, "percentVar") <- percentVar
attr(pca, "rotation") <- rotation
str(pca)

# Store PCA results
reducedDim(combined, "PCA") <- pca
```

     num [1:49665, 1:50] 5.5265 -0.8773 -0.0921 -3.2691 6.8074 ...
     - attr(*, "dimnames")=List of 2
      ..$ : chr [1:49665] "Control1_AAACAAGCAACAAGTTACTTTAGG-1" "Control1_AAACAAGCAACTAGTGACTTTAGG-1" "Control1_AAACAAGCAGTTATCCACTTTAGG-1" "Control1_AAACAAGCATAGCCGGACTTTAGG-1" ...
      ..$ : chr [1:50] "PC1" "PC2" "PC3" "PC4" ...
     - attr(*, "varExplained")= num [1:50] 158.3 89 67.1 52 32.7 ...
     - attr(*, "percentVar")= num [1:50] 8.64 4.85 3.66 2.84 1.78 ...
     - attr(*, "rotation")= num [1:3500, 1:50] 0.00311 0.10852 0.01264 0.07427 -0.01582 ...
      ..- attr(*, "dimnames")=List of 2
      .. ..$ : chr [1:3500] "GNLY" "CD74" "NKG7" "FOS" ...
      .. ..$ : chr [1:50] "PC1" "PC2" "PC3" "PC4" ...



```R
show <- 20 # Show first N PCs

# Scree plot
as.data.frame(percentVar) %>% rownames_to_column("PC") %>% mutate(PC = as.numeric(PC)) %>% head(show) %>%
    ggplot(aes(x = PC, y = percentVar)) + geom_point(size = 3) +
    scale_x_continuous(breaks = seq(0, 50, by = 5)) + theme_cowplot(20) + 
    guides(color = guide_legend(title = "Method")) + theme(legend.position = "top") + 
    labs(x = "PC", y = "Variance explained (%)")
```


    
![png](Integrated_files/Integrated_47_0.png)
    


### Run `runTSNE` and `runUMAP` to find low-dimensional representations of the *uncorrected* data

For `n_dimred`, better to include more PCs in integrated dataset than that used in a single dataset analysis. The default is 50, i.e. all dimensions.


```R
set.seed(12345)
combined <- runTSNE(combined, dimred = "PCA", name = "TSNE", n_dimred = 50, n_threads = nthreads, BPPARAM = bpp)
```


```R
set.seed(12345)
combined <- runUMAP(combined, dimred = "PCA", name = "UMAP", n_dimred = 50,
                    n_neighbors = 30, spread = 1, min_dist = 0.3, n_threads = nthreads, BPPARAM = bpp)
```


```R
# Added PCA, TSNE, UMAP reducedDims
combined
```


    class: SingleCellExperiment 
    dim: 18129 49665 
    metadata(24): Control1_Samples Control1_cyclone ... LungCancer_DoubletDensity LungCancer_runInfo
    assays(2): counts logcounts
    rownames(18129): SAMD11 NOC2L ... MT-ND6 MT-CYB
    rowData names(5): ID Symbol Type SEQNAME is_mito
    colnames(49665): Control1_AAACAAGCAACAAGTTACTTTAGG-1 Control1_AAACAAGCAACTAGTGACTTTAGG-1 ...
      LungCancer_TTTGTGAGTTGAGTCTAGCTGTGA-1 LungCancer_TTTGTGAGTTTACGACAGCTGTGA-1
    colData names(14): Sample Barcode ... label condition
    reducedDimNames(3): PCA TSNE UMAP
    mainExpName: Gene Expression
    altExpNames(1): Antibody Capture


### Visualise the *uncorrected* data in t-SNE and UMAP


```R
fig(width = 16, height = 8)
plotProjections(combined, "Sample", dimnames = c("TSNE", "UMAP"), feat_desc = "Sample", 
                feat_color = c_sample_col, guides_size = 4, point_size = 0.1, point_alpha = 0.1,
                legend_pos = "bottom", guides_nrow = 1, titles = c("TSNE, no correction","UMAP, no correction"))
reset.fig()
```


    
![png](Integrated_files/Integrated_53_0.png)
    



```R
fig(width = 16, height = 8)
plotProjections(combined, "CellType", dimnames = c("TSNE", "UMAP"), feat_desc = "Cell Type", 
                feat_color = c_celltype_col, guides_size = 4, point_size = 0.1, point_alpha = 0.1,
                legend_pos = "bottom", guides_nrow = 1, titles = c("TSNE, no correction","UMAP, no correction"))
reset.fig()
```


    
![png](Integrated_files/Integrated_54_0.png)
    


<div class="alert alert-info">
    <strong>About this dataset: </strong> (edits required)
    <ul>
        <li>Looks like the batch effect is very significant.</li>
        <li>We will perform correction and see how the connected data looks like.</li>
    </ul>
</div>

## Performing MNN correction

We use a fast version of the mutual nearest neighbors (MNN) method, `fastMNN`, in the `batchelor` package to correct for batch effects in single-cell expression data. The batch correction functions in batchelor are organized into three levels:

1. At the lowest level, we have the functions that implement a single type of correction, e.g., `fastMNN`.
2. At the next level, we have the `batchCorrect` function that wraps the single-correction function into a common interface.
3. At the highest level, we have the `correctExperiments` function that wraps the `batchCorrect` function.

Here, we run the `batchCorrect` function with the `FastMnnParam` flag to apply the `fastMNN` correction. 

The `multiBatchPCA` function is called when applying `fastMNN` to perform a PCA across multiple batches to project all cells to a common low-dimensional space.

### Controlling the merge order

By default, batches are merged in the user-supplied order in the input `SingleCellExperiment` objects or list. The merge order may occasionally be important as it determines the number of MNN pairs available at each merge step. MNN pairs results in greater stability of the batch vectors and increased likelihood of identifying shared subpopulations, which are important to the precision and accuracy of the MNN-based correction, respectively. The order can be changed by setting the `merge.order` or `auto.merge` arguments. See `?fastMNN` for more details.

**merge.order**: An integer vector containing the linear merge order of batches. Alternatively, a list of lists representing a tree structure specifying a hierarchical merge order. For example, a hierarchical merge to first merge together replicates with the same genotype, and then merge samples across different genotypes. See [Example](https://bioconductor.org/books/3.20/OSCA.multisample/chimeric-mouse-embryo-10x-genomics.html#merging)

**auto.merge**: Use `auto.merge=TRUE` to instruct `fastMNN` automatically identify the “best” merge order (i.e. largest number of MNNs). The aim is to improve the stability of the correction by first merging more similar batches with more MNN pairs. This can be somewhat time-consuming as MNN pairs need to be iteratively recomputed for all possible batch pairings.

Option 1: Use below code to use `auto.merge = TRUE`.


```R
set.seed(12345)
corrected <- batchCorrect(normed_l, subset.row = hvg_genes, 
                          PARAM = FastMnnParam(auto.merge = TRUE, 
                                               #weights = rep(1, length(normed_l)), # 1:1
                                               # use Randomized SVD or IRLBA for PCA
                                               #BSPARAM = BiocSingular::RandomParam(), # Use RandomParam
                                               BSPARAM = BiocSingular::IrlbaParam(),   # Or IrlbaParam
                                               BNPARAM = AnnoyParam(), BPPARAM = bpp))
```

Option 2: Use below code and use `merge.order` to set the merge order if your choice.


```R
#set.seed(12345)
#corrected <- batchCorrect(normed_l, subset.row = hvg_genes, 
#                          PARAM = FastMnnParam(merge.order = list(1, list(2, 3)), 
#                                               weights = rep(1, length(normed_l)), # 1:1
#                                               # use Randomized SVD or IRLBA for PCA
#                                               #BSPARAM = BiocSingular::RandomParam(), # Use RandomParam
#                                               BSPARAM = BiocSingular::IrlbaParam(),   # Or IrlbaParam
#                                               BNPARAM = AnnoyParam(), BPPARAM = bpp))
```

The returned object contains:

- A `batch` column in the `colData` slot, containing the batch of origin for each row (i.e., cell) in corrected.
- A `rotation` column the `rowData` slot, containing the rotation matrix used for the PCA.
- A `reconstructed` matrix in the `assays` slot, containing the low-rank reconstruction of the expression matrix. This can be interpreted as per-gene corrected log-expression values (after cosine normalization, if `cos.norm=TRUE` which is the default setting) but should not be used for quantitative analyses.
- A `corrected` matrix in the `reducedDims` slot, containing corrected low-dimensional coordinates for each cell. 

Additionally, the metadata of the output object contains:

- `merge.info`, a DataFrame of diagnostic information about each merge step.
- `pca.info`, a list of metadata produced by `multiBatchPCA`, such as the variance explained when `get.variance=TRUE`.


```R
corrected
```


    class: SingleCellExperiment 
    dim: 3500 49665 
    metadata(2): merge.info pca.info
    assays(1): reconstructed
    rownames(3500): GNLY CD74 ... PPP2R2A CD244
    rowData names(1): rotation
    colnames(49665): Control1_AAACAAGCAACAAGTTACTTTAGG-1 Control1_AAACAAGCAACTAGTGACTTTAGG-1 ...
      LungCancer_TTTGTGAGTTGAGTCTAGCTGTGA-1 LungCancer_TTTGTGAGTTTACGACAGCTGTGA-1
    colData names(1): batch
    reducedDimNames(1): corrected
    mainExpName: NULL
    altExpNames(0):


### Check the proportion of variance lost

For `fastMNN()`, one useful diagnostic is the proportion of variance within each batch that is lost during MNN correction. Specifically, this refers to the within-batch variance that is removed during orthogonalization with respect to the average correction vector at each merge step. This is returned via the `lost.var` field in the metadata of returned object, which contains a matrix of the variance lost in each batch (column) at each merge step (row).

Large proportions of lost variance (>10%) suggest that correction is removing genuine biological heterogeneity. This would occur due to violations of the assumption of orthogonality between the batch effect and the biological subspace. If the proportion of lost variance is small, it indicates that non-orthogonality is not a major concern.


```R
metadata(corrected)$merge.info
```


    DataFrame with 2 rows and 6 columns
                         left        right                               pairs batch.size   skipped
                       <List>       <List>                     <DataFrameList>  <numeric> <logical>
    1              LungCancer KidneyCancer 35629:20041,35630:29332,35630:31630   0.730185     FALSE
    2 LungCancer,KidneyCancer     Control1  35644:12864,35644:12006,35644:1220   0.744159     FALSE
                           lost.var
                           <matrix>
    1 0.0000000:0.0116568:0.0117827
    2 0.0112596:0.0181595:0.0164908



```R
data.frame(Order = 1:nrow(metadata(corrected)$merge.info),
           Left = stringi::stri_join_list(as.list(metadata(corrected)$merge.info$left), sep = ", "), 
           Right = stringi::stri_join_list(as.list(metadata(corrected)$merge.info$right)))
```


<table class="dataframe">
<caption>A data.frame: 2 × 3</caption>
<thead>
	<tr><th scope=col>Order</th><th scope=col>Left</th><th scope=col>Right</th></tr>
	<tr><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><td>1</td><td>LungCancer              </td><td>KidneyCancer</td></tr>
	<tr><td>2</td><td>LungCancer, KidneyCancer</td><td>Control1    </td></tr>
</tbody>
</table>




```R
metadata(corrected)$merge.info$lost.var
```


<table class="dataframe">
<caption>A matrix: 2 × 3 of type dbl</caption>
<thead>
	<tr><th scope=col>Control1</th><th scope=col>KidneyCancer</th><th scope=col>LungCancer</th></tr>
</thead>
<tbody>
	<tr><td>0.00000000</td><td>0.01165681</td><td>0.01178272</td></tr>
	<tr><td>0.01125956</td><td>0.01815954</td><td>0.01649082</td></tr>
</tbody>
</table>



## Save MNN results

<div class="alert alert-warning">
    <strong>Warning!</strong> The <code>reconstructed</code> matrix in the <code>assays</code> slot contains the corrected expression values for each gene in each cell, obtained by projecting the low-dimensional coordinates in corrected back into gene expression space. <strong>It is not recommended using this for anything other than visualisation.</strong> Use of the corrected values in any quantitative procedure should be treated with caution, and should be backed up by similar results from an analysis on the uncorrected values. See <a href="https://bioconductor.org/books/3.20/OSCA.multisample/using-corrected-values.html">OSCA reference</a>
</div>


```R
reducedDim(combined, "MNN") <- reducedDim(corrected, "corrected")
combined
```


    class: SingleCellExperiment 
    dim: 18129 49665 
    metadata(24): Control1_Samples Control1_cyclone ... LungCancer_DoubletDensity LungCancer_runInfo
    assays(2): counts logcounts
    rownames(18129): SAMD11 NOC2L ... MT-ND6 MT-CYB
    rowData names(5): ID Symbol Type SEQNAME is_mito
    colnames(49665): Control1_AAACAAGCAACAAGTTACTTTAGG-1 Control1_AAACAAGCAACTAGTGACTTTAGG-1 ...
      LungCancer_TTTGTGAGTTGAGTCTAGCTGTGA-1 LungCancer_TTTGTGAGTTTACGACAGCTGTGA-1
    colData names(14): Sample Barcode ... label condition
    reducedDimNames(4): PCA TSNE UMAP MNN
    mainExpName: Gene Expression
    altExpNames(1): Antibody Capture


## Visualising before & after MNN corrections with t-SNE


```R
set.seed(12345)
combined <- runTSNE(combined, dimred = "MNN", name = "MNN-TSNE", n_dimred = 50, 
                    n_threads = nthreads, BPPARAM = bpp)
```


```R
fig(width = 16, height = 8)
plotProjections(combined, "Sample", dimnames = c("TSNE", "MNN-TSNE"), feat_desc = "Sample", 
                feat_color = c_sample_col, guides_size = 4, point_size = 0.1, point_alpha = 0.1,
                legend_pos = "bottom", guides_nrow = 1, titles = c("TSNE, no correction","TSNE, MNN correction"))
reset.fig()
```


    
![png](Integrated_files/Integrated_71_0.png)
    



```R
fig(width = 16, height = 8)
plotProjections(combined, "CellType", dimnames = c("TSNE", "MNN-TSNE"), feat_desc = "Cell Type", 
                feat_color = c_celltype_col, guides_size = 4, point_size = 0.1, point_alpha = 0.1,
                legend_pos = "bottom", guides_nrow = 1, titles = c("TSNE, no correction","TSNE, MNN correction"))
reset.fig()
```


    
![png](Integrated_files/Integrated_72_0.png)
    


## Visualising before & after MNN corrections with UMAP


```R
set.seed(12345)
combined <- runUMAP(combined, dimred = "MNN", name = "MNN-UMAP", n_dimred = 50, 
                    n_neighbors = 30, spread = 1, min_dist = 0.3, n_threads = nthreads, BPPARAM = bpp)
```


```R
fig(width = 16, height = 8)
plotProjections(combined, "Sample", dimnames = c("UMAP", "MNN-UMAP"), feat_desc = "Sample", 
                feat_color = c_sample_col, guides_size = 4, point_size = 0.1, point_alpha = 0.1,
                legend_pos = "bottom", guides_nrow = 1, titles = c("UMAP, no correction","UMAP, MNN correction"))
reset.fig()
```


    
![png](Integrated_files/Integrated_75_0.png)
    



```R
fig(width = 16, height = 8)
plotProjections(combined, "CellType", dimnames = c("UMAP", "MNN-UMAP"), feat_desc = "Cell Type", 
                feat_color = c_celltype_col, guides_size = 4, point_size = 0.1, point_alpha = 0.1,
                legend_pos = "bottom", guides_nrow = 1, titles = c("UMAP, no correction","UMAP, MNN correction"))
reset.fig()
```


    
![png](Integrated_files/Integrated_76_0.png)
    


### Colour cells by... 

#### Library size `log10Sum`


```R
bk <- seq(min(combined$log10Sum), max(combined$log10Sum), max(combined$log10Sum)/20)
bk <- round(bk, 2)

fig(width = 16, height = 8)
plotProjections(combined, "log10Sum", dimnames = c("MNN-TSNE","MNN-UMAP"), feat_desc = "log10(Sum)", 
                feat_color = rev(rainbow(5)), color_breaks = bk, point_size = 0.1, point_alpha = 0.1, 
                guides_barheight = 20)
reset.fig()
```


    
![png](Integrated_files/Integrated_78_0.png)
    


#### Doublet score `DoubletDensity`


```R
combined$DoubletDensity <- unlist(lapply(all.common.normed, function(x) metadata(x)[['DoubletDensity']]))

fig(width = 16, height = 8)
plotProjections(combined, I(log1p(combined$DoubletDensity)), c("MNN-TSNE","MNN-UMAP"), 
                feat_desc = "Doublet Score (log1p)", feat_color = c_heatmap_col1, 
                point_size = 0.1, point_alpha = 0.1)
reset.fig()
```


    
![png](Integrated_files/Integrated_80_0.png)
    


# 5. Clustering of cells

Here, we recommend using any or all of graph-based methods, for example **Walktrap**, **Louvain** and **Leiden** algorithms, for cell clustering.

Also, to consider a <strong>2-pass clustering procedure</strong>. A one-time clustering on the full dataset might not yield the desired resolution or subclusters. For a two-step process, first, the cells are separated into smaller subsets based on a initial round of coarse-grain clustering. Then, a second round of clustering is performed on each subset. Finally the clustering results are combined and stored into the full object <code>combined</code>.

Below shows examples of uring `clusterRows()` function to perform the 3 clustering methods. Edit the workflow to perform clustering as one sees fit.

```
set.seed(12345)

# walktrap
clusterRows(mat, full = TRUE, 
            NNGraphParam(cluster.fun = "walktrap", k = k, type = "rank", # default in scran
                         BNPARAM = AnnoyParam(), num.threads = nthreads))

# louvain
clusterRows(mat, full = TRUE, 
            NNGraphParam(cluster.fun = "louvain", k = k, type = "jaccard", # default in Seurat
                         cluster.args = list(resolution = 1),
                         BNPARAM = AnnoyParam(), num.threads = nthreads))

# leiden
clusterRows(mat, full = TRUE, 
            NNGraphParam(cluster.fun = "leiden", k = k, 
                         cluster.args = list(resolution_parameter = 1),
                         BNPARAM = AnnoyParam(), num.threads = nthreads))
```


```R
# Stores clustering labels
my.clusters <- vector(mode = "list")
communities <- vector(mode = "list")
```

## First-pass clustering

Here, we use the faster **Leiden** algorithm to perform a quick coarse clustering.


```R
method <- "leiden"
dimname <- "MNN" # or "PCA" without batch correction
n_dimred <- 20 # number of dimensions to use; default is 50
k <- 15

mat <- reducedDim(combined, dimname)[, seq_len(n_dimred), drop = FALSE]

set.seed(12345)
communities[[method]][[dimname]] <- clusterRows(mat, full = TRUE,
                                                NNGraphParam(cluster.fun = method, k = k, 
                                                             cluster.args = list(resolution_parameter = 0.01),
                                                             BNPARAM = AnnoyParam(), num.threads = nthreads))
my.clusters[[method]][[dimname]] <- factor(communities[[method]][[dimname]]$clusters)
```


```R
print(sprintf("%s cluster assignments from %s:", stringr::str_to_title(method), dimname))
table(colData(combined)$Sample, my.clusters[[method]][[dimname]])
plotSilhouette(mat, my.clusters[[method]][[dimname]], printDiff = FALSE, plot = FALSE)

combined$first.pass <- my.clusters[[method]][[dimname]]
```

    [1] "Leiden cluster assignments from MNN:"



                  
                      1    2    3    4    5    6    7    8    9   10   11   12   13
      Control1      330 5458 2780 1611  496 1421 3430  467 2014  253   86  111    3
      KidneyCancer  209 1961 7380  655  155 2839  111 2408  279   17   23  594  537
      LungCancer    275  375 3474 2053   69 3070   79 2132 1281   12    8  465  744


    Silhouette width summary:


       Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    -0.4945  0.1119  0.2256  0.2267  0.3688  0.6199 


### (Optional) Remove randomly scattered cells that high `DoubletDensity` scores

If there are many cells with high doublet scores randomly scattered on the dimensionality reduction plots (i.e. TSNE or UMAP), we can choose to do a round of doublet cell cleaning now. *The aim is to NOT remove doublet cluster(s).*

<div style="width: 100%; height: 25px; border-bottom: 1px dashed black; text-align: center">
    <span style="font-size: 20px; background-color: yellow; padding: 0 10px;">Begin Optional Step (remove doublet cells)</span>
</div>


```R
fig(width = 16, height = 8)
plotProjections(combined, I(log1p(combined$DoubletDensity)), c("MNN-TSNE","MNN-UMAP"), 
                feat_desc = "Doublet Score (log1p)", feat_color = c_heatmap_col1, text_by = "first.pass", 
                point_size = 0.1, point_alpha = 0.1)
reset.fig()
```


    
![png](Integrated_files/Integrated_88_0.png)
    


#### Define per-sample per-cluster `DoubletDensity` score upperbound thresholds


```R
# Upper limit/outlier: 75% quantile + (IQR * 3)
upper <- as.data.frame(colData(combined)[,c("Sample","first.pass","DoubletDensity")]) %>% 
    group_by(Sample, first.pass) %>% 
    reframe(iqr = iqr(log1p(DoubletDensity)), 
            enframe(quantile(log1p(DoubletDensity), c(0.25, 0.5, 0.75)), "quantile", "Score")) %>% 
    filter(quantile == "75%") %>% mutate(upper = iqr * 3 + Score) %>% arrange(Sample, first.pass)
head(upper)
```


<table class="dataframe">
<caption>A tibble: 6 × 6</caption>
<thead>
	<tr><th scope=col>Sample</th><th scope=col>first.pass</th><th scope=col>iqr</th><th scope=col>quantile</th><th scope=col>Score</th><th scope=col>upper</th></tr>
	<tr><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td>Control1</td><td>1</td><td>0.3219096</td><td>75%</td><td>2.3887628</td><td>3.3544916</td></tr>
	<tr><td>Control1</td><td>2</td><td>0.2972515</td><td>75%</td><td>0.3364722</td><td>1.2282268</td></tr>
	<tr><td>Control1</td><td>3</td><td>0.5072478</td><td>75%</td><td>0.6205765</td><td>2.1423199</td></tr>
	<tr><td>Control1</td><td>4</td><td>0.1133287</td><td>75%</td><td>0.1133287</td><td>0.4533147</td></tr>
	<tr><td>Control1</td><td>5</td><td>0.3560794</td><td>75%</td><td>2.5790796</td><td>3.6473179</td></tr>
	<tr><td>Control1</td><td>6</td><td>0.2231436</td><td>75%</td><td>0.2623643</td><td>0.9317949</td></tr>
</tbody>
</table>




```R
doublets <- as.data.frame(colData(combined)[,c("Sample","first.pass","DoubletDensity")]) %>%
    rownames_to_column("Barcode") %>% left_join(upper[,c("Sample","first.pass","upper")]) %>%
    filter(log1p(DoubletDensity) > upper)
head(doublets)

table("Is doublet" = colnames(combined) %in% doublets$Barcode)
```

    [1m[22mJoining with `by = join_by(Sample, first.pass)`



<table class="dataframe">
<caption>A data.frame: 6 × 5</caption>
<thead>
	<tr><th></th><th scope=col>Barcode</th><th scope=col>Sample</th><th scope=col>first.pass</th><th scope=col>DoubletDensity</th><th scope=col>upper</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>1</th><td>Control1_AAAGGCTTCCAACATTACTTTAGG-1</td><td>Control1</td><td>2 </td><td>2.64</td><td>1.2282268</td></tr>
	<tr><th scope=row>2</th><td>Control1_AAATCCTTCATCAGCCACTTTAGG-1</td><td>Control1</td><td>12</td><td>0.72</td><td>0.3812407</td></tr>
	<tr><th scope=row>3</th><td>Control1_AACCAGGTCATGGAACACTTTAGG-1</td><td>Control1</td><td>2 </td><td>4.48</td><td>1.2282268</td></tr>
	<tr><th scope=row>4</th><td>Control1_AACCATAAGAGGCGGAACTTTAGG-1</td><td>Control1</td><td>3 </td><td>9.76</td><td>2.1423199</td></tr>
	<tr><th scope=row>5</th><td>Control1_AACCATAAGCTTATCGACTTTAGG-1</td><td>Control1</td><td>9 </td><td>2.72</td><td>0.8067847</td></tr>
	<tr><th scope=row>6</th><td>Control1_AACCATTTCGCGGATAACTTTAGG-1</td><td>Control1</td><td>9 </td><td>2.44</td><td>0.8067847</td></tr>
</tbody>
</table>




    Is doublet
    FALSE  TRUE 
    48494  1171 



```R
fig(width = 16, height = 8)
plotProjections(combined, colnames(combined) %in% doublets$Barcode, c("MNN-TSNE","MNN-UMAP"), 
                feat_desc = "Doublets to be removed", feat_color = c("yellow","red"), 
                point_size = 0.1, point_alpha = 0.1, guides_size = 4, rel_widths = c(9, 1))
reset.fig()
```


    
![png](Integrated_files/Integrated_92_0.png)
    


#### Remove marked doublet cells and re-do `runTSNE`and `runUMAP`


```R
combined <- combined[, !colnames(combined) %in% doublets$Barcode]
combined
```


    class: SingleCellExperiment 
    dim: 18129 48494 
    metadata(24): Control1_Samples Control1_cyclone ... LungCancer_DoubletDensity LungCancer_runInfo
    assays(2): counts logcounts
    rownames(18129): SAMD11 NOC2L ... MT-ND6 MT-CYB
    rowData names(5): ID Symbol Type SEQNAME is_mito
    colnames(48494): Control1_AAACAAGCAACAAGTTACTTTAGG-1 Control1_AAACAAGCAACTAGTGACTTTAGG-1 ...
      LungCancer_TTTGTGAGTTGAGTCTAGCTGTGA-1 LungCancer_TTTGTGAGTTTACGACAGCTGTGA-1
    colData names(16): Sample Barcode ... DoubletDensity first.pass
    reducedDimNames(6): PCA TSNE ... MNN-TSNE MNN-UMAP
    mainExpName: Gene Expression
    altExpNames(1): Antibody Capture



```R
set.seed(12345)
combined <- runTSNE(combined, dimred = "MNN", name = "MNN-TSNE", n_dimred = 50, 
                    n_threads = nthreads, BPPARAM = bpp)
```


```R
set.seed(12345)
combined <- runUMAP(combined, dimred = "MNN", name = "MNN-UMAP", n_dimred = 50, 
                    n_neighbors = 30, spread = 1, min_dist = 0.3, n_threads = nthreads, BPPARAM = bpp)
```

#### Re-run first-pass clustering using Leiden algorithm


```R
method <- "leiden"
dimname <- "MNN" # or "PCA" without batch correction
n_dimred <- 50 # number of dimensions to use; default is 50
k <- 10

mat <- reducedDim(combined, dimname)[, seq_len(n_dimred), drop = FALSE]

set.seed(12345)
communities[[method]][[dimname]] <- clusterRows(mat, full = TRUE,
                                                NNGraphParam(cluster.fun = method, k = k, 
                                                             cluster.args = list(resolution_parameter = 0.01),
                                                             BNPARAM = AnnoyParam(), num.threads = nthreads))
my.clusters[[method]][[dimname]] <- factor(communities[[method]][[dimname]]$clusters)
```


```R
print(sprintf("%s cluster assignments from %s:", stringr::str_to_title(method), dimname))
table(colData(combined)$Sample, my.clusters[[method]][[dimname]])
plotSilhouette(mat, my.clusters[[method]][[dimname]], printDiff = FALSE, plot = FALSE)

combined$first.pass <- my.clusters[[method]][[dimname]]
```

    [1] "Leiden cluster assignments from MNN:"



                  
                      1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18   19
      Control1      329 2164 2837 1007 3145  491 1348 3102  556  115 1997   28   54   79  285  307   95   95    4
      KidneyCancer  201  956  659  175  964  155 2679   81  457 1216  280   16 3738   23   24  497 3579  547  529
      LungCancer    269   59  670 1739  325   71 2591   67  261 1093 1254    8 1615    8    7   87 2363  458  735


    Silhouette width summary:


       Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    -0.4876  0.0592  0.1276  0.1317  0.1931  0.5620 


<div style="width: 100%; height: 25px; border-bottom: 1px dashed black; text-align: center">
    <span style="font-size: 20px; background-color: yellow; padding: 0 10px;">End Optional Step (remove doublet cells)</span>
</div>


```R
p1 <- plotProjections(combined, "CellType", dimnames = c("MNN-TSNE", "MNN-UMAP"), feat_desc = "Cell Type", 
                      feat_color = c_celltype_col, text_by = "first.pass", point_size = 0.1, point_alpha = 0.1, 
                      legend_pos = "none", add_void = TRUE)

p2 <- plotProjections(combined, "first.pass", dimnames = c("MNN-TSNE", "MNN-UMAP"), feat_desc = method, 
                      feat_color = c30(), text_by = "first.pass", point_size = 0.1, point_alpha = 0.1,
                      guides_size = 4)

fig(width = 16, height = 15)    
plot_grid(p1, p2, ncol = 1)
reset.fig()
```


    
![png](Integrated_files/Integrated_101_0.png)
    


## Create subsets

You may divide your cells into subsets containing fewer clusters.

In this example, we create 2 subsets:

<div class="alert alert-info">
    <ul>
        <li>Subset 1: 2, 7, 10, 16, 18</li>
        <li>Subset 2: 3, 5, 12, 13, 17</li>
        <li>Subset 3: 1, 4, 6, 8, 9, 11, 14, 15, 19</li>
    </ul>
</div>

### Create TSNE & UMAP (Subset 1)


```R
subset1 <- combined[, combined$first.pass %in% c(2, 7, 10, 16, 18)]
colData(subset1) <- droplevels(colData(subset1))
dim(subset1)

set.seed(12345)
subset1 <- runTSNE(subset1, dimred = "MNN", name = "MNN-TSNE", n_dimred = 50, n_threads = nthreads, BPPARAM = bpp)

set.seed(12345)
subset1 <- runUMAP(subset1, dimred = "MNN", name = "MNN-UMAP", n_dimred = 50, 
                   n_neighbors = 30, spread = 1, min_dist = 0.3, n_threads = nthreads, BPPARAM = bpp)

p1 <- plotProjections(subset1, "first.pass", dimnames = c("MNN-TSNE", "MNN-UMAP"), feat_desc = "Cluster", 
                      feat_color = c30(), text_by = "first.pass", point_size = 0.1, point_alpha = 0.1, 
                      guides_size = 4, rel_widths = c(8,1))

p2 <- plotProjections(subset1, "CellCycle", dimnames = c("MNN-TSNE", "MNN-UMAP"), feat_desc = "Cell Cycle", 
                      feat_color = c_phase_col, text_by = "first.pass", point_size = 0.1, point_alpha = 0.1,
                      guides_size = 4, rel_widths = c(8,1))

fig(width = 16, height = 15)
plot_grid(p1, p2, ncol = 1)
reset.fig()
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>18129</li><li>14212</li></ol>




    
![png](Integrated_files/Integrated_104_1.png)
    


### Create TSNE & UMAP (Subset 2)


```R
subset2 <- combined[, combined$first.pass %in% c(3, 5, 12, 13, 17)]
colData(subset2) <- droplevels(colData(subset2))
dim(subset2)

set.seed(12345)
subset2 <- runTSNE(subset2, dimred = "MNN", name = "MNN-TSNE", n_dimred = 50, n_threads = nthreads, BPPARAM = bpp)

set.seed(12345)
subset2 <- runUMAP(subset2, dimred = "MNN", name = "MNN-UMAP", n_dimred = 50, 
                   n_neighbors = 30, spread = 1, min_dist = 0.3, n_threads = nthreads, BPPARAM = bpp)

p1 <- plotProjections(subset2, "first.pass", dimnames = c("MNN-TSNE", "MNN-UMAP"), feat_desc = "Cluster", 
                      feat_color = c30(), text_by = "first.pass", point_size = 0.1, point_alpha = 0.1, 
                      guides_size = 4, rel_widths = c(8,1))

p2 <- plotProjections(subset2, "CellCycle", dimnames = c("MNN-TSNE", "MNN-UMAP"), feat_desc = "Cell Cycle", 
                      feat_color = c_phase_col, text_by = "first.pass", point_size = 0.1, point_alpha = 0.1, 
                      guides_size = 4, rel_widths = c(8,1))

fig(width = 16, height = 15)
plot_grid(p1, p2, ncol = 1)
reset.fig()
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>18129</li><li>20096</li></ol>




    
![png](Integrated_files/Integrated_106_1.png)
    


### Create TSNE & UMAP (Subset 3)


```R
subset3 <- combined[, combined$first.pass %in% c(1, 4, 6, 8, 9, 11, 14, 15, 19)]
colData(subset3) <- droplevels(colData(subset3))
dim(subset3)

set.seed(12345)
subset3 <- runTSNE(subset3, dimred = "MNN", name = "MNN-TSNE", n_dimred = 50, n_threads = nthreads, BPPARAM = bpp)

set.seed(12345)
subset3 <- runUMAP(subset3, dimred = "MNN", name = "MNN-UMAP", n_dimred = 50, 
                   n_neighbors = 30, spread = 1, min_dist = 0.3, n_threads = nthreads, BPPARAM = bpp)

p1 <- plotProjections(subset3, "first.pass", dimnames = c("MNN-TSNE", "MNN-UMAP"), feat_desc = "Cluster", 
                      feat_color = c30(), text_by = "first.pass", point_size = 0.1, point_alpha = 0.1, 
                      guides_size = 4, rel_widths = c(8,1))

p2 <- plotProjections(subset3, "CellCycle", dimnames = c("MNN-TSNE", "MNN-UMAP"), feat_desc = "Cell Cycle", 
                      feat_color = c_phase_col, text_by = "first.pass", point_size = 0.1, point_alpha = 0.1,
                      guides_size = 4, rel_widths = c(8,1))

fig(width = 16, height = 15)
plot_grid(p1, p2, ncol = 1)
reset.fig()
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>18129</li><li>14186</li></ol>




    
![png](Integrated_files/Integrated_108_1.png)
    


## Second-pass clustering

Here, we use the **Louvain** algorithm to perform clustering, but feel free to use different clustering algorithms that suit your dataset

**Subset 1**


```R
method <- "louvain"
dimname <- "MNN" # or "PCA" without batch correction
n_dimred <- 50 # number of dimensions to use; default is 50
k <- 15

mat <- reducedDim(subset1, dimname)[, seq_len(n_dimred), drop = FALSE]

set.seed(12345)
communities[[method]][[dimname]] <- clusterRows(mat, full = TRUE,
                                                NNGraphParam(cluster.fun = method, k = k, type = "jaccard", 
                                                             cluster.args = list(resolution = 2.5), 
                                                             BNPARAM = AnnoyParam(), num.threads = nthreads))
my.clusters[[method]][[dimname]] <- factor(communities[[method]][[dimname]]$clusters)
```


```R
print(sprintf("%s cluster assignments from %s:", stringr::str_to_title(method), dimname))
table(colData(subset1)$Sample, my.clusters[[method]][[dimname]])
plotSilhouette(mat, my.clusters[[method]][[dimname]], printDiff = FALSE, plot = FALSE)

subset1$merged.louvain <- my.clusters[[method]][[dimname]]
```

    [1] "Louvain cluster assignments from MNN:"



                  
                      1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18   19
      Control1     1342  817  774   87  454   10  204   95  114   13   25    4   68    2   13    3    1    1    1
      KidneyCancer   11  368   16  387   10  189    5  179  499 1287   11  157  405   21  729   11  705   39  187
      LungCancer      6   13   84   94   22  147    6  153   90  310  433  109  323  668   58  109  424  255  188
                  
                     20   21   22   23
      Control1        1    0    0    0
      KidneyCancer  615   61    3    0
      LungCancer     43  383  245  125


    Silhouette width summary:


        Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
    -0.36644  0.01334  0.07794  0.07275  0.13726  0.36216 


<div class="alert alert-warning">
  <strong>Note!</strong> After partially running the workflow, usually after preliminary manual curation of per-cluster cell types, you may want/need to come back to this stage to fine-tune the clusters.
</div>

As an example of changing pre-defined clusters, here we merge cluster 1 and 2, which are both **naive CD8 T cells**.


```R
levels(subset1$merged.louvain)[2] <- 1
levels(subset1$merged.louvain) <- seq(1:nlevels(subset1$merged.louvain))
table(subset1$condition, subset1$merged.louvain)
```


             
                 1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18   19   20
      Control 2159  774   87  454   10  204   95  114   13   25    4   68    2   13    3    1    1    1    1    0
      Cancer   398  100  481   32  336   11  332  589 1597  444  266  728  689  787  120 1129  294  375  658  444
             
                21   22
      Control    0    0
      Cancer   248  125



```R
p1 <- plotProjections(subset1, "merged.louvain", dimnames = c("MNN-TSNE", "MNN-UMAP"), feat_desc = "Cluster", 
                      feat_color = c30(), text_by = "merged.louvain", point_size = 0.3, point_alpha = 0.1,
                      guides_size = 4, rel_widths = c(8,1))

p2 <- plotProjections(subset1, "CellCycle", dimnames = c("MNN-TSNE", "MNN-UMAP"), feat_desc = "Cell Cycle", 
                      feat_color = c_phase_col, text_by = "merged.louvain", point_size = 0.3, point_alpha = 0.1,
                      guides_size = 4, rel_widths = c(8,1))

fig(width = 16, height = 15)
plot_grid(p1, p2, ncol = 1)
reset.fig()
```


    
![png](Integrated_files/Integrated_114_0.png)
    


**Subset 2**


```R
method <- "louvain"
dimname <- "MNN" # or "PCA" without batch correction
n_dimred <- 50 # number of dimensions to use; default is 50
k <- 10

mat <- reducedDim(subset2, dimname)[, seq_len(n_dimred), drop = FALSE]

set.seed(12345)
communities[[method]][[dimname]] <- clusterRows(mat, full = TRUE,
                                                NNGraphParam(cluster.fun = method, k = k, type = "jaccard", 
                                                             cluster.args = list(resolution = 0.6), 
                                                             BNPARAM = AnnoyParam(), num.threads = nthreads))
my.clusters[[method]][[dimname]] <- factor(communities[[method]][[dimname]]$clusters)
```


```R
print(sprintf("%s cluster assignments from %s:", stringr::str_to_title(method), dimname))
table(colData(subset2)$Sample, my.clusters[[method]][[dimname]])
plotSilhouette(mat, my.clusters[[method]][[dimname]], printDiff = FALSE, plot = FALSE)

subset2$merged.louvain <- my.clusters[[method]][[dimname]]
```

    [1] "Louvain cluster assignments from MNN:"



                  
                      1    2    3    4    5    6    7    8    9
      Control1     2161 3139  288  323   14  163   28   24   19
      KidneyCancer  292  942  246   28  839 2477   16 2978 1138
      LungCancer    217  307  315   21 1288 1683    8  337  805


    Silhouette width summary:


        Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
    -0.25800  0.03423  0.11327  0.12293  0.20316  0.41417 



```R
p1 <- plotProjections(subset2, "merged.louvain", dimnames = c("MNN-TSNE", "MNN-UMAP"), feat_desc = "Cluster", 
                      feat_color = c30(), text_by = "merged.louvain", point_size = 0.1, point_alpha = 0.1,
                      guides_size = 4, rel_widths = c(8,1))

p2 <- plotProjections(subset2, "CellCycle", dimnames = c("MNN-TSNE", "MNN-UMAP"), feat_desc = "Cell Cycle", 
                      feat_color = c_phase_col, text_by = "merged.louvain", point_size = 0.1, point_alpha = 0.1, 
                      guides_size = 4, rel_widths = c(8,1))

fig(width = 16, height = 15)
plot_grid(p1, p2, ncol = 1)
reset.fig()
```


    
![png](Integrated_files/Integrated_118_0.png)
    


**Subset 3**


```R
method <- "louvain"
dimname <- "MNN" # or "PCA" without batch correction
n_dimred <- 50 # number of dimensions to use; default is 50
k <- 10

mat <- reducedDim(subset3, dimname)[, seq_len(n_dimred), drop = FALSE]

set.seed(12345)
communities[[method]][[dimname]] <- clusterRows(mat, full = TRUE,
                                                NNGraphParam(cluster.fun = method, k = k, type = "jaccard", 
                                                             cluster.args = list(resolution = 1.1), 
                                                             BNPARAM = AnnoyParam(), num.threads = nthreads))
my.clusters[[method]][[dimname]] <- factor(communities[[method]][[dimname]]$clusters)
```


```R
print(sprintf("%s cluster assignments from %s:", stringr::str_to_title(method), dimname))
table(colData(subset3)$Sample, my.clusters[[method]][[dimname]])
plotSilhouette(mat, my.clusters[[method]][[dimname]], printDiff = FALSE, plot = FALSE)

subset3$merged.louvain <- my.clusters[[method]][[dimname]]
```

    [1] "Louvain cluster assignments from MNN:"



                  
                      1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18
      Control1      328  993  407 1174  317 1201 1541  729  141   80  181  286  248   82  134    3    5    0
      KidneyCancer  202    2  145    6  342   40   14   37    4   23   77   25  114   10  185  526  173    0
      LungCancer    272   67   57   22  186   25  104   22   11    8  195    7   72   14  791  733 1672  153


    Silhouette width summary:


       Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    -0.3929  0.1135  0.2045  0.1928  0.2816  0.5179 


<div class="alert alert-warning">
  <strong>Note!</strong> After partially running the workflow, usually after preliminary manual curation of per-cluster cell types, you may want/need to come back to this stage to fine-tune the clusters.
</div>

As an example of changing pre-defined clusters, here we merge clusters of **NK cells**.


```R
levels(subset3$merged.louvain)[c(9, 11, 15)] <- 7
levels(subset3$merged.louvain) <- seq(1:nlevels(subset3$merged.louvain))
table(subset3$condition, subset3$merged.louvain)
```


             
                 1    2    3    4    5    6    7    8    9   10   11   12   13   14   15
      Control  328  993  407 1174  317 1201 1997  729   80  286  248   82    3    5    0
      Cancer   474   69  202   28  528   65 1381   59   31   32  186   24 1259 1845  153



```R
p1 <- plotProjections(subset3, "merged.louvain", dimnames = c("MNN-TSNE", "MNN-UMAP"), feat_desc = "Cluster", 
                      feat_color = c30(), text_by = "merged.louvain", point_size = 0.1, point_alpha = 0.1,
                      guides_size = 4, rel_widths = c(8,1))

p2 <- plotProjections(subset3, "CellCycle", dimnames = c("MNN-TSNE", "MNN-UMAP"), feat_desc = "Cell Cycle", 
                      feat_color = c_phase_col, text_by = "merged.louvain", point_size = 0.1, point_alpha = 0.1, 
                      guides_size = 4, rel_widths = c(8,1))

fig(width = 16, height = 15)
plot_grid(p1, p2, ncol = 1)
reset.fig()
```


    
![png](Integrated_files/Integrated_124_0.png)
    


## Add assigned clusters to `sce` objects

Choose the desired clustering result and assign to the default "label" column using the `colLabels` function.


```R
# Combine cluster labels from all subsets
# Assign temporary IDs for now
combined$merged.louvain <- rbind(colData(subset1)[, c("Barcode","merged.louvain")] %>% as.data.frame %>% 
                                 mutate(merged.louvain = paste0("a", merged.louvain)), 
                                 colData(subset2)[, c("Barcode","merged.louvain")] %>% as.data.frame %>% 
                                 mutate(merged.louvain = paste0("b", merged.louvain)), 
                                 colData(subset3)[, c("Barcode","merged.louvain")] %>% as.data.frame %>% 
                                 mutate(merged.louvain = paste0("c", merged.louvain)))[colnames(combined),"merged.louvain"]
combined$merged.louvain <- factor(combined$merged.louvain, 
                                  levels = gtools::mixedsort(unique(combined$merged.louvain)))
colLabels(combined) <- combined$merged.louvain
```

### Rename and reorder clusters

<div class="alert alert-warning">
    <strong>Skip this step in the initial run!</strong>
    <br />Come back and run this step after cell clusters of the same "coarse cell type" and/or lineage are identified preliminarily. We do this so that the clusters of the same "<i>type</i>" are in sequential order by their cluster number.    
</div>


```R
levels(colLabels(combined)) <- c(1, 23, 13, 24, 21, 43, 41, 42, 10, 19, 17, 40, 12, 22, 11, 14, 18, 15, 2, 
                                 26, 20, 16, 8, 3, 9, 25, 7, 4, 29, 6, 5, 44, 30, 45, 34, 33, 35, 27, 36, 
                                 39, 38, 32, 46, 37, 31, 28)
colLabels(combined) <- factor(colLabels(combined), levels = gtools::mixedsort(levels(colLabels(combined))))
```


```R
# Set colours for cell clusters
c_clust_col <- choosePalette(colLabels(combined), 
                             pals::kovesi.rainbow_bgyrm_35_85_c69(nlevels(colLabels(combined)))) # more than 40 clusters
```

### Print number and percentage of cells


```R
table_samples_by_clusters <- table(Sample = combined$Sample, Cluster = combined$label)
table_samples_by_clusters
```


                  Cluster
    Sample            1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18   19
      Control1     2159    1 3139  163   19   24   14 2161  288   13    3    2   87    1    1    0    4    1   25
      KidneyCancer  379  615  942 2477 1138 2978  839  292  246 1287   11   21  387  705  187    0  157   39   11
      LungCancer     19   43  307 1683  805  337 1288  217  315  310  109  668   94  424  188  125  109  255  433
                  Cluster
    Sample           20   21   22   23   24   25   26   27   28   29   30   31   32   33   34   35   36   37   38
      Control1        0   10   13  774  454  323    0 1997    0   28  993    5  248  317 1174 1201  729    3  286
      KidneyCancer    3  189  729   16   10   28   61  280    0   16    2  173  114  342    6   40   37  526   25
      LungCancer    245  147   58   84   22   21  383 1101  153    8   67 1672   72  186   22   25   22  733    7
                  Cluster
    Sample           39   40   41   42   43   44   45   46
      Control1       80   68   95  114  204  328  407   82
      KidneyCancer   23  405  179  499    5  202  145   10
      LungCancer      8  323  153   90    6  272   57   14



```R
# No. of cells 
fig(width = 16, height = 5)
ggplot(data.frame(table_samples_by_clusters), aes(Sample, Freq, fill = Cluster)) + 
    geom_bar(position = "stack", stat = "identity", linewidth = 0.2, color = "black") + coord_flip() +
    scale_y_continuous("Number of cells", labels = comma) + guides(fill = guide_legend(ncol = 4)) +
    scale_fill_manual(values = c_clust_col) + theme_cowplot(20) + theme(axis.title.y = element_blank())
reset.fig()
```


    
![png](Integrated_files/Integrated_132_0.png)
    



```R
# Percenatge cells
fig(width = 16, height = 5)
ggplot(data.frame(table_samples_by_clusters), aes(Sample, Freq, fill = Cluster)) +
    geom_bar(position = "fill", stat = "identity", linewidth = 0.2, color = "black") + coord_flip() +
    scale_y_continuous("Percentage", labels = percent_format()) + guides(fill = guide_legend(ncol = 4)) +
    scale_fill_manual(values = c_clust_col) + theme_cowplot(20) + theme(axis.title.y = element_blank())
reset.fig()
```


    
![png](Integrated_files/Integrated_133_0.png)
    


## t-SNE and UMAP plots

<div class="alert alert-warning">
    <strong>Warning!</strong> change the <code>dimnames</code> names to view plot in the corresponding reduced dimension results.
</div>


```R
# Coloured by clusters
fig(width = 16, height = 9)
plotProjections(combined, "label", dimnames = c("MNN-TSNE", "MNN-UMAP"), feat_desc = "Cluster", 
                feat_color = c_clust_col, text_by = "label", point_size = 0.1, point_alpha = 0.1,
                text_size = 6, guides_nrow = 3, guides_size = 4, legend_pos = "bottom", rel_heights = c(8,1))
reset.fig()
```


    
![png](Integrated_files/Integrated_135_0.png)
    



```R
bk <- seq(min(combined$log10Sum), max(combined$log10Sum), max(combined$log10Sum)/20)
bk <- round(bk, 2)

# Coloured by log10Sum
fig(width = 16, height = 8)
plotProjections(combined, "log10Sum", dimnames = c("MNN-TSNE", "MNN-UMAP"), feat_desc = "log10(Sum)", 
                feat_color = rev(rainbow(5)), color_breaks = bk, point_size = 0.1, point_alpha = 0.1,
                text_by = "label", text_size = 6, guides_barheight = 20)
reset.fig()
```


    
![png](Integrated_files/Integrated_136_0.png)
    



```R
# Coloured by cell cycle phases
fig(width = 16, height = 9)
plotProjections(combined, "CellCycle", dimnames = c("MNN-TSNE", "MNN-UMAP"), feat_desc = "Cell Cycle Phases", 
                feat_color = c_phase_col, text_by = "label", point_size = 0.1, point_alpha = 0.1,
                text_size = 6, guides_nrow = 1, guides_size = 4, legend_pos = "bottom")
reset.fig()
```


    
![png](Integrated_files/Integrated_137_0.png)
    



```R
fig(width = 16, height = 5)
colData(combined) %>% as.data.frame %>% group_by(label, CellCycle) %>% summarise(counts = n()) %>%
    ggplot(aes(label, counts, fill = CellCycle)) + geom_col(position = "fill") + 
    scale_fill_manual(values = c_phase_col) + theme_cowplot(14) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
reset.fig()
```

    [1m[22m`summarise()` has grouped output by 'label'. You can override using the `.groups` argument.



    
![png](Integrated_files/Integrated_138_1.png)
    



```R
# Coloured by samples
fig(width = 16, height = 9)
plotProjections(combined, "Sample", dimnames = c("MNN-TSNE", "MNN-UMAP"), feat_desc = "Sample", 
                feat_color = c_sample_col, text_by = "label", point_size = 0.1, point_alpha = 0.1,
                text_size = 6, guides_nrow = 1, guides_size = 4, legend_pos = "bottom")
reset.fig()
```


    
![png](Integrated_files/Integrated_139_0.png)
    



```R
fig(width = 16, height = 5)
colData(combined) %>% as.data.frame %>% group_by(label, Sample) %>% summarise(counts = n()) %>%
    ggplot(aes(label, counts, fill = Sample)) + geom_col(position = "fill") + 
    scale_fill_manual(values = c_sample_col) + theme_cowplot(14) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
reset.fig()
```

    [1m[22m`summarise()` has grouped output by 'label'. You can override using the `.groups` argument.



    
![png](Integrated_files/Integrated_140_1.png)
    



```R
# Coloured by condition
fig(width = 16, height = 9)
plotProjections(combined, "condition", dimnames = c("MNN-TSNE", "MNN-UMAP"), feat_desc = "Stage", 
                feat_color = c_cond_col, text_by = "label", point_size = 0.1, point_alpha = 0.1,
                text_size = 6, guides_nrow = 1, guides_size = 4, legend_pos = "bottom")
reset.fig()
```


    
![png](Integrated_files/Integrated_141_0.png)
    



```R
fig(width = 16, height = 5)
colData(combined) %>% as.data.frame %>% group_by(label, condition) %>% summarise(counts = n()) %>%
    ggplot(aes(label, counts, fill = condition)) + geom_col(position = "fill") + 
    scale_fill_manual(values = c_cond_col) + theme_cowplot(14) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
reset.fig()
```

    [1m[22m`summarise()` has grouped output by 'label'. You can override using the `.groups` argument.



    
![png](Integrated_files/Integrated_142_1.png)
    



```R
# Coloured by cell types
fig(width = 16, height = 9)
plotProjections(combined, "CellType", dimnames = c("MNN-TSNE", "MNN-UMAP"), feat_desc = "Cell Type", 
                feat_color = c_celltype_col, text_by = "label", point_size = 0.1, point_alpha = 0.1,
                text_size = 6, guides_nrow = 1, guides_size = 4, legend_pos = "bottom")
reset.fig()
```


    
![png](Integrated_files/Integrated_143_0.png)
    



```R
table(Cluster = combined$label, CellType = combined$CellType)
```


           CellType
    Cluster B-cells CD4+ T-cells CD8+ T-cells   DC  HSC Monocytes Neutrophils NK cells Platelets T-cells
         1        0          908         1649    0    0         0           0        0         0       0
         2        0          488          171    0    0         0           0        0         0       0
         3        0         4349           39    0    0         0           0        0         0       0
         4        0         4138          185    0    0         0           0        0         0       0
         5        0         1427          534    0    0         0           0        0         0       1
         6        0         3058          281    0    0         0           0        0         0       0
         7        0          722         1418    0    0         0           0        0         0       1
         8        0         2069          599    0    0         0           0        0         0       2
         9        0          779           70    0    0         0           0        0         0       0
         10       0            0         1391    0    0         0           0      219         0       0
         11       0            0           78    0    0         0           0       45         0       0
         12       0            0          543    0    0         0           0      148         0       0
         13       0          106          462    0    0         0           0        0         0       0
         14       0           18         1112    0    0         0           0        0         0       0
         15       0            0          224    0    0         0           0      152         0       0
         16       0            0          105    0    0         0           0       20         0       0
         17       0            0          166    0    0         0           0      104         0       0
         18       0            0          245    0    0         0           0       50         0       0
         19       0            0            1    0    0         0           0      468         0       0
         20       0            0            7    0    0         0           0      241         0       0
         21       0            0           21    0    0         0           0      325         0       0
         22       0            0           34    0    0         0           0      766         0       0
         23       0            0          647    0    0         0           0      227         0       0
         24       0            0          337    0    0         0           0      149         0       0
         25       0           59          311    0    0         0           0        2         0       0
         26       0            0          383    0    0         0           0       61         0       0
         27       0            0           13    0    0         0           0     3365         0       0
         28       0            0            0    0    0         0           0      153         0       0
         29       0            9           43    0    0         0           0        0         0       0
         30    1062            0            0    0    0         0           0        0         0       0
         31    1850            0            0    0    0         0           0        0         0       0
         32     434            0            0    0    0         0           0        0         0       0
         33     845            0            0    0    0         0           0        0         0       0
         34       0            0            0    0    0      1202           0        0         0       0
         35       0            0            0    0    0      1266           0        0         0       0
         36       0            0            0    0    0       788           0        0         0       0
         37       0            0            0    0    0      1262           0        0         0       0
         38       0            0            0    1    0       317           0        0         0       0
         39      23            5            9    0   33        40           0        1         0       0
         40       0          523          265    0    0         0           0        7         1       0
         41       0            1          203    0    0         0           0      223         0       0
         42       0           31          657    0    0         0           0       15         0       0
         43       0            0          103    0    0         0           0      112         0       0
         44     301          197          166    0    0         0           0      138         0       0
         45       0            2           11    0    0       576           1       19         0       0
         46       2            0            0    0    0       104           0        0         0       0



```R
# Coloured by DoubletDensity scores
fig(width = 16, height = 8)
plotProjections(combined, I(log1p(combined$DoubletDensity)), c("MNN-TSNE","MNN-UMAP"), 
                feat_desc = "Doublet Score (log1p)", feat_color = c_heatmap_col1, 
                text_by = "label", text_size = 6, point_size = 0.1, point_alpha = 0.1)
reset.fig()
```


    
![png](Integrated_files/Integrated_145_0.png)
    



```R
# Plot DoubletDensity scores
fig(width = 16, height = 7)
as.data.frame(colData(combined)[,c("Sample","label","DoubletDensity")]) %>% 
ggplot(aes(label, log1p(DoubletDensity), color = Sample)) + 
    geom_boxplot(position = position_dodge(width = 0.5), width = 0.3, outlier.size = 0.5) + 
    scale_color_manual(values = c_sample_col) + theme_cowplot(16) + 
    theme(legend.position = "top", axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
    ylab("log Score")
reset.fig()
```


    
![png](Integrated_files/Integrated_146_0.png)
    


## Visualising gene expressions in cells

### Single gene expression

In the following plots we visualise the expression of `XIST` (expressed in female) and `DDX3Y` (expressed in male) in individual cells in human samples. The genes in mouse are `Xist` and `Ddx3y` respectively.

<div class="alert alert-warning">
    For <strong>Flex Gene Expression</strong> dataset, <code>XIST</code>/<code>Xist</code> is not available in the feature list, use other X-linked genes instead, such as <code>RLIM</code>/<code>Rlim</code>, <code>LAMP2</code>/<code>Lamp2</code>, and <code>ATRX</code>/<code>Atrx</code>.
</div>


```R
p1 <- plotProjection(combined, "RLIM", "MNN-TSNE", feat_color = c_heatmap_col1, other_fields = "Sample", 
                     point_size = 0.1, theme_size = 16) + facet_wrap(~ Sample, nrow = 1)
p2 <- plotProjection(combined, "DDX3Y", "MNN-TSNE", feat_color = c_heatmap_col1, other_fields = "Sample", 
                     point_size = 0.1, theme_size = 16) + facet_wrap(~ Sample, nrow = 1)

fig(width = 16, height = 12)
plot_grid(p1, p2, align = "vh", nrow = 2)
reset.fig()
```


    
![png](Integrated_files/Integrated_148_0.png)
    


## Output average expression (logcounts) 

Set outfile prefix


```R
# Set outfile ID
file_id <- "160k_All"
file_id
```


'160k_All'


Create a `dgCMatrix` logcounts matrix for faster computation in some steps.


```R
sce_l <- as(logcounts(combined, withDimnames = FALSE), "dgCMatrix")
rownames(sce_l) <- rownames(combined)
```

### Average `logcounts` expression across samples


```R
# Use logcounts from cdScAnnot
ave.expr.sample <- sumCountsAcrossCells(sce_l, average = TRUE, BPPARAM = bpp,
                                       ids = DataFrame(Sample = combined$Sample)) %>% 
    `colnames<-`(.$Sample) %>% assay %>% as.data.frame %>% rownames_to_column("Symbol")
head(ave.expr.sample)

outfile <- paste0(file_id, "_average_logcounts_in_samples.tsv")
print(paste("Write to file:", outfile))
write.table(ave.expr.sample, file = outfile, sep = "\t", quote = F, row.names = F, col.names = T)
```


<table class="dataframe">
<caption>A data.frame: 6 × 4</caption>
<thead>
	<tr><th></th><th scope=col>Symbol</th><th scope=col>Control1</th><th scope=col>KidneyCancer</th><th scope=col>LungCancer</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>1</th><td>SAMD11 </td><td>0.0000000000</td><td>0.0002540986</td><td>0.000347090</td></tr>
	<tr><th scope=row>2</th><td>NOC2L  </td><td>0.3079627448</td><td>0.1777618222</td><td>0.239002679</td></tr>
	<tr><th scope=row>3</th><td>KLHL17 </td><td>0.0767702381</td><td>0.0262869859</td><td>0.054636288</td></tr>
	<tr><th scope=row>4</th><td>PLEKHN1</td><td>0.0160913945</td><td>0.0118817931</td><td>0.009661303</td></tr>
	<tr><th scope=row>5</th><td>PERM1  </td><td>0.0001657382</td><td>0.0002445113</td><td>0.002616299</td></tr>
	<tr><th scope=row>6</th><td>HES4   </td><td>0.1119948700</td><td>0.0252295532</td><td>0.027066362</td></tr>
</tbody>
</table>



    [1] "Write to file: 160k_All_average_logcounts_in_samples.tsv"


### Average `logcounts` expression across clusters


```R
# Use logcounts from cdScAnnot
ave.expr.label <- sumCountsAcrossCells(sce_l, average = TRUE, BPPARAM = bpp,
                                       ids = DataFrame(cluster = combined$label)) %>% 
    `colnames<-`(paste0("Cluster", .$cluster)) %>% assay %>% as.data.frame %>% rownames_to_column("Symbol")
head(ave.expr.label)

outfile <- paste0(file_id, "_average_logcounts_in_clusters.tsv")
print(paste("Write to file:", outfile))
write.table(ave.expr.label, file = outfile, sep = "\t", quote = F, row.names = F, col.names = T)
```


<table class="dataframe">
<caption>A data.frame: 6 × 47</caption>
<thead>
	<tr><th></th><th scope=col>Symbol</th><th scope=col>Cluster1</th><th scope=col>Cluster2</th><th scope=col>Cluster3</th><th scope=col>Cluster4</th><th scope=col>Cluster5</th><th scope=col>Cluster6</th><th scope=col>Cluster7</th><th scope=col>Cluster8</th><th scope=col>Cluster9</th><th scope=col>⋯</th><th scope=col>Cluster37</th><th scope=col>Cluster38</th><th scope=col>Cluster39</th><th scope=col>Cluster40</th><th scope=col>Cluster41</th><th scope=col>Cluster42</th><th scope=col>Cluster43</th><th scope=col>Cluster44</th><th scope=col>Cluster45</th><th scope=col>Cluster46</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>⋯</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>1</th><td>SAMD11 </td><td>0.000000000</td><td>0.000000000</td><td>0.0000000000</td><td>0.0001851207</td><td>0.000281499</td><td>0.0000000000</td><td>0.0004142627</td><td>0.000000000</td><td>0.000000000</td><td>⋯</td><td>0.003930734</td><td>0.000000000</td><td>0.000000000</td><td>0.00000000</td><td>0.00000000</td><td>0.000000000</td><td>0.00000000</td><td>0.000000000</td><td>0.00000000</td><td>0.000000000</td></tr>
	<tr><th scope=row>2</th><td>NOC2L  </td><td>0.289073072</td><td>0.228131873</td><td>0.2552326726</td><td>0.2005263889</td><td>0.177815553</td><td>0.1446207341</td><td>0.1837185189</td><td>0.262239658</td><td>0.184388779</td><td>⋯</td><td>0.299139082</td><td>0.318064380</td><td>0.367755827</td><td>0.22306059</td><td>0.26079603</td><td>0.237045133</td><td>0.33632238</td><td>0.282463396</td><td>0.32322653</td><td>0.366947259</td></tr>
	<tr><th scope=row>3</th><td>KLHL17 </td><td>0.060047975</td><td>0.016069338</td><td>0.0575419993</td><td>0.0319562699</td><td>0.026365587</td><td>0.0207637371</td><td>0.0395390713</td><td>0.063304630</td><td>0.062215989</td><td>⋯</td><td>0.078312564</td><td>0.074023435</td><td>0.043037060</td><td>0.07034118</td><td>0.07672287</td><td>0.054938099</td><td>0.12193443</td><td>0.063856410</td><td>0.07120116</td><td>0.066152143</td></tr>
	<tr><th scope=row>4</th><td>PLEKHN1</td><td>0.005441157</td><td>0.002791956</td><td>0.0037047209</td><td>0.0067637874</td><td>0.015719263</td><td>0.0160182473</td><td>0.0152492059</td><td>0.030616690</td><td>0.015085388</td><td>⋯</td><td>0.036951024</td><td>0.016046641</td><td>0.005048371</td><td>0.02249621</td><td>0.01225767</td><td>0.014251745</td><td>0.02315386</td><td>0.008603543</td><td>0.04525137</td><td>0.022012800</td></tr>
	<tr><th scope=row>5</th><td>PERM1  </td><td>0.000000000</td><td>0.000000000</td><td>0.0008175283</td><td>0.0012098753</td><td>0.001723708</td><td>0.0006784237</td><td>0.0010308344</td><td>0.000000000</td><td>0.004149378</td><td>⋯</td><td>0.002706363</td><td>0.001836137</td><td>0.000000000</td><td>0.00000000</td><td>0.00000000</td><td>0.002219976</td><td>0.00000000</td><td>0.000000000</td><td>0.00000000</td><td>0.005395099</td></tr>
	<tr><th scope=row>6</th><td>HES4   </td><td>0.005462020</td><td>0.014186765</td><td>0.0180847976</td><td>0.0137838744</td><td>0.020078066</td><td>0.0138197940</td><td>0.0131753584</td><td>0.006752527</td><td>0.006932925</td><td>⋯</td><td>0.108856059</td><td>0.055838898</td><td>0.000000000</td><td>0.02814627</td><td>0.02894928</td><td>0.023659862</td><td>0.03694294</td><td>0.026656862</td><td>0.20506764</td><td>0.323875106</td></tr>
</tbody>
</table>



    [1] "Write to file: 160k_All_average_logcounts_in_clusters.tsv"


### Average `logcounts` expression across conditions


```R
ave.expr.cond <- sumCountsAcrossCells(sce_l, average = TRUE, BPPARAM = bpp,
                                      ids = DataFrame(condition = combined$condition)) %>% 
    `colnames<-`(.$condition) %>% assay %>% as.data.frame %>% rownames_to_column("Symbol")
head(ave.expr.cond)

outfile <- paste0(file_id, "_average_logcounts_in_conditions.tsv")
print(paste("Write to file:", outfile))

write.table(ave.expr.cond, file = outfile, sep = "\t", quote = F, row.names = F, col.names = T)
```


<table class="dataframe">
<caption>A data.frame: 6 × 3</caption>
<thead>
	<tr><th></th><th scope=col>Symbol</th><th scope=col>Control</th><th scope=col>Cancer</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>1</th><td>SAMD11 </td><td>0.0000000000</td><td>0.0002958678</td></tr>
	<tr><th scope=row>2</th><td>NOC2L  </td><td>0.3079627448</td><td>0.2052695355</td></tr>
	<tr><th scope=row>3</th><td>KLHL17 </td><td>0.0767702381</td><td>0.0390207151</td></tr>
	<tr><th scope=row>4</th><td>PLEKHN1</td><td>0.0160913945</td><td>0.0108844097</td></tr>
	<tr><th scope=row>5</th><td>PERM1  </td><td>0.0001657382</td><td>0.0013098532</td></tr>
	<tr><th scope=row>6</th><td>HES4   </td><td>0.1119948700</td><td>0.0260545975</td></tr>
</tbody>
</table>



    [1] "Write to file: 160k_All_average_logcounts_in_conditions.tsv"


### Average `logcounts` expression across condition and clusters


```R
ave.expr.gp <- sumCountsAcrossCells(sce_l, average = TRUE, BPPARAM = bpp,
                                    ids = DataFrame(group = paste(combined$label, combined$condition, sep = "_"))) %>% 
    `colnames<-`(paste0("Cluster", .$group)) 

ave.expr.gp$group <- factor(ave.expr.gp$group, levels = gtools::mixedsort(unique(ave.expr.gp$group)))
ave.expr.gp <- ave.expr.gp[,order(ave.expr.gp$group)]
ave.expr.gp <- ave.expr.gp %>% assay %>% as.data.frame %>% rownames_to_column("Symbol")
head(ave.expr.gp)

outfile <- paste0(file_id, "_average_logcounts_in_groups.tsv")
print(paste("Write to file:", outfile))

write.table(ave.expr.gp, file = outfile, sep = "\t", quote = F, row.names = F, col.names = T)
```


<table class="dataframe">
<caption>A data.frame: 6 × 89</caption>
<thead>
	<tr><th></th><th scope=col>Symbol</th><th scope=col>Cluster1_Cancer</th><th scope=col>Cluster1_Control</th><th scope=col>Cluster2_Cancer</th><th scope=col>Cluster2_Control</th><th scope=col>Cluster3_Cancer</th><th scope=col>Cluster3_Control</th><th scope=col>Cluster4_Cancer</th><th scope=col>Cluster4_Control</th><th scope=col>Cluster5_Cancer</th><th scope=col>⋯</th><th scope=col>Cluster42_Cancer</th><th scope=col>Cluster42_Control</th><th scope=col>Cluster43_Cancer</th><th scope=col>Cluster43_Control</th><th scope=col>Cluster44_Cancer</th><th scope=col>Cluster44_Control</th><th scope=col>Cluster45_Cancer</th><th scope=col>Cluster45_Control</th><th scope=col>Cluster46_Cancer</th><th scope=col>Cluster46_Control</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>⋯</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>1</th><td>SAMD11 </td><td>0.000000000</td><td>0.000000000</td><td>0.00000000</td><td>0</td><td>0.000000000</td><td>0.000000000</td><td>0.0001923742</td><td>0.00000000</td><td>0.0002842517</td><td>⋯</td><td>0.000000000</td><td>0.00000000</td><td>0.00000000</td><td>0.00000000</td><td>0.000000000</td><td>0.00000000</td><td>0.00000000</td><td>0.00000000</td><td>0.00000000</td><td>0.000000000</td></tr>
	<tr><th scope=row>2</th><td>NOC2L  </td><td>0.158008350</td><td>0.313234146</td><td>0.22847858</td><td>0</td><td>0.186600208</td><td>0.282541353</td><td>0.1961165309</td><td>0.31307246</td><td>0.1772989458</td><td>⋯</td><td>0.204633090</td><td>0.40450736</td><td>0.23687279</td><td>0.34168486</td><td>0.253968789</td><td>0.32364158</td><td>0.26773367</td><td>0.35076844</td><td>0.33753737</td><td>0.375555030</td></tr>
	<tr><th scope=row>3</th><td>KLHL17 </td><td>0.026161921</td><td>0.066294686</td><td>0.01609376</td><td>0</td><td>0.026093871</td><td>0.070055128</td><td>0.0301900701</td><td>0.07703229</td><td>0.0259657770</td><td>⋯</td><td>0.045589660</td><td>0.10323837</td><td>0.08412407</td><td>0.12397322</td><td>0.058159362</td><td>0.07208934</td><td>0.05320892</td><td>0.08013097</td><td>0.03701892</td><td>0.074678941</td></tr>
	<tr><th scope=row>4</th><td>PLEKHN1</td><td>0.002080362</td><td>0.006060701</td><td>0.00279620</td><td>0</td><td>0.004585689</td><td>0.003354186</td><td>0.0065048462</td><td>0.01337235</td><td>0.0151928461</td><td>⋯</td><td>0.014240085</td><td>0.01431199</td><td>0.07604164</td><td>0.02030206</td><td>0.007506229</td><td>0.01018929</td><td>0.04235259</td><td>0.04669008</td><td>0.00000000</td><td>0.028455571</td></tr>
	<tr><th scope=row>5</th><td>PERM1  </td><td>0.000000000</td><td>0.000000000</td><td>0.00000000</td><td>0</td><td>0.001917775</td><td>0.000379743</td><td>0.0012572815</td><td>0.00000000</td><td>0.0017405632</td><td>⋯</td><td>0.001560006</td><td>0.00562982</td><td>0.00000000</td><td>0.00000000</td><td>0.000000000</td><td>0.00000000</td><td>0.00000000</td><td>0.00000000</td><td>0.00000000</td><td>0.006974153</td></tr>
	<tr><th scope=row>6</th><td>HES4   </td><td>0.000000000</td><td>0.006468914</td><td>0.01420833</td><td>0</td><td>0.022829395</td><td>0.016196935</td><td>0.0139107089</td><td>0.01054687</td><td>0.0202744032</td><td>⋯</td><td>0.023182708</td><td>0.02612516</td><td>0.00000000</td><td>0.03893496</td><td>0.025618463</td><td>0.02815748</td><td>0.08773256</td><td>0.26330274</td><td>0.00000000</td><td>0.418667820</td></tr>
</tbody>
</table>



    [1] "Write to file: 160k_All_average_logcounts_in_groups.tsv"


# 6. Expression of manually selected genes

<div class="alert alert-info">
    <strong>Tip!</strong> Use known markers to help manual annotation of cell clusters
</div>

**Show expression profiles of manually selected genes.**

T cells: CD3D, CD3E, CD3G

- CD4+ T: CD4, CD28, IL7R
  - Naive: (see naive markers)
  - Memory: ICOS, CCR4, GPR183, SLC2A3, NR3C1
  - Th17: RORA, CCR6, RORC
  - Effector: CCL5, CXCR3, GZMK, PDE4B
  - Regulatory T (Treg): FOXP3, STAM, CTLA4
  - CD40LG

- CD8+ T: CD8A, KLRK1
  - Naive: (see naive markers)
  - Tumour-infiltrating (TIL) (with Naive markers): ZEB1, SETD1B, TMEM123, STAT3
  - Activated: DUSP4, IL21R, DUSP2, CXCR3
  - Terminally differentiated effector (TE): CX3CR1, KLRG1, ADRB2
  - gamma delta (γδ) T: PTPRC (CD45), TRDC, not KLRC1
    
- Double-negative (DN) T: CD3+ CD4- CD8-
- Double-positive (DP) T: CD4, CD8A, CXCR3, CCR5
- Naive: ATM, LEF1, CCR7, TCF7, NELL2, TGFBR2
- Cytotoxic : GZMH, CST7, GNLY, CTSW, FGFBP2, GZMB, PRF1, KLRC1, KLRC3, KLRC4
- Inflamed: CCL3, CCL4, IFIT2, IFIT3, ZC3HAV1

Natural Killer (NK) cells: CD7, MATK, IL2RB, KLRF1, NCR1

Innate lymphoid cells: GATA3, ICAM3
   - Group 2 (ILC2): RORA, MAF, PTGDR2, MBOAT2

B cells: CD79A, CD22, MS4A1 (CD20)
  - Naive: TCL1A, BACH2, IGHM, IGHD
  - Activated: CD83, CD55, BACH1, CXCR4, JUND
  - Unswitched memory (USM): ARHGAP25, CD24,PARP15, TCF4
  - Switched memory (SM): TNFRSF13B, IGHA1, IGHG1

Monocytes:
- CD14+ (classical): CD14, LYZ, VCAN, S100A8, S100A9, S100A12, TREM1
  - Tumour-associated: VEGFA, HIF1A, IER3, THBD, FOSL1
  - Cytokine stimulated: WARS1, IFITM3, GBP1
- CD16+ (non-classical): SIGLEC10, TCF7L2, FCGR3A (CD16a), CDKN1C, HES4, CX3CR1

Dendritic cells:
- CD1C+ DCs: CD1C, CD1D, CD33, CLEC10A, FCER1A
- Plasmacytoid dendritic cellS (pDCs): ITM2C, TCF4, LILRA4, CLEC4C, SERPINF1, MZB1


```R
geneNames <- c("CD3E","CD3D","CD3G","CD4","CD8A","TGFBR2","CCR7","LEF1","IL7R","TCF7","NELL2","SETD1B","STAT3",
               "NR3C1","ZEB1","TMEM123","CCR4","CD28","ICOS","CCR6","GPR183","SLC2A3","RORA","MAF","CD40LG",
               "STAM","FOXP3","CTLA4","GZMB","GNLY","GZMH","FGFBP2","KLRG1","ADRB2","CCL5","CST7","KLRK1","CTSW",
               "MATK","DUSP4","PDE4B","CXCR4","IL21R","DUSP2","GZMK","CXCR3","CCL3","CCL4","ZC3HAV1","IFIT2",
               "IFIT3","KLRC1","KLRC3","KLRC4","CCR5","KLRF1","NCR1","TRDC","ATM","CD7","PRF1","IL2RB","GATA3",
               "ICAM3","RORC","PTGDR2","MBOAT2","BACH2","CD24","CD22","CD79A","MS4A1","TCL1A","IGHM","IGHD",
               "CD83","BACH1","CD55","JUND","IGHA1","TNFRSF13B","IGHG1","ARHGAP25","PARP15","S100A12","TREM1",
               "VEGFA","LYZ","CD14","VCAN","S100A8","S100A9","GBP1","TCF7L2","SIGLEC10","WARS1","IFITM3",
               "HES4","CDKN1C","FCGR3A","CX3CR1","PTPRC","HIF1A","THBD","IER3","FOSL1","CD1D","CD33",
               "FCER1A","CD1C","CLEC10A","TCF4","MZB1","ITM2C","LILRA4","CLEC4C","SERPINF1")
length(geneNames)
```


117



```R
fig(width = 16, height = 28)
plotDots(combined, features = geneNames, group = "label", zlim = c(-3, 3),
              center = TRUE, scale = TRUE) + scale_size(limits = c(0, 1), range = c(0.1, 6)) + 
    guides(colour = guide_colourbar(title = "Row (Gene) Z-Score", barwidth = 10), 
           size = guide_legend(title = "Proportion Detected")) + theme_cowplot(16) + #coord_flip() + 
    theme(panel.grid.major = element_line(colour = "gray90"), 
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
          legend.position = "top", legend.justification = "center", legend.title.position = "top") + 
    labs(x = "Cluster", y = "Marker Genes")
reset.fig()
```

    [1m[22mScale for [32msize[39m is already present.
    Adding another scale for [32msize[39m, which will replace the existing scale.



    
![png](Integrated_files/Integrated_163_1.png)
    



```R
dimname <- "MNN-UMAP"

p <- as.data.frame(reducedDim(combined, dimname)) %>% 
    bind_cols(as_tibble(t(logcounts(combined[geneNames,])))) %>% rename(V1 = 1, V2 = 2) %>% 
    gather(., key = "Symbol", value ="Expression", -c(V1, V2) ) %>% 
    mutate_at(vars(Symbol), factor) %>% mutate(Symbol = factor(Symbol, levels = geneNames)) %>%
    ggplot(aes(x = V1, y = V2, color = Expression)) + geom_point(size = 0.1, alpha = 0.1) + 
    facet_wrap(~ Symbol, ncol = 5) + scale_color_viridis(option = "plasma", direction = -1) +
    theme_classic(base_size = 20) + labs(x = paste(dimname, "1"), y = paste(dimname, "2"))

fig(width = 16, height = 70)
p
reset.fig()
```


    
![png](Integrated_files/Integrated_164_0.png)
    


# 7. Define per-cluster cell types

Introduces 2 new cell type labels:

- `CellType_1`: Coarse cell type annotation (same as `ClusterCellType`)
- `CellType_2`: Fine (per-cluster) cell type annotation

## Define "Coarse Cell Type" annotation


```R
combined$CellType_1 <- combined$label
levels(combined$CellType_1) <- c("CD8+ T","CD8+ T","CD4+ T","CD4+ T","CD4+ T","CD4+ T","CD4+ T","CD4+ T","CD4+ T",
                                 "CD4+ T","CD4+ T","CD4+ T","CD8+ T","CD8+ T","CD8+ T","CD8+ T","CD8+ T","CD8+ T",
                                 "CD8+ T","CD8+ T","gdT","gdT","CD8+ T","CD8+ T","Other T","Other T","NK","NK",
                                 "ILC","B","B","B","B","Monocytes","Monocytes","Monocytes","Monocytes","DCs",
                                 "DCs","Unknown","Unknown","Doublets","Doublets","Doublets","Doublets","Doublets")
combined$ClusterCellType <- combined$CellType_1
```


```R
table("Coarse Cell Type" = combined$CellType_1, "Condition" = combined$condition)
```


                    Condition
    Coarse Cell Type Control Cancer
           CD8+ T       3507   4550
           CD4+ T       5826  16270
           gdT            23   1123
           Other T       323    493
           NK           1997   1534
           ILC            28     24
           B            1563   2628
           Monocytes    3107   1411
           DCs           366     63
           Unknown       163   1060
           Doublets     1135   1300



```R
table("Coarse Cell Type" = combined$CellType_1, "Clusters" = combined$label)
```


                    Clusters
    Coarse Cell Type    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18
           CD8+ T    2557  659    0    0    0    0    0    0    0    0    0    0  568 1130  376  125  270  295
           CD4+ T       0    0 4388 4323 1962 3339 2141 2670  849 1610  123  691    0    0    0    0    0    0
           gdT          0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
           Other T      0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
           NK           0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
           ILC          0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
           B            0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
           Monocytes    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
           DCs          0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
           Unknown      0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
           Doublets     0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
                    Clusters
    Coarse Cell Type   19   20   21   22   23   24   25   26   27   28   29   30   31   32   33   34   35   36
           CD8+ T     469  248    0    0  874  486    0    0    0    0    0    0    0    0    0    0    0    0
           CD4+ T       0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
           gdT          0    0  346  800    0    0    0    0    0    0    0    0    0    0    0    0    0    0
           Other T      0    0    0    0    0    0  372  444    0    0    0    0    0    0    0    0    0    0
           NK           0    0    0    0    0    0    0    0 3378  153    0    0    0    0    0    0    0    0
           ILC          0    0    0    0    0    0    0    0    0    0   52    0    0    0    0    0    0    0
           B            0    0    0    0    0    0    0    0    0    0    0 1062 1850  434  845    0    0    0
           Monocytes    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0 1202 1266  788
           DCs          0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
           Unknown      0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
           Doublets     0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
                    Clusters
    Coarse Cell Type   37   38   39   40   41   42   43   44   45   46
           CD8+ T       0    0    0    0    0    0    0    0    0    0
           CD4+ T       0    0    0    0    0    0    0    0    0    0
           gdT          0    0    0    0    0    0    0    0    0    0
           Other T      0    0    0    0    0    0    0    0    0    0
           NK           0    0    0    0    0    0    0    0    0    0
           ILC          0    0    0    0    0    0    0    0    0    0
           B            0    0    0    0    0    0    0    0    0    0
           Monocytes 1262    0    0    0    0    0    0    0    0    0
           DCs          0  318  111    0    0    0    0    0    0    0
           Unknown      0    0    0  796  427    0    0    0    0    0
           Doublets     0    0    0    0    0  703  215  802  609  106



```R
# Set colours for CellType_1
set.seed(10010)
c_celltype_col2 <- choosePalette(combined$CellType_1, sample(pals::cols25(nlevels(combined$CellType_1))))

# Coloured by CellType_1
fig(width = 16, height = 9)
plotProjections(combined, "CellType_1", dimnames = c("MNN-TSNE", "MNN-UMAP"), feat_desc = "Coarse Cell Type", 
                feat_color = c_celltype_col2, text_by = "label", point_size = 0.1, point_alpha = 0.1,
                text_size = 6, guides_nrow = 1, guides_size = 4, legend_pos = "bottom")
reset.fig()
```


    
![png](Integrated_files/Integrated_169_0.png)
    



```R
fig(width = 16, height = 6)
data.frame(table("ct" = combined$CellType_1, "label" = combined$label, "Condition" = combined$condition)) %>% 
    ggplot(aes(Condition, Freq, fill = ct)) + geom_col() + 
    facet_wrap(~ ct, nrow = 2, scales = "free_y") + scale_fill_manual(values = c_celltype_col2) + 
    theme_cowplot(16) + guides(fill = guide_legend("Coarse CT", ncol = 1)) +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1)) + ylab("Number of cells")
reset.fig()
```


    
![png](Integrated_files/Integrated_170_0.png)
    



```R
fig(width = 16, height = 6)
data.frame(prop.table(table("Condition" = combined$condition, "ct" = combined$CellType_1, "label" = combined$label), 1)) %>% 
    ggplot(aes(Condition, Freq, fill = ct)) + geom_col() + 
    facet_wrap(~ ct, nrow = 2, scales = "free_y") + scale_fill_manual(values = c_celltype_col2) + 
    scale_y_continuous(labels = percent) + theme_cowplot(16) + guides(fill = guide_legend("Coarse CT", ncol = 1)) + 
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1)) + ylab("% Cells in condition")
reset.fig()
```


    
![png](Integrated_files/Integrated_171_0.png)
    


## Define "Fine Cell Type" annotation


```R
combined$CellType_2 <- combined$label
levels(combined$CellType_2) <- paste(1:nlevels(combined$label), 
                                     c("Naive CD8+ T","CD8+ TILs","Naive CD4+ T","ICOS(hi) Momory CD4+ T",
                                       "ICOS(hi) Momory CD4+ T","Th17-like T","Effector CD4+ T","CD40LG+ CD4+ T",
                                       "Memory Tregs","CTSW- Cytotoxic CD4+ T","Cytotoxic CD4+ T",
                                       "Inflamed Cytotoxic CD4+ T","Naive-like Activated CD8+ T",
                                       "Activated CD8(hi) T","Activated CD8(lo) T","Inflamed Activated CD8+ T",
                                       "CD16- Cytotoxic CD8+ T","DUSP2/4+ Cytotoxic CD8+ T",
                                       "CD16+ Cytotoxic CD8+ T","CD16+ Inflamed Cytotoxic CD8+ T","CD45(lo) gdT",
                                       "CD45(hi) gdT","TE CD8(hi) T","TE CD8(lo) T","DP T","Activated DN T","NK",
                                       "Inflamed NK","ILC2","Naive B","Activated B","USM B","SM B",
                                       "CD14 Monocytes","Stimulated CD14 Monocytes","CD16 Monocytes","CD14 TAM",
                                       "CD1C+ DCs","pDCs","Cytotoxic-like","Momory-like","Naive T/TIL/NK",
                                       "Naive T/Effect T/NK","B/T","Monocytes/T","B/Monocytes"))
```


```R
table("Fine Cell Type" = combined$CellType_2, "Condition" = combined$condition)
```


                                        Condition
    Fine Cell Type                       Control Cancer
      1 Naive CD8+ T                        2159    398
      2 CD8+ TILs                              1    658
      3 Naive CD4+ T                        3139   1249
      4 ICOS(hi) Momory CD4+ T               163   4160
      5 ICOS(hi) Momory CD4+ T                19   1943
      6 Th17-like T                           24   3315
      7 Effector CD4+ T                       14   2127
      8 CD40LG+ CD4+ T                      2161    509
      9 Memory Tregs                         288    561
      10 CTSW- Cytotoxic CD4+ T               13   1597
      11 Cytotoxic CD4+ T                      3    120
      12 Inflamed Cytotoxic CD4+ T             2    689
      13 Naive-like Activated CD8+ T          87    481
      14 Activated CD8(hi) T                   1   1129
      15 Activated CD8(lo) T                   1    375
      16 Inflamed Activated CD8+ T             0    125
      17 CD16- Cytotoxic CD8+ T                4    266
      18 DUSP2/4+ Cytotoxic CD8+ T             1    294
      19 CD16+ Cytotoxic CD8+ T               25    444
      20 CD16+ Inflamed Cytotoxic CD8+ T       0    248
      21 CD45(lo) gdT                         10    336
      22 CD45(hi) gdT                         13    787
      23 TE CD8(hi) T                        774    100
      24 TE CD8(lo) T                        454     32
      25 DP T                                323     49
      26 Activated DN T                        0    444
      27 NK                                 1997   1381
      28 Inflamed NK                           0    153
      29 ILC2                                 28     24
      30 Naive B                             993     69
      31 Activated B                           5   1845
      32 USM B                               248    186
      33 SM B                                317    528
      34 CD14 Monocytes                     1174     28
      35 Stimulated CD14 Monocytes          1201     65
      36 CD16 Monocytes                      729     59
      37 CD14 TAM                              3   1259
      38 CD1C+ DCs                           286     32
      39 pDCs                                 80     31
      40 Cytotoxic-like                       68    728
      41 Momory-like                          95    332
      42 Naive T/TIL/NK                      114    589
      43 Naive T/Effect T/NK                 204     11
      44 B/T                                 328    474
      45 Monocytes/T                         407    202
      46 B/Monocytes                          82     24



```R
# Set colours for CellType_2
c_celltype_col3 <- choosePalette(combined$CellType_2, c_clust_col) # same colour as cluster colours

# Coloured by CellType_2
fig(width = 16, height = 12)
plotProjections(combined, "CellType_2", dimnames = c("MNN-TSNE", "MNN-UMAP"), feat_desc = "Fine Cell Type", 
                feat_color = c_celltype_col3, text_by = "label", point_size = 0.1, point_alpha = 0.1,
                text_size = 6, guides_nrow = 12, guides_size = 4, legend_pos = "bottom", rel_height = c(5, 2))
reset.fig()
```


    
![png](Integrated_files/Integrated_175_0.png)
    



```R
fig(width = 16, height = 22)
data.frame(table("ct" = combined$CellType_2, "label" = combined$label, "Condition" = combined$condition)) %>% 
    ggplot(aes(Condition, Freq, fill = ct)) + geom_col() + 
    facet_wrap(~ ct, ncol = 4, scales = "free_y") + scale_fill_manual(values = c_celltype_col3) + 
    theme_cowplot(16) + guides(fill = guide_legend("Fine CT", ncol = 1)) +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1)) + ylab("Number of cells")
reset.fig()
```


    
![png](Integrated_files/Integrated_176_0.png)
    



```R
fig(width = 16, height = 22)
data.frame(prop.table(table("Condition" = combined$condition, "ct" = combined$CellType_2, "label" = combined$label), 1)) %>% 
    ggplot(aes(Condition, Freq, fill = ct)) + geom_col() + 
    facet_wrap(~ ct, ncol = 4, scales = "free_y") + scale_fill_manual(values = c_celltype_col3) + 
    scale_y_continuous(labels = percent) + theme_cowplot(16) + guides(fill = guide_legend("Fine CT", ncol = 1)) + 
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1)) + ylab("% Cells in condition")
reset.fig()
```


    
![png](Integrated_files/Integrated_177_0.png)
    


# 8. Marker gene detection (per-cluster)

Potential marker genes are identified by taking the top set of DE genes from each pairwise comparison between clusters. The results are arranged into a single output table that allows a marker set to be easily defined for a user-specified size of the top set. For example, to construct a marker set from the top 10 genes of each comparison, one would filter `marker.set` to retain rows with Top less than or equal to 10.

<div class="alert alert-warning">
  <strong>Warning!</strong> Take note of any poor clusters in the silhouette plot and consider how trustworthy the marker genes actually are.
</div>

Here we are going to use `findMarkers` from `scran` package as this is much faster. 

## Decided the `pval.type` and `min.prop` settings

- Use `pval.type` to specify how p-values are to be combined across pairwise comparisons for a given group/cluster. Defaults `pval.type` is `"any"`.
- Use `min.prop` to specify the minimum proportion of significant comparisons per gene. Defaults `min.prop` is 0.5 when `pval.type="some"`, and zero for all other `pval.type` options.

The choice of `pval.type` determines whether the highly ranked genes are those that are DE between the current group and:

### 1. DE against any other cluster ("any")

If `pval.type="any"`, the null hypothesis is that the **gene is not DE in as least 1 contrasts**. This approach does not explicitly favour genes that are uniquely expressed in a cluster. Rather, it focuses on combinations of genes that - together - drive separation of a cluster from the others. This is more general and robust but tends to yield a less focused marker set compared to the other `pval.type` settings.

Using `pval.type="any"`, the result will contain a `Top` column that shows the minimum rank across all pairwise comparisons. For example, if we define a marker set with an `T` of 1 for a given cluster. The set of genes with `Top <= 1` will contain the top gene from each pairwise comparison to every other cluster. If `T` is `5`, the set will consist of the union of the top 5 genes from each pairwise comparison. This approach does not explicitly favour genes that are uniquely expressed in a cluster. Rather, it focuses on combinations of genes that - together - drive separation of a cluster from the others. This is more general and robust but tends to yield a less focused marker set compared to the other `pval.type` settings.

### 2. DE against all other clusters ("all")

If `pval.type="all"`, the null hypothesis is that the **gene is not DE in all contrasts**. This strategy is particularly effective when dealing with distinct clusters that have a unique expression profile. In such cases, it yields a highly focused marker set that concisely captures the differences between clusters. However, it can be too stringent if the cluster's separation is driven by combinations of gene expression. 

### 3. DE against some (50%) other clusters ("some")

The `pval.type="some"` setting serves as a compromise between "all" and "any". A combined p-value is calculated by taking the middlemost value of the Holm-corrected p-values for each gene. Here, the null hypothesis is that the **gene is not DE in at least half of the contrasts**.

### 4. DE against some other clusters (rank-style)

This is achieved by setting `pval.type="any"` with `min.prop` set to some positive value in (0, 1). Here, we are selecting high-ranked genes that are among the top-ranked (`Top`) genes in at least `min.prop` of the pairwise comparisons

For example, if `pval.type="any", min.prop=0.3`, any gene with a value of `Top` less than or equal to 5 will be in the top 5 DEGs of at least 30% of the comparisons. This method increases the stringency of the `"any"` setting in a safer manner than `pval.type="some"`.

More explanation can be found by `?findMarkers`, `?combineMarkers` and `?combinePValues`.

## Find marker genes for clusters

### Run `findMarkers` (both directions)

Considers both up- and downregulated genes to be potential markers.


```R
# Set pval.type (and min.prop if using "any")
pval.type <- "any"
min.prop <- 0.3

# Using logcounts from combined
marker.genes.cluster <- findMarkers(sce_l, groups = combined$label, pval.type = pval.type, min.prop = min.prop,
                                    block = combined$condition, #blocking on uninteresting factors
                                    BPPARAM = bpp)
marker.genes.cluster
```


    List of length 46
    names(46): 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 ... 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46


Print the number of markers that passed the FDR or `Top` threshold. This will be the number of genes as inut for `enrichR`.


```R
printMarkerStats(marker.genes.cluster, pval.type = pval.type, min.prop = min.prop)
```

    Number of selected markers (Top 200 genes of at least 30.0% comparisons):
    - Cluster1: 167; Up = 74; Down = 93; Max. P-value = 4.5e-138.
    - Cluster2: 164; Up = 95; Down = 69; Max. P-value = 3.6e-160.
    - Cluster3: 175; Up = 110; Down = 65; Max. P-value = 5.4e-221.
    - Cluster4: 152; Up = 106; Down = 46; Max. P-value = 0.
    - Cluster5: 159; Up = 114; Down = 45; Max. P-value = 2.3e-77.
    - Cluster6: 155; Up = 122; Down = 33; Max. P-value = 0.
    - Cluster7: 167; Up = 121; Down = 46; Max. P-value = 2e-310.
    - Cluster8: 145; Up = 79; Down = 66; Max. P-value = 1.4e-274.
    - Cluster9: 166; Up = 79; Down = 87; Max. P-value = 3.1e-119.
    - Cluster10: 140; Up = 90; Down = 50; Max. P-value = 0.
    - Cluster11: 138; Up = 74; Down = 64; Max. P-value = 1.1e-45.
    - Cluster12: 159; Up = 129; Down = 30; Max. P-value = 9.1e-197.
    - Cluster13: 134; Up = 50; Down = 84; Max. P-value = 1.6e-101.
    - Cluster14: 148; Up = 105; Down = 43; Max. P-value = 3.2e-267.
    - Cluster15: 150; Up = 99; Down = 51; Max. P-value = 3.6e-50.
    - Cluster16: 156; Up = 130; Down = 26; Max. P-value = 2.8e-20.
    - Cluster17: 118; Up = 56; Down = 62; Max. P-value = 4.9e-67.
    - Cluster18: 147; Up = 105; Down = 42; Max. P-value = 1.5e-65.
    - Cluster19: 158; Up = 113; Down = 45; Max. P-value = 5.3e-230.
    - Cluster20: 163; Up = 122; Down = 41; Max. P-value = 2.7e-35.
    - Cluster21: 139; Up = 66; Down = 73; Max. P-value = 2.7e-71.
    - Cluster22: 142; Up = 81; Down = 61; Max. P-value = 3.8e-157.
    - Cluster23: 121; Up = 71; Down = 50; Max. P-value = 2.6e-210.
    - Cluster24: 90; Up = 57; Down = 33; Max. P-value = 3.5e-207.
    - Cluster25: 118; Up = 50; Down = 68; Max. P-value = 4.2e-76.
    - Cluster26: 186; Up = 158; Down = 28; Max. P-value = 4.4e-46.
    - Cluster27: 181; Up = 142; Down = 39; Max. P-value = 0.
    - Cluster28: 175; Up = 111; Down = 64; Max. P-value = 1.4e-54.
    - Cluster29: 144; Up = 4; Down = 140; Max. P-value = 1.4e-27.
    - Cluster30: 193; Up = 76; Down = 117; Max. P-value = 5.9e-221.
    - Cluster31: 207; Up = 162; Down = 45; Max. P-value = 0.
    - Cluster32: 200; Up = 70; Down = 130; Max. P-value = 7.5e-95.
    - Cluster33: 208; Up = 94; Down = 114; Max. P-value = 1.7e-179.
    - Cluster34: 210; Up = 155; Down = 55; Max. P-value = 5.5e-227.
    - Cluster35: 217; Up = 152; Down = 65; Max. P-value = 2.4e-242.
    - Cluster36: 226; Up = 145; Down = 81; Max. P-value = 0.
    - Cluster37: 221; Up = 215; Down = 6; Max. P-value = 1.9e-200.
    - Cluster38: 195; Up = 93; Down = 102; Max. P-value = 8.5e-152.
    - Cluster39: 165; Up = 10; Down = 155; Max. P-value = 3.4e-35.
    - Cluster40: 171; Up = 33; Down = 138; Max. P-value = 2.9e-175.
    - Cluster41: 157; Up = 23; Down = 134; Max. P-value = 4.2e-70.
    - Cluster42: 121; Up = 62; Down = 59; Max. P-value = 3.1e-59.
    - Cluster43: 78; Up = 49; Down = 29; Max. P-value = 3.2e-26.
    - Cluster44: 185; Up = 125; Down = 60; Max. P-value = 1.3e-75.
    - Cluster45: 225; Up = 205; Down = 20; Max. P-value = 8e-176.
    - Cluster46: 198; Up = 102; Down = 96; Max. P-value = 4.4e-26.
    * Upregulated when logFC > 0.0 and downregulated when logFC < 0.0.



```R
for(ID in names(marker.genes.cluster)) {
    # Append Ensembl ID and Symbol from rowData
    marker.genes.cluster[[ID]] <- cbind(rowData(combined)[rownames(marker.genes.cluster[[ID]]),][,1:2], 
                                        marker.genes.cluster[[ID]])
}

exportResList(marker.genes.cluster, col_anno = c("ID","Symbol"), prefix = paste0(file_id, "_cluster"))
```

    Detecting findMarkers input.
    Creating file: 160k_All_cluster_findMarkers_Cluster1.tsv
    Creating file: 160k_All_cluster_findMarkers_Cluster2.tsv
    Creating file: 160k_All_cluster_findMarkers_Cluster3.tsv
    Creating file: 160k_All_cluster_findMarkers_Cluster4.tsv
    Creating file: 160k_All_cluster_findMarkers_Cluster5.tsv
    Creating file: 160k_All_cluster_findMarkers_Cluster6.tsv
    Creating file: 160k_All_cluster_findMarkers_Cluster7.tsv
    Creating file: 160k_All_cluster_findMarkers_Cluster8.tsv
    Creating file: 160k_All_cluster_findMarkers_Cluster9.tsv
    Creating file: 160k_All_cluster_findMarkers_Cluster10.tsv
    Creating file: 160k_All_cluster_findMarkers_Cluster11.tsv
    Creating file: 160k_All_cluster_findMarkers_Cluster12.tsv
    Creating file: 160k_All_cluster_findMarkers_Cluster13.tsv
    Creating file: 160k_All_cluster_findMarkers_Cluster14.tsv
    Creating file: 160k_All_cluster_findMarkers_Cluster15.tsv
    Creating file: 160k_All_cluster_findMarkers_Cluster16.tsv
    Creating file: 160k_All_cluster_findMarkers_Cluster17.tsv
    Creating file: 160k_All_cluster_findMarkers_Cluster18.tsv
    Creating file: 160k_All_cluster_findMarkers_Cluster19.tsv
    Creating file: 160k_All_cluster_findMarkers_Cluster20.tsv
    Creating file: 160k_All_cluster_findMarkers_Cluster21.tsv
    Creating file: 160k_All_cluster_findMarkers_Cluster22.tsv
    Creating file: 160k_All_cluster_findMarkers_Cluster23.tsv
    Creating file: 160k_All_cluster_findMarkers_Cluster24.tsv
    Creating file: 160k_All_cluster_findMarkers_Cluster25.tsv
    Creating file: 160k_All_cluster_findMarkers_Cluster26.tsv
    Creating file: 160k_All_cluster_findMarkers_Cluster27.tsv
    Creating file: 160k_All_cluster_findMarkers_Cluster28.tsv
    Creating file: 160k_All_cluster_findMarkers_Cluster29.tsv
    Creating file: 160k_All_cluster_findMarkers_Cluster30.tsv
    Creating file: 160k_All_cluster_findMarkers_Cluster31.tsv
    Creating file: 160k_All_cluster_findMarkers_Cluster32.tsv
    Creating file: 160k_All_cluster_findMarkers_Cluster33.tsv
    Creating file: 160k_All_cluster_findMarkers_Cluster34.tsv
    Creating file: 160k_All_cluster_findMarkers_Cluster35.tsv
    Creating file: 160k_All_cluster_findMarkers_Cluster36.tsv
    Creating file: 160k_All_cluster_findMarkers_Cluster37.tsv
    Creating file: 160k_All_cluster_findMarkers_Cluster38.tsv
    Creating file: 160k_All_cluster_findMarkers_Cluster39.tsv
    Creating file: 160k_All_cluster_findMarkers_Cluster40.tsv
    Creating file: 160k_All_cluster_findMarkers_Cluster41.tsv
    Creating file: 160k_All_cluster_findMarkers_Cluster42.tsv
    Creating file: 160k_All_cluster_findMarkers_Cluster43.tsv
    Creating file: 160k_All_cluster_findMarkers_Cluster44.tsv
    Creating file: 160k_All_cluster_findMarkers_Cluster45.tsv
    Creating file: 160k_All_cluster_findMarkers_Cluster46.tsv


### Run `findMarkers` (upregulated genes)

Set `direction='up'` to only consider upregulated genes as potential markers.


```R
# Set pval.type (and min.prop if using "any")
pval.type <- "any"
min.prop <- 0.3
direction <- "up"

# Using logcounts from combined
marker.genes.cluster.up <- findMarkers(sce_l, groups = combined$label, pval.type = pval.type, min.prop = min.prop,
                                       block = combined$condition, #blocking on uninteresting factors
                                       lfc = 0.5, direction = direction, BPPARAM = bpp)
marker.genes.cluster.up
```


    List of length 46
    names(46): 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 ... 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46


Print the number of markers that passed the FDR or `Top` threshold. This will be the number of genes as inut for `enrichR`.


```R
printMarkerStats(marker.genes.cluster.up, pval.type = pval.type, min.prop = min.prop)
```

    Number of selected markers (Top 200 genes of at least 30.0% comparisons):
    - Cluster1: 262; Up = 262; Down = 0; Max. P-value = 0.084.
    - Cluster2: 251; Up = 251; Down = 0; Max. P-value = 1.1e-10.
    - Cluster3: 250; Up = 250; Down = 0; Max. P-value = 9e-73.
    - Cluster4: 241; Up = 241; Down = 0; Max. P-value = 6.1e-34.
    - Cluster5: 269; Up = 269; Down = 0; Max. P-value = 9.1e-160.
    - Cluster6: 255; Up = 255; Down = 0; Max. P-value = 2.9e-28.
    - Cluster7: 254; Up = 254; Down = 0; Max. P-value = 0.0026.
    - Cluster8: 238; Up = 238; Down = 0; Max. P-value = 0.0058.
    - Cluster9: 247; Up = 247; Down = 0; Max. P-value = 1.
    - Cluster10: 250; Up = 250; Down = 0; Max. P-value = 7.4e-42.
    - Cluster11: 253; Up = 253; Down = 0; Max. P-value = 2e-31.
    - Cluster12: 257; Up = 257; Down = 0; Max. P-value = 9.8e-101.
    - Cluster13: 243; Up = 243; Down = 0; Max. P-value = 1.7e-68.
    - Cluster14: 254; Up = 254; Down = 0; Max. P-value = 1.3e-191.
    - Cluster15: 252; Up = 252; Down = 0; Max. P-value = 0.00014.
    - Cluster16: 256; Up = 256; Down = 0; Max. P-value = 3.9e-16.
    - Cluster17: 239; Up = 239; Down = 0; Max. P-value = 1.
    - Cluster18: 260; Up = 260; Down = 0; Max. P-value = 5.3e-16.
    - Cluster19: 254; Up = 254; Down = 0; Max. P-value = 2.2e-56.
    - Cluster20: 257; Up = 257; Down = 0; Max. P-value = 5.4e-09.
    - Cluster21: 248; Up = 248; Down = 0; Max. P-value = 0.093.
    - Cluster22: 249; Up = 249; Down = 0; Max. P-value = 1e-08.
    - Cluster23: 249; Up = 249; Down = 0; Max. P-value = 0.86.
    - Cluster24: 247; Up = 247; Down = 0; Max. P-value = 2e-73.
    - Cluster25: 247; Up = 247; Down = 0; Max. P-value = 9e-06.
    - Cluster26: 264; Up = 264; Down = 0; Max. P-value = 4.1e-74.
    - Cluster27: 239; Up = 239; Down = 0; Max. P-value = 3.2e-286.
    - Cluster28: 249; Up = 249; Down = 0; Max. P-value = 5.4e-11.
    - Cluster29: 232; Up = 232; Down = 0; Max. P-value = 0.06.
    - Cluster30: 242; Up = 242; Down = 0; Max. P-value = 6e-08.
    - Cluster31: 228; Up = 228; Down = 0; Max. P-value = 3.7e-44.
    - Cluster32: 243; Up = 243; Down = 0; Max. P-value = 4.2e-18.
    - Cluster33: 237; Up = 237; Down = 0; Max. P-value = 3.5e-21.
    - Cluster34: 241; Up = 241; Down = 0; Max. P-value = 1e-211.
    - Cluster35: 267; Up = 267; Down = 0; Max. P-value = 2.8e-264.
    - Cluster36: 253; Up = 253; Down = 0; Max. P-value = 0.
    - Cluster37: 246; Up = 246; Down = 0; Max. P-value = 1e-253.
    - Cluster38: 248; Up = 248; Down = 0; Max. P-value = 7.7e-37.
    - Cluster39: 244; Up = 244; Down = 0; Max. P-value = 0.01.
    - Cluster40: 249; Up = 249; Down = 0; Max. P-value = 1.
    - Cluster41: 247; Up = 247; Down = 0; Max. P-value = 1.
    - Cluster42: 215; Up = 215; Down = 0; Max. P-value = 1.1e-09.
    - Cluster43: 232; Up = 232; Down = 0; Max. P-value = 0.0044.
    - Cluster44: 213; Up = 213; Down = 0; Max. P-value = 2.4e-16.
    - Cluster45: 237; Up = 237; Down = 0; Max. P-value = 1.3e-89.
    - Cluster46: 247; Up = 247; Down = 0; Max. P-value = 5.5e-09.
    * Upregulated when logFC > 0.0 and downregulated when logFC < 0.0.



```R
# Append Ensembl ID and Symbol from rowData
for(ID in names(marker.genes.cluster.up)) {
    marker.genes.cluster.up[[ID]] <- cbind(rowData(combined)[rownames(marker.genes.cluster.up[[ID]]),][,1:2], 
                                           marker.genes.cluster.up[[ID]])
}

exportResList(marker.genes.cluster.up, col_anno = c("ID","Symbol"), prefix = paste0(file_id, "_cluster"), 
              direction = direction)
```

    Detecting findMarkers input.
    Creating file: 160k_All_cluster_findMarkers_upregulated_Cluster1.tsv
    Creating file: 160k_All_cluster_findMarkers_upregulated_Cluster2.tsv
    Creating file: 160k_All_cluster_findMarkers_upregulated_Cluster3.tsv
    Creating file: 160k_All_cluster_findMarkers_upregulated_Cluster4.tsv
    Creating file: 160k_All_cluster_findMarkers_upregulated_Cluster5.tsv
    Creating file: 160k_All_cluster_findMarkers_upregulated_Cluster6.tsv
    Creating file: 160k_All_cluster_findMarkers_upregulated_Cluster7.tsv
    Creating file: 160k_All_cluster_findMarkers_upregulated_Cluster8.tsv
    Creating file: 160k_All_cluster_findMarkers_upregulated_Cluster9.tsv
    Creating file: 160k_All_cluster_findMarkers_upregulated_Cluster10.tsv
    Creating file: 160k_All_cluster_findMarkers_upregulated_Cluster11.tsv
    Creating file: 160k_All_cluster_findMarkers_upregulated_Cluster12.tsv
    Creating file: 160k_All_cluster_findMarkers_upregulated_Cluster13.tsv
    Creating file: 160k_All_cluster_findMarkers_upregulated_Cluster14.tsv
    Creating file: 160k_All_cluster_findMarkers_upregulated_Cluster15.tsv
    Creating file: 160k_All_cluster_findMarkers_upregulated_Cluster16.tsv
    Creating file: 160k_All_cluster_findMarkers_upregulated_Cluster17.tsv
    Creating file: 160k_All_cluster_findMarkers_upregulated_Cluster18.tsv
    Creating file: 160k_All_cluster_findMarkers_upregulated_Cluster19.tsv
    Creating file: 160k_All_cluster_findMarkers_upregulated_Cluster20.tsv
    Creating file: 160k_All_cluster_findMarkers_upregulated_Cluster21.tsv
    Creating file: 160k_All_cluster_findMarkers_upregulated_Cluster22.tsv
    Creating file: 160k_All_cluster_findMarkers_upregulated_Cluster23.tsv
    Creating file: 160k_All_cluster_findMarkers_upregulated_Cluster24.tsv
    Creating file: 160k_All_cluster_findMarkers_upregulated_Cluster25.tsv
    Creating file: 160k_All_cluster_findMarkers_upregulated_Cluster26.tsv
    Creating file: 160k_All_cluster_findMarkers_upregulated_Cluster27.tsv
    Creating file: 160k_All_cluster_findMarkers_upregulated_Cluster28.tsv
    Creating file: 160k_All_cluster_findMarkers_upregulated_Cluster29.tsv
    Creating file: 160k_All_cluster_findMarkers_upregulated_Cluster30.tsv
    Creating file: 160k_All_cluster_findMarkers_upregulated_Cluster31.tsv
    Creating file: 160k_All_cluster_findMarkers_upregulated_Cluster32.tsv
    Creating file: 160k_All_cluster_findMarkers_upregulated_Cluster33.tsv
    Creating file: 160k_All_cluster_findMarkers_upregulated_Cluster34.tsv
    Creating file: 160k_All_cluster_findMarkers_upregulated_Cluster35.tsv
    Creating file: 160k_All_cluster_findMarkers_upregulated_Cluster36.tsv
    Creating file: 160k_All_cluster_findMarkers_upregulated_Cluster37.tsv
    Creating file: 160k_All_cluster_findMarkers_upregulated_Cluster38.tsv
    Creating file: 160k_All_cluster_findMarkers_upregulated_Cluster39.tsv
    Creating file: 160k_All_cluster_findMarkers_upregulated_Cluster40.tsv
    Creating file: 160k_All_cluster_findMarkers_upregulated_Cluster41.tsv
    Creating file: 160k_All_cluster_findMarkers_upregulated_Cluster42.tsv
    Creating file: 160k_All_cluster_findMarkers_upregulated_Cluster43.tsv
    Creating file: 160k_All_cluster_findMarkers_upregulated_Cluster44.tsv
    Creating file: 160k_All_cluster_findMarkers_upregulated_Cluster45.tsv
    Creating file: 160k_All_cluster_findMarkers_upregulated_Cluster46.tsv


### Run `findMarkers` (downregulated genes)

Set `direction='down'` to only consider downregulated genes as potential markers.


```R
# Set pval.type (and min.prop if using "any")
pval.type <- "any"
min.prop <- 0.3
direction <- "down"

# Using logcounts from combined
marker.genes.cluster.dn <- findMarkers(sce_l, groups = combined$label, pval.type = pval.type, min.prop = min.prop,
                                       block = combined$condition, #blocking on uninteresting factors
                                       lfc = 0.5, direction = direction, BPPARAM = bpp)
marker.genes.cluster.dn
```


    List of length 46
    names(46): 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 ... 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46


Print the number of markers that passed the FDR or `Top` threshold. This will be the number of genes as inut for `enrichR`.


```R
printMarkerStats(marker.genes.cluster.dn, pval.type = pval.type, min.prop = min.prop)
```

    Number of selected markers (Top 200 genes of at least 30.0% comparisons):
    - Cluster1: 189; Up = 0; Down = 189; Max. P-value = 1.4e-25.
    - Cluster2: 196; Up = 0; Down = 196; Max. P-value = 1.2e-37.
    - Cluster3: 190; Up = 0; Down = 190; Max. P-value = 2.2e-78.
    - Cluster4: 156; Up = 0; Down = 156; Max. P-value = 9.3e-159.
    - Cluster5: 132; Up = 0; Down = 132; Max. P-value = 0.0039.
    - Cluster6: 163; Up = 0; Down = 163; Max. P-value = 1.2e-24.
    - Cluster7: 102; Up = 0; Down = 102; Max. P-value = 1.9e-11.
    - Cluster8: 123; Up = 0; Down = 123; Max. P-value = 3.4e-116.
    - Cluster9: 155; Up = 0; Down = 155; Max. P-value = 4.6e-14.
    - Cluster10: 132; Up = 0; Down = 132; Max. P-value = 0.
    - Cluster11: 150; Up = 0; Down = 150; Max. P-value = 2.6e-18.
    - Cluster12: 133; Up = 0; Down = 133; Max. P-value = 5.1e-85.
    - Cluster13: 137; Up = 0; Down = 137; Max. P-value = 1.4e-26.
    - Cluster14: 119; Up = 0; Down = 119; Max. P-value = 3e-27.
    - Cluster15: 116; Up = 0; Down = 116; Max. P-value = 5.4e-86.
    - Cluster16: 124; Up = 0; Down = 124; Max. P-value = 9.2e-91.
    - Cluster17: 137; Up = 0; Down = 137; Max. P-value = 4.2e-53.
    - Cluster18: 96; Up = 0; Down = 96; Max. P-value = 2.7e-90.
    - Cluster19: 148; Up = 0; Down = 148; Max. P-value = 0.
    - Cluster20: 128; Up = 0; Down = 128; Max. P-value = 0.00012.
    - Cluster21: 162; Up = 0; Down = 162; Max. P-value = 4.6e-43.
    - Cluster22: 155; Up = 0; Down = 155; Max. P-value = 4.4e-49.
    - Cluster23: 135; Up = 0; Down = 135; Max. P-value = 0.31.
    - Cluster24: 97; Up = 0; Down = 97; Max. P-value = 1.
    - Cluster25: 139; Up = 0; Down = 139; Max. P-value = 2.4e-97.
    - Cluster26: 102; Up = 0; Down = 102; Max. P-value = 0.26.
    - Cluster27: 141; Up = 0; Down = 141; Max. P-value = 5.5e-113.
    - Cluster28: 175; Up = 0; Down = 175; Max. P-value = 7.7e-72.
    - Cluster29: 161; Up = 0; Down = 161; Max. P-value = 4.8e-21.
    - Cluster30: 205; Up = 0; Down = 205; Max. P-value = 4e-106.
    - Cluster31: 203; Up = 0; Down = 203; Max. P-value = 0.
    - Cluster32: 202; Up = 0; Down = 202; Max. P-value = 1.5e-72.
    - Cluster33: 207; Up = 0; Down = 207; Max. P-value = 1e-63.
    - Cluster34: 207; Up = 0; Down = 207; Max. P-value = 6.4e-306.
    - Cluster35: 214; Up = 0; Down = 214; Max. P-value = 6.3e-89.
    - Cluster36: 203; Up = 0; Down = 203; Max. P-value = 6.9e-27.
    - Cluster37: 218; Up = 0; Down = 218; Max. P-value = 4.5e-113.
    - Cluster38: 201; Up = 0; Down = 201; Max. P-value = 1.2e-27.
    - Cluster39: 210; Up = 0; Down = 210; Max. P-value = 1.1e-14.
    - Cluster40: 183; Up = 0; Down = 183; Max. P-value = 0.
    - Cluster41: 158; Up = 0; Down = 158; Max. P-value = 4.2e-48.
    - Cluster42: 101; Up = 0; Down = 101; Max. P-value = 6.4e-40.
    - Cluster43: 78; Up = 0; Down = 78; Max. P-value = 5.5e-13.
    - Cluster44: 142; Up = 0; Down = 142; Max. P-value = 1.2e-118.
    - Cluster45: 150; Up = 0; Down = 150; Max. P-value = 4.7e-35.
    - Cluster46: 193; Up = 0; Down = 193; Max. P-value = 2.8e-41.
    * Upregulated when logFC > 0.0 and downregulated when logFC < 0.0.



```R
# Append Ensembl ID and Symbol from rowData
for(ID in names(marker.genes.cluster.dn)) {
    marker.genes.cluster.dn[[ID]] <- cbind(rowData(combined)[rownames(marker.genes.cluster.dn[[ID]]),][,1:2], 
                                                marker.genes.cluster.dn[[ID]])
}

exportResList(marker.genes.cluster.dn, col_anno = c("ID","Symbol"), prefix = paste0(file_id, "_cluster"), 
              direction = direction)
```

    Detecting findMarkers input.
    Creating file: 160k_All_cluster_findMarkers_downregulated_Cluster1.tsv
    Creating file: 160k_All_cluster_findMarkers_downregulated_Cluster2.tsv
    Creating file: 160k_All_cluster_findMarkers_downregulated_Cluster3.tsv
    Creating file: 160k_All_cluster_findMarkers_downregulated_Cluster4.tsv
    Creating file: 160k_All_cluster_findMarkers_downregulated_Cluster5.tsv
    Creating file: 160k_All_cluster_findMarkers_downregulated_Cluster6.tsv
    Creating file: 160k_All_cluster_findMarkers_downregulated_Cluster7.tsv
    Creating file: 160k_All_cluster_findMarkers_downregulated_Cluster8.tsv
    Creating file: 160k_All_cluster_findMarkers_downregulated_Cluster9.tsv
    Creating file: 160k_All_cluster_findMarkers_downregulated_Cluster10.tsv
    Creating file: 160k_All_cluster_findMarkers_downregulated_Cluster11.tsv
    Creating file: 160k_All_cluster_findMarkers_downregulated_Cluster12.tsv
    Creating file: 160k_All_cluster_findMarkers_downregulated_Cluster13.tsv
    Creating file: 160k_All_cluster_findMarkers_downregulated_Cluster14.tsv
    Creating file: 160k_All_cluster_findMarkers_downregulated_Cluster15.tsv
    Creating file: 160k_All_cluster_findMarkers_downregulated_Cluster16.tsv
    Creating file: 160k_All_cluster_findMarkers_downregulated_Cluster17.tsv
    Creating file: 160k_All_cluster_findMarkers_downregulated_Cluster18.tsv
    Creating file: 160k_All_cluster_findMarkers_downregulated_Cluster19.tsv
    Creating file: 160k_All_cluster_findMarkers_downregulated_Cluster20.tsv
    Creating file: 160k_All_cluster_findMarkers_downregulated_Cluster21.tsv
    Creating file: 160k_All_cluster_findMarkers_downregulated_Cluster22.tsv
    Creating file: 160k_All_cluster_findMarkers_downregulated_Cluster23.tsv
    Creating file: 160k_All_cluster_findMarkers_downregulated_Cluster24.tsv
    Creating file: 160k_All_cluster_findMarkers_downregulated_Cluster25.tsv
    Creating file: 160k_All_cluster_findMarkers_downregulated_Cluster26.tsv
    Creating file: 160k_All_cluster_findMarkers_downregulated_Cluster27.tsv
    Creating file: 160k_All_cluster_findMarkers_downregulated_Cluster28.tsv
    Creating file: 160k_All_cluster_findMarkers_downregulated_Cluster29.tsv
    Creating file: 160k_All_cluster_findMarkers_downregulated_Cluster30.tsv
    Creating file: 160k_All_cluster_findMarkers_downregulated_Cluster31.tsv
    Creating file: 160k_All_cluster_findMarkers_downregulated_Cluster32.tsv
    Creating file: 160k_All_cluster_findMarkers_downregulated_Cluster33.tsv
    Creating file: 160k_All_cluster_findMarkers_downregulated_Cluster34.tsv
    Creating file: 160k_All_cluster_findMarkers_downregulated_Cluster35.tsv
    Creating file: 160k_All_cluster_findMarkers_downregulated_Cluster36.tsv
    Creating file: 160k_All_cluster_findMarkers_downregulated_Cluster37.tsv
    Creating file: 160k_All_cluster_findMarkers_downregulated_Cluster38.tsv
    Creating file: 160k_All_cluster_findMarkers_downregulated_Cluster39.tsv
    Creating file: 160k_All_cluster_findMarkers_downregulated_Cluster40.tsv
    Creating file: 160k_All_cluster_findMarkers_downregulated_Cluster41.tsv
    Creating file: 160k_All_cluster_findMarkers_downregulated_Cluster42.tsv
    Creating file: 160k_All_cluster_findMarkers_downregulated_Cluster43.tsv
    Creating file: 160k_All_cluster_findMarkers_downregulated_Cluster44.tsv
    Creating file: 160k_All_cluster_findMarkers_downregulated_Cluster45.tsv
    Creating file: 160k_All_cluster_findMarkers_downregulated_Cluster46.tsv


## Save `findMarkers` results to `metadata`

<div class="alert alert-info">
    <strong>Tip!</strong> In the accompanied Shiny App, it will look for list(s) named with "findMarkers_" in the prefix in <code>metadata()</code> and show their content under the <u>Gene markers</u> section of the website. In the example below, the content of the 3 lists will be displayed in their respective sub-menus under the following titles: 
    <ul>
        <li>Cluster</li>
        <li>Cluster: up</li>
        <li>Cluster: dn</li>
    </ul>
</div>


```R
metadata(combined)[['findMarkers_Cluster']] <- marker.genes.cluster
metadata(combined)[['findMarkers_Cluster_up']] <- marker.genes.cluster.up
metadata(combined)[['findMarkers_Cluster_dn']] <- marker.genes.cluster.dn
```

## Use cluster marker genes to show cluster similarities


```R
nGene <- 200
geneNames <- sapply(marker.genes.cluster, function(x) rownames(x[1:nGene,]))

geneNames <- unique(as.character(geneNames)) # Remove duplicated genes
print(paste("Number of genes to plot:", length(geneNames)))
```

    [1] "Number of genes to plot: 1694"



```R
fig(width = 16, height = 10)
plotGroupedHeatmap(combined, features = geneNames, group = "CellType_2", clustering_method = "ward.D2", 
                   border_color = "black", color = c_heatmap_col2, fontsize = 16, angle_col = 90,
                   center = TRUE, scale = TRUE, zlim = c(-3, 3), main = "Row-scaled", show_rownames = FALSE)
reset.fig()
```


    
![png](Integrated_files/Integrated_197_0.png)
    


## Visualise first N cluster marker genes

Use `findMarkers` result from both directions. We aimed to present between 50 to 100 genes in the heatmap.


```R
nGene <- 4
geneNames <- sapply(marker.genes.cluster, function(x) rownames(x[1:nGene,]))
t(geneNames) # A matrix

geneNames <- unique(as.character(geneNames)) # Remove duplicated genes
print(paste("Number of genes to plot:", length(geneNames)))
```


<table class="dataframe">
<caption>A matrix: 46 × 4 of type chr</caption>
<tbody>
	<tr><th scope=row>1</th><td>CCL5 </td><td>NELL2 </td><td>NKG7  </td><td>EFHD2    </td></tr>
	<tr><th scope=row>2</th><td>NELL2</td><td>CCL5  </td><td>EFHD2 </td><td>NKG7     </td></tr>
	<tr><th scope=row>3</th><td>CCR7 </td><td>MAL   </td><td>LEF1  </td><td>CCL5     </td></tr>
	<tr><th scope=row>4</th><td>CCR4 </td><td>MAL   </td><td>CD4   </td><td>CCR7     </td></tr>
	<tr><th scope=row>5</th><td>CCR4 </td><td>MAL   </td><td>CCL5  </td><td>NKG7     </td></tr>
	<tr><th scope=row>6</th><td>MAL  </td><td>CCR6  </td><td>NELL2 </td><td>CR1      </td></tr>
	<tr><th scope=row>7</th><td>CD4  </td><td>CCL5  </td><td>NKG7  </td><td>IL7R     </td></tr>
	<tr><th scope=row>8</th><td>NKG7 </td><td>CCL5  </td><td>EFHD2 </td><td>IL7R     </td></tr>
	<tr><th scope=row>9</th><td>CCL5 </td><td>NKG7  </td><td>GNLY  </td><td>EFHD2    </td></tr>
	<tr><th scope=row>10</th><td>GNLY </td><td>FGFBP2</td><td>NKG7  </td><td>CD4      </td></tr>
	<tr><th scope=row>11</th><td>GNLY </td><td>TRAC  </td><td>NKG7  </td><td>IL32     </td></tr>
	<tr><th scope=row>12</th><td>IFIT2</td><td>OASL  </td><td>GNLY  </td><td>HEATR9   </td></tr>
	<tr><th scope=row>13</th><td>GNLY </td><td>NKG7  </td><td>EFHD2 </td><td>IL7R     </td></tr>
	<tr><th scope=row>14</th><td>CD8A </td><td>CCL5  </td><td>NKG7  </td><td>NELL2    </td></tr>
	<tr><th scope=row>15</th><td>CCL5 </td><td>NKG7  </td><td>GNLY  </td><td>SLC7A5   </td></tr>
	<tr><th scope=row>16</th><td>CCL5 </td><td>PMAIP1</td><td>NFKBIZ</td><td>IFIT2    </td></tr>
	<tr><th scope=row>17</th><td>NKG7 </td><td>CCL5  </td><td>LTB   </td><td>GNLY     </td></tr>
	<tr><th scope=row>18</th><td>GNLY </td><td>CTSW  </td><td>CCL5  </td><td>NKG7     </td></tr>
	<tr><th scope=row>19</th><td>GNLY </td><td>CTSW  </td><td>GZMH  </td><td>ADGRG1   </td></tr>
	<tr><th scope=row>20</th><td>IFIT2</td><td>PMAIP1</td><td>GNLY  </td><td>RGS3     </td></tr>
	<tr><th scope=row>21</th><td>NKG7 </td><td>GNLY  </td><td>CST7  </td><td>IL7R     </td></tr>
	<tr><th scope=row>22</th><td>KLRC3</td><td>TRDC  </td><td>GNLY  </td><td>NKG7     </td></tr>
	<tr><th scope=row>23</th><td>CCL5 </td><td>GZMH  </td><td>NKG7  </td><td>IL32     </td></tr>
	<tr><th scope=row>24</th><td>CCL5 </td><td>KLRG1 </td><td>GNLY  </td><td>NKG7     </td></tr>
	<tr><th scope=row>25</th><td>MYBL1</td><td>PDZD4 </td><td>PFN1  </td><td>COTL1    </td></tr>
	<tr><th scope=row>26</th><td>DUSP2</td><td>KLRB1 </td><td>SLC7A5</td><td>CCL5     </td></tr>
	<tr><th scope=row>27</th><td>GNLY </td><td>IL2RB </td><td>NKG7  </td><td>CTSW     </td></tr>
	<tr><th scope=row>28</th><td>CD5  </td><td>CD3E  </td><td>IFIT2 </td><td>TRAC     </td></tr>
	<tr><th scope=row>29</th><td>CD5  </td><td>CD6   </td><td>PYHIN1</td><td>TBX21    </td></tr>
	<tr><th scope=row>30</th><td>CD74 </td><td>CD79A </td><td>IGHD  </td><td>IGHM     </td></tr>
	<tr><th scope=row>31</th><td>IGHM </td><td>IGHD  </td><td>CD79A </td><td>NIBAN3   </td></tr>
	<tr><th scope=row>32</th><td>CD74 </td><td>CD3E  </td><td>CCL5  </td><td>CD79A    </td></tr>
	<tr><th scope=row>33</th><td>CD74 </td><td>BLK   </td><td>CD79A </td><td>CD3E     </td></tr>
	<tr><th scope=row>34</th><td>CD3E </td><td>ZAP70 </td><td>CSF3R </td><td>VCAN     </td></tr>
	<tr><th scope=row>35</th><td>ZAP70</td><td>IL32  </td><td>CST7  </td><td>CD3E     </td></tr>
	<tr><th scope=row>36</th><td>CD3E </td><td>IL32  </td><td>CST7  </td><td>CSF1R    </td></tr>
	<tr><th scope=row>37</th><td>SPI1 </td><td>LYZ   </td><td>IFI30 </td><td>CST3     </td></tr>
	<tr><th scope=row>38</th><td>CD3E </td><td>ZAP70 </td><td>BCL11B</td><td>ETS1     </td></tr>
	<tr><th scope=row>39</th><td>CCL5 </td><td>CST7  </td><td>SAMD3 </td><td>TBX21    </td></tr>
	<tr><th scope=row>40</th><td>B2M  </td><td>NKG7  </td><td>UBA52 </td><td>GNLY     </td></tr>
	<tr><th scope=row>41</th><td>B2M  </td><td>LTB   </td><td>EEF1G </td><td>GNLY     </td></tr>
	<tr><th scope=row>42</th><td>CCL5 </td><td>NKG7  </td><td>GNLY  </td><td>EFHD2    </td></tr>
	<tr><th scope=row>43</th><td>NKG7 </td><td>GNLY  </td><td>HDLBP </td><td>TRAK1    </td></tr>
	<tr><th scope=row>44</th><td>CD74 </td><td>CD79A </td><td>MS4A1 </td><td>TNFRSF13C</td></tr>
	<tr><th scope=row>45</th><td>SPI1 </td><td>FTH1  </td><td>CLEC7A</td><td>IFI30    </td></tr>
	<tr><th scope=row>46</th><td>CST7 </td><td>CD3E  </td><td>SAMD3 </td><td>SPI1     </td></tr>
</tbody>
</table>



    [1] "Number of genes to plot: 70"


Alternatively, use genes identified from up- and downregulated `findMarkers` results.


```R
#nGene <- 5
#geneNames1 <- sapply(marker.genes.cluster.up, function(x) rownames(x[1:nGene,]))
#geneNames1 # A matrix

#geneNames2 <- sapply(marker.genes.cluster.dn, function(x) rownames(x[1:nGene,]))
#geneNames2 # A matrix

#geneNames <- unique(c(as.character(geneNames1), as.character(geneNames2))) # Remove duplicated genes
#print(paste("Number of genes to plot:", length(geneNames)))
```


```R
p1 <- plotGroupedHeatmap(combined, features = geneNames, group = "label", clustering_method = "ward.D2", 
                         border_color = "black", color = c_heatmap_col1, fontsize = 14, angle_col = 90, 
                         main = "Unscaled", silent = T)

p2 <- plotGroupedHeatmap(combined, features = geneNames, group = "label", clustering_method = "ward.D2", 
                         border_color = "black", color = c_heatmap_col2, fontsize = 14, angle_col = 90,
                         center = TRUE, scale = TRUE, zlim = c(-3, 3), main = "Row-scaled", silent = T)

fig(width = 16, height = 16)
plot(p1$gtable)
plot(p2$gtable)
reset.fig()
```


    
![png](Integrated_files/Integrated_202_0.png)
    



    
![png](Integrated_files/Integrated_202_1.png)
    



```R
fig(width = 16, height = 18)
plotDots(combined, features = geneNames[p2$tree_row$order], group = "label", zlim = c(-3, 3), 
         center = TRUE, scale = TRUE) + scale_size(limits = c(0, 1), range = c(0.1, 6)) + 
    scale_y_discrete(limits = geneNames[p2$tree_row$order]) + # order genes based on heatmap p2 above
    scale_x_discrete(limits = p2$tree_col$labels[p2$tree_col$order]) + # order clusters based on heatmap p2 above
    guides(colour = guide_colourbar(title = "Row (Gene) Z-Score", barwidth = 10), 
           size = guide_legend(title = "Proportion Detected")) + theme_cowplot(16) + #coord_flip() +
    theme(panel.grid.major = element_line(colour = "gray90"), 
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
          legend.position = "top", legend.justification = "center", legend.title.position = "top") + 
    labs(x = "Cluster", y = "Genes")
reset.fig()
```

    [1m[22mScale for [32msize[39m is already present.
    Adding another scale for [32msize[39m, which will replace the existing scale.



    
![png](Integrated_files/Integrated_203_1.png)
    


## Visualise first upregulated marker genes from each cluster


```R
geneNames <- sapply(marker.genes.cluster.up, function(x) rownames(x[1,]))
#geneNames

df <- data.frame(cluster = names(geneNames), gene = geneNames) %>% group_by(gene) %>% 
    summarise(cluster = paste(cluster, collapse = ",")) %>% 
    mutate(order = as.numeric(str_replace(cluster, ",.*", ""))) %>% arrange(order)
df

print(paste("Number of genes to plot:", nrow(df)))
```


<table class="dataframe">
<caption>A tibble: 21 × 3</caption>
<thead>
	<tr><th scope=col>gene</th><th scope=col>cluster</th><th scope=col>order</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td>NELL2   </td><td>1,2                 </td><td> 1</td></tr>
	<tr><td>TCF7    </td><td>3                   </td><td> 3</td></tr>
	<tr><td>MAL     </td><td>4,6                 </td><td> 4</td></tr>
	<tr><td>CCR4    </td><td>5                   </td><td> 5</td></tr>
	<tr><td>IL7R    </td><td>7,13                </td><td> 7</td></tr>
	<tr><td>LTB     </td><td>8                   </td><td> 8</td></tr>
	<tr><td>TRAC    </td><td>9                   </td><td> 9</td></tr>
	<tr><td>FGFBP2  </td><td>10                  </td><td>10</td></tr>
	<tr><td>CCL5    </td><td>11,15,16,20,23      </td><td>11</td></tr>
	<tr><td>GNLY    </td><td>12,18,19,22,27,41,42</td><td>12</td></tr>
	<tr><td>CD8A    </td><td>14                  </td><td>14</td></tr>
	<tr><td>NKG7    </td><td>17,21,24,28,43      </td><td>17</td></tr>
	<tr><td>PFN1    </td><td>25                  </td><td>25</td></tr>
	<tr><td>DUSP2   </td><td>26                  </td><td>26</td></tr>
	<tr><td>KLRB1   </td><td>29                  </td><td>29</td></tr>
	<tr><td>CD74    </td><td>30,32,33,38,44,46   </td><td>30</td></tr>
	<tr><td>IGHM    </td><td>31                  </td><td>31</td></tr>
	<tr><td>CSF3R   </td><td>34                  </td><td>34</td></tr>
	<tr><td>SPI1    </td><td>35,36,37,45         </td><td>35</td></tr>
	<tr><td>ITM2C   </td><td>39                  </td><td>39</td></tr>
	<tr><td>TNFRSF25</td><td>40                  </td><td>40</td></tr>
</tbody>
</table>



    [1] "Number of genes to plot: 21"



```R
p1 <- plotGroupedHeatmap(combined, features = df$gene, group = "label", clustering_method = "ward.D2", 
                         border_color = "black", color = c_heatmap_col1, fontsize = 14, angle_col = 90, 
                         main = "Unscaled", silent = T)

p2 <- plotGroupedHeatmap(combined, features = df$gene, group = "label", clustering_method = "ward.D2", 
                         border_color = "black", color = c_heatmap_col2, fontsize = 14, angle_col = 90,
                         center = TRUE, scale = TRUE, zlim = c(-3, 3), main = "Row-scaled", silent = T)

fig(width = 16, height = 7)
plot(p1$gtable)
plot(p2$gtable)
reset.fig()
```


    
![png](Integrated_files/Integrated_206_0.png)
    



    
![png](Integrated_files/Integrated_206_1.png)
    



```R
fig(width = 16, height = 7)
plotDots(combined, features = df$gene[p2$tree_row$order], group = "label", zlim = c(-3, 3), 
         center = TRUE, scale = TRUE) + scale_size(limits = c(0, 1), range = c(0.1, 6)) + 
    scale_y_discrete(limits = df$gene[p2$tree_row$order]) + # order genes based on heatmap p2 above
    scale_x_discrete(limits = p2$tree_col$labels[p2$tree_col$order]) + # order clusters based on heatmap p2 above
    guides(colour = guide_colourbar(title = "Row (Gene) Z-Score", barwidth = 10), 
           size = guide_legend(title = "Proportion Detected")) + theme_cowplot(16) + #coord_flip() +
    theme(panel.grid.major = element_line(colour = "gray90"), 
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
          legend.position = "top", legend.justification = "center", legend.title.position = "top") + 
    labs(x = "Cluster", y = "Genes")
reset.fig()
```

    [1m[22mScale for [32msize[39m is already present.
    Adding another scale for [32msize[39m, which will replace the existing scale.



    
![png](Integrated_files/Integrated_207_1.png)
    



```R
dimname <- "MNN-TSNE"

fig(width = 16, height = 15)
as.data.frame(reducedDim(combined, dimname)) %>% `colnames<-`(c("V1","V2")) %>%
    bind_cols(as_tibble(t(logcounts(combined[df$gene,])))) %>%
    gather(., key = "Symbol", value = "Expression", -c(V1, V2) ) %>% 
    mutate_at(vars(Symbol), factor) %>% mutate(Symbol = factor(Symbol, levels = df$gene)) %>%
    ggplot(aes(x = V1, y = V2, color = Expression)) + geom_point(size = 0.3, alpha = 0.3) + 
    facet_wrap(~ Symbol, ncol = 5, 
               labeller = as_labeller(function(x) paste0(df$cluster,": ", x))) + # add cluster id
    scale_color_viridis(option = "plasma", direction = -1) +
    theme_classic(base_size = 20) + labs(x = paste(dimname, "1"), y = paste(dimname, "2"))
reset.fig()
```


    
![png](Integrated_files/Integrated_208_0.png)
    


# 9. DE analysis between conditions

## Prepare input

**Excluded celltypes and/or clusters**

*Assuming excluding 'Other T', 'ILC', 'DC', 'Unknown' and 'Doublets'*

- Cluster 25 DP T
- Cluster 26 Activated DN T
- Cluster 29 ILC2
- Cluster 38 CD1C+ DCs
- Cluster 39 pDCs
- Cluster 40 Cytotoxic-like
- Cluster 41 Momory-like
- Clusters 42-46 Doublets


```R
table(combined$CellType_1, combined$Sample)
```


               
                Control1 KidneyCancer LungCancer
      CD8+ T        3507         2509       2041
      CD4+ T        5826        10231       6039
      gdT             23          918        205
      Other T        323           89        404
      NK            1997          280       1254
      ILC             28           16          8
      B             1563          631       1997
      Monocytes     3107          609        802
      DCs            366           48         15
      Unknown        163          584        476
      Doublets      1135          861        439



```R
table("Cells kept" = !combined$CellType_1 %in% c("Other T","ILC","DC","Unknown","Doublets"))

kept <- combined[,!combined$CellType_1 %in% c("Other T","ILC","DC","Unknown","Doublets")]
colData(kept) <- droplevels(colData(kept))

# remove whitespaces and '+' sign from CellType_1, which we'll use later
levels(kept$CellType_1) <- gsub("\\s*|\\+","", levels(kept$CellType_1))

# Remove genes not expressed at all
is.exp <- rowSums(counts(kept) > 0) > 1
table("Is expressed" = is.exp)

kept <- kept[is.exp,]
kept
```


    Cells kept
    FALSE  TRUE 
     4526 43968 



    Is expressed
    FALSE  TRUE 
     2201 15928 



    class: SingleCellExperiment 
    dim: 15928 43968 
    metadata(27): Control1_Samples Control1_cyclone ... findMarkers_Cluster_up findMarkers_Cluster_dn
    assays(2): counts logcounts
    rownames(15928): SAMD11 NOC2L ... MT-ND6 MT-CYB
    rowData names(5): ID Symbol Type SEQNAME is_mito
    colnames(43968): Control1_AAACAAGCAACTAGTGACTTTAGG-1 Control1_AAACAAGCAGTTATCCACTTTAGG-1 ...
      LungCancer_TTTGTGAGTTGAGTCTAGCTGTGA-1 LungCancer_TTTGTGAGTTTACGACAGCTGTGA-1
    colData names(20): Sample Barcode ... ClusterCellType CellType_2
    reducedDimNames(6): PCA TSNE ... MNN-TSNE MNN-UMAP
    mainExpName: Gene Expression
    altExpNames(1): Antibody Capture


Create a `dgCMatrix` counts matrix for faster computation in some steps.


```R
kept_c <- as(counts(kept, withDimnames = FALSE), "dgCMatrix")
rownames(kept_c) <- rownames(kept)
```

## Creating pseudo-bulk samples

We sum counts together from cells with the same combination of selected features. See [OSCA reference](https://bioconductor.org/books/3.20/OSCA.multisample/multi-sample-comparisons.html)

<div class="alert alert-warning">
    <strong>Warning!</strong> Change the <code>"ids" DataFrame</code> to include the factors where the unique combination of levels is used to define a group.
</div>


```R
summed <- aggregateAcrossCells(kept, id = colData(kept)[,c("Sample","condition","CellType_1")], BPPARAM = bpp)
colData(summed) <- droplevels(colData(summed))
sizeFactors(summed) <- NULL
summed <- logNormCounts(summed) # Add logcounts
summed
```


    class: SingleCellExperiment 
    dim: 15928 21 
    metadata(27): Control1_Samples Control1_cyclone ... findMarkers_Cluster_up findMarkers_Cluster_dn
    assays(2): counts logcounts
    rownames(15928): SAMD11 NOC2L ... MT-ND6 MT-CYB
    rowData names(5): ID Symbol Type SEQNAME is_mito
    colnames: NULL
    colData names(24): Sample Barcode ... ncells sizeFactor
    reducedDimNames(6): PCA TSNE ... MNN-TSNE MNN-UMAP
    mainExpName: Gene Expression
    altExpNames(0):


### Add column names to `summed`

<div class="alert alert-warning">
    <strong>Warning!</strong> Very important to define the column name of the <code>summed</code> object, i.e. a unique ID for each pseudo-bulk sample in <code>colData</code>. This is so that the <code>DGEList</code> object in the following edgeR analysis is created with the approprate row names for the <code>samples</code> slot.
</div>


```R
summed$PseudoSample <- paste(summed$Sample, summed$CellType_1, sep = "_")
summed$PseudoSample <- as.factor(summed$PseudoSample)
summed$PseudoSample <- factor(summed$PseudoSample, levels = gtools::mixedsort(levels(summed$PseudoSample)))
colnames(summed) <- summed$PseudoSample
```


```R
fig(width = 16, height = 8)
as.data.frame(colData(summed)[,c("PseudoSample","condition","ncells")]) %>% 
    ggplot(aes(PseudoSample, ncells, fill = condition)) + geom_col() + 
    geom_text(aes(label = ncells), hjust = -0.1, size = 4) + coord_flip() + 
    scale_y_continuous(limits = c(0, 11000), expand = c(0, 0)) + scale_fill_manual(values = c_cond_col) + 
    theme_cowplot(16) + theme(legend.position = "top") + ylab("Number of cells")
reset.fig()
```


    
![png](Integrated_files/Integrated_218_0.png)
    


## Perform DE analysis between conditions

[DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) or [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html): For data with a high level of multimodality (data heterogeneity) methods that consider the behavior of each individual gene, such as DESeq2 show better true positive rates (sensitivities). If the level of multimodality is low, however, edgeR can provide higher precision. [Reference](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-2599-6)

Here we use edgeR to perform DE analysis.

### Create DGEList data object

<div class="alert alert-warning">
    <strong>Warning!</strong> Change the information to be added to the <code>samples</code> slot.
</div>


```R
y <- DGEList(counts(summed), samples = colData(summed)[,c("PseudoSample","ncells","condition","CellType_1")], 
             genes = rowData(summed)[,c("ID","Symbol")])

# Removing sample that have very few or lowly-sequenced cells
# In this case, we remove combinations containing fewer than 10 cells
discarded <- y$samples$ncells < 10
y <- y[,!discarded]
summary(discarded)
```


       Mode   FALSE 
    logical      21 


### Set design matrix

Read more about design matrices with and without intercept term [here](https://bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/designmatrices.html#design-matrices-with-and-without-intercept-term), and ordering of terms in the `model.matrix` when using `0+` in the model formula [here](https://support.bioconductor.org/p/110327/#110379).


```R
my.design <- model.matrix(~ 0 + y$samples$condition + y$samples$CellType_1, y$samples)

colnames(my.design) <- gsub("y\\$samples\\$", "", colnames(my.design)) # simplify term name
data.frame(coef = 1:ncol(my.design), term = colnames(my.design))
```


<table class="dataframe">
<caption>A data.frame: 8 × 2</caption>
<thead>
	<tr><th scope=col>coef</th><th scope=col>term</th></tr>
	<tr><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><td>1</td><td>conditionControl   </td></tr>
	<tr><td>2</td><td>conditionCancer    </td></tr>
	<tr><td>3</td><td>CellType_1CD4T     </td></tr>
	<tr><td>4</td><td>CellType_1gdT      </td></tr>
	<tr><td>5</td><td>CellType_1NK       </td></tr>
	<tr><td>6</td><td>CellType_1B        </td></tr>
	<tr><td>7</td><td>CellType_1Monocytes</td></tr>
	<tr><td>8</td><td>CellType_1DCs      </td></tr>
</tbody>
</table>



### Construct contrast vector using `makeContrasts` to make comparison easier


```R
my.contrasts <- makeContrasts(
    Cancer_Control = conditionCancer-conditionControl,
    levels = my.design
)
```

### Remove lowly expressed genes


```R
# Use the filterByExpr() function from edgeR
keep <- filterByExpr(y, design = my.design, group = NULL)
y <- y[keep, , keep.lib.sizes = FALSE]
summary(keep)

# Correct for composition biases by computing normalization factors
y <- calcNormFactors(y)
DataFrame(y$samples)
```


       Mode   FALSE    TRUE 
    logical    3299   12629 



    DataFrame with 21 rows and 7 columns
                            group  lib.size norm.factors         PseudoSample    ncells condition CellType_1
                         <factor> <numeric>    <numeric>             <factor> <integer>  <factor>   <factor>
    Control1_CD8T               1  15189582      1.22579        Control1_CD8T      3507   Control       CD8T
    Control1_CD4T               1  24092369      1.21040        Control1_CD4T      5826   Control       CD4T
    Control1_gdT                1     91934      1.33602        Control1_gdT         23   Control       gdT 
    Control1_NK                 1   9944841      1.12732        Control1_NK        1997   Control       NK  
    Control1_B                  1   6331726      1.14387        Control1_B         1563   Control       B   
    ...                       ...       ...          ...                  ...       ...       ...        ...
    LungCancer_gdT              1   1217279     1.113303 LungCancer_gdT             205    Cancer  gdT      
    LungCancer_NK               1   7880633     0.992904 LungCancer_NK             1254    Cancer  NK       
    LungCancer_B                1  10328638     0.890531 LungCancer_B              1997    Cancer  B        
    LungCancer_Monocytes        1   6996744     0.641612 LungCancer_Monocytes       802    Cancer  Monocytes
    LungCancer_DCs              1    168084     1.110995 LungCancer_DCs              15    Cancer  DCs      


### Create multi-dimensional scaling (MDS) plot


```R
fig(width = 8, height = 8)
limma::plotMDS(cpm(y, log = TRUE), col = str_replace_all(y$samples$condition, c_cond_col),
        main = "PCoA", cex.main = 2, cex.lab = 1.5, cex.axis = 1.5, cex = 1)
reset.fig()
```


    
![png](Integrated_files/Integrated_228_0.png)
    


### Estimating the dispersions


```R
y <- estimateDisp(y, my.design)
summary(y$trended.dispersion)

# Visualise dispersion estimates
fig(width = 7, height = 7)
plotBCV(y, cex.lab = 1.5, cex.axis = 1.5, cex = 0.5, main = "Show common, trended and genewise BCV estimates")
reset.fig()
```


       Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    0.08356 0.08728 0.10067 0.14687 0.18655 0.34303 



    
![png](Integrated_files/Integrated_230_1.png)
    


### Estimating the quasi-likelihood (QL) dispersions


```R
fit <- glmQLFit(y, my.design, robust = TRUE)

# Visualise QL dispersion estimates
fig(width = 7, height = 7)
plotQLDisp(fit, cex.lab = 1.5, cex.axis = 1.5, cex = 0.5, 
           main = "Show common, trended and genewise QL dispersion estimates")
reset.fig()
```


    
![png](Integrated_files/Integrated_232_0.png)
    


### Performs QL F-test to determine DE genes

#### DE: Cancer vs. Control


```R
res <- glmQLFTest(fit, contrast = my.contrasts[,"Cancer_Control"])

# Print the total number of differentially expressed genes at 5% FDR
summary(decideTests(res))

# Print top 10 DE genes
topTags(res)

# Output DE results to text file
outfile <- paste0(file_id, "_edgeR_DE_results_Cancer_vs_Ctrl.tsv")
print(paste("Write to file:", outfile))
write.table(topTags(res, n = nrow(res), sort.by = "PValue", adjust.method = "BH"), file = outfile, 
            sep = "\t", quote = F, row.names = F, col.names = T)
```


           -1*conditionControl 1*conditionCancer
    Down                                    3707
    NotSig                                  5414
    Up                                      3508



<dl>
	<dt>$table</dt>
		<dd><table class="dataframe">
<caption>A data.frame: 10 × 7</caption>
<thead>
	<tr><th></th><th scope=col>ID</th><th scope=col>Symbol</th><th scope=col>logFC</th><th scope=col>logCPM</th><th scope=col>F</th><th scope=col>PValue</th><th scope=col>FDR</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>ETF1</th><td>ENSG00000120705</td><td>ETF1     </td><td> 2.541868</td><td>6.252779</td><td>393.2064</td><td>2.811069e-14</td><td>2.762753e-10</td></tr>
	<tr><th scope=row>MRM1</th><td>ENSG00000278619</td><td>MRM1     </td><td>-3.221633</td><td>2.967170</td><td>297.1798</td><td>8.377542e-14</td><td>2.762753e-10</td></tr>
	<tr><th scope=row>TNFAIP8L2</th><td>ENSG00000163154</td><td>TNFAIP8L2</td><td>-3.753401</td><td>4.225871</td><td>276.5543</td><td>8.654676e-14</td><td>2.762753e-10</td></tr>
	<tr><th scope=row>ZNF644</th><td>ENSG00000122482</td><td>ZNF644   </td><td> 2.103136</td><td>9.074763</td><td>391.4011</td><td>9.598485e-14</td><td>2.762753e-10</td></tr>
	<tr><th scope=row>ATP5PO</th><td>ENSG00000241837</td><td>ATP5PO   </td><td> 7.191873</td><td>4.490992</td><td>292.1230</td><td>1.093813e-13</td><td>2.762753e-10</td></tr>
	<tr><th scope=row>YTHDF3</th><td>ENSG00000185728</td><td>YTHDF3   </td><td> 2.313850</td><td>7.910533</td><td>367.1672</td><td>1.691840e-13</td><td>3.561041e-10</td></tr>
	<tr><th scope=row>BTG3</th><td>ENSG00000154640</td><td>BTG3     </td><td> 3.289743</td><td>6.314007</td><td>316.3214</td><td>2.145984e-13</td><td>3.871662e-10</td></tr>
	<tr><th scope=row>GBP1</th><td>ENSG00000117228</td><td>GBP1     </td><td>-3.875099</td><td>6.764013</td><td>334.6406</td><td>2.710475e-13</td><td>4.278823e-10</td></tr>
	<tr><th scope=row>RLF</th><td>ENSG00000117000</td><td>RLF      </td><td> 2.485302</td><td>6.924447</td><td>335.4624</td><td>3.100979e-13</td><td>4.351363e-10</td></tr>
	<tr><th scope=row>GIMAP6</th><td>ENSG00000133561</td><td>GIMAP6   </td><td>-3.317655</td><td>6.081295</td><td>227.1101</td><td>6.545696e-13</td><td>8.266560e-10</td></tr>
</tbody>
</table>
</dd>
	<dt>$adjust.method</dt>
		<dd>'BH'</dd>
	<dt>$comparison</dt>
		<dd>'-1*conditionControl 1*conditionCancer'</dd>
	<dt>$test</dt>
		<dd>'glm'</dd>
</dl>



    [1] "Write to file: 160k_All_edgeR_DE_results_Cancer_vs_Ctrl.tsv"


### Save DE results to `metadata`

<div class="alert alert-info">
    <strong>Tip!</strong> In the accompanied Shiny App, it will look for list(s) named with "edgeR_" in the prefix in <code>metadata()</code> and show their content under the <u>DEA (edgeR)</u> section of the website. In the example below, the content of the 3 <code>TopTags</code> objects will be displayed in the <u>Condition</u> sub-menu.
</div>


```R
metadata(combined)[['edgeR_Cancer_Control']] <- list("Cancer_Control" = res)
```

## Visualise first N condition DE genes

#### Vis: Cancer vs. Control

Assuming we are interested in looking at cells grouped by both Sample ID and Coarse Cell Type. Let's define these groupings here.


```R
kept$Group <- as.factor(paste0(kept$Sample, kept$CellType_1))
kept$Group <- factor(kept$Group, levels = as.vector(outer(levels(kept$Sample), levels(kept$CellType_1), paste0)))
#levels(kept$Group)

summed$Group <- as.factor(paste0(summed$Sample, summed$CellType_1))
summed$Group <- factor(summed$Group, levels = levels(kept$Group))
#levels(summed$Group)
```


```R
nGene <- 50
geneNames <- rownames(topTags(res, nGene))
geneNames
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'ETF1'</li><li>'MRM1'</li><li>'TNFAIP8L2'</li><li>'ZNF644'</li><li>'ATP5PO'</li><li>'YTHDF3'</li><li>'BTG3'</li><li>'GBP1'</li><li>'RLF'</li><li>'GIMAP6'</li><li>'GNA13'</li><li>'RPRD1B'</li><li>'PYCR2'</li><li>'CDK17'</li><li>'SUCO'</li><li>'PSMB10'</li><li>'NAP1L2'</li><li>'LRRC61'</li><li>'CENPC'</li><li>'MED11'</li><li>'GIMAP4'</li><li>'DNAJC2'</li><li>'SOCS5'</li><li>'GIMAP5'</li><li>'UBXN2A'</li><li>'METTL13'</li><li>'RSL24D1'</li><li>'ARL11'</li><li>'ELL2'</li><li>'DCP1A'</li><li>'GIMAP8'</li><li>'MAPK6'</li><li>'MYOF'</li><li>'RNF139'</li><li>'RAB1A'</li><li>'RAB2A'</li><li>'CSF1R'</li><li>'PARS2'</li><li>'ADNP2'</li><li>'PLA2G4C'</li><li>'SNIP1'</li><li>'ZNF14'</li><li>'HNRNPK'</li><li>'ZNF639'</li><li>'THNSL2'</li><li>'WDR48'</li><li>'AARS1'</li><li>'UBQLN1'</li><li>'MFSD14A'</li><li>'ZNF331'</li></ol>



Using `summed` pseudo-bulk samples.


```R
logexp <- logcounts(summed)[geneNames,] # logcounts

expr <- data.frame(Group = summed$Group, t(logexp)) %>% group_by(Group) %>% data.frame %>% 
    column_to_rownames("Group") %>% t()

ann_col <- data.frame(Sample = summed$Sample, Condition = summed$condition, row.names = colnames(expr))
condition_colors <- list(Sample = c_sample_col, Condition = c_cond_col)

p <- pheatmap(expr, scale = "row", clustering_method = "ward.D2", cluster_cols = TRUE, 
              border_color = "black", fontsize = 12, angle_col = 45, color = c_heatmap_col2, breaks = breaks,
              annotation_col = ann_col, annotation_colors = condition_colors, 
              main = "Row-scaled 'summed' logcounts", silent = T)

fig(width = 14, height = 15)
plot(p$gtable)
reset.fig()
```


    
![png](Integrated_files/Integrated_241_0.png)
    


Using `kept` filtered single cells.


```R
fig(width = 12, height = 16)
plotDots(kept, features = geneNames[p$tree_row$order], group = "Group", 
         center = TRUE, scale = TRUE, zlim = c(-3, 3)) + scale_size(limits = c(0, 1), range = c(0.1, 6)) + 
    scale_x_discrete(limits = p$tree_col$labels) + # order groups based on heatmap p above
    scale_y_discrete(limits = geneNames[p$tree_row$order]) + # order genes based on heatmap p above
    guides(colour = guide_colourbar(title = "Row (Gene) Z-Score", barwidth = 10), 
           size = guide_legend(title = "Proportion Detected")) + theme_cowplot(16) + #coord_flip() + 
    theme(panel.grid.major = element_line(colour = "gray90"), 
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
          legend.position = "top", legend.justification = "center", legend.title.position = "top") + 
    labs(title = "Cancer vs. Control", x = "Group", y ="Genes")
reset.fig()
```

    [1m[22mScale for [32msize[39m is already present.
    Adding another scale for [32msize[39m, which will replace the existing scale.



    
![png](Integrated_files/Integrated_243_1.png)
    


# 10. DE analysis between conditions in each cluster/cell type

For this analysis, each cell is consider a sample. [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) and [glmGamPoi](https://bioconductor.org/packages/release/bioc/html/glmGamPoi.html) are used to fit the gene-wise dispersion, its trend and calculate the MAP based on the quasi-likelihood framework for single-cell data. See [recommendations for single-cell analysis](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#recommendations-for-single-cell-analysis) for more information.

> If **glmGamPoi** is used in published research, please cite:
    Ahlmann-Eltze, C., Huber, W. (2020) glmGamPoi: Fitting Gamma-Poisson
    Generalized Linear Models on Single Cell Count Data. Bioinformatics.
    https://doi.org/10.1093/bioinformatics/btaa1009


```R
# Use previously prepared input
table(kept$CellType_1, kept$Sample)
```


               
                Control1 KidneyCancer LungCancer
      CD8T          3507         2509       2041
      CD4T          5826        10231       6039
      gdT             23          918        205
      NK            1997          280       1254
      B             1563          631       1997
      Monocytes     3107          609        802
      DCs            366           48         15


## Performs likelihood ratio test (LRT) using `DESeq2`

LRT test is testing whether the term(s) removed in the 'reduced' model explains a significant amount of variation in the data. The LRT p-values are determined solely by the difference in deviance between the 'full' and 'reduced' model formula, i.e. p-values is identical in all contrasts. The test statistic and p-values may involve the testing of one or more log2 fold changes.

If the adjusted p-value (`padj`) is small, then for the set of genes with those small adjusted p-values, the additional coefficient in 'full' and not in 'reduced' increased the log likelihood more than would be expected if their true value was zero.

**Define thresholds, variable of interest and include unwanted variation present in the data**


```R
minCell <- 10 # Ignore clusters that have fewer than 10 cells in test and reference conditions
padj_cutoff <- 0.05 # FDR threshold

test <- "condition"

full <- ~ condition # full model
reduced <- ~ 1 # reduced model
```


```R
outputRes <- function(dds, test = "condition", ref = NULL, old_ref = NULL) {
    # Build results tables
    res <- list()
    freq <- data.frame(table(dds$condition))
    ref <- if(is.null(ref)) levels(dds$condition)[1] else ref
    t <- NULL
    for(j in resultsNames(dds)) {
        if(grepl(test, j)) {
            t <- sub(test, "", j)
            if(!is.null(old_ref) && t == old_ref) next
            
            if(freq$Freq[freq$Var1 == ref] >= minCell & freq$Freq[freq$Var1 == t] >= minCell) {
                resultsname <- j
            } else next
        } else next
        
        id <- paste0(t, "_vs_", ref)
        title <- paste0(t, "_vs_", ref, "_", i)
        message(paste(t, "vs.", ref, "in Cluster", i))

        res[[id]] <- results(dds, name = resultsname, alpha = padj_cutoff)
        res[[id]]$ID <- mcols(dds)$ID
        res[[id]]$Symbol <- mcols(dds)$Symbol
        slot(mcols(res[[id]]), 
             "listData")$type <- c(slot(mcols(res[[id]]), "listData")$type[1:6], 
                                   "annotation", "annotation")
        slot(mcols(res[[id]]), 
             "listData")$description <- c(slot(mcols(res[[id]]), "listData")$description[1:6], 
                                          "ID", "Symbol")
        res[[id]]<- res[[id]][,c(7,8,1:6)]

        # Print the total number of DEGs at defined FDR threshold
        print(summary(res[[id]]))
    
        # Print top 10 genes, ordered by FDR then p-values
        df <- as.data.frame(res[[id]]) %>% arrange(padj, pvalue)
        print(head(df, 10))

        # Output DE results to text file
        outfile <- paste0(file_id, "_DESeq2_DE_results_", title, ".tsv")
        print(paste("Write to file:", outfile))
        write.table(df, file = outfile, sep = "\t", quote = F, row.names = F, col.names = T)
        flush.console()
    }
    return(res)
}
```

The code below currently only work with a `test` term with 3 or less levels.

<div class="alert alert-warning">
    In the example below, we group cells and compare between conditions in each <b>coarse cell type</b>, i.e. <code>CellType_1</code>.
</div>


```R
dds <- list()
res <- list()

# Loop through available cell groups
for(i in levels(kept$CellType_1)) {
    tmp <- kept[,kept$CellType_1 == i]
    # Build DESeqDataSet
    dds[[i]] <- DESeqDataSetFromMatrix(countData = kept_c[, kept$CellType_1 == i], 
                                       colData = droplevels(colData(tmp)), 
                                       rowData = rowData(tmp), design = full)
    
    # Use sizeFactors from SingleCellExperiment
    sizeFactors(dds[[i]]) <- sizeFactors(kept)[kept$CellType_1 == i]
    
    # Run LRT test
    dds[[i]] <- DESeq(dds[[i]], test = "LRT", useT = TRUE, minmu = 1e-6, minReplicatesForReplace = Inf, 
                      fitType = "glmGamPoi", sfType = "poscounts", reduced = reduced, quiet = TRUE)

    # Build results tables
    res[[i]] <- outputRes(dds[[i]], test = test, ref = levels(colData(dds[[i]])[,test])[1])
    
    # Relevel to get final comparison
    if(nlevels(colData(dds[[i]])[,test]) == 3) {
        colData(dds[[i]])[,test] <- relevel(colData(dds[[i]])[,test], ref = levels(colData(dds[[i]])[,test])[2])
        dds[[i]] <- nbinomLRT(dds[[i]], minmu = 1e-6, type = "glmGamPoi", reduced = reduced, quiet = TRUE)
        res[[i]] <- c(res[[i]], outputRes(dds[[i]], test = test,
                                          ref = levels(colData(dds[[i]])[,test])[1], 
                                          old_ref = levels(colData(tmp[[i]])[,test])[1]))
    }
    flush.console()
}
```

    converting counts to integer mode
    
    Cancer vs. Control in Cluster CD8T
    


    
    out of 15003 with nonzero total read count
    adjusted p-value < 0.05
    LFC > 0 (up)       : 3602, 24%
    LFC < 0 (down)     : 6974, 46%
    outliers [1]       : 0, 0%
    low counts [2]     : 290, 1.9%
    (mean count < 0)
    [1] see 'cooksCutoff' argument of ?results
    [2] see 'independentFiltering' argument of ?results
    
    NULL
                         ID  Symbol  baseMean log2FoldChange lfcSE     stat pvalue padj
    SKI     ENSG00000157933     SKI 3.4592466       2.160663    NA 4168.760      0    0
    KCNAB2  ENSG00000069424  KCNAB2 1.9573693      -1.315292    NA 1982.481      0    0
    PIK3CD  ENSG00000171608  PIK3CD 2.7159255      -1.287634    NA 2387.496      0    0
    PRDM2   ENSG00000116731   PRDM2 1.6966803       1.247248    NA 1798.102      0    0
    CAPZB   ENSG00000077549   CAPZB 1.9914869      -1.357994    NA 2566.519      0    0
    RUNX3   ENSG00000020633   RUNX3 2.3073798       1.589740    NA 2292.349      0    0
    LDLRAP1 ENSG00000157978 LDLRAP1 0.9145429      -1.964889    NA 2115.266      0    0
    SYTL1   ENSG00000142765   SYTL1 1.6139938      -1.473589    NA 2166.530      0    0
    LCK     ENSG00000182866     LCK 1.6376228      -1.941380    NA 4723.661      0    0
    ZC3H12A ENSG00000163874 ZC3H12A 1.1459845       2.760562    NA 2348.123      0    0
    [1] "Write to file: 160k_All_DESeq2_DE_results_Cancer_vs_Control_CD8T.tsv"


    Warning message in asMethod(object):
    “sparse->dense coercion: allocating vector of size 2.6 GiB”
    converting counts to integer mode
    
    Cancer vs. Control in Cluster CD4T
    


    
    out of 15576 with nonzero total read count
    adjusted p-value < 0.05
    LFC > 0 (up)       : 3859, 25%
    LFC < 0 (down)     : 7516, 48%
    outliers [1]       : 0, 0%
    low counts [2]     : 0, 0%
    (mean count < 0)
    [1] see 'cooksCutoff' argument of ?results
    [2] see 'independentFiltering' argument of ?results
    
    NULL
                          ID   Symbol  baseMean log2FoldChange lfcSE     stat pvalue padj
    SKI      ENSG00000157933      SKI 3.7293479       1.506668    NA 4835.401      0    0
    TNFRSF14 ENSG00000157873 TNFRSF14 0.8062638      -1.210239    NA 2362.961      0    0
    ENO1     ENSG00000074800     ENO1 0.3854514      -1.650642    NA 2407.095      0    0
    PIK3CD   ENSG00000171608   PIK3CD 2.2587332      -1.113516    NA 4244.321      0    0
    PRDM2    ENSG00000116731    PRDM2 1.6696914       1.113826    NA 2665.380      0    0
    CAPZB    ENSG00000077549    CAPZB 1.6121652      -1.205107    NA 3981.370      0    0
    NIPAL3   ENSG00000001461   NIPAL3 0.2300899      -1.858950    NA 1539.109      0    0
    RUNX3    ENSG00000020633    RUNX3 1.7586888       2.065901    NA 4311.249      0    0
    LDLRAP1  ENSG00000157978  LDLRAP1 0.7523601      -1.838951    NA 4581.246      0    0
    MAN1C1   ENSG00000117643   MAN1C1 0.2849186      -2.194859    NA 2387.884      0    0
    [1] "Write to file: 160k_All_DESeq2_DE_results_Cancer_vs_Control_CD4T.tsv"


    converting counts to integer mode
    
    Cancer vs. Control in Cluster gdT
    


    
    out of 12961 with nonzero total read count
    adjusted p-value < 0.05
    LFC > 0 (up)       : 379, 2.9%
    LFC < 0 (down)     : 880, 6.8%
    outliers [1]       : 0, 0%
    low counts [2]     : 1986, 15%
    (mean count < 0)
    [1] see 'cooksCutoff' argument of ?results
    [2] see 'independentFiltering' argument of ?results
    
    NULL
                        ID Symbol   baseMean log2FoldChange lfcSE     stat       pvalue         padj
    STAT1  ENSG00000115415  STAT1  0.4494064      -3.868866    NA 224.8718 9.702889e-47 1.064892e-42
    IL32   ENSG00000008517   IL32  4.1113872      -2.361776    NA 190.5925 2.329404e-40 1.278260e-36
    CD300A ENSG00000167851 CD300A  0.6893034      -3.457464    NA 171.6694 9.104694e-37 3.330801e-33
    GBP5   ENSG00000154451   GBP5  0.8341454      -2.713382    NA 125.0078 1.104422e-27 3.030258e-24
    CORO1A ENSG00000102879 CORO1A  2.6483719      -1.855612    NA 122.0071 4.354021e-27 9.557075e-24
    GIMAP4 ENSG00000133574 GIMAP4  0.5301034      -3.278169    NA 121.3703 5.827676e-27 1.065979e-23
    MT-ND3 ENSG00000198840 MT-ND3 11.2251075       2.195788    NA 115.2769 9.559431e-26 1.498782e-22
    NLRC5  ENSG00000140853  NLRC5  2.0013631      -1.969602    NA 114.7273 1.231158e-25 1.559191e-22
    MT-CO3 ENSG00000198938 MT-CO3 23.0284929       2.015318    NA 114.6452 1.278608e-25 1.559191e-22
    RIPOR2 ENSG00000111913 RIPOR2  0.7834515      -2.878145    NA 114.0462 1.684886e-25 1.849163e-22
    [1] "Write to file: 160k_All_DESeq2_DE_results_Cancer_vs_Control_gdT.tsv"


    converting counts to integer mode
    
    Cancer vs. Control in Cluster NK
    


    
    out of 14397 with nonzero total read count
    adjusted p-value < 0.05
    LFC > 0 (up)       : 2803, 19%
    LFC < 0 (down)     : 6112, 42%
    outliers [1]       : 0, 0%
    low counts [2]     : 0, 0%
    (mean count < 0)
    [1] see 'cooksCutoff' argument of ?results
    [2] see 'independentFiltering' argument of ?results
    
    NULL
                          ID   Symbol baseMean log2FoldChange lfcSE     stat pvalue padj
    KCNAB2   ENSG00000069424   KCNAB2 3.991316      -1.765893    NA 2204.794      0    0
    FGR      ENSG00000000938      FGR 3.746266      -1.988989    NA 2650.883      0    0
    ZC3H12A  ENSG00000163874  ZC3H12A 1.334464       3.030274    NA 2440.874      0    0
    DENND2D  ENSG00000162777  DENND2D 1.939116      -2.309743    NA 2257.058      0    0
    TENT5C   ENSG00000183508   TENT5C 1.013045       3.223403    NA 1976.294      0    0
    ARHGAP30 ENSG00000186517 ARHGAP30 2.595008      -1.581996    NA 1963.546      0    0
    IER5     ENSG00000162783     IER5 3.698717       2.914985    NA 3354.113      0    0
    RGS2     ENSG00000116741     RGS2 1.193566       3.379973    NA 2010.045      0    0
    YPEL5    ENSG00000119801    YPEL5 2.962236       1.892579    NA 2036.240      0    0
    PLEK     ENSG00000115956     PLEK 2.628495      -2.013504    NA 2293.591      0    0
    [1] "Write to file: 160k_All_DESeq2_DE_results_Cancer_vs_Control_NK.tsv"


    converting counts to integer mode
    
    Cancer vs. Control in Cluster B
    


    
    out of 14591 with nonzero total read count
    adjusted p-value < 0.05
    LFC > 0 (up)       : 3092, 21%
    LFC < 0 (down)     : 6209, 43%
    outliers [1]       : 0, 0%
    low counts [2]     : 282, 1.9%
    (mean count < 0)
    [1] see 'cooksCutoff' argument of ?results
    [2] see 'independentFiltering' argument of ?results
    
    NULL
                          ID   Symbol   baseMean log2FoldChange lfcSE     stat pvalue padj
    CASZ1    ENSG00000130940    CASZ1  1.3347142       3.245481    NA 1866.982      0    0
    PRDM2    ENSG00000116731    PRDM2  4.4159749       1.722118    NA 2561.127      0    0
    LAPTM5   ENSG00000162511   LAPTM5 38.9600479       1.659145    NA 4725.619      0    0
    SMAP2    ENSG00000084070    SMAP2  7.1755398       1.597823    NA 2115.373      0    0
    ZNF644   ENSG00000122482   ZNF644  3.5515181       1.878388    NA 2248.015      0    0
    CD53     ENSG00000143119     CD53  1.9004307      -1.783876    NA 2049.031      0    0
    IFI16    ENSG00000163565    IFI16  1.4631673      -2.107295    NA 1813.024      0    0
    SLAMF6   ENSG00000162739   SLAMF6  0.8808965      -3.064470    NA 2311.689      0    0
    ARHGAP30 ENSG00000186517 ARHGAP30  1.0749031      -2.243309    NA 1953.135      0    0
    SELL     ENSG00000188404     SELL  2.3867811      -2.758018    NA 2638.676      0    0
    [1] "Write to file: 160k_All_DESeq2_DE_results_Cancer_vs_Control_B.tsv"


    converting counts to integer mode
    
    Cancer vs. Control in Cluster Monocytes
    


    
    out of 14692 with nonzero total read count
    adjusted p-value < 0.05
    LFC > 0 (up)       : 3578, 24%
    LFC < 0 (down)     : 6167, 42%
    outliers [1]       : 0, 0%
    low counts [2]     : 0, 0%
    (mean count < 0)
    [1] see 'cooksCutoff' argument of ?results
    [2] see 'independentFiltering' argument of ?results
    
    NULL
                         ID  Symbol   baseMean log2FoldChange lfcSE     stat pvalue padj
    NADK    ENSG00000008130    NADK  2.3882142      -2.506941    NA 3366.123      0    0
    DHRS3   ENSG00000162496   DHRS3  0.5079112       5.595458    NA 2984.811      0    0
    PLEKHM2 ENSG00000116786 PLEKHM2  2.0904786       1.676018    NA 2379.320      0    0
    RSRP1   ENSG00000117616   RSRP1  6.7344313      -1.185579    NA 1830.319      0    0
    LAPTM5  ENSG00000162511  LAPTM5 19.3224569       1.561981    NA 8324.010      0    0
    GBP1    ENSG00000117228    GBP1  2.7424973      -5.586310    NA 3732.804      0    0
    GBP5    ENSG00000154451    GBP5  2.4623455      -4.919611    NA 2880.421      0    0
    ZNF644  ENSG00000122482  ZNF644  1.2415164       1.918292    NA 1937.897      0    0
    PLEKHO1 ENSG00000023902 PLEKHO1  7.1295468      -1.172884    NA 1911.966      0    0
    S100A10 ENSG00000197747 S100A10  6.3500356       1.493481    NA 2361.434      0    0
    [1] "Write to file: 160k_All_DESeq2_DE_results_Cancer_vs_Control_Monocytes.tsv"


    converting counts to integer mode
    
    Cancer vs. Control in Cluster DCs
    


    
    out of 13627 with nonzero total read count
    adjusted p-value < 0.05
    LFC > 0 (up)       : 1691, 12%
    LFC < 0 (down)     : 3155, 23%
    outliers [1]       : 0, 0%
    low counts [2]     : 1833, 13%
    (mean count < 0)
    [1] see 'cooksCutoff' argument of ?results
    [2] see 'independentFiltering' argument of ?results
    
    NULL
                         ID  Symbol  baseMean log2FoldChange lfcSE     stat       pvalue         padj
    MT-CYB  ENSG00000198727  MT-CYB  3.422211       1.914089    NA 349.6834 1.423456e-58 1.678824e-54
    MT-CO3  ENSG00000198938  MT-CO3 12.952293       1.383054    NA 290.0922 6.937726e-51 4.091177e-47
    NAGK    ENSG00000124357    NAGK  3.277954      -3.037566    NA 251.7807 1.301644e-45 5.117196e-42
    PSMB9   ENSG00000240065   PSMB9  2.058332      -3.526640    NA 246.7615 6.704739e-45 1.976892e-41
    MT-ATP6 ENSG00000198899 MT-ATP6 13.032460       1.244259    NA 243.8406 1.749934e-44 3.834807e-41
    MT-ND3  ENSG00000198840  MT-ND3  6.242289       1.331496    NA 243.5103 1.950894e-44 3.834807e-41
    JAML    ENSG00000160593    JAML  6.382330      -2.811850    NA 242.0323 3.175517e-44 5.350293e-41
    YTHDF3  ENSG00000185728  YTHDF3  0.587758       2.522687    NA 235.9006 2.423251e-43 3.572478e-40
    FGD2    ENSG00000146192    FGD2  3.833157      -2.935036    NA 234.9986 3.272532e-43 4.288472e-40
    MT-ND1  ENSG00000198888  MT-ND1  2.572526       1.614478    NA 211.8737 8.315503e-40 9.807304e-37
    [1] "Write to file: 160k_All_DESeq2_DE_results_Cancer_vs_Control_DCs.tsv"


### Save DESeq2 DE results to `metadata`

<div class="alert alert-info">
    <strong>Tip!</strong> In the accompanied Shiny App, it will look for list(s) named with "DESeq2_" in the prefix in <code>metadata()</code> and show their content under the <u>DEA (DESeq2)</u> section of the website.
</div>


```R
# Check list names
names(res)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'CD8T'</li><li>'CD4T'</li><li>'gdT'</li><li>'NK'</li><li>'B'</li><li>'Monocytes'</li><li>'DCs'</li></ol>




```R
metadata(combined)[["DESeq2_Cancer_Control"]] <- purrr::map(res, "Cancer_vs_Control")
```

## Visualise MA-plot

The function `plotMA` shows the log2 fold changes attributable to a given variable over the mean of normalized counts for all the samples in the comparison. Points will be colored blue if the adjusted p value is less than `padj_cutoff`. Points which fall out of the window are plotted as open triangles pointing either up or down.

<div class="alert alert-warning">
    <strong>Warning!</strong> Usually, log2 fold change (LFC) shrinkage is performed with the <code>lfcShrink</code> function to generate more accurate estimates. However, shrinkage cannot be performed on results fitted by <code>glmGamPoiM</code>. Therefore, you may observed large fold changes from genes with low information, including genes that have low counts or high dispersion values.
</div>


```R
tmp <- purrr::map(res, "Cancer_vs_Control") # all levels have same pvalue/padj in LRT
```


```R
minexp <- 0.1 # red line at baseMean = 0.1
plist <- list()

j <- 1
for(i in 1:length(tmp)) {
    plist[[j]] <- ggplotify::as.ggplot(~{
        DESeq2::plotMA(tmp[[i]], ylim = c(-4, 4))
        title(paste("Cancer vs. Control, ", names(tmp)[i]), line = 0.5)
        abline(v = minexp, lwd = 2, lty = 2, col = "red")
    })
    j <- j + 1
}

fig(width = 18, height = 10)
plot_grid(plotlist = plist, ncol = 4, align = "vh")
reset.fig()

```


    
![png](Integrated_files/Integrated_256_0.png)
    


## Visualise gene expression as boxplots

Show top N genes with `baseMean` > `minexp`.

<div class="alert alert-warning">
    The <code>baseMean</code> represents the mean of normalized counts for all samples in the dataset, not the subset of samples specified by 'contrast'.
</div>

<div class="alert alert-warning">
    In the example below, we group cells and compare between conditions in each <b>coarse cell type</b>, i.e. <code>CellType_1</code>.
</div>


```R
minexp <- 0.1 # red line at baseMean = 0.1
nTop <- 10

fig(width = 16, height = 6)
for(i in names(tmp)) {
    g <- if(is.numeric(i)) paste("Cluster", i) else i
    df <- as.data.frame(tmp[[i]]) %>% arrange(padj, pvalue) %>% filter(padj < padj_cutoff & baseMean > minexp)
    if(nrow(df) < 1) {
        message(sprintf("There is no DE genes that passed the defined thresholds in '%s'.", g))
        next
    }
    
    geneNames <- head(rownames(df), nTop)
    p <- plotBox(kept[,kept$CellType_1 == i], features = geneNames, group_by = c("CellType_1","condition"), 
                 color_by = "Detected", box_colors = c_heatmap_col2, guides_barheight = 10, 
                 facet_ncol = 10, x.text_angle = 90, x.text_size = 14, theme_size = 16) + 
    geom_boxplot(color = "black", linewidth = 0.1, outlier.size = 0.2) +
    ggtitle(g, subtitle = sprintf("Total %d genes with FDR < %.2f and baseMean > %.1f", 
                                                    nrow(df), padj_cutoff, minexp))
    
    print(p)
    flush.console()
}
reset.fig()
```


    
![png](Integrated_files/Integrated_258_0.png)
    



    
![png](Integrated_files/Integrated_258_1.png)
    



    
![png](Integrated_files/Integrated_258_2.png)
    



    
![png](Integrated_files/Integrated_258_3.png)
    



    
![png](Integrated_files/Integrated_258_4.png)
    



    
![png](Integrated_files/Integrated_258_5.png)
    



    
![png](Integrated_files/Integrated_258_6.png)
    


# 11. Functional analysis using `enrichR`

## Enriched pathways

We use the `enrichR` package to access the [Enrichr](https://maayanlab.cloud/Enrichr/) website to perform gene set enrichment analysis.

### Set up `enrichR`

All the available gene-set libraries are listed in `dbs`


```R
# This function generates the whole list of database for the enrichR. 
dbs <- listEnrichrDbs()
print(paste("Number of available databases from Enrichr:", nrow(dbs)))
#head(dbs)
```

    [1] "Number of available databases from Enrichr: 240"


Change `dbsSel` to remove or include more gene-set libraries in the enrichment analysis.


```R
# Human
dbsSel <- c("GO_Biological_Process_2023", # Ontologies
#            "GO_Molecular_Function_2023", # Ontologies
#            "GO_Cellular_Component_2023", # Ontologies
            "Reactome_Pathways_2024",     # Pathways
            "WikiPathways_2024_Human",    # Pathways
            "CellMarker_2024")            # Cell types

# Mouse
#dbsSel <- c("GO_Biological_Process_2023", # Ontologies
#            "GO_Molecular_Function_2023", # Ontologies
#            "GO_Cellular_Component_2023", # Ontologies
#            "Reactome_Pathways_2024",     # Pathways
#            "WikiPathways_2024_Mouse",    # Pathways
#            "CellMarker_2024")            # Cell types
```

### Run enrichR

**On upregulated 'Cluster' marker genes (from `findMarkers`)**


```R
input <- metadata(combined)[['findMarkers_Cluster_up']]
metadata(combined)[['enrichR_findMarkers_Cluster_up']] <- runEnrichR(input, dbs = dbsSel, direction = "up", 
                                                                     column_by = "Symbol")
```

    Detecting findMarkers input.
    Connection changed to https://maayanlab.cloud/Enrichr/
    
    Connection is Live!
    
    Running enrichR on 'Cluster1' with 262 up-regulated genes.


    Uploading data to Enrichr... Done.
      Querying GO_Biological_Process_2023... Done.
      Querying Reactome_Pathways_2024... Done.
      Querying WikiPathways_2024_Human... Done.
      Querying CellMarker_2024... Done.
    Parsing results... Done.


    Running enrichR on 'Cluster2' with 251 up-regulated genes.


    Uploading data to Enrichr... Done.
      Querying GO_Biological_Process_2023... Done.
      Querying Reactome_Pathways_2024... Done.
      Querying WikiPathways_2024_Human... Done.
      Querying CellMarker_2024... Done.
    Parsing results... Done.


    Running enrichR on 'Cluster3' with 250 up-regulated genes.


    Uploading data to Enrichr... Done.
      Querying GO_Biological_Process_2023... Done.
      Querying Reactome_Pathways_2024... Done.
      Querying WikiPathways_2024_Human... Done.
      Querying CellMarker_2024... Done.
    Parsing results... Done.


    Running enrichR on 'Cluster4' with 241 up-regulated genes.


    Uploading data to Enrichr... Done.
      Querying GO_Biological_Process_2023... Done.
      Querying Reactome_Pathways_2024... Done.
      Querying WikiPathways_2024_Human... Done.
      Querying CellMarker_2024... Done.
    Parsing results... Done.


    Running enrichR on 'Cluster5' with 269 up-regulated genes.


    Uploading data to Enrichr... Done.
      Querying GO_Biological_Process_2023... Done.
      Querying Reactome_Pathways_2024... Done.
      Querying WikiPathways_2024_Human... Done.
      Querying CellMarker_2024... Done.
    Parsing results... Done.


    Running enrichR on 'Cluster6' with 255 up-regulated genes.


    Uploading data to Enrichr... Done.
      Querying GO_Biological_Process_2023... Done.
      Querying Reactome_Pathways_2024... Done.
      Querying WikiPathways_2024_Human... Done.
      Querying CellMarker_2024... Done.
    Parsing results... Done.


    Running enrichR on 'Cluster7' with 254 up-regulated genes.


    Uploading data to Enrichr... Done.
      Querying GO_Biological_Process_2023... Done.
      Querying Reactome_Pathways_2024... Done.
      Querying WikiPathways_2024_Human... Done.
      Querying CellMarker_2024... Done.
    Parsing results... Done.


    Running enrichR on 'Cluster8' with 238 up-regulated genes.


    Uploading data to Enrichr... Done.
      Querying GO_Biological_Process_2023... Done.
      Querying Reactome_Pathways_2024... Done.
      Querying WikiPathways_2024_Human... Done.
      Querying CellMarker_2024... Done.
    Parsing results... Done.


    Running enrichR on 'Cluster9' with 247 up-regulated genes.


    Uploading data to Enrichr... Done.
      Querying GO_Biological_Process_2023... Done.
      Querying Reactome_Pathways_2024... Done.
      Querying WikiPathways_2024_Human... Done.
      Querying CellMarker_2024... Done.
    Parsing results... Done.


    Running enrichR on 'Cluster10' with 250 up-regulated genes.


    Uploading data to Enrichr... Done.
      Querying GO_Biological_Process_2023... Done.
      Querying Reactome_Pathways_2024... Done.
      Querying WikiPathways_2024_Human... Done.
      Querying CellMarker_2024... Done.
    Parsing results... Done.


    Running enrichR on 'Cluster11' with 253 up-regulated genes.


    Uploading data to Enrichr... Done.
      Querying GO_Biological_Process_2023... Done.
      Querying Reactome_Pathways_2024... Done.
      Querying WikiPathways_2024_Human... Done.
      Querying CellMarker_2024... Done.
    Parsing results... Done.


    Running enrichR on 'Cluster12' with 257 up-regulated genes.


    Uploading data to Enrichr... Done.
      Querying GO_Biological_Process_2023... Done.
      Querying Reactome_Pathways_2024... Done.
      Querying WikiPathways_2024_Human... Done.
      Querying CellMarker_2024... Done.
    Parsing results... Done.


    Running enrichR on 'Cluster13' with 243 up-regulated genes.


    Uploading data to Enrichr... Done.
      Querying GO_Biological_Process_2023... Done.
      Querying Reactome_Pathways_2024... Done.
      Querying WikiPathways_2024_Human... Done.
      Querying CellMarker_2024... Done.
    Parsing results... Done.


    Running enrichR on 'Cluster14' with 254 up-regulated genes.


    Uploading data to Enrichr... Done.
      Querying GO_Biological_Process_2023... Done.
      Querying Reactome_Pathways_2024... Done.
      Querying WikiPathways_2024_Human... Done.
      Querying CellMarker_2024... Done.
    Parsing results... Done.


    Running enrichR on 'Cluster15' with 252 up-regulated genes.


    Uploading data to Enrichr... Done.
      Querying GO_Biological_Process_2023... Done.
      Querying Reactome_Pathways_2024... Done.
      Querying WikiPathways_2024_Human... Done.
      Querying CellMarker_2024... Done.
    Parsing results... Done.


    Running enrichR on 'Cluster16' with 256 up-regulated genes.


    Uploading data to Enrichr... Done.
      Querying GO_Biological_Process_2023... Done.
      Querying Reactome_Pathways_2024... Done.
      Querying WikiPathways_2024_Human... Done.
      Querying CellMarker_2024... Done.
    Parsing results... Done.


    Running enrichR on 'Cluster17' with 239 up-regulated genes.


    Uploading data to Enrichr... Done.
      Querying GO_Biological_Process_2023... Done.
      Querying Reactome_Pathways_2024... Done.
      Querying WikiPathways_2024_Human... Done.
      Querying CellMarker_2024... Done.
    Parsing results... Done.


    Running enrichR on 'Cluster18' with 260 up-regulated genes.


    Uploading data to Enrichr... Done.
      Querying GO_Biological_Process_2023... Done.
      Querying Reactome_Pathways_2024... Done.
      Querying WikiPathways_2024_Human... Done.
      Querying CellMarker_2024... Done.
    Parsing results... Done.


    Running enrichR on 'Cluster19' with 254 up-regulated genes.


    Uploading data to Enrichr... Done.
      Querying GO_Biological_Process_2023... Done.
      Querying Reactome_Pathways_2024... Done.
      Querying WikiPathways_2024_Human... Done.
      Querying CellMarker_2024... Done.
    Parsing results... Done.


    Running enrichR on 'Cluster20' with 257 up-regulated genes.


    Uploading data to Enrichr... Done.
      Querying GO_Biological_Process_2023... Done.
      Querying Reactome_Pathways_2024... Done.
      Querying WikiPathways_2024_Human... Done.
      Querying CellMarker_2024... Done.
    Parsing results... Done.


    Running enrichR on 'Cluster21' with 248 up-regulated genes.


    Uploading data to Enrichr... Done.
      Querying GO_Biological_Process_2023... Done.
      Querying Reactome_Pathways_2024... Done.
      Querying WikiPathways_2024_Human... Done.
      Querying CellMarker_2024... Done.
    Parsing results... Done.


    Running enrichR on 'Cluster22' with 249 up-regulated genes.


    Uploading data to Enrichr... Done.
      Querying GO_Biological_Process_2023... Done.
      Querying Reactome_Pathways_2024... Done.
      Querying WikiPathways_2024_Human... Done.
      Querying CellMarker_2024... Done.
    Parsing results... Done.


    Running enrichR on 'Cluster23' with 249 up-regulated genes.


    Uploading data to Enrichr... Done.
      Querying GO_Biological_Process_2023... Done.
      Querying Reactome_Pathways_2024... Done.
      Querying WikiPathways_2024_Human... Done.
      Querying CellMarker_2024... Done.
    Parsing results... Done.


    Running enrichR on 'Cluster24' with 247 up-regulated genes.


    Uploading data to Enrichr... Done.
      Querying GO_Biological_Process_2023... Done.
      Querying Reactome_Pathways_2024... Done.
      Querying WikiPathways_2024_Human... Done.
      Querying CellMarker_2024... Done.
    Parsing results... Done.


    Running enrichR on 'Cluster25' with 247 up-regulated genes.


    Uploading data to Enrichr... Done.
      Querying GO_Biological_Process_2023... Done.
      Querying Reactome_Pathways_2024... Done.
      Querying WikiPathways_2024_Human... Done.
      Querying CellMarker_2024... Done.
    Parsing results... Done.


    Running enrichR on 'Cluster26' with 264 up-regulated genes.


    Uploading data to Enrichr... Done.
      Querying GO_Biological_Process_2023... Done.
      Querying Reactome_Pathways_2024... Done.
      Querying WikiPathways_2024_Human... Done.
      Querying CellMarker_2024... Done.
    Parsing results... Done.


    Running enrichR on 'Cluster27' with 239 up-regulated genes.


    Uploading data to Enrichr... Done.
      Querying GO_Biological_Process_2023... Done.
      Querying Reactome_Pathways_2024... Done.
      Querying WikiPathways_2024_Human... Done.
      Querying CellMarker_2024... Done.
    Parsing results... Done.


    Running enrichR on 'Cluster28' with 249 up-regulated genes.


    Uploading data to Enrichr... Done.
      Querying GO_Biological_Process_2023... Done.
      Querying Reactome_Pathways_2024... Done.
      Querying WikiPathways_2024_Human... Done.
      Querying CellMarker_2024... Done.
    Parsing results... Done.


    Running enrichR on 'Cluster29' with 232 up-regulated genes.


    Uploading data to Enrichr... Done.
      Querying GO_Biological_Process_2023... Done.
      Querying Reactome_Pathways_2024... Done.
      Querying WikiPathways_2024_Human... Done.
      Querying CellMarker_2024... Done.
    Parsing results... Done.


    Running enrichR on 'Cluster30' with 242 up-regulated genes.


    Uploading data to Enrichr... Done.
      Querying GO_Biological_Process_2023... Done.
      Querying Reactome_Pathways_2024... Done.
      Querying WikiPathways_2024_Human... Done.
      Querying CellMarker_2024... Done.
    Parsing results... Done.


    Running enrichR on 'Cluster31' with 228 up-regulated genes.


    Uploading data to Enrichr... Done.
      Querying GO_Biological_Process_2023... Done.
      Querying Reactome_Pathways_2024... Done.
      Querying WikiPathways_2024_Human... Done.
      Querying CellMarker_2024... Done.
    Parsing results... Done.


    Running enrichR on 'Cluster32' with 243 up-regulated genes.


    Uploading data to Enrichr... Done.
      Querying GO_Biological_Process_2023... Done.
      Querying Reactome_Pathways_2024... Done.
      Querying WikiPathways_2024_Human... Done.
      Querying CellMarker_2024... Done.
    Parsing results... Done.


    Running enrichR on 'Cluster33' with 237 up-regulated genes.


    Uploading data to Enrichr... Done.
      Querying GO_Biological_Process_2023... Done.
      Querying Reactome_Pathways_2024... Done.
      Querying WikiPathways_2024_Human... Done.
      Querying CellMarker_2024... Done.
    Parsing results... Done.


    Running enrichR on 'Cluster34' with 241 up-regulated genes.


    Uploading data to Enrichr... Done.
      Querying GO_Biological_Process_2023... Done.
      Querying Reactome_Pathways_2024... Done.
      Querying WikiPathways_2024_Human... Done.
      Querying CellMarker_2024... Done.
    Parsing results... Done.


    Running enrichR on 'Cluster35' with 267 up-regulated genes.


    Uploading data to Enrichr... Done.
      Querying GO_Biological_Process_2023... Done.
      Querying Reactome_Pathways_2024... Done.
      Querying WikiPathways_2024_Human... Done.
      Querying CellMarker_2024... Done.
    Parsing results... Done.


    Running enrichR on 'Cluster36' with 253 up-regulated genes.


    Uploading data to Enrichr... Done.
      Querying GO_Biological_Process_2023... Done.
      Querying Reactome_Pathways_2024... Done.
      Querying WikiPathways_2024_Human... Done.
      Querying CellMarker_2024... Done.
    Parsing results... Done.


    Running enrichR on 'Cluster37' with 246 up-regulated genes.


    Uploading data to Enrichr... Done.
      Querying GO_Biological_Process_2023... Done.
      Querying Reactome_Pathways_2024... Done.
      Querying WikiPathways_2024_Human... Done.
      Querying CellMarker_2024... Done.
    Parsing results... Done.


    Running enrichR on 'Cluster38' with 248 up-regulated genes.


    Uploading data to Enrichr... Done.
      Querying GO_Biological_Process_2023... Done.
      Querying Reactome_Pathways_2024... Done.
      Querying WikiPathways_2024_Human... Done.
      Querying CellMarker_2024... Done.
    Parsing results... Done.


    Running enrichR on 'Cluster39' with 244 up-regulated genes.


    Uploading data to Enrichr... Done.
      Querying GO_Biological_Process_2023... Done.
      Querying Reactome_Pathways_2024... Done.
      Querying WikiPathways_2024_Human... Done.
      Querying CellMarker_2024... Done.
    Parsing results... Done.


    Running enrichR on 'Cluster40' with 249 up-regulated genes.


    Uploading data to Enrichr... Done.
      Querying GO_Biological_Process_2023... Done.
      Querying Reactome_Pathways_2024... Done.
      Querying WikiPathways_2024_Human... Done.
      Querying CellMarker_2024... Done.
    Parsing results... Done.


    Running enrichR on 'Cluster41' with 247 up-regulated genes.


    Uploading data to Enrichr... Done.
      Querying GO_Biological_Process_2023... Done.
      Querying Reactome_Pathways_2024... Done.
      Querying WikiPathways_2024_Human... Done.
      Querying CellMarker_2024... Done.
    Parsing results... Done.


    Running enrichR on 'Cluster42' with 215 up-regulated genes.


    Uploading data to Enrichr... Done.
      Querying GO_Biological_Process_2023... Done.
      Querying Reactome_Pathways_2024... Done.
      Querying WikiPathways_2024_Human... Done.
      Querying CellMarker_2024... Done.
    Parsing results... Done.


    Running enrichR on 'Cluster43' with 232 up-regulated genes.


    Uploading data to Enrichr... Done.
      Querying GO_Biological_Process_2023... Done.
      Querying Reactome_Pathways_2024... Done.
      Querying WikiPathways_2024_Human... Done.
      Querying CellMarker_2024... Done.
    Parsing results... Done.


    Running enrichR on 'Cluster44' with 213 up-regulated genes.


    Uploading data to Enrichr... Done.
      Querying GO_Biological_Process_2023... Done.
      Querying Reactome_Pathways_2024... Done.
      Querying WikiPathways_2024_Human... Done.
      Querying CellMarker_2024... Done.
    Parsing results... Done.


    Running enrichR on 'Cluster45' with 237 up-regulated genes.


    Uploading data to Enrichr... Done.
      Querying GO_Biological_Process_2023... Done.
      Querying Reactome_Pathways_2024... Done.
      Querying WikiPathways_2024_Human... Done.
      Querying CellMarker_2024... Done.
    Parsing results... Done.


    Running enrichR on 'Cluster46' with 247 up-regulated genes.


    Uploading data to Enrichr... Done.
      Querying GO_Biological_Process_2023... Done.
      Querying Reactome_Pathways_2024... Done.
      Querying WikiPathways_2024_Human... Done.
      Querying CellMarker_2024... Done.
    Parsing results... Done.


**On upregulated and downregulated 'Cancer vs. Control' DE genes (from `edgeR`)**


```R
input <- metadata(combined)[['edgeR_Cancer_Control']]
metadata(combined)[['enrichR_edgeR_Cancer_Control_up']] <- runEnrichR(input, dbs = dbsSel, direction = "up", 
                                                                 column_by = "Symbol")
```

    Detecting edgeR input.
    Connection changed to https://maayanlab.cloud/Enrichr/
    
    Connection is Live!
    
    Running enrichR on 'Cancer_Control' with 1000 up-regulated genes.


    Uploading data to Enrichr... Done.
      Querying GO_Biological_Process_2023... Done.
      Querying Reactome_Pathways_2024... Done.
      Querying WikiPathways_2024_Human... Done.
      Querying CellMarker_2024... Done.
    Parsing results... Done.



```R
input <- metadata(combined)[['edgeR_Cancer_Control']]
metadata(combined)[['enrichR_edgeR_Cancer_Control_dn']] <- runEnrichR(input, dbs = dbsSel, direction = "down", 
                                                                 column_by = "Symbol")
```

    Detecting edgeR input.
    Connection changed to https://maayanlab.cloud/Enrichr/
    
    Connection is Live!
    
    Running enrichR on 'Cancer_Control' with 1000 down-regulated genes.


    Uploading data to Enrichr... Done.
      Querying GO_Biological_Process_2023... Done.
      Querying Reactome_Pathways_2024... Done.
      Querying WikiPathways_2024_Human... Done.
      Querying CellMarker_2024... Done.
    Parsing results... Done.


**On upregulated and downregulated 'Cancer vs. Control' DE genes (from `DESeq2`)**


```R
input <- metadata(combined)[['DESeq2_Cancer_Control']]
metadata(combined)[['enrichR_DESeq2_Cancer_Control_up']] <- runEnrichR(input, dbs = dbsSel, direction = "up", 
                                                                       column_by = "Symbol")
```

    Detecting DESeq2 input.
    Connection changed to https://maayanlab.cloud/Enrichr/
    
    Connection is Live!
    
    Running enrichR on 'CD8T' with 1000 up-regulated genes.


    Uploading data to Enrichr... Done.
      Querying GO_Biological_Process_2023... Done.
      Querying Reactome_Pathways_2024... Done.
      Querying WikiPathways_2024_Human... Done.
      Querying CellMarker_2024... Done.
    Parsing results... Done.


    Running enrichR on 'CD4T' with 1000 up-regulated genes.


    Uploading data to Enrichr... Done.
      Querying GO_Biological_Process_2023... Done.
      Querying Reactome_Pathways_2024... Done.
      Querying WikiPathways_2024_Human... Done.
      Querying CellMarker_2024... Done.
    Parsing results... Done.


    Running enrichR on 'gdT' with 379 up-regulated genes.


    Uploading data to Enrichr... Done.
      Querying GO_Biological_Process_2023... Done.
      Querying Reactome_Pathways_2024... Done.
      Querying WikiPathways_2024_Human... Done.
      Querying CellMarker_2024... Done.
    Parsing results... Done.


    Running enrichR on 'NK' with 1000 up-regulated genes.


    Uploading data to Enrichr... Done.
      Querying GO_Biological_Process_2023... Done.
      Querying Reactome_Pathways_2024... Done.
      Querying WikiPathways_2024_Human... Done.
      Querying CellMarker_2024... Done.
    Parsing results... Done.


    Running enrichR on 'B' with 1000 up-regulated genes.


    Uploading data to Enrichr... Done.
      Querying GO_Biological_Process_2023... Done.
      Querying Reactome_Pathways_2024... Done.
      Querying WikiPathways_2024_Human... Done.
      Querying CellMarker_2024... Done.
    Parsing results... Done.


    Running enrichR on 'Monocytes' with 1000 up-regulated genes.


    Uploading data to Enrichr... Done.
      Querying GO_Biological_Process_2023... Done.
      Querying Reactome_Pathways_2024... Done.
      Querying WikiPathways_2024_Human... Done.
      Querying CellMarker_2024... Done.
    Parsing results... Done.


    Running enrichR on 'DCs' with 1000 up-regulated genes.


    Uploading data to Enrichr... Done.
      Querying GO_Biological_Process_2023... Done.
      Querying Reactome_Pathways_2024... Done.
      Querying WikiPathways_2024_Human... Done.
      Querying CellMarker_2024... Done.
    Parsing results... Done.



```R
input <- metadata(combined)[['DESeq2_Cancer_Control']]
metadata(combined)[['enrichR_DESeq2_Cancer_Control_dn']] <- runEnrichR(input, dbs = dbsSel, direction = "down", 
                                                                       column_by = "Symbol")
```

    Detecting DESeq2 input.
    Connection changed to https://maayanlab.cloud/Enrichr/
    
    Connection is Live!
    
    Running enrichR on 'CD8T' with 1000 down-regulated genes.


    Uploading data to Enrichr... Done.
      Querying GO_Biological_Process_2023... Done.
      Querying Reactome_Pathways_2024... Done.
      Querying WikiPathways_2024_Human... Done.
      Querying CellMarker_2024... Done.
    Parsing results... Done.


    Running enrichR on 'CD4T' with 1000 down-regulated genes.


    Uploading data to Enrichr... Done.
      Querying GO_Biological_Process_2023... Done.
      Querying Reactome_Pathways_2024... Done.
      Querying WikiPathways_2024_Human... Done.
      Querying CellMarker_2024... Done.
    Parsing results... Done.


    Running enrichR on 'gdT' with 880 down-regulated genes.


    Uploading data to Enrichr... Done.
      Querying GO_Biological_Process_2023... Done.
      Querying Reactome_Pathways_2024... Done.
      Querying WikiPathways_2024_Human... Done.
      Querying CellMarker_2024... Done.
    Parsing results... Done.


    Running enrichR on 'NK' with 1000 down-regulated genes.


    Uploading data to Enrichr... Done.
      Querying GO_Biological_Process_2023... Done.
      Querying Reactome_Pathways_2024... Done.
      Querying WikiPathways_2024_Human... Done.
      Querying CellMarker_2024... Done.
    Parsing results... Done.


    Running enrichR on 'B' with 1000 down-regulated genes.


    Uploading data to Enrichr... Done.
      Querying GO_Biological_Process_2023... Done.
      Querying Reactome_Pathways_2024... Done.
      Querying WikiPathways_2024_Human... Done.
      Querying CellMarker_2024... Done.
    Parsing results... Done.


    Running enrichR on 'Monocytes' with 1000 down-regulated genes.


    Uploading data to Enrichr... Done.
      Querying GO_Biological_Process_2023... Done.
      Querying Reactome_Pathways_2024... Done.
      Querying WikiPathways_2024_Human... Done.
      Querying CellMarker_2024... Done.
    Parsing results... Done.


    Running enrichR on 'DCs' with 1000 down-regulated genes.


    Uploading data to Enrichr... Done.
      Querying GO_Biological_Process_2023... Done.
      Querying Reactome_Pathways_2024... Done.
      Querying WikiPathways_2024_Human... Done.
      Querying CellMarker_2024... Done.
    Parsing results... Done.


### Plot enrichR results

Using `GO_Biological_Process_2023` as example.

**On upregulated 'Cluster' marker genes**


```R
fig(width = 16, height = 5)
plotEnrichR(metadata(combined)[['enrichR_findMarkers_Cluster_up']], db = "GO_Biological_Process_2023")
reset.fig()
```


    
![png](Integrated_files/Integrated_272_0.png)
    



    
![png](Integrated_files/Integrated_272_1.png)
    



    
![png](Integrated_files/Integrated_272_2.png)
    



    
![png](Integrated_files/Integrated_272_3.png)
    



    
![png](Integrated_files/Integrated_272_4.png)
    



    
![png](Integrated_files/Integrated_272_5.png)
    



    
![png](Integrated_files/Integrated_272_6.png)
    



    
![png](Integrated_files/Integrated_272_7.png)
    



    
![png](Integrated_files/Integrated_272_8.png)
    



    
![png](Integrated_files/Integrated_272_9.png)
    



    
![png](Integrated_files/Integrated_272_10.png)
    



    
![png](Integrated_files/Integrated_272_11.png)
    



    
![png](Integrated_files/Integrated_272_12.png)
    



    
![png](Integrated_files/Integrated_272_13.png)
    



    
![png](Integrated_files/Integrated_272_14.png)
    



    
![png](Integrated_files/Integrated_272_15.png)
    



    
![png](Integrated_files/Integrated_272_16.png)
    



    
![png](Integrated_files/Integrated_272_17.png)
    



    
![png](Integrated_files/Integrated_272_18.png)
    



    
![png](Integrated_files/Integrated_272_19.png)
    



    
![png](Integrated_files/Integrated_272_20.png)
    



    
![png](Integrated_files/Integrated_272_21.png)
    


    Warning message in plotEnrich(object[[group]][[db]], showTerms = showTerms, numChar = numChar, :
    “There are duplicated trimmed names in the plot, consider increasing the 'numChar' setting.”



    
![png](Integrated_files/Integrated_272_23.png)
    



    
![png](Integrated_files/Integrated_272_24.png)
    



    
![png](Integrated_files/Integrated_272_25.png)
    



    
![png](Integrated_files/Integrated_272_26.png)
    



    
![png](Integrated_files/Integrated_272_27.png)
    



    
![png](Integrated_files/Integrated_272_28.png)
    



    
![png](Integrated_files/Integrated_272_29.png)
    



    
![png](Integrated_files/Integrated_272_30.png)
    



    
![png](Integrated_files/Integrated_272_31.png)
    



    
![png](Integrated_files/Integrated_272_32.png)
    



    
![png](Integrated_files/Integrated_272_33.png)
    



    
![png](Integrated_files/Integrated_272_34.png)
    



    
![png](Integrated_files/Integrated_272_35.png)
    



    
![png](Integrated_files/Integrated_272_36.png)
    



    
![png](Integrated_files/Integrated_272_37.png)
    



    
![png](Integrated_files/Integrated_272_38.png)
    



    
![png](Integrated_files/Integrated_272_39.png)
    


    Warning message in plotEnrich(object[[group]][[db]], showTerms = showTerms, numChar = numChar, :
    “There are duplicated trimmed names in the plot, consider increasing the 'numChar' setting.”



    
![png](Integrated_files/Integrated_272_41.png)
    



    
![png](Integrated_files/Integrated_272_42.png)
    



    
![png](Integrated_files/Integrated_272_43.png)
    



    
![png](Integrated_files/Integrated_272_44.png)
    



    
![png](Integrated_files/Integrated_272_45.png)
    



    
![png](Integrated_files/Integrated_272_46.png)
    



    
![png](Integrated_files/Integrated_272_47.png)
    


**On upregulated 'Cancer vs. Control' DE genes (from `edgeR`)**


```R
fig(width = 16, height = 5)
plotEnrichR(metadata(combined)[['enrichR_edgeR_Cancer_Control_up']], db = "GO_Biological_Process_2023")
reset.fig()
```

    Warning message in plotEnrich(object[[group]][[db]], showTerms = showTerms, numChar = numChar, :
    “There are duplicated trimmed names in the plot, consider increasing the 'numChar' setting.”



    
![png](Integrated_files/Integrated_274_1.png)
    


**On downregulated 'Cancer vs. Control' DE genes (from `edgeR`)**


```R
fig(width = 16, height = 5)
plotEnrichR(metadata(combined)[['enrichR_edgeR_Cancer_Control_dn']], db = "GO_Biological_Process_2023")
reset.fig()
```


    
![png](Integrated_files/Integrated_276_0.png)
    


**On upregulated 'Cancer vs. Control' DE genes (from `DESeq2`)**


```R
fig(width = 16, height = 5)
plotEnrichR(metadata(combined)[['enrichR_DESeq2_Cancer_Control_up']], db = "GO_Biological_Process_2023")
reset.fig()
```


    
![png](Integrated_files/Integrated_278_0.png)
    



    
![png](Integrated_files/Integrated_278_1.png)
    



    
![png](Integrated_files/Integrated_278_2.png)
    



    
![png](Integrated_files/Integrated_278_3.png)
    



    
![png](Integrated_files/Integrated_278_4.png)
    



    
![png](Integrated_files/Integrated_278_5.png)
    



    
![png](Integrated_files/Integrated_278_6.png)
    


**On downregulated 'Cancer vs. Control' DE genes (from `DESeq2`)**


```R
fig(width = 16, height = 5)
plotEnrichR(metadata(combined)[['enrichR_DESeq2_Cancer_Control_dn']], db = "GO_Biological_Process_2023")
reset.fig()
```


    
![png](Integrated_files/Integrated_280_0.png)
    



    
![png](Integrated_files/Integrated_280_1.png)
    



    
![png](Integrated_files/Integrated_280_2.png)
    



    
![png](Integrated_files/Integrated_280_3.png)
    



    
![png](Integrated_files/Integrated_280_4.png)
    



    
![png](Integrated_files/Integrated_280_5.png)
    



    
![png](Integrated_files/Integrated_280_6.png)
    


### Print enrichR results to files

#### Create a sub-folder `Enrichr` to store enrichR results


```R
dir.create(file.path("Enrichr"), showWarnings = FALSE)
```

**On upregulated 'Cluster' marker genes**


```R
printEnrichR(metadata(combined)[["enrichR_findMarkers_Cluster_up"]], 
             prefix = file.path("Enrichr", paste0(file_id, "_findMarkers_upregulated")))
```

    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster1_GO_Biological_Process_2023.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster1_Reactome_Pathways_2024.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster1_WikiPathways_2024_Human.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster1_CellMarker_2024.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster2_GO_Biological_Process_2023.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster2_Reactome_Pathways_2024.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster2_WikiPathways_2024_Human.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster2_CellMarker_2024.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster3_GO_Biological_Process_2023.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster3_Reactome_Pathways_2024.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster3_WikiPathways_2024_Human.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster3_CellMarker_2024.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster4_GO_Biological_Process_2023.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster4_Reactome_Pathways_2024.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster4_WikiPathways_2024_Human.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster4_CellMarker_2024.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster5_GO_Biological_Process_2023.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster5_Reactome_Pathways_2024.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster5_WikiPathways_2024_Human.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster5_CellMarker_2024.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster6_GO_Biological_Process_2023.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster6_Reactome_Pathways_2024.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster6_WikiPathways_2024_Human.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster6_CellMarker_2024.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster7_GO_Biological_Process_2023.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster7_Reactome_Pathways_2024.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster7_WikiPathways_2024_Human.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster7_CellMarker_2024.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster8_GO_Biological_Process_2023.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster8_Reactome_Pathways_2024.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster8_WikiPathways_2024_Human.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster8_CellMarker_2024.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster9_GO_Biological_Process_2023.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster9_Reactome_Pathways_2024.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster9_WikiPathways_2024_Human.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster9_CellMarker_2024.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster10_GO_Biological_Process_2023.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster10_Reactome_Pathways_2024.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster10_WikiPathways_2024_Human.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster10_CellMarker_2024.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster11_GO_Biological_Process_2023.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster11_Reactome_Pathways_2024.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster11_WikiPathways_2024_Human.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster11_CellMarker_2024.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster12_GO_Biological_Process_2023.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster12_Reactome_Pathways_2024.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster12_WikiPathways_2024_Human.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster12_CellMarker_2024.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster13_GO_Biological_Process_2023.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster13_Reactome_Pathways_2024.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster13_WikiPathways_2024_Human.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster13_CellMarker_2024.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster14_GO_Biological_Process_2023.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster14_Reactome_Pathways_2024.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster14_WikiPathways_2024_Human.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster14_CellMarker_2024.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster15_GO_Biological_Process_2023.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster15_Reactome_Pathways_2024.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster15_WikiPathways_2024_Human.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster15_CellMarker_2024.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster16_GO_Biological_Process_2023.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster16_Reactome_Pathways_2024.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster16_WikiPathways_2024_Human.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster16_CellMarker_2024.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster17_GO_Biological_Process_2023.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster17_Reactome_Pathways_2024.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster17_WikiPathways_2024_Human.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster17_CellMarker_2024.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster18_GO_Biological_Process_2023.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster18_Reactome_Pathways_2024.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster18_WikiPathways_2024_Human.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster18_CellMarker_2024.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster19_GO_Biological_Process_2023.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster19_Reactome_Pathways_2024.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster19_WikiPathways_2024_Human.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster19_CellMarker_2024.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster20_GO_Biological_Process_2023.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster20_Reactome_Pathways_2024.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster20_WikiPathways_2024_Human.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster20_CellMarker_2024.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster21_GO_Biological_Process_2023.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster21_Reactome_Pathways_2024.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster21_WikiPathways_2024_Human.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster21_CellMarker_2024.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster22_GO_Biological_Process_2023.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster22_Reactome_Pathways_2024.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster22_WikiPathways_2024_Human.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster22_CellMarker_2024.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster23_GO_Biological_Process_2023.tsv


    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster23_Reactome_Pathways_2024.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster23_WikiPathways_2024_Human.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster23_CellMarker_2024.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster24_GO_Biological_Process_2023.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster24_Reactome_Pathways_2024.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster24_WikiPathways_2024_Human.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster24_CellMarker_2024.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster25_GO_Biological_Process_2023.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster25_Reactome_Pathways_2024.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster25_WikiPathways_2024_Human.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster25_CellMarker_2024.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster26_GO_Biological_Process_2023.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster26_Reactome_Pathways_2024.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster26_WikiPathways_2024_Human.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster26_CellMarker_2024.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster27_GO_Biological_Process_2023.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster27_Reactome_Pathways_2024.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster27_WikiPathways_2024_Human.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster27_CellMarker_2024.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster28_GO_Biological_Process_2023.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster28_Reactome_Pathways_2024.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster28_WikiPathways_2024_Human.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster28_CellMarker_2024.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster29_GO_Biological_Process_2023.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster29_Reactome_Pathways_2024.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster29_WikiPathways_2024_Human.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster29_CellMarker_2024.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster30_GO_Biological_Process_2023.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster30_Reactome_Pathways_2024.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster30_WikiPathways_2024_Human.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster30_CellMarker_2024.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster31_GO_Biological_Process_2023.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster31_Reactome_Pathways_2024.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster31_WikiPathways_2024_Human.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster31_CellMarker_2024.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster32_GO_Biological_Process_2023.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster32_Reactome_Pathways_2024.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster32_WikiPathways_2024_Human.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster32_CellMarker_2024.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster33_GO_Biological_Process_2023.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster33_Reactome_Pathways_2024.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster33_WikiPathways_2024_Human.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster33_CellMarker_2024.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster34_GO_Biological_Process_2023.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster34_Reactome_Pathways_2024.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster34_WikiPathways_2024_Human.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster34_CellMarker_2024.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster35_GO_Biological_Process_2023.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster35_Reactome_Pathways_2024.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster35_WikiPathways_2024_Human.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster35_CellMarker_2024.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster36_GO_Biological_Process_2023.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster36_Reactome_Pathways_2024.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster36_WikiPathways_2024_Human.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster36_CellMarker_2024.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster37_GO_Biological_Process_2023.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster37_Reactome_Pathways_2024.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster37_WikiPathways_2024_Human.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster37_CellMarker_2024.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster38_GO_Biological_Process_2023.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster38_Reactome_Pathways_2024.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster38_WikiPathways_2024_Human.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster38_CellMarker_2024.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster39_GO_Biological_Process_2023.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster39_Reactome_Pathways_2024.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster39_WikiPathways_2024_Human.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster39_CellMarker_2024.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster40_GO_Biological_Process_2023.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster40_Reactome_Pathways_2024.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster40_WikiPathways_2024_Human.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster40_CellMarker_2024.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster41_GO_Biological_Process_2023.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster41_Reactome_Pathways_2024.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster41_WikiPathways_2024_Human.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster41_CellMarker_2024.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster42_GO_Biological_Process_2023.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster42_Reactome_Pathways_2024.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster42_WikiPathways_2024_Human.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster42_CellMarker_2024.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster43_GO_Biological_Process_2023.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster43_Reactome_Pathways_2024.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster43_WikiPathways_2024_Human.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster43_CellMarker_2024.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster44_GO_Biological_Process_2023.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster44_Reactome_Pathways_2024.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster44_WikiPathways_2024_Human.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster44_CellMarker_2024.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster45_GO_Biological_Process_2023.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster45_Reactome_Pathways_2024.tsv


    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster45_WikiPathways_2024_Human.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster45_CellMarker_2024.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster46_GO_Biological_Process_2023.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster46_Reactome_Pathways_2024.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster46_WikiPathways_2024_Human.tsv
    Creating file: Enrichr/160k_All_findMarkers_upregulated_Cluster46_CellMarker_2024.tsv


**On upregulated and downregulated 'Cancer vs. Control' DE genes (from `edgeR`)**


```R
printEnrichR(metadata(combined)[["enrichR_edgeR_Cancer_Control_up"]], 
             prefix = file.path("Enrichr", paste0(file_id, "_Cancer_Control_edgeR_upregulated")))
```

    Creating file: Enrichr/160k_All_Cancer_Control_edgeR_upregulated_Cancer_Control_GO_Biological_Process_2023.tsv
    Creating file: Enrichr/160k_All_Cancer_Control_edgeR_upregulated_Cancer_Control_Reactome_Pathways_2024.tsv
    Creating file: Enrichr/160k_All_Cancer_Control_edgeR_upregulated_Cancer_Control_WikiPathways_2024_Human.tsv
    Creating file: Enrichr/160k_All_Cancer_Control_edgeR_upregulated_Cancer_Control_CellMarker_2024.tsv



```R
printEnrichR(metadata(combined)[["enrichR_edgeR_Cancer_Control_dn"]], 
             prefix = file.path("Enrichr", paste0(file_id, "_Cancer_Control_edgeR_downregulated")))
```

    Creating file: Enrichr/160k_All_Cancer_Control_edgeR_downregulated_Cancer_Control_GO_Biological_Process_2023.tsv
    Creating file: Enrichr/160k_All_Cancer_Control_edgeR_downregulated_Cancer_Control_Reactome_Pathways_2024.tsv
    Creating file: Enrichr/160k_All_Cancer_Control_edgeR_downregulated_Cancer_Control_WikiPathways_2024_Human.tsv
    Creating file: Enrichr/160k_All_Cancer_Control_edgeR_downregulated_Cancer_Control_CellMarker_2024.tsv


**On upregulated and downregulated 'Cancer vs. Control' DE genes (from `DESeq2`)**


```R
printEnrichR(metadata(combined)[["enrichR_DESeq2_Cancer_Control_up"]], 
             prefix = file.path("Enrichr", paste0(file_id, "_Cancer_Control_DESeq2_upregulated")))
```

    Creating file: Enrichr/160k_All_Cancer_Control_DESeq2_upregulated_CD8T_GO_Biological_Process_2023.tsv
    Creating file: Enrichr/160k_All_Cancer_Control_DESeq2_upregulated_CD8T_Reactome_Pathways_2024.tsv
    Creating file: Enrichr/160k_All_Cancer_Control_DESeq2_upregulated_CD8T_WikiPathways_2024_Human.tsv
    Creating file: Enrichr/160k_All_Cancer_Control_DESeq2_upregulated_CD8T_CellMarker_2024.tsv
    Creating file: Enrichr/160k_All_Cancer_Control_DESeq2_upregulated_CD4T_GO_Biological_Process_2023.tsv
    Creating file: Enrichr/160k_All_Cancer_Control_DESeq2_upregulated_CD4T_Reactome_Pathways_2024.tsv
    Creating file: Enrichr/160k_All_Cancer_Control_DESeq2_upregulated_CD4T_WikiPathways_2024_Human.tsv
    Creating file: Enrichr/160k_All_Cancer_Control_DESeq2_upregulated_CD4T_CellMarker_2024.tsv
    Creating file: Enrichr/160k_All_Cancer_Control_DESeq2_upregulated_gdT_GO_Biological_Process_2023.tsv
    Creating file: Enrichr/160k_All_Cancer_Control_DESeq2_upregulated_gdT_Reactome_Pathways_2024.tsv
    Creating file: Enrichr/160k_All_Cancer_Control_DESeq2_upregulated_gdT_WikiPathways_2024_Human.tsv
    Creating file: Enrichr/160k_All_Cancer_Control_DESeq2_upregulated_gdT_CellMarker_2024.tsv
    Creating file: Enrichr/160k_All_Cancer_Control_DESeq2_upregulated_NK_GO_Biological_Process_2023.tsv
    Creating file: Enrichr/160k_All_Cancer_Control_DESeq2_upregulated_NK_Reactome_Pathways_2024.tsv
    Creating file: Enrichr/160k_All_Cancer_Control_DESeq2_upregulated_NK_WikiPathways_2024_Human.tsv
    Creating file: Enrichr/160k_All_Cancer_Control_DESeq2_upregulated_NK_CellMarker_2024.tsv
    Creating file: Enrichr/160k_All_Cancer_Control_DESeq2_upregulated_B_GO_Biological_Process_2023.tsv
    Creating file: Enrichr/160k_All_Cancer_Control_DESeq2_upregulated_B_Reactome_Pathways_2024.tsv
    Creating file: Enrichr/160k_All_Cancer_Control_DESeq2_upregulated_B_WikiPathways_2024_Human.tsv
    Creating file: Enrichr/160k_All_Cancer_Control_DESeq2_upregulated_B_CellMarker_2024.tsv
    Creating file: Enrichr/160k_All_Cancer_Control_DESeq2_upregulated_Monocytes_GO_Biological_Process_2023.tsv
    Creating file: Enrichr/160k_All_Cancer_Control_DESeq2_upregulated_Monocytes_Reactome_Pathways_2024.tsv
    Creating file: Enrichr/160k_All_Cancer_Control_DESeq2_upregulated_Monocytes_WikiPathways_2024_Human.tsv
    Creating file: Enrichr/160k_All_Cancer_Control_DESeq2_upregulated_Monocytes_CellMarker_2024.tsv
    Creating file: Enrichr/160k_All_Cancer_Control_DESeq2_upregulated_DCs_GO_Biological_Process_2023.tsv
    Creating file: Enrichr/160k_All_Cancer_Control_DESeq2_upregulated_DCs_Reactome_Pathways_2024.tsv
    Creating file: Enrichr/160k_All_Cancer_Control_DESeq2_upregulated_DCs_WikiPathways_2024_Human.tsv
    Creating file: Enrichr/160k_All_Cancer_Control_DESeq2_upregulated_DCs_CellMarker_2024.tsv



```R
printEnrichR(metadata(combined)[["enrichR_DESeq2_Cancer_Control_dn"]], 
             prefix = file.path("Enrichr", paste0(file_id, "_Cancer_Control_DESeq2_downregulated")))
```

    Creating file: Enrichr/160k_All_Cancer_Control_DESeq2_downregulated_CD8T_GO_Biological_Process_2023.tsv
    Creating file: Enrichr/160k_All_Cancer_Control_DESeq2_downregulated_CD8T_Reactome_Pathways_2024.tsv
    Creating file: Enrichr/160k_All_Cancer_Control_DESeq2_downregulated_CD8T_WikiPathways_2024_Human.tsv
    Creating file: Enrichr/160k_All_Cancer_Control_DESeq2_downregulated_CD8T_CellMarker_2024.tsv
    Creating file: Enrichr/160k_All_Cancer_Control_DESeq2_downregulated_CD4T_GO_Biological_Process_2023.tsv
    Creating file: Enrichr/160k_All_Cancer_Control_DESeq2_downregulated_CD4T_Reactome_Pathways_2024.tsv
    Creating file: Enrichr/160k_All_Cancer_Control_DESeq2_downregulated_CD4T_WikiPathways_2024_Human.tsv
    Creating file: Enrichr/160k_All_Cancer_Control_DESeq2_downregulated_CD4T_CellMarker_2024.tsv
    Creating file: Enrichr/160k_All_Cancer_Control_DESeq2_downregulated_gdT_GO_Biological_Process_2023.tsv
    Creating file: Enrichr/160k_All_Cancer_Control_DESeq2_downregulated_gdT_Reactome_Pathways_2024.tsv
    Creating file: Enrichr/160k_All_Cancer_Control_DESeq2_downregulated_gdT_WikiPathways_2024_Human.tsv
    Creating file: Enrichr/160k_All_Cancer_Control_DESeq2_downregulated_gdT_CellMarker_2024.tsv
    Creating file: Enrichr/160k_All_Cancer_Control_DESeq2_downregulated_NK_GO_Biological_Process_2023.tsv
    Creating file: Enrichr/160k_All_Cancer_Control_DESeq2_downregulated_NK_Reactome_Pathways_2024.tsv
    Creating file: Enrichr/160k_All_Cancer_Control_DESeq2_downregulated_NK_WikiPathways_2024_Human.tsv
    Creating file: Enrichr/160k_All_Cancer_Control_DESeq2_downregulated_NK_CellMarker_2024.tsv
    Creating file: Enrichr/160k_All_Cancer_Control_DESeq2_downregulated_B_GO_Biological_Process_2023.tsv
    Creating file: Enrichr/160k_All_Cancer_Control_DESeq2_downregulated_B_Reactome_Pathways_2024.tsv
    Creating file: Enrichr/160k_All_Cancer_Control_DESeq2_downregulated_B_WikiPathways_2024_Human.tsv
    Creating file: Enrichr/160k_All_Cancer_Control_DESeq2_downregulated_B_CellMarker_2024.tsv
    Creating file: Enrichr/160k_All_Cancer_Control_DESeq2_downregulated_Monocytes_GO_Biological_Process_2023.tsv
    Creating file: Enrichr/160k_All_Cancer_Control_DESeq2_downregulated_Monocytes_Reactome_Pathways_2024.tsv
    Creating file: Enrichr/160k_All_Cancer_Control_DESeq2_downregulated_Monocytes_WikiPathways_2024_Human.tsv
    Creating file: Enrichr/160k_All_Cancer_Control_DESeq2_downregulated_Monocytes_CellMarker_2024.tsv
    Creating file: Enrichr/160k_All_Cancer_Control_DESeq2_downregulated_DCs_GO_Biological_Process_2023.tsv
    Creating file: Enrichr/160k_All_Cancer_Control_DESeq2_downregulated_DCs_Reactome_Pathways_2024.tsv
    Creating file: Enrichr/160k_All_Cancer_Control_DESeq2_downregulated_DCs_WikiPathways_2024_Human.tsv
    Creating file: Enrichr/160k_All_Cancer_Control_DESeq2_downregulated_DCs_CellMarker_2024.tsv


## Ingenuity Pathway Analysis (IPA)

Here we generate a file that can be imported directly into IPA for downstream analysis.

### Concatenate `findMarkers` 'Cluster' marker gene results


```R
exportResList(metadata(combined)[['findMarkers_Cluster']], concatenate = TRUE, col_anno = c("ID","Symbol"), 
              prefix = paste0(file_id, "_cluster"))
```

    Detecting findMarkers input.
    Creating a concatenated file: 160k_All_cluster_findMarkers_concatenated.tsv


### Concatenate `edgeR` 'Cancer vs. Control' DE gene results


```R
exportResList(metadata(combined)[['edgeR_Cancer_Control']], concatenate = TRUE, col_anno = c("ID","Symbol"), 
              prefix = paste0(file_id, "_Cancer_Control"))
```

    Detecting edgeR input.
    Creating a concatenated file: 160k_All_Cancer_Control_edgeR_concatenated.tsv


### Concatenate `DESeq2` 'Cancer vs. Control' DE gene results


```R
exportResList(metadata(combined)[['DESeq2_Cancer_Control']], concatenate = TRUE, col_anno = c("ID","Symbol"), 
              prefix = paste0(file_id, "_Cancer_Control"))
```

    Detecting DESeq2 input.
    Creating a concatenated file: 160k_All_Cancer_Control_DESeq2_concatenated.tsv


# Save `runInfo` to `metadata`

<div class="alert alert-info">
    <strong>Tip!</strong> In the accompanied Shiny App, the Run Information will be displayed under the <u>Overview</u>.
</div>


```R
metadata(combined)[['runInfo']] <- runInfo
combined
```


    class: SingleCellExperiment 
    dim: 18129 48494 
    metadata(35): Control1_Samples Control1_cyclone ... enrichR_DESeq2_Cancer_Control_dn runInfo
    assays(2): counts logcounts
    rownames(18129): SAMD11 NOC2L ... MT-ND6 MT-CYB
    rowData names(5): ID Symbol Type SEQNAME is_mito
    colnames(48494): Control1_AAACAAGCAACAAGTTACTTTAGG-1 Control1_AAACAAGCAACTAGTGACTTTAGG-1 ...
      LungCancer_TTTGTGAGTTGAGTCTAGCTGTGA-1 LungCancer_TTTGTGAGTTTACGACAGCTGTGA-1
    colData names(20): Sample Barcode ... ClusterCellType CellType_2
    reducedDimNames(6): PCA TSNE ... MNN-TSNE MNN-UMAP
    mainExpName: Gene Expression
    altExpNames(1): Antibody Capture


# Save objects


```R
# Remove unwanted information stored in the SingleCellExperiment before saving the object
# For example, this removed "PCA" from the reducedDims slot
# reducedDim(combined, "PCA") <- NULL
# For example, this removed "MNN" from the reducedDims slot
# reducedDim(combined, "MNN") <- NULL
# For example, this removed "reconstructed" from the assays slot
# assay(combined, "reconstructed") <- NULL
```


```R
# For HDF5-based SummarizedExperiment object
HDF5Array::saveHDF5SummarizedExperiment(combined, dir = "combined_h5_sce", replace = TRUE, verbose = FALSE)

# Print file size
"Folder: combined_h5_sce"
paste("Size:", utils:::format.object_size(file.info("combined_h5_sce/assays.h5")$size + 
                                          file.info("combined_h5_sce/se.rds")$size, "auto"))
```


'Folder: combined_h5_sce'



'Size: 2.4 Gb'


# Session Info


```R
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
sessionInfo()
```


    R version 4.4.3 (2025-02-28)
    Platform: x86_64-conda-linux-gnu
    Running under: Ubuntu 20.04.6 LTS
    
    Matrix products: default
    BLAS/LAPACK: /home/ihsuan/miniconda3/envs/github_sc/lib/libopenblasp-r0.3.29.so;  LAPACK version 3.12.0
    
    locale:
     [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C               LC_TIME=en_GB.UTF-8       
     [4] LC_COLLATE=en_GB.UTF-8     LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8   
     [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
    [10] LC_TELEPHONE=C             LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       
    
    time zone: Europe/London
    tzcode source: system (glibc)
    
    attached base packages:
    [1] grid      stats4    stats     graphics  grDevices utils     datasets  methods   base     
    
    other attached packages:
     [1] kableExtra_1.4.0            lubridate_1.9.4             forcats_1.0.0              
     [4] stringr_1.5.1               dplyr_1.1.4                 purrr_1.0.4                
     [7] readr_2.1.5                 tidyr_1.3.1                 tibble_3.2.1               
    [10] tidyverse_2.0.0             scRUtils_0.3.8              viridis_0.6.5              
    [13] viridisLite_0.4.2           scran_1.34.0                scater_1.34.1              
    [16] scuttle_1.16.0              scales_1.3.0                pheatmap_1.0.12            
    [19] ggforce_0.4.2               ggplot2_3.5.2               enrichR_3.4                
    [22] edgeR_4.4.2                 limma_3.62.2                DESeq2_1.46.0              
    [25] cowplot_1.1.3               bluster_1.16.0              BiocParallel_1.40.2        
    [28] BiocNeighbors_2.0.1         batchelor_1.22.0            SingleCellExperiment_1.28.1
    [31] SummarizedExperiment_1.36.0 Biobase_2.66.0              GenomicRanges_1.58.0       
    [34] GenomeInfoDb_1.42.3         IRanges_2.40.1              S4Vectors_0.44.0           
    [37] BiocGenerics_0.52.0         MatrixGenerics_1.18.1       matrixStats_1.5.0          
    
    loaded via a namespace (and not attached):
      [1] RColorBrewer_1.1-3        rstudioapi_0.17.1         jsonlite_2.0.0            shape_1.4.6.1            
      [5] magrittr_2.0.3            ggbeeswarm_0.7.2          rmarkdown_2.29            farver_2.1.2             
      [9] fs_1.6.6                  GlobalOptions_0.1.2       zlibbioc_1.52.0           vctrs_0.6.5              
     [13] Cairo_1.6-2               DelayedMatrixStats_1.28.1 base64enc_0.1-3           htmltools_0.5.8.1        
     [17] S4Arrays_1.6.0            curl_6.2.2                Rhdf5lib_1.28.0           gridGraphics_0.5-1       
     [21] rhdf5_2.50.2              SparseArray_1.6.2         ResidualMatrix_1.16.0     uuid_1.2-1               
     [25] igraph_2.1.4              lifecycle_1.0.4           iterators_1.0.14          pkgconfig_2.0.3          
     [29] rsvd_1.0.5                Matrix_1.7-3              R6_2.6.1                  fastmap_1.2.0            
     [33] GenomeInfoDbData_1.2.13   clue_0.3-66               digest_0.6.37             colorspace_2.1-1         
     [37] ggnewscale_0.5.1          dqrng_0.4.1               irlba_2.3.5.1             beachmat_2.22.0          
     [41] labeling_0.4.3            WriteXLS_6.7.0            timechange_0.3.0          httr_1.4.7               
     [45] polyclip_1.10-7           abind_1.4-8               compiler_4.4.3            withr_3.0.2              
     [49] doParallel_1.0.17         maps_3.4.2.1              HDF5Array_1.34.0          MASS_7.3-65              
     [53] DelayedArray_0.32.0       rjson_0.2.23              gtools_3.9.5              tools_4.4.3              
     [57] vipor_0.4.7               beeswarm_0.4.0            glmGamPoi_1.18.0          glue_1.8.0               
     [61] rhdf5filters_1.18.1       Rtsne_0.17                pbdZMQ_0.3-14             cluster_2.1.8.1          
     [65] generics_0.1.3            gtable_0.3.6              tzdb_0.5.0                hms_1.1.3                
     [69] xml2_1.3.8                BiocSingular_1.22.0       ScaledMatrix_1.14.0       metapod_1.14.0           
     [73] XVector_0.46.0            RcppAnnoy_0.0.22          ggrepel_0.9.6             foreach_1.5.2            
     [77] pillar_1.10.2             yulab.utils_0.2.0         pals_1.10                 IRdisplay_1.1            
     [81] circlize_0.4.16           tweenr_2.0.3              lattice_0.22-7            tidyselect_1.2.1         
     [85] ComplexHeatmap_2.22.0     locfit_1.5-9.12           knitr_1.50                gridExtra_2.3            
     [89] svglite_2.1.3             xfun_0.52                 statmod_1.5.0             stringi_1.8.7            
     [93] UCSC.utils_1.2.0          evaluate_1.0.3            codetools_0.2-20          ggplotify_0.1.2          
     [97] cli_3.6.4                 uwot_0.2.3                IRkernel_1.3.2            systemfonts_1.2.2        
    [101] repr_1.1.7                munsell_0.5.1             dichromat_2.0-0.1         Rcpp_1.0.14              
    [105] mapproj_1.2.11            png_0.1-8                 parallel_4.4.3            sparseMatrixStats_1.18.0 
    [109] crayon_1.5.3              GetoptLong_1.0.5          rlang_1.1.6              

