# BCF Single Cell Gene Expression Shiny App

This repository holds the Shiny application tailored to work with `SingleCellExperiment` objects (in HDF5 format) created after running the UoM-BCF single-cell gene expression data analysis notebooks.

## Prerequisite

1. The App requires R version >=4.4 and Bioconductor version 3.20.
2. Required R packages (other dependencies not shown):

| Package | Source |
| --- | --- |
| cowplot | CRAN |
| dplyr | CRAN |
| DT | CRAN |
| edgeR | Bioconductor |
| fontawesome | CRAN |
| gtools | CRAN |
| HDF5Array | Bioconductor |
| htmlwidgets | CRAN |
| pals | CRAN |
| pheatmap | CRAN |
| plotly | CRAN |
| scater | Bioconductor |
| shiny | CRAN |
| shinycssloaders | CRAN |
| shinydashboard | CRAN |
| shinyWidgets | CRAN |
| tidyr | CRAN |

3. Optional: DESeq2 (if DESeq2 is used and `DESeqResults` objects are stored as part of the analysis)

Set up a Shiny environment using `mamba` (or `conda`):

```
mamba create -n shiny -c conda-forge -c bioconda r-base=4.4 r-shiny r-shinycssloaders \
r-shinywidgets r-shinydashboard r-dt r-dplyr r-plotly bioconductor-scater r-cowplot \
bioconductor-edger bioconductor-hdf5array r-htmlwidgets r-pals r-tidyr r-pheatmap  \
r-gtools r-fontawesome
```

To also install DESeq2:

```
mamba install -c bioconda bioconductor-deseq2
```

## How to configured the App

Open the file `app.R` in an editor and change the following variables to that applicable to your dataset.

The `filepaths` variable is pre-configured with paths of the demo dataset 
that will be generated if running the Jupyter Notebooks inside the `notebook` folder.

``` r
#--------- Adjustment required; START -------- #
  
# Set the paths to the directories where the HDF5-based object were saved
filepaths <- c(Combined         = "160k_Human_Integrated/combined_h5_sce",
               Control1         = "160k_Human_Control1_PBMC/Control1_h5_sce",
               "Kidney Cancer"  = "160k_Human_Kidney_Cancer_PBMC/KidneyCancer_h5_sce",
               "Lung Cancer"    = "160k_Human_Lung_Cancer_PBMC/LungCancer_h5_sce")

# Set the default selected reducedDimName in the dropdown menu, accessible from `reducedDimNames(sce)`
# If not found, TSNE is selected
default.dimred <- "MNN-TSNE"

# Add custom cell features (and their descriptions) for colouring, accessible from `colData(sce)`
# Built-in recognised feature names: "Sample", "sum", "detected", "subsets_Mt_percent", "CellCycle",
# "DoubletDensity", "DoubletDensity_log1p", "label", "CellType", "ClusterCellType"
my.color_by <- c("condition")
my.color_by.desc <- c("Condition")

# Add custom cell features (and their descriptions) for grouped presentation
# (typically in the dotplots, boxplots and heatmap), accessible from `colData(sce)`
# Built-in recognised feature names: "Sample","label","CellType","ClusterCellType"
my.group_by <- c("condition")
my.group_by.desc <- c("Condition")

# Add custom clustering methods (and their descriptions), accessible from `colData(sce)`
# Built-in recognised methods: "hclust", "walktrap", "louvain", "leiden", 
# or "merged.hclust", "merged.walktrap", "merged.louvain", "merged.leiden" when the "merged.*" is found in `colData(sce)`
my.cluster.methods <- c("merged.custom")
my.cluster.methods.desc <- c("Custom method")

# Set the maximum allowable genes in multi-gene plots
# The more genes used, more time the App takes to create the plots
max.gene <- 100

# Define a vector containing gene symbols to be used as example in multi-gene plots, for example highly variable genes (HVGs)
# The symbols has to be identifiable in `rownames(sce)`, only the first `max.gene` number of genes are used in the App
# If none is given, the App will first use the HVGs stored in `metadata(sce)[["runInfo"]][["HVG"]][["Genes"]]`
# If this data structure cannot be found in `metadata(sce)`, the App will randomly select `max.gene` number of genes from `rownames(sce)`
my.example.genes <- c()

# Set the maximum number of cell features to allow App users to choose and combine them into a new feature to group cells in multi-gene plots
# Beware, by allowing App users to select more cell features, the more unique groups there will be in the final plots
# For example, when `multi_max_options` is 2 and App user chooses "Sample" (N=2) and "label" (N=10), there will be 20 groups
multi_max_options <- 2

#--------- Adjustment required; END -------- #
```

## How to run the App

Run a Shiny App using the command `runApp`. Stop the App by interrupt R, i.e. pressing Ctrl+C or Esc)

- Use `port=` to specify the TCP port that the application should listen on.
- Use `host =` to specify the IPv4 address that the application should listen on.

To run it in R console:

```r
shiny::runApp('.')

shiny::runApp('.', host = '127.0.0.1', port = 7010)
```

To run it in terminal/console:

```r
R -e "shiny::runApp('.')"

R -e "shiny::runApp('.', host = '127.0.0.1', port = 7010)"
```

