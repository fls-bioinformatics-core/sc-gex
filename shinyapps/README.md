# BCF Single Cell Gene Expression Shiny App

This repository holds the Shiny application tailored to work with `SingleCellExperiment` objects (in HDF5 format) created after running the BCF single-cell gene expression single-sample and integrated notebooks.

## Prerequisite

1. The App requires R version >=4.1 and Bioconductor version 3.13.
2. Required R packages (other dependencies not shown):

| Package | Source |
| --- | --- |
| cowplot | CRAN |
| dplyr | CRAN |
| DT | CRAN |
| edgeR | Bioconductor |
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

Set up a Shiny environment using `conda`:

```
conda create -n shiny -c conda-forge -c bioconda r-base=4.1 r-shiny r-shinycssloaders r-shinywidgets r-shinydashboard r-dt r-dplyr r-plotly bioconductor-scater r-cowplot bioconductor-edger bioconductor-hdf5array r-htmlwidgets r-pals r-tidyr r-pheatmap r-gtools
```

## Demo site showing six 10x Genomics snRNA-seq datasets

*Accessible only within University of Manchester's network (or with UoM VPN)*

- 1k Mouse Kidney
- 5k Adult Mouse Brain
- 5k Adult Mouse Heart
- 5k Adult Mouse Liver
- 5k Adult Mouse Lung
- 5k Human Jejunum

http://doublecell-1.ls.manchester.ac.uk:7169/

## How to run an App

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

## Opening ports

To open port allowing access to Shiny App on BCF's doublecell VMs:

```
sudo ufw allow proto tcp to any port PORT_NUM
``
