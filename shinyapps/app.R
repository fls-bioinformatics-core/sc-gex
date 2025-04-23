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

####################
# App info and settings
####################

app.version <- "v0.7.8"
app.header <- "BCF Single Cell GEX"
app.title <- "BCF Single Cell Gene Expression Shiny App"
app.author <- "I-Hsuan Lin [Author, Creator], Syed Murtuza baker [Contributor]"
app.license <- " GPL-3.0"
app.description <- paste0("The <em>", app.title, "</em> is a user friendly interface allowing users to explore processed single-cell RNA-seq datasets. It is designed to work with <code>SingleCellExperiment</code> objects that have been processed with the BCF in-house single cell gene expression workflow.")

sidebar_width <- 210
base_size <- 14
print_log <- TRUE

####################
# Additional manipulations (run once when app is launched)
####################
# Build continuous palette list
continuous.n <- 25
continuous <- list(
  yellowGreen   = pals::brewer.ylgn(continuous.n),
  yellowBlue    = pals::brewer.ylgnbu(continuous.n),
  yellowRed     = pals::brewer.ylorrd(continuous.n),
  blues		= pals::brewer.blues(continuous.n),
  greens	= pals::brewer.greens(continuous.n),
  oranges	= pals::brewer.oranges(continuous.n),
  purples	= pals::brewer.purples(continuous.n),
  reds		= pals::brewer.reds(continuous.n),
  viridis       = pals::viridis(continuous.n),
  coolwarm      = pals::coolwarm(continuous.n),
  warmcool      = pals::warmcool(continuous.n),
  cubehelix	= pals::cubehelix(continuous.n),
  gnuplot	= pals::gnuplot(continuous.n), 
  jet		= pals::jet(continuous.n), 
  parula	= pals::parula(continuous.n), 
  cividis	= pals::cividis(continuous.n), 
  tol.rainbow   = pals::tol.rainbow(continuous.n), 
  turbo 	= pals::turbo(continuous.n)
)

discrete <- c("Accent","Dark2","Paired","Pastel1","Pastel2","Set1","Set2","Set3")

# Build discrete palette names
discrete_func <- function(name, n) {
  nc <- ifelse(n > 2, n, 3) # minimal value for n is 3
  if(name == "Accent") colors <- pals::brewer.accent(nc) 
  else if(name == "Paired") colors <- pals::brewer.paired(nc) 
  else if(name == "Pastel1") colors <- pals::brewer.pastel1(nc) 
  else if(name == "Pastel2") colors <- pals::brewer.pastel2(nc) 
  else if(name == "Set1") colors <- pals::brewer.set1(nc) 
  else if(name == "Set2") colors <- pals::brewer.set2(nc) 
  else if(name == "Set3") colors <- pals::brewer.set3(nc)
  else colors <-  pals::brewer.dark2(nc) # default palette

  head(colors, n) # returns the correct number of colours
}

# Build a vector of "color_by" elements
default.color_by <- c(
  "Sample", "sum", "detected", "subsets_Mt_percent", "CellCycle", "DoubletDensity", 
  "DoubletDensity_log1p", "label", "CellType", "ClusterCellType"
)

default.color_by.desc <- c(
  "Sample", "Library size", "Detected genes", "Percent Mito", "Cell cycle", "Doublet density score", 
  "Doublet density score (log1p)", "Cluster label", "Cell type", "Cluster-specific cell type"
)
names(default.color_by) <- default.color_by.desc

# Append custom cell features for colouring
if(length(my.color_by > 0)) {
  names(my.color_by) <- my.color_by.desc
  default.color_by <- c(default.color_by[!default.color_by %in% my.color_by], my.color_by)
}
default.color_by <- default.color_by[order(names(default.color_by))]

# Build a vector of "cluster.methods" elements
default.cluster.methods <- c("hclust", "walktrap", "louvain", "leiden")
default.cluster.methods.desc <- paste(c("Hclust", "Walktrap", "Louvain", "Leiden"), "method")

# Build a vector of "group_by" elements (typically in the dot/bar/box/heatmap)
default.group_by <- c("Sample","label","CellType","ClusterCellType")
default.group_by <- default.color_by[default.color_by %in% default.group_by]

# Append custom cell features for grouping
if(length(my.group_by > 0)) {
  names(my.group_by) <- my.group_by.desc
  default.group_by <- c(default.group_by[!default.group_by %in% my.group_by], my.group_by)
}
default.group_by <- default.group_by[order(names(default.group_by))]

# plotly legend for Cell features projection
pl.legend <- list(font = list(size = 12, color = "#000"), bordercolor = "#000", borderwidth = 1, itemsizing = "constant")

####################
# Load R libraries
####################
suppressPackageStartupMessages({
  # Load shiny packages
  library(shiny)
  library(shinycssloaders) # withSpinner
  library(shinyWidgets)
  library(shinydashboard)
  library(fontawesome)

  # Load additional packages
  library(DT)
  library(dplyr)
  library(plotly)
  library(scater)

  # Not loaded, use :: to call required functions
  #library(cowplot)     # theme_cowplot, plot_grid
  #library(edgeR)       # topTags
  #library(HDF5Array)   # loadHDF5SummarizedExperiment
  #library(htmlwidgets) # JS
  #library(pals)        # continuous color palettes
  #library(tidyr)       # gather
})

####################
# Load sce files
####################
# Check files/paths
fileExists <- file.exists(filepaths)
if(length(filepaths[!fileExists]) > 0) warning(sprintf("File/path(s) not exists: %s", filepaths[!fileExists]))

# Subset to available files/paths
filepaths <- filepaths[fileExists]
if(length(filepaths) == 0) stop("None of the file/path(s) in 'filepaths' is valid.")

sce.list <- list()
for(file in names(filepaths)) {
  filepath <- filepaths[file]
  if(file.info(filepath)$isdir) {
    if(file.exists(file.path(filepath, "se.rds")) & file.exists(file.path(filepath, "assays.h5"))) {
      print(sprintf("Loading '%s' from '%s'.", file, filepath))
      sce.list[[file]] <- HDF5Array::loadHDF5SummarizedExperiment(dir = filepath)
    } else warning(sprintf("HDF5 file '%s' skipped, 'se.rds' and 'assays.h5' is missing.", filepath))
  } else {
    if(grepl("\\.rds|\\.RDS", filepath)) {
      print(sprintf("Loading '%s' from '%s'.", file, filepath))
      sce.list[[file]] <- readRDS(filepath)
    } else warning(sprintf("RDS file '%s' skipped, requires RDS file with a '.rds' or '.RDS' extension.", filepath))
  }
}
init.sce <- sce.list[[1]]

####################
# Build user interface (ui) object
####################
ui <- dashboardPage(
  header = dashboardHeader(title = app.header),
  sidebar = dashboardSidebar(
    # Remove the sidebar toggle element
    tags$script(htmlwidgets::JS("document.getElementsByClassName('sidebar-toggle')[0].style.visibility = 'hidden';")),
    width = sidebar_width,
    sidebarMenuOutput("menuItms")
  ),
  body = dashboardBody(
    tags$head(tags$style(HTML("
/* logo */
.skin-blue .main-header .logo{
  font-size:20px;
  color:#ffd500;
  background-color: #6b2c91
}

/* logo when hovered */
.skin-blue .main-header .logo:hover{
  background-color: #6b2c91
}

/* navbar (rest of the header) */
.skin-blue .main-header .navbar{
  background-color: #6b2c91
}

/* other links in the sidebarmenu when hovered */
.skin-blue .main-sidebar .sidebar .sidebar-menu a:hover{
  background-color: #6b2c91
}

/* toggle button when hovered  */
.skin-blue .main-header .navbar .sidebar-toggle:hover{
  background-color: #6b2c91
}

/* sliderInput font size */
.irs-grid-text{
  font-size: 10px
}


/* Change primary status box */
.box.box-primary{
  border-top-color: #6b2c91
}

.box.box-solid.box-primary{
  border: 1px solid #6b2c91
}

.box.box-solid.box-primary>.box-header{
  color: #ffd500;
  background: #6b2c91;
  background-color: #6b2c91
}

.box.box-solid.box-primary>.box-header .btn,.box.box-solid.box-primary>.box-header a{
  color: #ffd500
}

.box.box-primary h3{
  font-size:20px;
  font-weight:500;
  line-height:1;
  display:inline-block;
  margin:0
}

/* show header text in single line */
table.dataTable thead th {
  white-space: nowrap
}

/* allow multi-line verbatimTextOutput() */
.not_found pre {
  white-space: pre-wrap;
  overflow-y:scroll;
  max-height: 250px;
}

.not_expr pre {
  white-space: pre-wrap;
  overflow-y:scroll;
  max-height: 150px;
}
"))),
    uiOutput("dynamicItems")
  ),
#  footer = dashboardFooter(
#    left = footer.left,
#    right = footer.right
#  ),
  title = paste0(app.header, " (",  basename(getwd()), ")")
)

####################
# Build server function
####################
server <- function(input, output, session) {
  # Empty ggplot
  empty <- ggplot() + theme_void()

  # Calculate dot size
  dot_size <- function(n) {
    size <- ceiling((-1.303 * log(n) + 18.097) / 0.5) * 0.5
    ifelse(size < 1, 1, ifelse(size > 15, 15, size)) # min 1, max 15
  }

  ####################
  # Reset submit counter when switch sce files
  ####################
  reset.cell.submit <- reset.gene.submit <- reset.multi.submit <- reactiveVal(1)
  reset.ORAfm.submit <- reset.ORAer.submit <- reset.ORAde.submit <- reactiveVal(1)

  observeEvent(input$cell.submit, { reset.cell.submit(0) })
  observeEvent(input$gene.submit, { reset.gene.submit(0) })
  observeEvent(input$multi.submit, { reset.multi.submit(0) })
  observeEvent(input$ORAfm.submit, { reset.ORAfm.submit(0) })
  observeEvent(input$ORAer.submit, { reset.ORAer.submit(0) })
  observeEvent(input$ORAde.submit, { reset.ORAde.submit(0) })

  observeEvent(input$sce.file, {
    reset.cell.submit(1)
    reset.gene.submit(1)
    reset.multi.submit(1)
    reset.ORAfm.submit(1)
    reset.ORAer.submit(1)
    reset.ORAde.submit(1)
    updateTextAreaInput(session, "multi.ganenames", value = "")
  })

  ####################
  # Reset submit counter when viewing ORAfm, ORAer, ORAde sub-menu
  ####################
  observeEvent(input$menuItems, {
    if(input$menuItems %in% c("ORAfm","ORAer","ORAde")) {
      reset.ORAfm.submit(1)
      reset.ORAer.submit(1)
      reset.ORAde.submit(1)
    }
  })

  ####################
  # Choose sce file
  ####################
  # initialize reactive value as first sce
  choose_sce <- reactiveVal(names(filepaths)[1])

  # Stop using the default
  observeEvent(input$sce.file, ignoreInit = F, {
    choose_sce(input$sce.file)
  })

  # Print access logs
  observeEvent(input$menuItems, {
    if(print_log) {
      clint_ip <- ifelse(!is.null(session$request$HTTP_X_FORWARDED_FOR), session$request$HTTP_X_FORWARDED_FOR, session$request$REMOTE_ADDR)
      print(sprintf("'%s' is viewing '%s' in '%s' on %s", clint_ip, input$menuItems, input$sce.file, format(Sys.time(), usetz = TRUE)))
    }
  })

  ####################
  # Load data
  ####################
  sce <- reactive({
    sce <- sce.list[[choose_sce()]]

    # 10X vs. Parse Biosciences names
    if("sample" %in% colnames(colData(sce)) & !"Sample" %in% colnames(colData(sce))) sce$Sample <- sce$sample
    if(!"Barcode" %in% colnames(colData(sce))) sce$Barcode <- colnames(sce)

    # Encode vector as factor
    colLabels(sce) <- if(is.factor(colLabels(sce))) droplevels(colLabels(sce)) else as.factor(colLabels(sce))
    sce$Sample <- if(is.factor(sce$Sample)) droplevels(sce$Sample) else as.factor(sce$Sample)
    sce$CellType <- if("CellType" %in% colnames(colData(sce))) {
      if(is.factor(sce$CellType)) droplevels(sce$CellType) else as.factor(sce$CellType)
    } else colLabels(sce)
    sce$ClusterCellType <- if("ClusterCellType" %in% colnames(colData(sce))) { 
      if(is.factor(sce$ClusterCellType)) droplevels(sce$ClusterCellType) else as.factor(sce$ClusterCellType)
    } else colLabels(sce)

    # Add log1p converted scores
    if("DoubletDensity" %in% colnames(colData(sce))) sce$DoubletDensity_log1p <- log1p(sce$DoubletDensity)
    sce
  })

  ####################
  # Build sidebarMenu
  ####################
  output$menuItms <- renderMenu({
    sce <- if(is.null(input$sce.file)) init.sce else sce()

    # findMarkers
    fm.listnames <- names(metadata(sce))[grep("^findMarkers", names(metadata(sce)))]
    if(length(fm.listnames) > 0) {
      fm.menu <- menuItem("Gene markers", tabName = "findMarkers", icon = fa_i("thumbtack", fill = "ffd500"), startExpanded = TRUE,
        lapply(fm.listnames, function(listname) {
          tab.name <- gsub("_", ": ", gsub("findMarkers_", "", listname))
          menuSubItem(tab.name, tabName = listname)
        })
      )
    } else fm.menu <- NULL

    # edgeR
    edger.listnames <- names(metadata(sce))[grep("^edgeR", names(metadata(sce)))]
    if(length(edger.listnames) > 0) {
      edger.menu <- menuItem("DEA (edgeR)", tabName = "edgeR", icon = fa_i("magnifying-glass-chart", fill = "ffd500"), startExpanded = TRUE,
        lapply(edger.listnames, function(listname) {
          tab.name <- gsub("_", " vs. ", gsub("edgeR_", "", listname))
          menuSubItem(tab.name, tabName = listname)
        })
      )
    } else edger.menu <- NULL

    # DESeq2
    deseq.listnames <- names(metadata(sce))[grep("^DESeq2", names(metadata(sce)))]
    if(length(deseq.listnames) > 0) {
      deseq.menu <- menuItem("DEA (DESeq2)", tabName = "DESeq2", icon = fa_i("magnifying-glass-chart", fill = "ffd500"), startExpanded = TRUE,
        lapply(deseq.listnames, function(listname) {
          tab.name <- gsub("_", " vs. ", gsub("DESeq2_", "", listname))
          menuSubItem(tab.name, tabName = listname)
        })
      )
    } else deseq.menu <- NULL

    # enrichR (findMarkers/edgeR/DESeq2)
    enrichr.listnames <- names(metadata(sce))[grep("^enrichR", names(metadata(sce)))]
    if(length(enrichr.listnames) > 0) {
      enrichr.menu <- menuItem("Enrichment analysis", tabName = "enrichR", icon = fa_i("magnifying-glass-chart", fill = "ffd500"), startExpanded = TRUE,
                               if(sum(unlist(lapply(metadata(sce)[enrichr.listnames[(grepl("findMarkers", enrichr.listnames))]], length))) > 0)
                                       do.call(tagList, list(menuSubItem("Gene markers", tabName = "ORAfm"))),
                               if(sum(unlist(lapply(metadata(sce)[enrichr.listnames[(grepl("edgeR", enrichr.listnames))]], length))) > 0)
                                       do.call(tagList, list(menuSubItem("DEA (edgeR)", tabName = "ORAer"))),
                               if(sum(unlist(lapply(metadata(sce)[enrichr.listnames[(grepl("DESeq2", enrichr.listnames))]], length))) > 0)
                                       do.call(tagList, list(menuSubItem("DEA (DESeq2)", tabName = "ORAde")))
      )
    } else enrichr.menu <- NULL

    menu <- c(
      list(id = "menuItems"),
      list(menuItem(paste0("Overview: [", choose_sce(), "]"), tabName = "overview", icon = fa_i("house", fill = "ffd500"))),
      list(menuItem("Cell feature", tabName = "cellFeatures", icon = fa_i("dna", fill = "ffd500"))),
      list(menuItem("Gene expression", tabName = "geneExpression", icon = fa_i("signal", fill = "ffd500"), startExpanded = TRUE,
                    menuSubItem("Single gene", tabName = "onegeneExpression"),
                    menuSubItem("Multiple genes", tabName = "multigeneExpression"))),
      list(fm.menu),
      list(edger.menu),
      list(deseq.menu),
      list(enrichr.menu),
      list(menuItem("About", tabName = "about", icon = fa_i("circle-info", fill = "ffd500")))
    )
    do.call(sidebarMenu, menu)
  })

  ####################
  # Build tabItems
  ####################
  output$dynamicItems <- renderUI({
    items <- c(
      ####################
      # Overview panel
      ####################
      list(tabItem(tabName = "overview",
        fluidRow(
          column(width = 12, radioGroupButtons(inputId = "sce.file", label = "Showing:", direction = "horizontal", choices = names(filepaths), selected = isolate(input$sce.file)))
        ),
        fluidRow(
          uiOutput("overview.ui"),
          uiOutput("overview.plot.ui"),
          uiOutput("runinfo.ui")
        )
      )),
      ####################
      # Cell features panel
      ####################
      list(tabItem(tabName = "cellFeatures",
        fluidRow(
          column(width = 3, uiOutput("cell.menu.ui")),
          column(width = 9, uiOutput("cell.plot.ui"))
        )
      )),
      ####################
      # Single-gene expression panel
      ####################
      list(tabItem(tabName = "onegeneExpression",
        fluidRow(
          column(width = 3, uiOutput("gene.menu.ui")),
          column(width = 9, uiOutput("gene.plot.ui"), uiOutput("gene.morePlots.ui"))
        )
      )),
      ####################
      # Multi-gene expression panel
      ####################
      list(tabItem(tabName = "multigeneExpression",
        fluidRow(
          box(title = "Input", width = 12, status = "primary", solidHeader = TRUE, collapsible = TRUE,
            fluidRow(
               column(width = 3,
                 textAreaInput("multi.ganenames", "Add gene symbols:", height = "325px", resize = "none",
                               placeholder = sprintf("Paste a set of valid Entrez gene symbols, one gene in a row. Maximum %d genes used in plots.", max.gene)),
                 actionButton("multi.submit", "Submit", width = "75px", class = "btn-primary", style = "color: #fff;"),
                 actionButton("multi.example", "Use example", width = "100px", class = "btn-info", style = "color: #fff;"),
                 actionButton("multi.reset", "Reset input", width = "100px")
               ),
               column(width = 4,
                 radioGroupButtons(inputId = "multi.plot_type", label = "Plot type:", direction = "horizontal",
                                   choices = c("Dot", "Heatmap", "Boxplot", "Projection"), selected = "Dot"),
                 uiOutput("multi.menu.ui")
               ),
               uiOutput("multi.stats")
            )
          )
        ),
          fluidRow(uiOutput("multi.plot.ui")
        )
      )),
      ####################
      # findMarkers panel
      ####################
      fm_items$items,
      ####################
      # edgeR panel
      ####################
      edger_items$items,
      ####################
      # DEseq2 panel
      ####################
      deseq_items$items,
      ####################
      # enrichR: findMarkers panel
      ####################
      list(tabItem(tabName = "ORAfm",
        fluidRow(
          column(width = 12, uiOutput("ORAfm.menu.ui")),
          column(width = 12, uiOutput("ORAfm.output.ui"))
        )
      )),
      ####################
      # enrichR: edgeR panel
      ####################
      list(tabItem(tabName = "ORAer",
        fluidRow(
          column(width = 12, uiOutput("ORAer.menu.ui")),
          column(width = 12, uiOutput("ORAer.output.ui"))
        )
      )),
      ####################
      # enrichR: DESeq2 panel
      ####################
      list(tabItem(tabName = "ORAde",
        fluidRow(
          column(width = 12, uiOutput("ORAde.menu.ui")),
          column(width = 12, uiOutput("ORAde.output.ui"))
        )
      )),
      ####################
      # About panel
      #################### 
      list(tabItem(tabName = "about",
        fluidRow(
          box(title = "About the App", width = 6, status = "primary", solidHeader = TRUE,
            column(width = 12, h3(app.title), h4(app.version),
            p(paste("Author:", app.author)),
	    p(paste("License:", app.license)),
            HTML(paste0("<p>", app.description, "</p>")))
          )
        )
      ))
    )
    do.call(tabItems, items)
  })

  ####################
  # Build findMarkers tabItems
  ####################
  fm_items <- reactiveValues(items = NULL)

  observeEvent(input$sce.file, {
    sce <- if(is.null(input$sce.file)) init.sce else sce()
    fm.listnames <- names(metadata(sce))[grep("^findMarkers", names(metadata(sce)))]
    if(length(fm.listnames) > 0) {
      fm_items$items <- NULL # reset and re-build
      for(listname in fm.listnames) {
        tab.name <- gsub("_", ": ", gsub("findMarkers_", "", listname))
        fm.sublistnames <- names(metadata(sce)[[listname]])

        fm_items$items[[length(fm_items$items)+1]] <- tabItem(tabName = listname,
          fluidRow(box(title = tab.name, width = 12, status = "primary", solidHeader = TRUE, collapsible = FALSE,
            lapply(fm.sublistnames, function(sublistname) {
              fm.df <- as.data.frame(metadata(sce)[[listname]][[sublistname]]) %>% dplyr::arrange(Top, FDR)
              fluidRow(box(title = span(fa("caret-right", fill = "purple"), sublistname, fa("caret-left", fill = "purple")),
                           width = 12, status = "primary", solidHeader = FALSE, collapsible = TRUE,
                           renderDT(datatable(fm.df, options = list(searching = TRUE, pageLength = 10, scrollX = TRUE, lengthChange = FALSE), 
					      rownames = FALSE, selection = "none", class = "white-space: nowrap") %>% 
			   formatRound(columns = colnames(fm.df)[grep("logFC", colnames(fm.df))], digits = 4) %>% 
			   formatSignif(columns = c("p.value", "FDR"), digits = 4))))
        }))))
      }
    updateTabItems(session, "menuItems", "about")
    updateTabItems(session, "menuItems", "overview")
    }
  })

  ####################
  # Build edgeR tabItems
  ####################
  edger_items <- reactiveValues(items = NULL)

  observeEvent(input$sce.file, {
    sce <- if(is.null(input$sce.file)) init.sce else sce()
    edger.listnames <- names(metadata(sce))[grep("^edgeR", names(metadata(sce)))]
    if(length(edger.listnames) > 0) {
      edger_items$items <- NULL # reset and re-build
      for(listname in edger.listnames) {
        tab.name <- gsub("_", " vs. ", gsub("edgeR_", "", listname))
        edger.sublistnames <- names(metadata(sce)[[listname]])

        edger_items$items[[length(edger_items$items)+1]] <- tabItem(tabName = listname,
          fluidRow(box(title = tab.name, width = 12, status = "primary", solidHeader = TRUE, collapsible = FALSE,
            lapply(edger.sublistnames, function(sublistname) {
              res <- metadata(sce)[[listname]][[sublistname]]
              edger.df <- edgeR::topTags(res, n = nrow(res), sort.by = "PValue", adjust.method = "BH")$table %>% as.data.frame %>% dplyr::arrange(FDR, Symbol)
              if("F" %in% colnames(edger.df)) { # QLFTest
                fluidRow(box(title = span(fa("caret-right", fill = "purple"), paste(sublistname, "(QLFTest)"),
                                          fa("caret-left", fill = "purple")), width = 12, status = "primary", solidHeader = FALSE, collapsible = TRUE,
                             renderDT(datatable(edger.df, options = list(searching = TRUE, pageLength = 10, scrollX = TRUE, lengthChange = FALSE),
                                      rownames = FALSE, selection = "none", class = "white-space: nowrap") %>%
                                      formatRound(columns = c("logFC","logCPM","F"), digits = 4) %>% formatSignif(columns = c("PValue", "FDR"), digits = 4))))
              } else { # LRT
                logFC_cols <- colnames(edger.df)[grep("logFC", colnames(edger.df))]
                fluidRow(box(title = span(fa("caret-right", fill = "purple"), paste(sublistname, "(LRT)"),
                                          fa("caret-left", fill = "purple")), width = 12, status = "primary", solidHeader = FALSE, collapsible = TRUE,
                             renderDT(datatable(edger.df, options = list(searching = TRUE, pageLength = 10, scrollX = TRUE, lengthChange = FALSE),
                                      rownames = FALSE, selection = "none", class = "white-space: nowrap") %>%
                                      formatRound(columns = c(logFC_cols,"logCPM","LR"), digits = 4) %>% formatSignif(columns = c("PValue", "FDR"), digits = 4))))
              }
        }))))
      }
#    updateTabItems(session, "menuItems", "about")
#    updateTabItems(session, "menuItems", "overview")
    }
  })

  ####################
  # Build DESeq2 tabItems
  ####################
  deseq_items <- reactiveValues(items = NULL)

  observeEvent(input$sce.file, {
    sce <- if(is.null(input$sce.file)) init.sce else sce()
    deseq.listnames <- names(metadata(sce))[grep("^DESeq2", names(metadata(sce)))]
    if(length(deseq.listnames) > 0) {
      deseq_items$items <- NULL # reset and re-build
      for(listname in deseq.listnames) {
        tab.name <- gsub("_", " vs. ", gsub("DESeq2_", "", listname))
        deseq.sublistnames <- names(metadata(sce)[[listname]])

        deseq_items$items[[length(deseq_items$items)+1]] <- tabItem(tabName = listname,
        # Check if DESeq2 is installed
        if(class(try(find.package("DESeq2"), silent = TRUE)) %in% c("try-error", "NULL")) {
          fluidRow(box(title = "The DESeq2 package is required to view stored results",
                       width = 6, status = "primary", solidHeader = TRUE, collapsible = FALSE,
                       column(width = 6, "Please install",
                              a("DESeq2", href = "https://bioconductor.org/packages/release/bioc/html/DESeq2.html", target = "_blank"),
                              "and restart the Shiny App.")))
        } else {
            fluidRow(box(title = tab.name, width = 12, status = "primary", solidHeader = TRUE, collapsible = FALSE,
              lapply(deseq.sublistnames, function(sublistname) {
                res <- metadata(sce)[[listname]][[sublistname]]
                deseq.df <- res %>% as.data.frame %>% select(-lfcSE) %>% dplyr::arrange(padj, pvalue, Symbol)
                fluidRow(box(title = span(fa("caret-right", fill = "purple"), sublistname, fa("caret-left", fill = "purple")),
                             width = 12, status = "primary", solidHeader = FALSE, collapsible = TRUE,
                             renderDT(datatable(deseq.df, options = list(searching = TRUE, pageLength = 10, scrollX = TRUE, lengthChange = FALSE),
                                                rownames = FALSE, selection = "none", class = "white-space: nowrap") %>%
                                      formatRound(columns = c("baseMean","log2FoldChange","stat"), digits = 4) %>% formatSignif(columns = c("pvalue", "padj"), digits = 4))))
          })))
        })
      }
    }
  })

  sce.vars <- reactive({
    sce <- sce()
    sce.ncells <- ncol(sce)
    sce.ngenes <- nrow(sce)

    # Unique names
    sce.samples <- levels(sce$Sample)
    sce.labels <- levels(colLabels(sce))
    sce.celltypes <- levels(sce$CellType)
    sce.dimreds <- reducedDimNames(sce)

    # Default reducedDim to show
    default.dimred <- if(default.dimred %in% sce.dimreds) default.dimred else if("TSNE" %in% sce.dimreds) "TSNE" else tail(sce.dimreds, 1)

    # Remove "color_by" elements not found in colData()
    sce.color_by <- default.color_by[default.color_by %in% colnames(colData(sce))]

    if(any(grepl("merged", colnames(colData(sce))))) { # if integrated sce, show only "merged" clustering results
      default.cluster.methods <- paste0("merged.", default.cluster.methods)
    }
    names(default.cluster.methods) <- default.cluster.methods.desc

    # Append custom clustering methods
    if(length(my.cluster.methods > 0)) {
      names(my.cluster.methods) <- my.cluster.methods.desc
      default.cluster.methods <- c(default.cluster.methods[!default.cluster.methods %in% my.cluster.methods], my.cluster.methods)
    }

    # Remove "cluster.methods" elements not found in colData()
    cluster.methods <- default.cluster.methods[default.cluster.methods %in% colnames(colData(sce))]

    # Determine the clustering method stored in colLabels
    default.cluster.method <- NULL
    for(m in cluster.methods) {
      if(identical(colData(sce)[, m], colLabels(sce))) default.cluster.method <- setNames(m, names(which(cluster.methods == m)))
    }

    # Combine cell features and clustering methods for cell colouring
    sce.color_by <- c(sce.color_by, cluster.methods)

    # Add "default" to default clustering method in drop-down list
    names(sce.color_by)[sce.color_by == default.cluster.method] <- paste(names(sce.color_by)[sce.color_by == default.cluster.method], "(default)")

    # Remove "group_by" elements not found in colData()
    sce.group_by <- default.group_by[default.group_by %in% colnames(colData(sce))]

    # Remove ClusterCellType from "group_by" if N/A
    sce.clustercelltypes <- levels(sce$ClusterCellType) # N/A will be NULL
    if(is.null(sce.clustercelltypes)) sce.group_by <- sce.group_by[sce.group_by != "ClusterCellType"]

    list(sce.ncells = sce.ncells, sce.ngenes = sce.ngenes, sce.samples = sce.samples, sce.labels = sce.labels, sce.celltypes = sce.celltypes, 
	 sce.clustercelltypes = sce.clustercelltypes, sce.dimreds = sce.dimreds, default.dimred = default.dimred, sce.color_by = sce.color_by, 
	 cluster.methods = cluster.methods, default.cluster.method = default.cluster.method, sce.group_by = sce.group_by)
  })

  sce.examples <- reactive({
    sce <- sce()
    sce.genenames <- setNames(sort(rownames(sce)), sort(rownames(sce)))

    # Build a vector containing gene symbols to be used as example in multi-gene plots
    if(!is.null(my.example.genes)) {
      example <- sce.genenames[sce.genenames %in% my.example.genes]
      if(length(example) == 0) {
        set.seed(12345)
        example <- sample(sce.genenames, min(c(max.gene, length(sce.genenames))))
      }
    } else {
      if(length(metadata(sce)[["runInfo"]][["HVG"]][["Genes"]]) > 0) {
        example <- head(metadata(sce)[["runInfo"]][["HVG"]][["Genes"]], max.gene)
      } else {
        set.seed(12345)
        example <- sample(sce.genenames, min(c(max.gene, length(sce.genenames))))
      }
    }
    example
  })

  ####################
  # Overview panel
  ####################
  output$overview.ui <- renderUI({
    sce.vars <- sce.vars()
    sce.ncells <- sce.vars[["sce.ncells"]]
    sce.samples <- sce.vars[["sce.samples"]]
    sce.labels <- sce.vars[["sce.labels"]]
    default.cluster.method <- sce.vars[["default.cluster.method"]]

    box(title = "Overview", width = 12, status = "primary", solidHeader = TRUE, 
      column(width = 4, align = "center", 
        box(title = h3("No. of cells", style = "display:inline-block;font-size:20px;margin:0;font-weight:bold;"), 
	    width = 12, status = "info", solidHeader = TRUE, style = "font-size: 30px; font-weight:bold;", sce.ncells)),
      column(width = 4, align = "center", 
        box(title = h3("No. of samples", style = "display:inline-block;font-size:20px;margin:0;font-weight:bold;"), 
	    width = 12, status = "info", solidHeader = TRUE, style = "font-size: 30px; font-weight:bold;", length(sce.samples))),
      column(width = 4, align = "center", 
        box(title = h3(if(!is.null(default.cluster.method)) paste0("No. of Clusters (based on ", names(default.cluster.method), ")") else paste0("No. of Clusters"), 
		       style = "display:inline-block;font-size:20px;margin:0;font-weight:bold;"), 
	    width = 12, status = "info", solidHeader = TRUE, style = "font-size: 30px; font-weight:bold;", length(sce.labels)))
    )
  })

  output$plotOverview <- renderPlotly({
    sce <- sce()
    sce.vars <- sce.vars()
    sce.samples <- sce.vars[["sce.samples"]]
    sce.labels <- sce.vars[["sce.labels"]]

    # Prepare DataFrame
    df1 <- as.data.frame(table(sce$Sample))
    df2 <- as.data.frame(table(sce$label))

    widths <- round(length(sce.samples)/(length(sce.samples)+length(sce.labels)), 2)
    if(widths > 0.7) widths <- 0.7 else if(widths < 0.2) widths <- 0.2
    widths <- c(widths, 1-widths)

    # Create plot
    fig1 <- plot_ly(df1, x = ~Var1, y = ~Freq) %>%
      add_bars(cliponaxis = FALSE, name = "", hovertemplate = "%{y:,} cells<br />in Sample %{x}") %>%
      layout(xaxis = list(title = "Samples"), yaxis = list(title = "No. of cells", range = c(0, max(df1$Freq)*1.15)))

    fig2 <- plot_ly(df2, x = ~Var1, y = ~Freq) %>%
      add_bars(cliponaxis = FALSE, name = "", hovertemplate = "%{y:,} cells<br />in Cluster %{x}") %>%
      layout(xaxis = list(title = "Cluster"), yaxis = list(title = "No. of cells" , range = c(0, max(df2$Freq)*1.15)))

    subplot(fig1, fig2, titleX = TRUE, titleY = TRUE, margin = c(0.06, 0.00, 0.0, 0.0), # c(left, right, top, bottom)
            nrows = 1, widths = widths) %>%
      layout(showlegend = FALSE, annotations = list(
        list(
          text = "<b>Cells in each sample</b>", x = -0.04, y = 1.02, xanchor = "left", yanchor = "bottom",
          showarrow = FALSE, xref = "paper", yref = "paper", font = list(size = 16)
        ),
        list(
          text = "<b>Cells in each cluster</b>", x = widths[1]+0.01, y = 1.02, xanchor = "left", yanchor = "bottom",
          showarrow = FALSE, xref = "paper", yref = "paper", font = list(size = 16)
        )
      ))
  })

  output$runinfo.ui <- renderUI({
    sce <- sce()
    runInfo <- metadata(sce)[["runInfo"]]

    if(all(c("Sample","Data") %in% names(runInfo))) {
      if(class(runInfo[["Sample"]]) == "character") { # Single-sample
        df.sample <- data.frame(t(runInfo[["Sample"]]), check.names = FALSE)
        df.data <- data.frame(t(runInfo[["Data"]]), check.names = FALSE)
      } else { # Integrated
        df.sample <- runInfo[["Sample"]]
        df.data <- runInfo[["Data"]]
      }
      box(title = "Run Info", width = 12, status = "primary", solidHeader = TRUE, collapsible = TRUE,
        box(title = "Sample Summary", width = 12, status = "primary", solidHeader = FALSE, collapsible = FALSE,
            column(width = 12, div(renderDT(datatable(df.sample, options = list(searching = FALSE, pageLength = 20, scrollX = TRUE, lengthChange = FALSE),
                                                      selection = "none")), style = "font-size:105%"))),
        box(title = "Data Analysis Summary", width = 12, status = "primary", solidHeader = FALSE, collapsible = FALSE,
            column(width = 12, div(renderDT(datatable(df.data, options = list(searching = FALSE, pageLength = 20, scrollX = TRUE, lengthChange = FALSE),
                                                      selection = "none")), style = "font-size:105%"))))
    } else { # pre-v2.0.0
      if(all(c("10X","BCF") %in% names(runInfo))) {
        if(class(runInfo[["10X"]]) == "character") { # 10X single-sample
          df.10x <- data.frame(t(runInfo[["10X"]]), check.names = FALSE)
          df.bcf <- data.frame(t(runInfo[["BCF"]]), check.names = FALSE)
        } else { # 10X integrated
          df.10x <- runInfo[["10X"]]
          df.bcf <- runInfo[["BCF"]]
        }
        box(title = "Run Info", width = 12, status = "primary", solidHeader = TRUE, collapsible = TRUE,
          box(title = "Sample Summary", width = 12, status = "primary", solidHeader = FALSE, collapsible = FALSE,
              column(width = 12, div(renderDT(datatable(df.10x, options = list(searching = FALSE, pageLength = 20, scrollX = TRUE, lengthChange = FALSE),
                                                        selection = "none")), style = "font-size:105%"))),
          box(title = "Data Analysis Summary", width = 12, status = "primary", solidHeader = FALSE, collapsible = FALSE,
              column(width = 12, div(renderDT(datatable(df.bcf, options = list(searching = FALSE, pageLength = 20, scrollX = TRUE, lengthChange = FALSE),
                                                        selection = "none")), style = "font-size:105%"))))
      } else if(all(c("Parse","Wells","BCF") %in% names(runInfo))) {
        df.parse <- data.frame("Sample Summary" = runInfo[["Parse"]], check.names = FALSE)
        df.bcf <- data.frame("Data Analysis Summary" = runInfo[["BCF"]], check.names = FALSE)
        df.wells <- data.frame("Well Number" = runInfo[["Wells"]], check.names = FALSE)
        box(title = "Run Info", width = 12, status = "primary", solidHeader = TRUE, collapsible = TRUE,
          column(width = 4, div(renderDT(datatable(df.parse, options = list(searching = FALSE, pageLength = 20, scrollX = TRUE, lengthChange = FALSE),
                                                   rownames = FALSE, selection = "none")), style = "font-size:105%")),
          column(width = 5, div(renderDT(datatable(df.bcf, options = list(searching = FALSE, pageLength = 20, scrollX = TRUE, lengthChange = FALSE),
                                                   rownames = FALSE, selection = "none")), style = "font-size:105%")),
          column(width = 3, div(renderDT(datatable(df.wells, options = list(searching = FALSE, pageLength = 20, scrollX = TRUE, lengthChange = FALSE),
                                                   rownames = FALSE, selection = "none")), style = "font-size:105%")))
      } else if(all(c("Parse Biosciences","BCF") %in% names(runInfo))) {
        df.parse <- runInfo[["Parse Biosciences"]]
        df.bcf <- runInfo[["BCF"]]
        box(title = "Run Info", width = 12, status = "primary", solidHeader = TRUE, collapsible = TRUE,
          box(title = "Sample Summary", width = 12, status = "primary", solidHeader = FALSE, collapsible = FALSE,
              column(width = 12, div(renderDT(datatable(df.parse, options = list(searching = FALSE, pageLength = 20, scrollX = TRUE, lengthChange = FALSE),
                                                        selection = "none")), style = "font-size:105%"))),
          box(title = "Data Analysis Summary", width = 12, status = "primary", solidHeader = FALSE, collapsible = FALSE,
              column(width = 12, div(renderDT(datatable(df.bcf, options = list(searching = FALSE, pageLength = 20, scrollX = TRUE, lengthChange = FALSE),
                                                        selection = "none")), style = "font-size:105%"))))
      }
    }
  })

  output$table_sample_label <- renderDT({
    sce <- sce()
    datatable(as.data.frame.matrix(table(Sample = sce$Sample, Cluster = sce$label)), options = list(searching = FALSE, pageLength = 10, scrollX = TRUE, lengthChange = FALSE), 
	      selection = "none")
  })

  output$overview.plot.ui <- renderUI({
    box(title = "Samples and Clusters", width = 12, status = "primary", solidHeader = TRUE, collapsible = TRUE,
      column(width = 12, plotlyOutput("plotOverview", width = "100%", height = "350px")),
      column(width = 12, DTOutput("table_sample_label"))
    )
  })

  ####################
  # Cell features & Gene expression projection (TSNE/UMAP)
  ####################
  observeEvent(input$menuItems, {
    # Server-side selectize
    if(input$menuItems == "onegeneExpression") {
      sce <- sce()
      sce.genenames <- setNames(sort(rownames(sce)), sort(rownames(sce)))
      if(reset.gene.submit() == 1) 
        updateSelectizeInput(session, "gene.color_by", choices = sce.genenames, selected = character(0), server = TRUE)
      else
        updateSelectizeInput(session, "gene.color_by", choices = sce.genenames, selected = isolate(input$gene.color_by), server = TRUE)
    }
  })

  output$cell.menu.ui <- renderUI({
    sce.vars <- sce.vars()
    sce.ncells <- sce.vars[["sce.ncells"]]
    sce.color_by <- sce.vars[["sce.color_by"]]
    sce.dimreds <- sce.vars[["sce.dimreds"]]
    default.dimred <- sce.vars[["default.dimred"]]
    sce.samples <- sce.vars[["sce.samples"]]
    default.cluster.method <- sce.vars[["default.cluster.method"]]
    sce.labels <- sce.vars[["sce.labels"]]
    sce.celltypes <- sce.vars[["sce.celltypes"]]

    if(length(sce.dimreds) > 0) {
      box(title = "Input", width = 12, status = "primary", solidHeader = TRUE,
        selectInput("cell.color_by", "Colour cells by:", choices = names(sce.color_by), multiple = FALSE, selectize = TRUE),
        selectInput("cell.dimred", "Select a projection:", choices = sce.dimreds, selected = default.dimred, multiple = FALSE, selectize = TRUE),
        actionButton("cell.submit", "Submit", width = "75px", class = "btn-primary", style = "color: #fff;"),
        br(), br(),
        selectInput("cell.palette_dis", "Select a palette for discrete data:", choices = discrete, selected = "Set1", multiple = FALSE, selectize = TRUE),
        selectInput("cell.palette_con", "Select a palette for continuous data:", choices = names(continuous), selected = "viridis", multiple = FALSE, selectize = TRUE),
        sliderTextInput("cell.dot_size","Dot size:", choices = seq(1, 15, by = 0.5), selected = dot_size(sce.ncells), grid = T),
        sliderInput("cell.dot_opacity", "Dot opacity:", 0, 1, 0.8, 0.1),
        pickerInput("cell.select_sample", "Select samples:", choices = sce.samples, selected = sce.samples,
          options = list(`actions-box` = TRUE, `selected-text-format` = paste0("count == ", length(sce.samples))), multiple = TRUE),
        pickerInput("cell.select_cluster",
          if(!is.null(default.cluster.method)) paste0("No. of Clusters (based on ", names(default.cluster.method), ")") else paste0("No. of Clusters"),
          choices = sce.labels, selected = sce.labels,
          options = list(`actions-box` = TRUE, `selected-text-format` = paste0("count == ", length(sce.labels))), multiple = TRUE),
        pickerInput("cell.select_celltype", "Select cell types:", choices = sce.celltypes, selected = sce.celltypes,
          options = list(`actions-box` = TRUE, `selected-text-format` = paste0("count == ", length(sce.celltypes))), multiple = TRUE)
      )
    } else {
      box(title = "Cell feature", width = 12, status = "primary", solidHeader = TRUE, column(width = 12, p("There is no dimensional reduction available on this data.")))
    }
  })

  output$gene.menu.ui <- renderUI({
    sce.vars <- sce.vars()
    sce.ncells <- sce.vars[["sce.ncells"]]
    sce.dimreds <- sce.vars[["sce.dimreds"]]
    default.dimred <- sce.vars[["default.dimred"]]
    sce.samples <- sce.vars[["sce.samples"]]
    default.cluster.method <- sce.vars[["default.cluster.method"]]
    sce.labels <- sce.vars[["sce.labels"]]
    sce.celltypes <- sce.vars[["sce.celltypes"]]
    sce.group_by <- sce.vars[["sce.group_by"]]

    if(length(sce.dimreds) > 0) {
      box(title = "Input", width = 12, status = "primary", solidHeader = TRUE,
        selectizeInput("gene.color_by", "Colour cells by:", choices = NULL, selected = NULL, multiple = FALSE,
                       options = list(placeholder = "Enter a gene symbol.", plugins = list('restore_on_backspace'))),
        actionButton("gene.submit", "Submit", width = "75px", class = "btn-primary", style = "color: #fff;"),
        tags$p(tags$div(style = "background-color: #6b2c91;", tags$b("Projection settings"), style = "text-align: center; color: #ffd500;")),
        selectInput("gene.dimred", "Select a projection:", choices = sce.dimreds, selected = default.dimred, multiple = FALSE, selectize = TRUE),
        selectInput("gene.palette_con", "Select a colour palette:", choices = names(continuous), selected = "yellowRed", multiple = FALSE, selectize = TRUE),
        sliderTextInput("gene.dot_size","Dot size:", choices = seq(1, 15, by = 0.5), selected = dot_size(sce.ncells), grid = T),
        sliderInput("gene.dot_opacity", "Dot opacity:", 0, 1, 0.8, 0.1),
        pickerInput("gene.select_sample", "Select samples:", choices = sce.samples, selected = sce.samples,
          options = list(`actions-box` = TRUE, `selected-text-format` = paste0("count == ", length(sce.samples))), multiple = TRUE),
        pickerInput("gene.select_cluster",
          if(!is.null(default.cluster.method)) paste0("No. of Clusters (based on ", names(default.cluster.method), ")") else paste0("No. of Clusters"),
          choices = sce.labels, selected = sce.labels,
          options = list(`actions-box` = TRUE, `selected-text-format` = paste0("count == ", length(sce.labels))), multiple = TRUE),
        pickerInput("gene.select_celltype", "Select cell types:", choices = sce.celltypes, selected = sce.celltypes,
          options = list(`actions-box` = TRUE, `selected-text-format` = paste0("count == ", length(sce.celltypes))), multiple = TRUE),
        tags$p(tags$div(style = "background-color: #6b2c91;", tags$b("Bar/boxplot settings"), style = "text-align: center; color: #ffd500;")),
        checkboxGroupInput("gene.group_by", "Group cells by:", sce.group_by, selected = "label"),
        selectInput("gene.palette_dis", "Select a colour palette:", choices = discrete, selected = "Set1", multiple = FALSE, selectize = TRUE)
      )
    } else {
      box(title = "Gene expression", width = 12, status = "primary", solidHeader = TRUE, column(width = 12, p("There is no dimensional reduction available on this data.")))
    }
  })

  cell.prepare <- eventReactive(input$cell.submit, {
    isolate({
      dimred <- input$cell.dimred
      color_by <- input$cell.color_by
      palette_con <- input$cell.palette_con
      palette_dis <- input$cell.palette_dis
      dot_size <- input$cell.dot_size
      dot_opacity <- input$cell.dot_opacity
      select_sample <- input$cell.select_sample
      select_cluster <- input$cell.select_cluster
      select_celltype <- input$cell.select_celltype
    })

    list(dimred = dimred, color_by = color_by, palette_con = palette_con, palette_dis = palette_dis, 
	 dot_size = dot_size, dot_opacity = dot_opacity, select_sample = select_sample, 
	 select_cluster = select_cluster, select_celltype = select_celltype)
  })

  gene.prepare <- eventReactive(input$gene.submit, {
    isolate({
      dimred <- input$gene.dimred
      color_by <- input$gene.color_by
      palette_con <- input$gene.palette_con
      palette_dis <- input$gene.palette_dis
      dot_size <- input$gene.dot_size
      dot_opacity <- input$gene.dot_opacity
      select_sample <- input$gene.select_sample
      select_cluster <- input$gene.select_cluster
      select_celltype <- input$gene.select_celltype
      group_by <- input$gene.group_by
    })

    list(dimred = dimred, color_by = color_by, palette_con = palette_con, palette_dis = palette_dis, 
	 dot_size = dot_size, dot_opacity = dot_opacity, select_sample = select_sample, 
	 select_cluster = select_cluster, select_celltype = select_celltype, group_by = group_by)
  })

  output$plotReducedDim1 <- output$plotReducedDim2 <- renderPlotly({
    whichMenu <- isolate(input$menuItems)
    if(whichMenu == "onegeneExpression") {
      if(reset.gene.submit() > 0) return()
      inputList <- gene.prepare()
      palette <- inputList[["palette_con"]]
    } else {
      if(reset.cell.submit() > 0) return()
      inputList <- cell.prepare()
      palette_dis <- inputList[["palette_dis"]]
      palette_con <- inputList[["palette_con"]]
    }

    dimred <- inputList[["dimred"]]
    color_by <- inputList[["color_by"]]
    dot_size <- inputList[["dot_size"]]
    dot_opacity <- inputList[["dot_opacity"]]
    select_sample <- inputList[["select_sample"]]
    select_cluster <- inputList[["select_cluster"]]
    select_celltype <- inputList[["select_celltype"]]

    sce <- sce()
    sce.vars <- sce.vars()
    sce.color_by <- sce.vars[["sce.color_by"]]

    # Prepare DataFrame
    df <- as.data.frame(reducedDim(sce, dimred)[,1:2])
    colnames(df) <- c("X","Y")
    df <- cbind(df, data.frame(colData(sce)[, c("Barcode", sce.color_by)]))

    if(whichMenu == "onegeneExpression") {
      df <- cbind(df, data.frame(expr = logcounts(sce)[color_by,]))
      expr.upper <- if(quantile(df$expr, 0.999) == 0) max(df$expr) else quantile(df$expr, 0.999)
    }

    # Calculate xaxis and yaxis
    xaxis <- c(-max(abs(df$X)), max(abs(df$X)))
    yaxis <- c(-max(abs(df$Y)), max(abs(df$Y)))

    # Subset DataFrame
    df <- df %>%
      filter(Sample %in% select_sample) %>%
      filter(label %in% select_cluster) %>%
      filter(CellType %in% select_celltype)

    # Plot title
    title <- paste0(dimred, " (coloured by ", color_by, ")")

    # Create plot
    fig <- plot_ly(df, x = ~X, y = ~Y, hoverinfo = "text", hoverlabel = list(bgcolor = "white"))

    if(whichMenu == "onegeneExpression") { # continuous data
      fig <- fig %>%
        add_trace(mode = "markers", type = "scatter",
		  color = ~expr, colors = continuous[[palette]], 
		  marker = list(size = dot_size, opacity = dot_opacity), 
		  text = ~paste0("<b>Cell:</b> ", Barcode, "<br />", 
				 "<b>Sample:</b> ", Sample, "<br />", 
				 "<b>Cell type:</b> ", CellType, "<br />", 
				 "<b>Cluster:</b> ", label, " (", ClusterCellType, ")<br />", 
				 "<b>Expression:</b> ", sprintf("%.4f", expr))) %>% 
        colorbar(title = NULL, limits = c(0, expr.upper))
    } else {
      color_type <- scale_type(df[, sce.color_by[color_by]])
      colors <- if(color_type == "continuous") continuous[[palette_con]] else discrete_func(palette_dis, length(unique(df[, sce.color_by[color_by]])))

      fig <- fig %>% 
        add_trace(mode = "markers", type = "scatter", 
		  color = as.formula(paste0("~", sce.color_by[color_by])), colors = colors, 
		  marker = list(size = dot_size, opacity = dot_opacity), 
		  text = ~paste0("<b>Cell:</b> ", Barcode, "<br />", 
				 "<b>Sample:</b> ", Sample, "<br />", 
				 "<b>Cell type:</b> ", CellType, "<br />", 
				 "<b>Cluster:</b> ", label, " (", ClusterCellType, ")"))

      if(color_type == "continuous") { # continuous data
        fig <- fig %>% colorbar(title = NULL) 
      } else { # discrete data
        fig <- fig %>% layout(legend = pl.legend)
      }
    }

    fig %>% layout(
      showlegend = TRUE,
      xaxis = list(title = paste(dimred, "1"), range = xaxis),
      yaxis = list(title = paste(dimred, "2"), range = yaxis),
      title = list(text = paste("<b>", title, "</b>"), x = 0, y = 1.01, xanchor = "left", yanchor = "bottom")
    )
  })

  output$cell.plot.ui <- renderUI({
    if(reset.cell.submit() == 0) {
      box(status = "primary", width = 12, 
        plotlyOutput("plotReducedDim1", width = "100%", height = 850) %>% withSpinner(type = getOption("spinner.type", default = 8))
      )
    }
  })

  output$gene.plot.ui <- renderUI({
    if(reset.gene.submit() == 0) {
      gene.color_by <- input$gene.color_by
      if(!is.null(gene.color_by) && nchar(gene.color_by) > 0) {
        box(status = "primary", width = 12,
            plotlyOutput("plotReducedDim2", width = "100%", height = 850) %>% withSpinner(type = getOption("spinner.type", default = 8))
        )
      }
    }
  })

  ####################
  # Gene expression panel (more plots)
  ####################
  # Determine number of levels in group_by
  gene.group_by_n <- eventReactive(input$gene.submit, {
    sce <- sce()
    inputList <- gene.prepare()
    groups <- inputList[["group_by"]]
    n <- c()
    for(i in 1:length(groups)) {
      group_by <- groups[i]
      n <- c(n, length(unique(colData(sce)[, group_by])))
    }
    n
  })

  for(i in 1:length(default.group_by)) { # set up maximum number of plots
    plotname <- paste0("gene.morePlots.", i)
    local({
      j <- i
      output[[plotname]] <- renderPlotly({
        if(reset.gene.submit() > 0) return()
        sce <- sce()
	sce.vars <- sce.vars()
	sce.color_by <- sce.vars[["sce.color_by"]]
        inputList <- gene.prepare()
        color_by <- inputList[["color_by"]] # Gene name
	group_by <- inputList[["group_by"]][j]
        
        group_by.desc <- names(which(sce.color_by == group_by))
	palette <- inputList[["palette_dis"]]

        # Prepare DataFrame
        df <- data.frame(colData(sce)[, c("Barcode", sce.color_by)])
        df <- cbind(df, data.frame(expr = logcounts(sce)[color_by,])) 

	# Calculate per-group sums and append to x-axis labels
	group_cells <- group_by(df, across(all_of(group_by)), .drop = FALSE) %>% summarize(num = n())

        # Create plot
	dodge <- position_dodge(width = 1)
        fig <- ggplot(df, aes(x = .data[[group_by]], y = .data[["expr"]], fill = .data[[group_by]])) +
		geom_violin(aes(text = paste("Group:", levels(pull(group_cells, 1))[after_stat(x)], "<br /># Cells:", after_stat(n)))) +
		geom_boxplot(color = "grey", alpha = 0.3) + 
		scale_fill_manual(values = discrete_func(palette, nrow(group_cells))) +
		scale_x_discrete(drop = FALSE) + coord_flip() + theme(legend.position = "none") +
                labs(title = paste(color_by, "expression"), x = group_by.desc, y = "logcounts", fill = group_by.desc)

        ggplotly(fig, tooltip = "text")
      })
    })
  }

  # Determine plot height
  gene.plot_height <- eventReactive(input$gene.submit, {
    n <- gene.group_by_n()
    heights <- c()
    if(!is.null(n)) {
      for(i in 1:length(n)) {
        height <- n[i] * (400/8) # 400 px for 8 features
        height <- if(height < 300) 300 else height # min 300 px
	heights <- c(heights, round(height))
      }
    }
    heights
  })

  output$gene.morePlots.ui <- renderUI({
    if(reset.gene.submit() == 0) {
      gene.color_by <- input$gene.color_by
      n.gene.group_by <- length(isolate(input$gene.group_by))
      if(!is.null(gene.color_by) && nchar(gene.color_by) > 0 && n.gene.group_by > 0) {
        heights <- gene.plot_height()
        plot_output_list <- lapply(1:n.gene.group_by, function(i) {
          box(status = "primary", width = 12,
            plotlyOutput(paste0("gene.morePlots.", i), height = heights[i], width = "100%") %>% withSpinner(type = getOption("spinner.type", default = 8))
          )
        })
        plot_output_list
      }
    }
  })

  ####################
  # Multi-gene expression panel
  ####################
  output$multi.menu.ui <- renderUI({
    sce.vars <- sce.vars()
    sce.dimreds <- sce.vars[["sce.dimreds"]]
    default.dimred <- sce.vars[["default.dimred"]]
    sce.group_by <- sce.vars[["sce.group_by"]]

    plot_group_text <- sprintf("Group cells by (max. %d; selection order recorded):", multi_max_options)
    none_selected_text <- sprintf("Select between 1 to %d features.", multi_max_options)

    if(input$multi.plot_type == "Projection") {
      tags$div(
        pickerInput("multi.dimred", "Select a projection:", choices = sce.dimreds, selected = default.dimred, multiple = FALSE),
        pickerInput("multi.palette", "Select a colour palette:", choices = names(continuous), selected = "yellowRed", multiple = FALSE),
        sliderTextInput("multi.dot_size","Dot size:", choices = seq(0.2, 3, by = 0.2), selected = 0.8, grid = T),
        sliderInput("multi.dot_opacity", "Dot opacity:", 0, 1, 0.8, 0.1)
      )
    } else if(input$multi.plot_type == "Boxplot") {
      tags$div(
        pickerInput("multi.plot_group", plot_group_text, choices = names(sce.group_by), selected = names(which(sce.group_by == "label")),
                    multiple = TRUE, options =  list("max-options" = multi_max_options, "none-selected-text" = none_selected_text)),
        fluidRow(
          column(width = 6, radioGroupButtons(inputId = "multi.color_by", label = "Colour cells by:", direction = "horizontal",
                                              choices = c("Detected", "Group"), selected = "Detected")),
          column(width = 6, radioGroupButtons(inputId = "multi.rotate_x", label = "Rotate X-axis labels:", direction = "horizontal",
                                              choices = c("No", "45°", "90°"), selected = "90°"))
        ),
        sliderInput("multi.facet_ncol", "Genes per row:", 1, 15, 3, 1),
        sliderInput("multi.fontsize.x", "X-axis label size:", 6, 14, 10, 1)
      )
    } else if(input$multi.plot_type == "Dot") {
      tags$div(
        pickerInput("multi.plot_group", plot_group_text, choices = names(sce.group_by), selected = names(which(sce.group_by == "label")),
                    multiple = TRUE, options =  list("max-options" = multi_max_options, "none-selected-text" = none_selected_text)),
        fluidRow(
          column(width = 6, radioGroupButtons(inputId = "multi.plot_scale", label = "Expression scaling:", direction = "horizontal",
                                              choices = c("None", "Centered & scaled"), selected = "Centered & scaled")),
          column(width = 6, radioGroupButtons(inputId = "multi.rotate_x", label = "Rotate X-axis labels:", direction = "horizontal",
                                              choices = c("No", "45°", "90°"), selected = "90°"))
        ),
        radioGroupButtons(inputId = "multi.plot_cluster", label = "Hierarchical clustering:", direction = "horizontal",
                          choices = c("None", "Group (X-axis)", "Gene (Y-axis)", "Both"), selected = "None"),
        fluidRow(
          column(width = 6, sliderInput("multi.dots_size", "Maximun dot size:", 3, 10, 6, 1)),
          column(width = 6, sliderInput("multi.base_size", "Axis label size:", 8, 22, base_size, 1))
        )
      )
    } else if(input$multi.plot_type == "Heatmap") {
      tags$div(
        pickerInput("multi.plot_group", plot_group_text, choices = names(sce.group_by), selected = names(which(sce.group_by == "label")),
                    multiple = TRUE, options =  list("max-options" = multi_max_options, "none-selected-text" = none_selected_text)),
        fluidRow(
          column(width = 6, radioGroupButtons(inputId = "multi.plot_scale", label = "Expression scaling:", direction = "horizontal",
                                              choices = c("None", "Centered & scaled"), selected = "Centered & scaled")),
          column(width = 6, radioGroupButtons(inputId = "multi.rotate_x", label = "Rotate X-axis labels:", direction = "horizontal",
                                              choices = c("No", "45°", "90°"), selected = "90°"))
        ),
        sliderInput("multi.base_size", "Axis label size:", 8, 22, base_size, 1)
      )
    }
  })

  # Create a reactiveValues to store the selection order of 'multi.plot_group'
  multi.plot_group_order <- reactiveValues(values = NULL)

  # Use reactive to get and sort the selected terms in the order of selection
  ordered_plot_group <- reactive({
    multi.plot_group_order$values <- if (length(multi.plot_group_order$values) > length(input$multi.plot_group))
            multi.plot_group_order$values[multi.plot_group_order$values %in% input$multi.plot_group]
    else c(multi.plot_group_order$values, input$multi.plot_group[!input$multi.plot_group %in% multi.plot_group_order$values])
    multi.plot_group_order$values
  })

  # Use observe to update the reactive function above
  observe({ ordered_plot_group() })

  # Show example genes
  observeEvent(input$multi.example, {
    example <- sce.examples()
    example <- paste(example, collapse = "\r\n")
    updateTextAreaInput(session, "multi.ganenames", value = example)
  })

  observeEvent(input$multi.reset, {
    sce.vars <- sce.vars()
    sce.color_by <- sce.vars[["sce.color_by"]]
    default.dimred <- sce.vars[["default.dimred"]]

    updateTextAreaInput(session, "multi.ganenames", value = "")
    updateRadioGroupButtons(session, "multi.plot_type", selected = "Dot")
    updateRadioGroupButtons(session, "multi.plot_scale", selected = "Centered & scaled")
    updateRadioGroupButtons(session, "multi.plot_group", selected = names(which(sce.color_by == "label")))
    updateRadioGroupButtons(session, "multi.color_by", selected = "Detected")
    updateRadioGroupButtons(session, "multi.rotate_x", selected = "90°")
    updateRadioGroupButtons(session, "multi.plot_cluster", selected = "None")
    updateSelectInput(session, "multi.dimred", selected = default.dimred)
    updateSelectInput(session, "multi.palette", selected = "yellowRed")
    updateSliderInput(session, "multi.dot_size", value = 0.8)
    updateSliderInput(session, "multi.dot_opacity", value = 0.8)
    updateSliderInput(session, "multi.facet_ncol", value = 3)
    updateSliderInput(session, "multi.fontsize.x", value = 10)
    updateSliderInput(session, "multi.dots_size", value = 6)
    updateSliderInput(session, "multi.base_size", value = base_size)
  })

  # Compare user input genes with sce
  multi.compare.name <- eventReactive(input$multi.submit, {
    sce <- sce()
    sce.genenames <- setNames(sort(rownames(sce)), sort(rownames(sce)))
    multi.input <- unlist(strsplit(x = isolate(input$multi.ganenames), split = "[\r\n]"))
    multi.input <- trimws(multi.input) # remove leading and trailing whitespace
    multi.input <- unique(multi.input) # remove duplicated names
    multi.input <- multi.input[multi.input != ""] # remove empty element
    res <- multi.input %in% sce.genenames # TRUE/FALSE
    setNames(res, multi.input)
  })

  # Check if user input genes are expressed
  multi.compare.expr <- eventReactive(input$multi.submit, {
    sce <- sce()
    res <- multi.compare.name()
    multi.input <- names(res[res == TRUE])

    # TRUE/FALSE
    if(length(multi.input) == 1) setNames(sum(logcounts(sce)[multi.input,]) > 0, multi.input) else setNames(rowSums(logcounts(sce)[multi.input,]) > 0, multi.input)
  })

  # Print 'multi.plot_group' selection order
  output$multi.ordered_plot_group <- renderPrint({
    cat(sprintf("Group cells by: %s", paste(ordered_plot_group(), collapse = ", ")))
  })

  # Print input gene stats to verbatimTextOutput()
  output$multi.n_matched <- renderPrint({ 
    res <- multi.compare.name()
    if(reset.multi.submit() == 0) cat(sprintf("Of the %d entries provided, %d genes found in dataset.", length(res), sum(res))) else cat("Calculating...")
  })

  # Print unmatched genes to verbatimTextOutput()
  output$multi.not_found <- renderPrint({
    res <- multi.compare.name()
    if(reset.multi.submit() == 0) {
      if(sum(res == FALSE) > 0) cat("Genes not found:", shQuote(names(res[res == FALSE]))) else cat("Genes not found: none")
    } else cat("Calculating...")
  })

  # Print unexpressed genes to verbatimTextOutput()
  output$multi.not_expr <- renderPrint({
    res <- multi.compare.expr()
    if(reset.multi.submit() == 0) {
      if(sum(res == FALSE) > 0) cat("Genes not expressed:", shQuote(names(res[res == FALSE])), "(not included in scaled plot)") else cat("Genes not expressed: none")
    } else cat("Calculating...")
  })

  output$multi.stats <- renderUI({
    if(reset.multi.submit() == 0)
      column(width = 5, verbatimTextOutput("multi.ordered_plot_group"), verbatimTextOutput("multi.n_matched"), div(class = "not_found", verbatimTextOutput("multi.not_found")),
             div(class = "not_expr", verbatimTextOutput("multi.not_expr")))
  })

  # Validate and subset genes
  multi.features <- eventReactive(input$multi.submit, {
    res <- multi.compare.name()
    if(sum(res) > 0) {
      features <- names(res[res])
      features <- head(features, max.gene)
    } else {
      features <- NULL
    }
    features
  })

  # Determine plot height
  multi.plot_height <- eventReactive(input$multi.submit, {
    features <- multi.features()
    n_group_by <- length(ordered_plot_group())
    height <- 0
    if(!is.null(features)) {
      if(input$multi.plot_type == "Projection") {
        # 250 px for 5 sets of genes
	height <- ceiling(length(features)/5) * 275
      } else if(input$multi.plot_type == "Boxplot") {
        height <- length(features) * (2000/50) # 2000px for 50 genes
        height <- if(height < 350) 350 else height # min 350 px

        # Adjust when facet_ncol is 1
        if(input$multi.facet_ncol == 1 & length(features) > 6) {
          height <- if(length(features) < 25) height * 1.5 else if(length(features) < 75) height * 1.75 else height * 2
        }
      } else {
        height <- length(features) * (900/50) # 800px for 50 genes
        height <- if(height < 300) 300 else height # min 300 px

        # Adjust for larger base_size
        if(input$multi.base_size > base_size) {
	  height <- height + (15 * (input$multi.base_size - base_size))
	}
      }
      # Adjust for rotated cell groups
      if(input$multi.rotate_x != "No") height <- height + 200
      # Adjust for appended cell groups
      if(n_group_by > 1) height <- height + (n_group_by * 50)
    }
    height
  })

  # Create plot
  multi.prepare_plot <- eventReactive(input$multi.submit, {
    features <- multi.features()
    n_group_by <- length(ordered_plot_group())

    if(!is.null(features) & n_group_by > 0) {
      sce <- sce()
      sce.vars <- sce.vars()
      sce.group_by <- sce.vars[["sce.group_by"]]
      keep <- multi.compare.expr()

      group_by <- sce.group_by[ordered_plot_group()]
      xlab <- paste(names(group_by), collapse = " + ")
      rotate_x_angle <- ifelse(input$multi.rotate_x == "90°", 90, ifelse(input$multi.rotate_x == "45°", 45, 0))

      # Create a column with unique name to store group_by results
      randStr <- paste0("BCF", as.numeric(Sys.time()))
      colData(sce)[, randStr] <- if(n_group_by == 1) if(is.factor(colData(sce)[, group_by])) droplevels(colData(sce)[, group_by]) else as.factor(colData(sce)[, group_by])
              else colData(sce)[, group_by] %>% as.data.frame() %>% tidyr::unite(Group, sep = " - ") %>% pull(Group) %>% as.factor()
      colData(sce)[, randStr] <- factor(colData(sce)[, randStr], levels = gtools::mixedsort(levels(colData(sce)[, randStr])))

      if(input$multi.plot_type == "Dot") {
        if(input$multi.plot_cluster != "None") {
          # Do clustering (through plotGroupedHeatmap())
          # clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", clustering_method = "ward.D2"
          heat <- if(input$multi.plot_scale == "None") 
                  plotGroupedHeatmap(sce, features = features, group = randStr, clustering_method = "ward.D2", silent = T)
                  else plotGroupedHeatmap(sce, features = features[keep], group = randStr, clustering_method = "ward.D2",
                                          center = TRUE, scale = TRUE, silent = T)
        }

        fig <- if(input$multi.plot_scale == "None") {
                plotDots(sce, features = features, group = randStr) +
                        guides(colour = guide_colorbar(title = "Average Expression", barwidth = 20))
        } else {
                plotDots(sce, features = features[keep], group = randStr, center = TRUE, scale = TRUE) +
                        guides(colour = guide_colorbar(title = "Row (Gene) Z-Score", barwidth = 20))
        }

        if(input$multi.plot_cluster %in% c("None","Group (X-axis)")) {
          fig <- if(input$multi.plot_scale == "None") fig + scale_y_discrete(limits = features) 
                  else fig + scale_y_discrete(limits = features[keep])
        }

        if(input$multi.plot_cluster == "Group (X-axis)") {
          fig <- fig + scale_x_discrete(limits = heat$tree_col$labels[heat$tree_col$order])
        } else if(input$multi.plot_cluster == "Gene (Y-axis)") {
          fig <- fig + scale_y_discrete(limits = heat$tree_row$labels[heat$tree_row$order])
        } else if(input$multi.plot_cluster == "Both") {
          fig <- fig + scale_x_discrete(limits = heat$tree_col$labels[heat$tree_col$order]) + 
                  scale_y_discrete(limits = heat$tree_row$labels[heat$tree_row$order])
        }

        fig <- fig + scale_size(range = c(1, input$multi.dots_size), limits = c(0, 1)) +
                guides(size = guide_legend(title = "Proportion Detected")) +
                theme_minimal(base_size = input$multi.base_size) + theme(legend.position = "top") + xlab(xlab) + ylab("Genes")
      } else if(input$multi.plot_type == "Heatmap") {
        fig <- if(input$multi.plot_scale == "None")
                plotGroupedHeatmap(sce, features = features, group = randStr, clustering_method = "ward.D2", border_color = "black",
                                   fontsize = input$multi.base_size, angle_col = rotate_x_angle)
                else plotGroupedHeatmap(sce, features = features[keep], group = randStr, clustering_method = "ward.D2", border_color = "black",
                                        fontsize = input$multi.base_size, angle_col = rotate_x_angle, center = TRUE, scale = TRUE)
      } else {
        expr <- t(logcounts(sce[features,]))
        expr <- if(class(expr) == "dgCMatrix") as.data.frame(as.matrix(expr)) else as.data.frame(expr)
        colnames(expr) <- features

        if(input$multi.plot_type == "Boxplot") {
          # Prepare DataFrame
          df <- cbind(expr, data.frame(Group = colData(sce)[, randStr]))

          # Compute summary statistics for groups of cells
          summarized <- scuttle::summarizeAssayByGroup(assay(sce, "logcounts")[as.character(features), , drop = FALSE],
                                                       ids = df$Group, statistics = c("prop.detected"), threshold = 0)
          num <- data.frame(Symbol = rownames(summarized), assay(summarized, "prop.detected"), check.names = FALSE) %>%
                  tidyr::gather(key = "Group", value = "prop.detected", -Symbol)

          # Change wide to long format
          df <- tidyr::gather(df, key = "Symbol", value = "Expression", -Group) %>%
                  mutate_at(vars(Symbol), factor) %>% mutate(Symbol = factor(Symbol, levels = features))

          # Combind expr and prop.detected
          df <- merge(df, num, by = c("Group","Symbol"))

          # Add number of cells to group labels
          levels(df$Group) <- paste0(levels(df$Group), " (", dplyr::count(df, Group)$n / length(features), ")")

	  # Create plot
          box.color.fill <- if(input$multi.color_by == "Detected") "prop.detected" else "Group"
          my.aes <- aes(x = .data[["Group"]], y = .data[["Expression"]], color = .data[[box.color.fill]], fill = .data[[box.color.fill]])
          fig <- ggplot(df, my.aes) + geom_boxplot(outlier.size = 0.5, alpha = 0.3) + 
                  facet_wrap(~ Symbol, scales = "free_y", ncol = input$multi.facet_ncol) +
                  cowplot::theme_cowplot() + xlab(xlab) + ylab("logcounts") + 
                  theme(axis.text.x = element_text(size = input$multi.fontsize.x))

          if(input$multi.color_by == "Detected") fig <- fig + scale_color_gradientn(colours = continuous[["turbo"]]) +
                  scale_fill_gradientn(colours = continuous[["turbo"]]) +
                  guides(color = guide_colorbar(title = "Proportion\nDetected", barheight = 15),
                         fill = guide_colorbar(title = "Proportion\nDetected", barheight = 15))
        } else if(input$multi.plot_type == "Projection") {
          # Prepare DataFrame
          df <- as.data.frame(reducedDim(sce, input$multi.dimred)[,1:2])
          colnames(df) <- c("X","Y")
          df <- cbind(df, expr)

          # Change wide to long format
          df <- tidyr::gather(df, key = "Symbol", value = "Expression", -c(X, Y) ) %>% 
		  mutate_at(vars(Symbol), factor) %>% mutate(Symbol = factor(Symbol, levels = features))

          # Create plot
          fig <- ggplot(df, aes(x = X, y = Y, color = Expression)) + 
		  geom_point(size = input$multi.dot_size, alpha = input$multi.dot_opacity) + 
		  facet_wrap(~ Symbol, ncol = 5) + cowplot::theme_cowplot() + 
		  scale_color_gradientn(colours =  continuous[[input$multi.palette]]) + 
		  xlab(paste(input$multi.dimred, "1")) + ylab(paste(input$multi.dimred, "2"))

          if(length(features) < 5) {
            pl <- list(fig, empty)
            if(length(features) == 1) fig <- cowplot::plot_grid(plotlist = pl, nrow = 1, rel_widths = c(0.27, 0.73)) 
	    else if(length(features) == 2) fig <- cowplot::plot_grid(plotlist = pl, nrow = 1, rel_widths = c(0.45, 0.55)) 
	    else if(length(features) == 3) fig <- cowplot:: plot_grid(plotlist = pl, nrow = 1, rel_widths = c(0.63, 0.37)) 
	    else fig <- cowplot::plot_grid(plotlist = pl, nrow = 1, rel_widths = c(0.83, 0.17))
	  }
	}
      }

      if(input$multi.plot_type %in% c("Dot", "Boxplot")) {
        if(rotate_x_angle == "90") fig <- fig + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
        else if(rotate_x_angle == "45") fig <- fig + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
      }
      fig
    }
  })

  # Determine plot type and show plot with dynamic height
  output$multi.plot.ui <- renderUI({
    sce.vars <- sce.vars()
    sce.dimreds <- sce.vars[["sce.dimreds"]]
    plot_type <- isolate(input$multi.plot_type)

    if(multi.plot_height() > 0 & reset.multi.submit() == 0) {
      if((plot_type != "Projection") | (plot_type == "Projection" & length(sce.dimreds) > 0)) {
        output$multi.plot <- renderPlot({
          multi.prepare_plot()
        })

        box(title = plot_type, status = "primary", width = 12,
            plotOutput("multi.plot", width = "99%", height = multi.plot_height()) %>% withSpinner(type = getOption("spinner.type", default = 8))
        )
      } else {
        box(title = "Gene expression: Projection", width = 12, status = "primary", solidHeader = TRUE,
            column(width = 12, p("There is no dimensional reduction available on this data.")))
      }
    }
  })

  ####################
  # enrichR: findMarkers/edgeR/DESeq2 panel
  ####################
  prep.ora.listnames <- function(sce, menuItems) {
    if(menuItems == "ORAfm") {
      listnames <- names(metadata(sce))[grep("^enrichR_findMarkers", names(metadata(sce)))]
      names(listnames) <- gsub("enrichR_findMarkers_", "", listnames)
      type <- "markers"
    }

    if(menuItems == "ORAer") {
      listnames <- names(metadata(sce))[grep("^enrichR_edgeR", names(metadata(sce)))]
      names(listnames) <- gsub("enrichR_edgeR_", "", listnames)
      type <- "DEGs"
    }

    if(menuItems == "ORAde") {
      listnames <- names(metadata(sce))[grep("^enrichR_DESeq2", names(metadata(sce)))]
      names(listnames) <- gsub("enrichR_DESeq2_", "", listnames)
      type <- "DEGs"
    }

    names(listnames) <- gsub("_up$", paste0(" (up-regulated ", type, ")"), names(listnames))
    names(listnames) <- gsub("_dn$|_down$", paste0(" (down-regulated ", type, ")"), names(listnames))
    listnames
  }

  observeEvent(input$ORAfm.listname, {
    sce <- sce()
    listnames <- prep.ora.listnames(sce, "ORAfm")
    listname <- listnames[input$ORAfm.listname]
    groups <- names(metadata(sce)[[listname]])
    updateSelectInput(session, "ORAfm.group", choices = groups, selected = groups[1])
  })

  observeEvent(input$ORAfm.group, {
    sce <- sce()
    listnames <- prep.ora.listnames(sce, "ORAfm")
    listname <- listnames[input$ORAfm.listname]
    dbs <- names(metadata(sce)[[listname]][[input$ORAfm.group]])
    updateSelectInput(session, "ORAfm.db", choices = dbs, selected = dbs[1])
  })

  observeEvent(input$ORAer.listname, {
    sce <- sce()
    listnames <- prep.ora.listnames(sce, "ORAer")
    listname <- listnames[input$ORAer.listname]
    groups <- names(metadata(sce)[[listname]])
    updateSelectInput(session, "ORAer.group", choices = groups, selected = groups[1])
  })

  observeEvent(input$ORAer.group, {
    sce <- sce()
    listnames <- prep.ora.listnames(sce, "ORAer")
    listname <- listnames[input$ORAer.listname]
    dbs <- names(metadata(sce)[[listname]][[input$ORAer.group]])
    updateSelectInput(session, "ORAer.db", choices = dbs, selected = dbs[1])
  })

  observeEvent(input$ORAde.listname, {
    sce <- sce()
    listnames <- prep.ora.listnames(sce, "ORAde")
    listname <- listnames[input$ORAde.listname]
    groups <- names(metadata(sce)[[listname]])
    updateSelectInput(session, "ORAde.group", choices = groups, selected = groups[1])
  })

  observeEvent(input$ORAde.group, {
    sce <- sce()
    listnames <- prep.ora.listnames(sce, "ORAde")
    listname <- listnames[input$ORAde.listname]
    dbs <- names(metadata(sce)[[listname]][[input$ORAde.group]])
    updateSelectInput(session, "ORAde.db", choices = dbs, selected = dbs[1])
  })

  output$ORAfm.menu.ui <- renderUI({
    sce <- sce()
    listnames <- prep.ora.listnames(sce, "ORAfm")
    # Default, first element
    groups <- names(metadata(sce)[[listnames[1]]])
    dbs <- names(metadata(sce)[[listnames[1]]][[groups[1]]])

    box(title = "Input", width = 12, status = "primary", solidHeader = TRUE,
      fluidRow(
        column(width = 3, selectInput("ORAfm.listname", "Select comparison:", choices = names(listnames), selected = names(listnames)[1], multiple = FALSE, selectize = FALSE)),
        column(width = 3, selectInput("ORAfm.group", "Select cluster/condition:", choices = groups, selected = groups[1], multiple = FALSE, selectize = FALSE)),
        column(width = 3, selectInput("ORAfm.db", "Select gene-set library:", choices = dbs, selected = dbs[1], multiple = FALSE, selectize = FALSE)),
        column(width = 3, sliderInput("ORAfm.nterms", "Number of terms:", 10, 30, 20, 5), style = "height: 70px;")
      ),
      fluidRow(
        column(width = 3, radioGroupButtons("ORAfm.xaxis", label = "Bar length by:", direction = "horizontal", choices = c("Gene count", "Gene ratio"),
                                            selected = "Gene count")),
        column(width = 3, radioGroupButtons("ORAfm.order_by", label = "Order bars by:", direction = "horizontal", choices = c("P-value", "FDR", "Combined score"),
                                            selected = "FDR")),
        column(width = 3, radioGroupButtons("ORAfm.color_by", label = "Colour bars by:", direction = "horizontal", choices = c("P-value", "FDR", "Combined score"),
                                            selected = "P-value")),
        column(width = 3, actionButton("ORAfm.submit", "Submit", width = "75px", class = "btn-primary", style = "color: #fff; margin-top: 25px;"))
      )
    )
  })

  output$ORAer.menu.ui <- renderUI({
    sce <- sce()
    listnames <- prep.ora.listnames(sce, "ORAer")
    # Default, first element
    groups <- names(metadata(sce)[[listnames[1]]])
    dbs <- names(metadata(sce)[[listnames[1]]][[groups[1]]])

    box(title = "Input", width = 12, status = "primary", solidHeader = TRUE,
      fluidRow(
        column(width = 3, selectInput("ORAer.listname", "Select DE condition:", choices = names(listnames), selected = names(listnames)[1], multiple = FALSE, selectize = FALSE)),
        column(width = 3, selectInput("ORAer.group", "Select DE result:", choices = groups, selected = groups[1], multiple = FALSE, selectize = FALSE)),
        column(width = 3, selectInput("ORAer.db", "Select gene-set library:", choices = dbs, selected = dbs[1], multiple = FALSE, selectize = FALSE)),
        column(width = 3, sliderInput("ORAer.nterms", "Number of terms:", 10, 30, 20, 5), style = "height: 70px;")
      ),
      fluidRow(
        column(width = 3, radioGroupButtons("ORAer.xaxis", label = "Bar length by:", direction = "horizontal", choices = c("Gene count", "Gene ratio"),
                                            selected = "Gene count")),
        column(width = 3, radioGroupButtons("ORAer.order_by", label = "Order bars by:", direction = "horizontal", choices = c("P-value", "FDR", "Combined score"),
                                            selected = "FDR")),
        column(width = 3, radioGroupButtons("ORAer.color_by", label = "Colour bars by:", direction = "horizontal", choices = c("P-value", "FDR", "Combined score"),
                                            selected = "P-value")),
        column(width = 3, actionButton("ORAer.submit", "Submit", width = "75px", class = "btn-primary", style = "color: #fff; margin-top: 25px;"))
      )
    )
  })

  output$ORAde.menu.ui <- renderUI({
    sce <- sce()
    listnames <- prep.ora.listnames(sce, "ORAde")
    # Default, first element
    groups <- names(metadata(sce)[[listnames[1]]])
    dbs <- names(metadata(sce)[[listnames[1]]][[groups[1]]])

    box(title = "Input", width = 12, status = "primary", solidHeader = TRUE,
      fluidRow(
        column(width = 3, selectInput("ORAde.listname", "Select DE condition:", choices = names(listnames), selected = names(listnames)[1], multiple = FALSE, selectize = FALSE)),
        column(width = 3, selectInput("ORAde.group", "Select DE result:", choices = groups, selected = groups[1], multiple = FALSE, selectize = FALSE)),
        column(width = 3, selectInput("ORAde.db", "Select gene-set library:", choices = dbs, selected = dbs[1], multiple = FALSE, selectize = FALSE)),
        column(width = 3, sliderInput("ORAde.nterms", "Number of terms:", 10, 30, 20, 5), style = "height: 70px;")
      ),
      fluidRow(
        column(width = 3, radioGroupButtons("ORAde.xaxis", label = "Bar length by:", direction = "horizontal", choices = c("Gene count", "Gene ratio"),
                                            selected = "Gene count")),
        column(width = 3, radioGroupButtons("ORAde.order_by", label = "Order bars by:", direction = "horizontal", choices = c("P-value", "FDR", "Combined score"),
                                            selected = "FDR")),
        column(width = 3, radioGroupButtons("ORAde.color_by", label = "Colour bars by:", direction = "horizontal", choices = c("P-value", "FDR", "Combined score"),
                                            selected = "P-value")),
        column(width = 3, actionButton("ORAde.submit", "Submit", width = "75px", class = "btn-primary", style = "color: #fff; margin-top: 25px;"))
      )
    )
  })

  ora.fm <- eventReactive(input$ORAfm.submit, {
    isolate({
      listname <- input$ORAfm.listname
      group <- input$ORAfm.group
      db <- input$ORAfm.db
      nterms <- input$ORAfm.nterms
      xaxis <- input$ORAfm.xaxis
      order_by <- input$ORAfm.order_by
      color_by <- input$ORAfm.color_by
    })
    list(listname = listname, group = group, db = db, nterms = nterms, xaxis = xaxis, order_by = order_by, color_by = color_by)
  })

  ora.er <- eventReactive(input$ORAer.submit, {
    isolate({
      listname <- input$ORAer.listname
      group <- input$ORAer.group
      db <- input$ORAer.db
      nterms <- input$ORAer.nterms
      xaxis <- input$ORAer.xaxis
      order_by <- input$ORAer.order_by
      color_by <- input$ORAer.color_by
    })
    list(listname = listname, group = group, db = db, nterms = nterms, xaxis = xaxis, order_by = order_by, color_by = color_by)
  })

  ora.de <- eventReactive(input$ORAde.submit, {
    isolate({
      listname <- input$ORAde.listname
      group <- input$ORAde.group
      db <- input$ORAde.db
      nterms <- input$ORAde.nterms
      xaxis <- input$ORAde.xaxis
      order_by <- input$ORAde.order_by
      color_by <- input$ORAde.color_by
    })
    list(listname = listname, group = group, db = db, nterms = nterms, xaxis = xaxis, order_by = order_by, color_by = color_by)
  })

  prep.ora.df <- function(sce, menuItems) {
    ora.listnames <- prep.ora.listnames(sce, menuItems)
    if(menuItems == "ORAfm") {
      ora.vars <- ora.fm()
    } else if(menuItems == "ORAer") {
      ora.vars <- ora.er()
    } else {
      ora.vars <- ora.de()
    }
    listname <- ora.vars[["listname"]]
    group <- ora.vars[["group"]]
    db <- ora.vars[["db"]]
    order_by <- ora.vars[["order_by"]]

    ora.listnames <- prep.ora.listnames(sce, menuItems)
    listname <- ora.listnames[listname]
    df <- metadata(sce)[[listname]][[group]][[db]]
    df$nGenes <- as.numeric(gsub("([0-9]+)/[0-9]+", "\\1", df$Overlap))
    df$rGenes <- df$nGenes/as.numeric(gsub("[0-9]+/([0-9]+)", "\\1", df$Overlap))
    df$Overlap <- paste0(sprintf("%.4f", df$rGenes), " (", df$Overlap, ")")
    if(order_by == "P-value") {
      df <- df %>% arrange(P.value)
    } else if(order_by == "FDR") {
      df <- df %>% arrange(Adjusted.P.value)
    } else {
      df <- df %>% arrange(desc(Combined.Score))
    }
    df
  }

  prep.ora.plot <- function(df, menuItems) {
    if(menuItems == "ORAfm") {
      ora.vars <- ora.fm()
    } else if(menuItems == "ORAer") {
      ora.vars <- ora.er()
    } else {
      ora.vars <- ora.de()
    }
    nterms <- ora.vars[["nterms"]]
    xaxis <- ora.vars[["xaxis"]]
    color_by <- ora.vars[["color_by"]]

    df <- head(df, nterms)
    df$Term <- factor(df$Term, levels = rev(df$Term))
    bar.x <- ifelse(xaxis == "Gene count", "nGenes", "rGenes")
    bar.fill <- ifelse(color_by == "P-value", "P.value", ifelse(color_by == "FDR", "Adjusted.P.value", "Combined.Score"))
    fig <- ggplot(df, aes(.data[[bar.x]], .data[["Term"]], fill = .data[[bar.fill]])) + geom_col() + cowplot::theme_cowplot() + xlab(xaxis)
    if(color_by == "Combined score") {
      fig + scale_fill_continuous(low = "blue", high = "red") + guides(fill = guide_colorbar(title = color_by, reverse = FALSE))
    } else {
      fig + scale_fill_continuous(low = "red", high = "blue") + guides(fill = guide_colorbar(title = color_by, reverse = TRUE))
    }
  }

  prep.ora.table <- function(df) {
    df <- df[, colnames(df) %in% c("Term","Overlap","P.value","Adjusted.P.value","Odds.Ratio","Combined.Score","Genes")]
    df$Genes <- gsub(";", ", ", df$Genes)
    datatable(df, options = list(searching = TRUE, pageLength = 10, scrollX = TRUE, lengthChange = FALSE),
              rownames = FALSE, selection = "none", class = "white-space: nowrap") %>%
    formatRound(columns = c("Odds.Ratio","Combined.Score"), digits = 4) %>%
    formatSignif(columns = c("P.value","Adjusted.P.value"), digits = 4)
  }

  output$ORAfm.output.ui <- renderUI({
    if(reset.ORAfm.submit() == 0) {
      sce <- sce()
      df <- prep.ora.df(sce, "ORAfm")

      output$ORAfm.plot <- renderPlot({
        prep.ora.plot(df, "ORAfm")
      })

      box(title = "Result", status = "primary", width = 12, solidHeader = TRUE, collapsible = FALSE,
          plotOutput("ORAfm.plot", width = "99%", height = "400px") %>% withSpinner(type = getOption("spinner.type", default = 8)),
          renderDT(prep.ora.table(df)))
    }
  })

  output$ORAer.output.ui <- renderUI({
    if(reset.ORAer.submit() == 0) {
      sce <- sce()
      df <- prep.ora.df(sce, "ORAer")

      output$ORAer.plot <- renderPlot({
        prep.ora.plot(df, "ORAer")
      })

      box(title = "Result", status = "primary", width = 12, solidHeader = TRUE, collapsible = FALSE,
          plotOutput("ORAer.plot", width = "99%", height = "400px") %>% withSpinner(type = getOption("spinner.type", default = 8)),
          renderDT(prep.ora.table(df)))
    }
  })

  output$ORAde.output.ui <- renderUI({
    if(reset.ORAde.submit() == 0) {
      sce <- sce()
      df <- prep.ora.df(sce, "ORAde")

      output$ORAde.plot <- renderPlot({
        prep.ora.plot(df, "ORAde")
      })

      box(title = "Result", status = "primary", width = 12, solidHeader = TRUE, collapsible = FALSE,
          plotOutput("ORAde.plot", width = "99%", height = "400px") %>% withSpinner(type = getOption("spinner.type", default = 8)),
          renderDT(prep.ora.table(df)))
    }
  })
}

####################
# Call to the shinyApp function
####################
shinyApp(ui = ui, server = server)
