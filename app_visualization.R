library(DT)
library(ape)
library(dplyr)
library(shiny)
library(ggtree)
library(ggplot2)
library(phytools)
library(phangorn)
library(tidytree)
library(TreeDist) 
library(TreeTools)
library(shinydashboard)


options(shiny.maxRequestSize = 100 * 1024^2)


robust_clean_name <- function(x) {
  if(is.null(x)) return(x)
  x <- as.character(x)
  x <- gsub("['\"]", "", x)
  x <- gsub("[ _]+", " ", x)
  x <- trimws(x)
  return(x)
}

# Clear ape::phylo
force_clean_phylo <- function(tr) {
  tr$tip.label <- robust_clean_name(tr$tip.label)
  tr$node.label <- NULL
  tr$edge.length <- NULL
  tr$boot <- NULL
  txt <- write.tree(tr)
  tr_clean <- read.tree(text = txt)
  tr_clean <- multi2di(tr_clean) # Binary tree
  return(tr_clean)
}


ln_double_fact <- function(n) {
  if (n < 3) return(0)
  lfactorial(2 * n - 5) - ((n - 3) * log(2) + lfactorial(n - 3))
}


split_information <- function(k, n) {
  if(k <= 0 || k >= n) return(0)
  log_prob <- (ln_double_fact(k) + ln_double_fact(n - k)) - ln_double_fact(n)
  return(max(0, -(log_prob / log(2))))
}


calculate_spi_metrics <- function(ref_tree, comp_tree) {
  t1 <- force_clean_phylo(ref_tree)
  t2 <- force_clean_phylo(comp_tree)
  
  common <- sort(intersect(t1$tip.label, t2$tip.label))
  if(length(common) < 4) return(list(raw=0, norm=0, jaccard=0))
  
  t1 <- keep.tip(t1, common)
  t2 <- keep.tip(t2, common)
  n <- length(common)
  t1 <- force_clean_phylo(t1)
  t2 <- force_clean_phylo(t2)
  bp1 <- prop.part(t1)
  bp2 <- prop.part(t2)
  
  get_splits <- function(bp, lab) {
    sapply(bp, function(x) {
      cl <- lab[x]
      if (length(cl) > n/2) cl <- setdiff(lab, cl)
      paste(sort(cl), collapse="|")
    })
  }
  
  s1 <- get_splits(bp1, attr(bp1, "labels"))
  s2 <- get_splits(bp2, attr(bp2, "labels"))
  
  # Info Ref
  info_ref <- sum(sapply(s1, function(s) {
    if(s=="") return(0)
    split_information(stringr::str_count(s, "\\|")+1, n)
  }))
  
  # Info Shared
  shared <- intersect(s1, s2)
  info_shared <- sum(sapply(shared, function(s) {
    if(s=="") return(0)
    split_information(stringr::str_count(s, "\\|")+1, n)
  }))
  
  # Jaccard Index
  u <- unique(c(s1, s2))
  jac <- if(length(u)>0) length(shared)/length(u) else 0
  
  norm_val <- if(info_ref>0) info_shared/info_ref else 0
  
  return(list(
    raw = round(info_shared, 2),
    norm = round(norm_val, 4),
    jaccard = round(jac, 4)
  ))
}

# UI
ui <- dashboardPage(
  dashboardHeader(title = "SAPP-Vizualizer"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("1. Upload & Settings", tabName = "setup_tab", icon = icon("upload")),
      menuItem("2. Analysis", tabName = "analysis_tab", icon = icon("cogs")),
      menuItem("3. Visualizations", tabName = "viz_tab", icon = icon("eye")),
      menuItem("4. Report", tabName = "report_tab", icon = icon("file-alt"))
    ),
    hr(),
    fileInput("ref_tree_file", "Reference Tree (.nwk)", accept = c(".nwk", ".newick", ".treefile")),
    fileInput("meta_file", "Metadata (.csv)", accept = ".csv"),
    fileInput("comp_files", "Comparison Trees", multiple = TRUE, accept = c(".nwk", ".newick", ".treefile"))
  ),
  dashboardBody(
    tabItems(
      # TAB 1
      tabItem(tabName = "setup_tab",
              fluidRow(
                box(title = "Metadata Preview", status = "primary", width = 12, DTOutput("meta_preview")),
                box(title = "Reference Tree Preview", status = "info", width = 12,
                    plotOutput("ref_plot_raw", height = "500px"))
              )
      ),

      tabItem(tabName = "analysis_tab",
              fluidRow(
                box(title = "Configuration", status = "warning", solidHeader = TRUE, width = 12,
                    h4("Rooting & Unification"),
                    verbatimTextOutput("rooting_info"),
                    fluidRow(
                      column(6, selectInput("map_id_col", "ID Column:", choices = NULL)),
                      column(6, selectInput("map_target_col", "Species Column:", choices = NULL))
                    ),
                    checkboxInput("clean_underscores", "Convert '_' to spaces", value = TRUE),
                    actionButton("run_unification", "Run Analysis Pipeline", class = "btn-success", icon = icon("play")),
                    hr(),
                    verbatimTextOutput("unification_status")
                )
              )
      ),

      tabItem(tabName = "viz_tab",
              tabsetPanel(
                tabPanel("Reference Tree", 
                         uiOutput("ref_color_selector"),
                         plotOutput("ref_plot_final", height = "800px"),
                         downloadButton("dl_ref_png", "PNG")),
                tabPanel("Tanglegram", 
                         uiOutput("comp_tree_selector_tang"),
                         uiOutput("tang_color_selector"),
                         plotOutput("tanglegram_plot", height = "900px"),
                         downloadButton("dl_tang_png", "PNG")),
                tabPanel("Detailed Phylogram", 
                         uiOutput("detail_tree_selector"),
                         uiOutput("detail_color_selector"),
                         plotOutput("detail_plot", height = "800px"),
                         downloadButton("dl_det_png", "PNG")),
                tabPanel("SPI Info (TreeDist)", 
                         p("Visualization using TreeDist::VisualizeMatching"),
                         p("Green = Matched Splits, Red = Conflict."),
                         uiOutput("spi_tree_selector"),
                         plotOutput("spi_plot", height = "1000px"),
                         downloadButton("dl_spi_png", "PNG"))
              )
      ),

      tabItem(tabName = "report_tab",
              box(width = 12, title = "Results", status = "success",
                  p("Metrics included: Robinson-Foulds (RF), Normalized RF, SPI (Shared Phylogenetic Info), Jaccard (Splits), Path Distance."),
                  downloadButton("download_txt", "Download CSV"),
                  DTOutput("metrics_table"))
      )
    )
  )
)

server <- function(input, output, session) {
  
  # Data
  meta_data <- reactive({
    req(input$meta_file)
    tryCatch({ read.csv(input$meta_file$datapath, stringsAsFactors = FALSE) }, error = function(e) NULL)
  })
  
  observeEvent(meta_data(), {
    cols <- colnames(meta_data())
    id_guess <- cols[grepl("id", tolower(cols))][1]
    sp_guess <- cols[grepl("species|name", tolower(cols))][1]
    updateSelectInput(session, "map_id_col", choices = cols, selected = ifelse(!is.na(id_guess), id_guess, cols[1]))
    updateSelectInput(session, "map_target_col", choices = cols, selected = ifelse(!is.na(sp_guess), sp_guess, cols[2]))
  })
  
  output$rooting_info <- renderText({
    req(meta_data(), input$map_target_col)
    s <- robust_clean_name(meta_data()[1, input$map_target_col])
    if(input$clean_underscores) s <- gsub("_", " ", s)
    paste("Outgroup:", s)
  })
  
  # Raw Reference
  raw_ref_tree <- reactive({ 
    req(input$ref_tree_file)
    tryCatch({
      tr <- read.tree(input$ref_tree_file$datapath)
      tr$tip.label <- robust_clean_name(tr$tip.label)
      tr$node.label <- NULL # Clean garbage
      return(tr)
    }, error = function(e) return(NULL))
  })
  
  output$ref_plot_raw <- renderPlot({
    req(raw_ref_tree())
    ggtree(raw_ref_tree(), layout="rectangular", linewidth=0.5) + geom_tiplab() + theme_tree2()
  })
  
  # Pipeline SAPP-Vizualizer
  pipeline_res <- eventReactive(input$run_unification, {
    req(raw_ref_tree(), input$comp_files, meta_data())
    
    ref <- raw_ref_tree()
    meta <- meta_data()
    comp_files <- input$comp_files
    
    id_col <- input$map_id_col
    tgt_col <- input$map_target_col
    do_clean <- input$clean_underscores
    
    fmt_name <- function(x) {
      x <- robust_clean_name(x)
      if(do_clean) gsub("_", " ", x) else x
    }
    
    meta[[id_col]] <- robust_clean_name(meta[[id_col]])
    meta[[tgt_col]] <- fmt_name(meta[[tgt_col]])
    dict <- setNames(meta[[tgt_col]], meta[[id_col]])
    outgroup <- meta[1, tgt_col]
    
    ref$tip.label <- fmt_name(ref$tip.label)
    if(outgroup %in% ref$tip.label) {
      ref <- tryCatch(root(ref, outgroup=outgroup, resolve.root=TRUE), error=function(e) ref)
    }

    meta_clean <- meta
    meta_clean$label <- meta_clean[[tgt_col]]
    
    metrics_list <- list()
    comps_list <- list()
    trees_list <- list()
    
    withProgress(message = "Analyzing...", value = 0, {
      for(i in 1:nrow(comp_files)) {
        incProgress(1/nrow(comp_files), detail = comp_files$name[i])
        
        tr <- tryCatch(read.tree(comp_files$datapath[i]), error=function(e) NULL)
        if(is.null(tr)) next
        
        # Unification
        tr$tip.label <- robust_clean_name(tr$tip.label)
        matches <- dict[tr$tip.label]
        tr$tip.label[!is.na(matches)] <- matches[!is.na(matches)]
        tr$tip.label <- fmt_name(tr$tip.label)
        
        # Root
        if(outgroup %in% tr$tip.label) {
          tr <- tryCatch(root(tr, outgroup=outgroup, resolve.root=TRUE), error=function(e) tr)
        }
        
        trees_list[[comp_files$name[i]]] <- tr
        
        # Metrics
        common <- intersect(ref$tip.label, tr$tip.label)
        if(length(common) > 3) {
          t1 <- keep.tip(ref, common)
          t2 <- keep.tip(tr, common)
          
          # RF Metrics
          rf <- RF.dist(t1, t2)
          max_rf <- 2*(length(common)-3)
          
          # SPI Metrics
          spi <- calculate_spi_metrics(ref, tr)
          
          # Path Distance
          path_d <- path.dist(t1, t2)
          
          metrics_list[[i]] <- data.frame(
            File = comp_files$name[i],
            Taxa = length(common),
            RF_Dist = rf,
            Norm_RF = round(rf/max_rf, 4),
            SPI_Bits = spi$raw,
            SPI_Norm = spi$norm,
            Jaccard = spi$jaccard,
            Path_Dist = round(path_d, 2),
            stringsAsFactors = FALSE
          )
          
          comps_list[[comp_files$name[i]]] <- cophylo(t1, t2, rotate=TRUE)
        }
      }
    })
    
    final_metrics <- if(length(metrics_list)>0) do.call(rbind, metrics_list) else data.frame(Error="No Data")
    
    list(ref=ref, trees=trees_list, cophylos=comps_list, metrics=final_metrics, meta=meta_clean, msg="Analysis Complete.")
  })
  
  # Selectors
  output$ref_color_selector <- renderUI({ req(pipeline_res()); selectInput("col_ref", "Color By:", choices=names(pipeline_res()$meta)) })
  output$detail_color_selector <- renderUI({ req(pipeline_res()); selectInput("col_det", "Color By:", choices=names(pipeline_res()$meta)) })
  output$tang_color_selector <- renderUI({ req(pipeline_res()); selectInput("col_tang", "Color Lines By:", choices=names(pipeline_res()$meta)) })
  output$comp_tree_selector_tang <- renderUI({ req(pipeline_res()); selectInput("sel_tang", "Comparison:", choices=names(pipeline_res()$cophylos)) })
  output$detail_tree_selector <- renderUI({ req(pipeline_res()); selectInput("sel_det", "Tree:", choices=names(pipeline_res()$trees)) })
  output$spi_tree_selector <- renderUI({ req(pipeline_res()); selectInput("sel_spi", "Comparison:", choices=names(pipeline_res()$trees)) })
  output$unification_status <- renderText({ req(pipeline_res()); pipeline_res()$msg })
  
  # PLOTS
  
  # Reference
  output$ref_plot_final <- renderPlot({
    req(pipeline_res(), input$col_ref)
    dat <- pipeline_res()
    p <- ggtree(dat$ref, layout="rectangular", linewidth=0.8)
    p$data <- dplyr::left_join(p$data, dat$meta, by="label")
    p + geom_tiplab(align=TRUE) + 
      geom_tippoint(aes(color=.data[[input$col_ref]]), size=4) +
      theme_tree2() + scale_x_continuous(expand=expansion(mult=0.3))
  })
  
  # Tanglegram
  output$tanglegram_plot <- renderPlot({
    req(pipeline_res(), input$sel_tang, input$col_tang)
    dat <- pipeline_res()
    obj <- dat$cophylos[[input$sel_tang]]
    
    if(!is.null(obj)) {
      map_col <- setNames(dat$meta[[input$col_tang]], robust_clean_name(dat$meta$label))
      cats <- unique(na.omit(map_col))
      pal <- tryCatch(scales::hue_pal()(length(cats)), error=function(e) rainbow(length(cats)))
      names(pal) <- cats
      left_tips <- robust_clean_name(obj$trees[[1]]$tip.label)
      cols <- sapply(1:nrow(obj$assoc), function(i) {
        val <- map_col[left_tips[obj$assoc[i,1]]]
        if(!is.na(val) && val %in% names(pal)) pal[val] else "grey"
      })
      par(mar=c(5,5,5,5))
      plot(obj, link.type="curved", link.lwd=2, link.col=make.transparent(cols, 0.6), fsize=0.8)
      legend("bottomleft", legend=names(pal), fill=pal, title=input$col_tang, bty="n")
    }
  })
  
  # Shared phylogenetic information visualization
  output$spi_plot <- renderPlot({
    req(pipeline_res(), input$sel_spi)
    dat <- pipeline_res()
    
    ref_full <- dat$ref
    comp_full <- dat$trees[[input$sel_spi]]
    
    common_taxa <- intersect(ref_full$tip.label, comp_full$tip.label)
    if(length(common_taxa) < 4) { plot.new(); text(0.5,0.5,"<4 common taxa"); return() }
    
    t1 <- keep.tip(ref_full, common_taxa)
    t2 <- keep.tip(comp_full, common_taxa)
    
    # Cleaning
    t1 <- force_clean_phylo(t1)
    t2 <- force_clean_phylo(t2)
    
    tryCatch({
      TreeDist::VisualizeMatching(TreeDist::SharedPhylogeneticInfo, t1, t2, 
                                  matchZeros = FALSE,
                                  prediction = 16, # Mniej cyfr
                                  setPar = TRUE)
    }, error=function(e) {
      plot.new(); text(0.5, 0.5, paste("SPI Plot Error:", e$message), col="red")
    })
  })

  output$detail_plot <- renderPlot({
    req(pipeline_res(), input$sel_det, input$col_det)
    dat <- pipeline_res()
    tr <- dat$trees[[input$sel_det]]
    p <- ggtree(tr, layout="rectangular", linewidth=0.8)
    p$data <- dplyr::left_join(p$data, dat$meta, by="label")
    p + geom_tiplab(align=TRUE) + 
      geom_tippoint(aes(color=.data[[input$col_det]]), size=4) +
      theme_tree2() + scale_x_continuous(expand=expansion(mult=0.3))
  })
  
  # Exports
  output$meta_preview <- renderDT(meta_data(), options=list(scrollX=TRUE))
  output$metrics_table <- renderDT(pipeline_res()$metrics, options=list(scrollX=TRUE))
  output$download_txt <- downloadHandler("report.csv", function(file) write.csv(pipeline_res()$metrics, file))
  output$dl_ref_png <- downloadHandler("ref.png", function(file) ggsave(file, plot_ref_reactive(), width=12, height=10))
  output$dl_tang_png <- downloadHandler("tang.png", function(file) { png(file, width=2000, height=2000, res=150); draw_tanglegram(); dev.off() })
  output$dl_spi_png <- downloadHandler("spi.png", function(file) { png(file, width=2000, height=2000, res=150); draw_spi_plot(); dev.off() })
}

shinyApp(ui, server)
