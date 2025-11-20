library(shiny)
library(dplyr)
library(ggplot2)
library(patchwork)
library(data.table)
library(DT)
library(devtools)

# Load your locally stored NIRS package
devtools::load_all("/Users/michaelwolfe/Documents/GitHub/nirsdexata")
# --- Utility: drop an entire wavelength from rawdata before analysis ---
# --- Utility: drop one wavelength robustly ---
remove_wavelength <- function(rawdata, wave) {
  # rawdata$info$cols is a data.frame with one row per channel
  cols_df <- rawdata$info$cols
  if (!wave %in% cols_df$col_name) {
    stop("Wavelength not found: ", wave)
  }
  keep <- cols_df$col_name != wave
  
  # subset the data columns (columns align 1:1 with rows of cols_df)
  rawdata$data <- as.data.frame(rawdata$data)[ , keep, drop = FALSE]
  
  # drop the entire row from info$cols
  rawdata$info$cols <- cols_df[keep, , drop = FALSE]
  
  return(rawdata)
}


# --- Asymptotic Regression Helpers ---
asym_regression <- function(x, Asym, R0, lrc) {
  Asym - (Asym - R0) * exp(-exp(lrc) * x)
}
asym_regression_normalized <- function(x, Asym, R0, lrc) {
  curve <- asym_regression(x, Asym, R0, lrc)
  curve - curve[1]
}

# --- OOK Plot Helpers ---
prepare_for_ook_plot <- function(
    df, sfreq, meas_start, meas_end,
    bads, removed, col_types, col_names, type = NULL
) {
  prepare_shiny_plot(
    df, sfreq, meas_start, meas_end,
    bads, removed, col_types, col_names, type
  )
}
# --- OOK Plot Helpers (updated) ---
# --- OOK Plot Helpers (updated) ---
# --- OOK Plot Helper (with TSI as ratio) ---
ggook <- function(data, x = NULL, y = NULL) {
  required <- c("ZeroedTime", "value", "col_name", "col_type")
  if (!all(required %in% names(data))) stop("Data not compatible with ggook()")
  
  # Build a label: 
  #  - for TSI channels: "TSI (Ratio)"
  #  - for HBO/HBR:   "Wavelength (Oxy)" or "Wavelength (Deoxy)"
  data$WavLabel <- ifelse(
    data$col_type == "TSI",
    paste0(data$col_name, " (Ratio)"),
    paste0(
      data$col_name,
      " (",
      ifelse(data$col_type == "HBO", "Oxy", "Deoxy"),
      ")"
    )
  )
  data$WavLabel <- factor(data$WavLabel, levels = unique(data$WavLabel))
  
  plt <- ggplot(data, aes(x = ZeroedTime, y = value, color = WavLabel)) +
    geom_line(na.rm = TRUE) +
    theme_bw() +
    labs(
      x     = "Time (sec)",
      y     = "",
      color = "Wavelength\n(Type)"
    )
  
  major_breaks <- seq(0, 960, 120)
  plt <- plt + scale_x_continuous(breaks = major_breaks)
  if (!is.null(x) || !is.null(y)) {
    plt <- plt + coord_cartesian(xlim = x, ylim = y, expand = FALSE)
  }
  for (v in major_breaks) plt <- plt + geom_vline(xintercept = v)
  plt
}

# --- Regression Plot Helper (with TSI as its own facet row) ---
plot_regressions <- function(models) {
  df <- models$data %>% filter(!Section %in% c("BegRest", "WarmUp"))
  
  # Same label logic as above
  df$WavLabel <- ifelse(
    df$col_type == "TSI",
    paste0(df$col_name, " (Ratio)"),
    paste0(
      df$col_name,
      " (",
      ifelse(df$col_type == "HBO", "Oxy", "Deoxy"),
      ")"
    )
  )
  df$WavLabel <- factor(df$WavLabel, levels = unique(df$WavLabel))
  
  # Facet‐row labels
  type_labs <- c(
    HBO = "Oxyhemoglobin",
    HBR = "Deoxyhemoglobin",
    TSI = "TSI (Oxy/Deoxy Ratio)"
  )
  
  ggplot(df, aes(x = SectionZeroedTime, group = WavLabel)) +
    facet_grid(
      col_type ~ Section,
      scales   = "free_y",
      labeller = labeller(col_type = type_labs)
    ) +
    geom_line(aes(y = value, color = WavLabel), size = 1, alpha = 0.6) +
    geom_line(aes(y = predicted, color = WavLabel),
              size = 1, linetype = "dashed", na.rm = TRUE) +
    labs(
      x     = "Time Since Start of Section (sec)",
      y     = "",
      color = "Wavelength\n(Type)"
    ) +
    scale_x_continuous(breaks = c(0, 30, 60, 90, 120), limits = c(0, 120)) +
    theme_light()
}

plot_rawdata <- function(rawdata) {
  pd <- prepare_for_ook_plot(
    rawdata$data,
    rawdata$info$sfreq,
    rawdata$info$bounds$meas_start,
    rawdata$info$bounds$meas_end,
    rawdata$info$bads,
    rawdata$info$bounds$removed,
    rawdata$info$cols$col_type,
    rawdata$info$cols$col_name
  )
  ggook(pd)
}


# --- Cycle File Processor ---
processCycleFile <- function(path) {
  rd <- tryCatch(
    read_nirs(path),
    error = function(e) stop("Failed to read cycle file: ", e$message)
  )
  if (is.atomic(rd)) {
    df <- tryCatch(
      read.table(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE),
      error = function(e) read.csv(path, stringsAsFactors = FALSE)
    )
    infoObj <- create_info(col_names = names(df))
    rd <- new_rawdata(data = df, info = infoObj)
  }
  rd
}

# --- Cleaning Module UI & Server ---
interact_rawdata_ui <- function(id) {
  ns <- NS(id)
  tagList(
    h4("Interactive Data Cleaning"),
    fluidRow(
      column(6,
             plotOutput(ns("ook_plot"), height = "300px",
                        brush = ns("brush"), dblclick = ns("dblclick"))
      ),
      column(6,
             plotOutput(ns("ook_plot_zoom"), height = "300px",
                        brush = ns("brush2"))
      )
    ),
    fluidRow(
      column(4, actionButton(ns("reset"), "Reset")),
      column(4, actionButton(ns("remove"), "Remove")),
      column(4, actionButton(ns("undo"),   "Undo"))
    ),
    fluidRow(
      column(12, actionButton(ns("done"), "Done", class = "btn-primary"))
    ),
    fluidRow(
      column(12, verbatimTextOutput(ns("rollingTau")))
    )
  )
}

interact_rawdata_server <- function(input, output, session,
                                    rawdata, subjID, cycleNumber) {
  vals <- reactiveValues(
    meas_start = NULL,
    meas_end   = NULL,
    bads       = NULL,
    removed    = NULL,
    removal_history = list()
  )
  
  observeEvent(rawdata(), ignoreNULL = TRUE, {
    rd <- rawdata()
    vals$meas_start <- meas_start(rd)
    vals$meas_end   <- meas_end(rd)
    vals$bads       <- bads(rd)
    vals$removed    <- rd$info$bounds$removed
  })
  
  ranges <- reactiveValues(xmin = NULL, xmax = NULL, ymin = NULL, ymax = NULL)
  observeEvent(input$dblclick, {
    b <- input$brush
    if (!is.null(b)) {
      ranges$xmin <- b$xmin; ranges$xmax <- b$xmax
      ranges$ymin <- b$ymin; ranges$ymax <- b$ymax
    } else {
      ranges[] <- NULL
    }
  })
  
  plotdata <- reactive({
    prepare_shiny_plot(
      rawdata()$data,
      sfreq       = rawdata()$info$sfreq,
      meas_start  = vals$meas_start,
      meas_end    = vals$meas_end,
      bads        = vals$bads,
      removed     = vals$removed,
      col_types   = rawdata()$info$cols$col_type,
      col_names   = rawdata()$info$cols$col_name,
      type        = NULL
    )
  })
  plotdata_zoom <- reactive({
    # just return the full dataset; we’ll apply zoom in the plotting step
    plotdata()
  })
  
  
  output$ook_plot      <- renderPlot(ggook(plotdata()), res = 96)
  output$ook_plot_zoom <- renderPlot({
    data <- plotdata()
    if (is.null(ranges$xmin)) {
      ggook(data)
    } else {
      ggook(
        data,
        x = c(ranges$xmin, ranges$xmax),
        y = c(ranges$ymin, ranges$ymax)
      )
    }
  }, res = 96)
  
  
  # Rolling Tau updates reactively
  output$rollingTau <- renderText({
    rd2 <- rawdata()
    rd2$info$bounds$meas_start <- vals$meas_start
    rd2$info$bounds$meas_end   <- vals$meas_end
    rd2$info$bounds$removed    <- vals$removed
    rd2$info$bads              <- vals$bads
    
    mod <- tryCatch(analysis(rd2), error = function(e) NULL)
    if (is.null(mod)) return("Rolling Tau: n/a")
    coefs <- mod$coefs %>%
      filter(col_type == "HBO", Section %in% c("Rest1","Rest2","Rest3"))
    paste0("Rolling Tau: ", round(mean(coefs$Tau, na.rm = TRUE), 2), " sec")
  })
  
  observeEvent(input$reset, {
    rd <- rawdata()
    vals$meas_start <- meas_start(rd)
    vals$meas_end   <- meas_end(rd)
    vals$bads       <- bads(rd)
    vals$removed    <- rd$info$bounds$removed
    vals$removal_history <- list()
  })
  observeEvent(input$remove, {
    pts    <- brushedPoints(plotdata_zoom(), input$brush2)
    newids <- unique(pts$samp_num)
    vals$removal_history <- c(vals$removal_history, list(vals$removed))
    vals$removed <- union(vals$removed, newids)
  })
  observeEvent(input$undo, {
    if (length(vals$removal_history) > 0) {
      last <- tail(vals$removal_history,1)[[1]]
      vals$removal_history <- head(vals$removal_history,-1)
      vals$removed <- last
    }
  })
  
  observeEvent(input$done, ignoreInit = TRUE, {
    rd <- rawdata()
    rd$info$bounds$meas_start <- vals$meas_start
    rd$info$bounds$meas_end   <- vals$meas_end
    rd$info$bounds$removed    <- vals$removed
    rd$info$bads              <- vals$bads
    session$userData$cleaned  <- rd
    
    # append removed
    newRemoved <- data.frame(
      ID           = as.numeric(subjID()),
      Cycle_number = as.numeric(cycleNumber()),
      Removed      = I(list(vals$removed)),
      Name         = paste0("id",subjID(),"c",cycleNumber()),
      stringsAsFactors = FALSE
    )
    globalRemovedData(rbind(globalRemovedData(), newRemoved))
    
    # append summary
    mod <- tryCatch(analysis(rd), error = function(e) {
      showModal(modalDialog(
        title = "Analysis Error",
        paste("Failed to analyze cycle:", e$message),
        easyClose = TRUE
      ))
      NULL
    })
    if (!is.null(mod)) {
      coefs <- mod$coefs %>%
        filter(col_type=="HBO", Section %in% c("Rest1","Rest2","Rest3")) %>%
        summarise(
          Asym        = mean(Asym, na.rm=TRUE),
          R0          = mean(R0,   na.rm=TRUE),
          lrc         = mean(lrc,  na.rm=TRUE),
          Tau         = mean(Tau,  na.rm=TRUE),
          steadystate = mean(steadystate, na.rm=TRUE)
        )
      newSum <- data.frame(
        ID           = as.numeric(subjID()),
        Cycle_number = as.numeric(cycleNumber()),
        Asym         = coefs$Asym,
        R0           = coefs$R0,
        lrc          = coefs$lrc,
        Tau          = coefs$Tau,
        steadystate  = coefs$steadystate,
        Name         = paste0("id",subjID(),"c",cycleNumber()),
        stringsAsFactors = FALSE
      )
      globalSummaryData(rbind(globalSummaryData(), newSum))
    }
    
    showNotification("Done: cleaned & summary updated", type="message")
  })
  
  reactive({ session$userData$cleaned })
}

# --- Global Storage ---
globalSummaryData <- reactiveVal(data.frame(
  ID           = numeric(),
  Cycle_number = numeric(),
  Asym         = numeric(),
  R0           = numeric(),
  lrc          = numeric(),
  Tau          = numeric(),
  steadystate  = numeric(),
  Name         = character(),
  stringsAsFactors = FALSE
))
globalRemovedData <- reactiveVal(data.frame(
  ID           = numeric(),
  Cycle_number = numeric(),
  Removed      = I(list()),
  Name         = character(),
  stringsAsFactors = FALSE
))

# --- Patient Plot Helper (with τ on curves) ---
generate_patient_plot <- function(pid, summary_df) {
  patient_df <- summary_df %>% filter(ID == pid)
  if (nrow(patient_df) == 0) return(NULL)
  
  curves_list <- lapply(seq_len(nrow(patient_df)), function(i) {
    row   <- patient_df[i, ]
    x_vals <- seq(0, 120, length.out = 500)
    y_vals <- asym_regression_normalized(x_vals, row$Asym, row$R0, row$lrc)
    data.frame(x = x_vals, y = y_vals, Cycle = row$Cycle_number)
  })
  curve_data <- bind_rows(curves_list)
  
  tau_pts <- patient_df %>%
    transmute(
      Cycle = Cycle_number,
      x     = Tau,
      y     = asym_regression(Tau, Asym, R0, lrc) -
        asym_regression(0, Asym, R0, lrc)
    )
  
  
  curve_plot <- ggplot(curve_data,
                       aes(x = x, y = y, group = factor(Cycle), color = factor(Cycle))) +
    geom_line(size = 1) +
    geom_point(data = tau_pts,
               aes(x = x, y = y, color = factor(Cycle)),
               size = 3) +
    labs(title = paste("Recovery Curves for Patient", pid),
         x = "Recovery Time (sec)", y = "",
         color = "Cycle") +
    theme_minimal()
  
  tau_plot <- ggplot(patient_df,
                     aes(x = Cycle_number, y = Tau)) +
    geom_line(size = 1, color = "grey40") +
    geom_point(size = 3, aes(color = factor(Cycle_number))) +
    scale_color_discrete(guide = FALSE) +
    labs(title = paste("Tau Progression for Patient", pid),
         x = "Cycle Number", y = "Tau (sec)") +
    theme_minimal()
  
  curve_plot / tau_plot
}

# --- UI ---
ui <- fluidPage(
  titlePanel("NIRS Analysis & Cycle Processing"),
  tabsetPanel(
    # 1) Cycle Processing tab
    tabPanel("Cycle Processing",
             sidebarLayout(
               sidebarPanel(
                 fileInput("cycleFile",   "Upload Single Cycle File (TXT)", accept = ".txt"),
                 textInput("cycleSubjID",  "Subject ID:"),
                 textInput("cycleNumber",  "Cycle Number:"),
                 numericInput("startMin",  "Start Minute:",  4,  min = 0, max = 60),
                 numericInput("startSec",  "Start Second:", 30, min = 0, max = 59),
                 actionButton("applyCycleSettings", "Apply Cycle Settings"),
                 selectInput(
                   "dropWav",
                   "Drop Wavelength:",
                   choices  = NULL,
                   selected = "None"
                 )
               ),
               mainPanel(
                 h4("Raw Data Plot"),
                 plotOutput("rawCyclePlot", height = "300px"),
                 h4("Regression Plot"),
                 plotOutput("regCyclePlot", height = "300px"),
                 h4("Cycle Info"),
                 verbatimTextOutput("cycleInfo")
               )
             )
    ),
    
    # 2) Data Cleaning Module tab
    tabPanel("Data Cleaning Module",
             interact_rawdata_ui("cleanModule")
    ),
    
    # 3) Summary tab
    tabPanel("Summary",
             sidebarLayout(
               sidebarPanel(
                 downloadButton("downloadSummary",     "Export Summary CSV"),
                 actionButton( "showGraphs",           "Show Graphs"),
                 downloadButton("downloadPatientPNGs", "Download Patient Plots (PNG ZIP)")
               ),
               mainPanel(
                 DT::dataTableOutput("summaryTable"),
                 uiOutput("exportedGraphsUI")
               )
             )
    ),
    
    # 4) Removed Values tab
    tabPanel("Removed Values",
             sidebarLayout(
               sidebarPanel(
                 downloadButton("downloadRemoved", "Export Removed CSV")
               ),
               mainPanel(
                 DT::dataTableOutput("removedTable")
               )
             )
    )
  )
)

# --- Server ---
server <- function(input, output, session) {
  cycleRaw <- reactive({
    req(input$cycleFile)
    tryCatch(
      processCycleFile(input$cycleFile$datapath),
      error = function(e) {
        showModal(modalDialog(
          title = "File Error", e$message, easyClose=TRUE
        ))
        NULL
      }
    )
  })
  
  cycleProcessed <- eventReactive(input$applyCycleSettings, {
    tryCatch({
      validate(
        need(input$cycleFile,   "Please upload a cycle file first."),
        need(input$cycleSubjID, "Please enter a Subject ID."),
        need(input$cycleNumber,"Please enter a Cycle Number.")
      )
      rd <- cycleRaw()
      if (is.null(rd)) stop("Invalid data")
      subj_id(rd)    <- input$cycleSubjID
      meas_id(rd)    <- input$cycleNumber
      st <- input$startMin * 60 + input$startSec
      meas_start(rd) <- st
      meas_end(rd)   <- st + 16*60
      session$userData$cleaned <- NULL
      rd
    }, error = function(e) {
      showModal(modalDialog(
        title = "Settings Error", e$message, easyClose=TRUE
      ))
      NULL
    })
  })
  
  
  
  # 2a) Populate the dropdown choices once a file is processed:
  observeEvent(cycleProcessed(), {
    wavs <- cycleProcessed()$info$cols$col_name
    updateSelectInput(
      session, 
      "dropWav",
      choices  = c("None", wavs),
      selected = "None"
    )
  })
  
  # 2b) Build a reactive version of the cycle data that removes the chosen channel:
  filteredCycle <- reactive({
    rd  <- cycleProcessed()
    sel <- input$dropWav
    
    # only drop if sel exactly matches one of the col_names
    if (!is.null(sel) &&
        nzchar(sel) &&
        sel != "None" &&
        sel %in% rd$info$cols$col_name
    ) {
      rd <- remove_wavelength(rd, sel)
    }
    
    rd
  })
  
  # after you define filteredCycle <- reactive({ ... })
  cleanedData <- callModule(
    interact_rawdata_server, "cleanModule",
    rawdata     = filteredCycle,
    subjID      = reactive(input$cycleSubjID),
    cycleNumber = reactive(input$cycleNumber)
  )
  # --- automatic replots whenever filteredCycle() changes ---
  output$rawCyclePlot <- renderPlot({
    req(filteredCycle())
    plot_rawdata(filteredCycle())
  }, height = 300)
  
  output$regCyclePlot <- renderPlot({
    req(filteredCycle())
    plot_regressions(analysis(filteredCycle()))
  }, height = 300)
  
  output$cycleInfo <- renderPrint({
    rd <- filteredCycle()
    list(
      Info           = info(rd),
      Bad_Cols       = bads(rd),
      Subject        = subj_id(rd),
      Measurement_ID = meas_id(rd),
      Start          = meas_start(rd),
      End            = meas_end(rd)
    )
  })
  
  
  
  output$summaryTable <- DT::renderDataTable({
    df <- globalSummaryData()
    df$Delete <- sprintf('<button class="delete-btn" data-row="%d">✖</button>',
                         seq_len(nrow(df)))
    df
  }, escape=FALSE, selection='none', server=FALSE,
  callback=JS(
    "table.on('click', '.delete-btn', function() {",
    "  var row=$(this).data('row');",
    "  Shiny.setInputValue('delete_row', row, {priority:'event'});",
    "});"
  ), options=list(pageLength=10))
  
  observeEvent(input$delete_row, {
    idx <- as.numeric(input$delete_row); df <- globalSummaryData()
    if (!is.na(idx) && idx>=1 && idx<=nrow(df)) globalSummaryData(df[-idx,])
  })
  
  output$downloadSummary <- downloadHandler(
    filename = function() sprintf("cycle_summary_%s.csv",Sys.Date()),
    content = function(f) write.csv(globalSummaryData(),f,row.names=FALSE)
  )
  output$removedTable <- DT::renderDataTable(globalRemovedData(), options=list(pageLength=10))
  output$downloadRemoved <- downloadHandler(
    filename = function() sprintf("removed_values_%s.csv",Sys.Date()),
    content = function(f){
      df<-globalRemovedData()
      if(nrow(df)>0)df$Removed<-sapply(df$Removed,paste,collapse=", ")
      write.csv(df,f,row.names=FALSE)
    }
  )
  output$downloadPatientPNGs <- downloadHandler(
    filename = function() {
      paste0("patient_plots_", Sys.Date(), ".zip")
    },
    content = function(zipfile) {
      # temp dir for the PNGs
      tmp <- tempdir()
      png_files <- vector("character", length = 0)
      
      # loop over each patient ID
      ids <- unique(globalSummaryData()$ID)
      for (id in ids) {
        png_path <- file.path(tmp, paste0("patient_", id, ".png"))
        png(png_path, width = 800, height = 600)
        print(generate_patient_plot(id, globalSummaryData()))
        dev.off()
        png_files <- c(png_files, png_path)
      }
      
      # zip them up
      # on most systems utils::zip will work
      old_wd <- setwd(tmp)
      utils::zip(zipfile, basename(png_files))
      setwd(old_wd)
    },
    contentType = "application/zip"
  )
  
  
  
  
  observeEvent(input$showGraphs, {
    df <- globalSummaryData(); req(nrow(df)>0)
    ids <- unique(df$ID)
    output$exportedGraphsUI <- renderUI({
      tagList(lapply(ids,function(id) plotOutput(paste0("patientPlot_",id),height="400px")))
    })
    for(id in ids) local({
      this_id <- id
      output[[paste0("patientPlot_",this_id)]] <- renderPlot({
        generate_patient_plot(this_id, globalSummaryData())
      })
    })
  })
}

# Run the app
shinyApp(ui, server)
