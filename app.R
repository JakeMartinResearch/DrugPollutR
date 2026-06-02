library(shiny)
library(bslib)
library(shinyWidgets)
library(dplyr)
library(ggplot2)
library(data.table)
library(DT)
library(scales)
library(patchwork)

# Source dose function
source("DoseSelectR.R")

# Load pre-processed data
env_data <- readRDS("env_data.rds")
compound_key <- readRDS("compound_key.rds")

setDT(env_data)
setDT(compound_key)

# Ensure key columns are character
char_cols <- c("compound_name", "compound_cas", "compound_cid", "parent_compound")

for (col in intersect(char_cols, names(env_data))) {
  set(env_data, j = col, value = as.character(env_data[[col]]))
}

for (col in intersect(char_cols, names(compound_key))) {
  set(compound_key, j = col, value = as.character(compound_key[[col]]))
}

# Ensure concentration is numeric and keep only positive detections
env_data[, value := suppressWarnings(as.numeric(value))]
env_data <- env_data[!is.na(value) & value > 0]

# -----------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------
safe_name <- function(x) {
  x <- paste(x, collapse = "_")
  x <- trimws(x)
  x <- gsub("[^A-Za-z0-9_\\-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  if (!nzchar(x)) x <- "chemical"
  x
}

make_choice_vector <- function(search_type, key_dt) {
  vals <- unique(key_dt[[search_type]])
  vals <- vals[!is.na(vals) & nzchar(vals)]
  sort(vals)
}

dose_input_status <- function(dat, max_doses, spacing) {
  messages <- character()
  level <- "success"
  
  n_values <- nrow(dat)
  n_unique <- uniqueN(dat$value)
  
  if (n_values < 3) {
    level <- "danger"
    messages <- c(
      messages,
      paste0(
        "Dose selection needs at least 3 positive concentration values; this result has ",
        n_values,
        "."
      )
    )
  } else if (n_values <= 10) {
    level <- "warning"
    messages <- c(
      messages,
      paste0(
        "Only ",
        n_values,
        " positive concentration values are available, so dose selection may be unstable."
      )
    )
  }
  
  if (n_unique < 2) {
    level <- "danger"
    messages <- c(messages, "Dose selection needs at least 2 unique positive values to build the density plot.")
  }
  
  if (is.na(spacing) || spacing <= 1) {
    level <- "danger"
    messages <- c(messages, "Dose spacing factor must be greater than 1.")
  }
  
  if (is.na(max_doses) || max_doses < 3) {
    level <- "danger"
    messages <- c(messages, "Number of doses must be at least 3.")
  }
  
  list(level = level, messages = messages)
}

export_citation_note <- paste(
  "This database was compiled for Martin et al (2025) Environ. Sci. Technol. Lett. 2025, 12, 10, 1308-1313",
  "(https://doi.org/10.1021/acs.estlett.5c00665). It is based on a filtered synthesis of three publicly available",
  "datasets: (1) the NORMAN EMPODAT database for chemical occurrence (accessed 18/03/2025), (2) the",
  "Umweltbundesamt Pharmaceuticals in the Environment database (PHARMS-UBA; accessed 19/12/2024), and (3)",
  "Wilkinson et al. (2022) Pharmaceutical Pollution of the World's Rivers database. Data were restricted to",
  "entries reported in mass per volume of water (e.g., \u00b5g/L) and relevant to surface water and wastewater",
  "matrices (for details of the filtering process, please refer to Martin et al. (2025)). Users of this data",
  "should cite the original sources of the data: NORMAN EMPODAT, PHARMS-UBA, and Wilkinson et al. (2022),",
  "in conjunction with Martin et al (2025) https://doi.org/10.1021/acs.estlett.5c00665"
)

write_export_with_citation <- function(dat, file) {
  dat <- as.data.table(dat)
  n_cols <- ncol(dat)
  note_row <- as.list(c(export_citation_note, rep("", max(0, n_cols - 1L))))
  blank_row <- as.list(rep("", n_cols))
  names(note_row) <- names(dat)
  names(blank_row) <- names(dat)
  
  fwrite(as.data.table(note_row), file, col.names = FALSE)
  fwrite(as.data.table(blank_row), file, col.names = FALSE, append = TRUE)
  fwrite(dat, file, append = TRUE, col.names = TRUE)
}

# -----------------------------------------------------------------------
# UI
# -----------------------------------------------------------------------
ui <- navbarPage(
  title = "PollutionScopeR",
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  
  tabPanel(
    "ChemicalExploreR",
    sidebarLayout(
      sidebarPanel(
        h4("Find compounds"),
        
        radioButtons(
          "search_type",
          "Search by",
          choices = c(
            "Name" = "compound_name",
            "CAS"  = "compound_cas",
            "CID"  = "compound_cid"
          ),
          selected = "compound_name"
        ),
        
        uiOutput("compound_search_ui"),
        
        checkboxInput(
          "group_parents",
          label = tags$span(
            "Group by parent compound",
            tags$small(
              class = "text-muted d-block",
              "Combines related forms under parent_compound"
            )
          ),
          value = TRUE
        ),
        
        hr(),
        
        selectInput(
          "matrix",
          "Matrix type",
          choices = c("Combined", "effluent", "surfacewater"),
          selected = "Combined"
        ),
        
        actionButton(
          "search",
          "Search",
          class = "btn-primary w-100 mt-2"
        ),
        
        hr(),
        
        uiOutput("result_picker_ui")
      ),
      
      mainPanel(
        h5("Note: statistics are calculated from positive environmental detections only."),
        br(),
        
        h3("Matched compounds"),
        DTOutput("matched_compounds"),
        br(),
        
        h3("Summary statistics (µg/L)"),
        DTOutput("summary_table"),
        br(),
        div(
          class = "d-flex gap-2 flex-wrap",
          downloadButton("download_summary", "Download summary statistics"),
          downloadButton("download_all_summary", "Download all search summaries")
        ),
        br(), br(),
        
        h3("Concentration distribution"),
        uiOutput("density_warning_ui"),
        plotOutput("density_plot", height = "500px"),
        br(),
        
        h3("Filtered data"),
        DTOutput("data_table"),
        br(),
        div(
          class = "d-flex gap-2 flex-wrap",
          downloadButton("download_filtered", "Download filtered data"),
          downloadButton("download_all_filtered", "Download all filtered search data")
        )
      )
    )
  ),
  
  tabPanel(
    "DoseSelectR",
    sidebarLayout(
      sidebarPanel(
        selectInput(
          "central",
          "Central tendency",
          choices = c("median", "mean", "mode")
        ),
        numericInput(
          "max_doses",
          "Number of doses",
          5,
          min = 3,
          max = 10
        ),
        numericInput(
          "spacing",
          "Dose spacing factor",
          3.2
        ),
        hr(),
        
        uiOutput("dose_result_picker_ui"),
        uiOutput("dose_warning_ui"),
        
        hr(),
        
        actionButton("run_doses", "Generate doses")
      ),
      mainPanel(
        h3("Recommended doses (µg/L)"),
        DTOutput("dose_table"),
        br(),
        downloadButton("download_doses", "Download dose table"),
        br(), br(),
        plotOutput("dose_plot", height = "600px")
      )
    )
  )
)

# -----------------------------------------------------------------------
# Server
# -----------------------------------------------------------------------
server <- function(input, output, session) {
  
  # Dynamic autocomplete dropdown based on selected search type
  output$compound_search_ui <- renderUI({
    choices <- make_choice_vector(input$search_type, compound_key)
    
    selectizeInput(
      "search_values",
      "Search terms",
      choices = choices,
      selected = NULL,
      multiple = TRUE,
      options = list(
        placeholder = "Start typing to search...",
        maxOptions = 2000,
        plugins = list("remove_button")
      )
    )
  })
  
  # Matched compounds from compound_key
  matched_key <- eventReactive(input$search, {
    req(input$search_type)
    req(length(input$search_values) > 0)
    
    search_col <- input$search_type
    vals <- as.character(input$search_values)
    
    compound_key[get(search_col) %in% vals]
  })
  
  output$matched_compounds <- renderDT({
    dat <- matched_key()
    req(nrow(dat) > 0)
    
    show_cols <- c("compound_name", "compound_cas", "compound_cid", "parent_compound")
    show_cols <- show_cols[show_cols %in% names(dat)]
    
    datatable(
      dat[, ..show_cols],
      options = list(pageLength = 10),
      rownames = FALSE
    )
  })
  
  # Toggle through selected compounds one at a time
  output$result_picker_ui <- renderUI({
    dat <- matched_key()
    req(nrow(dat) > 0)
    
    if (isTRUE(input$group_parents)) {
      vals <- unique(dat$parent_compound)
      vals <- vals[!is.na(vals) & nzchar(vals)]
      vals <- sort(vals)
      selectInput("active_result", "View result", choices = vals, selected = vals[1])
    } else {
      vals <- unique(dat$compound_name)
      vals <- vals[!is.na(vals) & nzchar(vals)]
      vals <- sort(vals)
      selectInput("active_result", "View result", choices = vals, selected = vals[1])
    }
  })
  
  output$dose_result_picker_ui <- renderUI({
    dat <- matched_key()
    req(nrow(dat) > 0)
    
    if (isTRUE(input$group_parents)) {
      vals <- unique(dat$parent_compound)
    } else {
      vals <- unique(dat$compound_name)
    }
    
    vals <- vals[!is.na(vals) & nzchar(vals)]
    vals <- sort(vals)
    req(length(vals) > 0)
    
    selected <- if (isTruthy(input$active_result) && input$active_result %in% vals) {
      input$active_result
    } else {
      vals[1]
    }
    
    selectInput("dose_active_result", "View result", choices = vals, selected = selected)
  })
  
  # Filter env_data for selected result
  filtered_data <- reactive({
    req(input$active_result)
    
    dat <- copy(env_data)
    
    if (isTRUE(input$group_parents)) {
      dat <- dat[parent_compound == input$active_result]
    } else {
      dat <- dat[compound_name == input$active_result]
    }
    
    if (input$matrix != "Combined") {
      dat <- dat[matrix_group == input$matrix]
    }
    
    dat <- dat[!is.na(value) & value > 0]
    dat
  })
  
  dose_filtered_data <- reactive({
    req(input$dose_active_result)
    
    dat <- copy(env_data)
    
    if (isTRUE(input$group_parents)) {
      dat <- dat[parent_compound == input$dose_active_result]
    } else {
      dat <- dat[compound_name == input$dose_active_result]
    }
    
    if (input$matrix != "Combined") {
      dat <- dat[matrix_group == input$matrix]
    }
    
    dat <- dat[!is.na(value) & value > 0]
    dat
  })
  
  # Filter env_data for all matched results
  all_filtered_data <- reactive({
    dat <- copy(env_data)
    key_matches <- matched_key()
    req(nrow(key_matches) > 0)
    
    if (isTRUE(input$group_parents)) {
      keep_vals <- unique(key_matches$parent_compound)
      keep_vals <- keep_vals[!is.na(keep_vals) & nzchar(keep_vals)]
      dat <- dat[parent_compound %in% keep_vals]
    } else {
      keep_vals <- unique(key_matches$compound_name)
      keep_vals <- keep_vals[!is.na(keep_vals) & nzchar(keep_vals)]
      dat <- dat[compound_name %in% keep_vals]
    }
    
    if (input$matrix != "Combined") {
      dat <- dat[matrix_group == input$matrix]
    }
    
    dat <- dat[!is.na(value) & value > 0]
    dat
  })
  
  # Summary statistics for selected result
  summary_stats <- reactive({
    dat <- filtered_data()
    req(nrow(dat) > 0)
    
    data.table(
      N      = nrow(dat),
      Mean   = round(mean(dat$value, na.rm = TRUE), 3),
      Median = round(median(dat$value, na.rm = TRUE), 3),
      Min    = round(min(dat$value, na.rm = TRUE), 3),
      Max    = round(max(dat$value, na.rm = TRUE), 3),
      SD     = round(sd(dat$value, na.rm = TRUE), 3),
      Q05    = round(quantile(dat$value, 0.05, na.rm = TRUE), 3),
      Q95    = round(quantile(dat$value, 0.95, na.rm = TRUE), 3)
    )
  })
  
  # Summary statistics for all matched results
  all_summary_stats <- reactive({
    dat <- all_filtered_data()
    req(nrow(dat) > 0)
    
    group_var <- if (isTRUE(input$group_parents)) "parent_compound" else "compound_name"
    
    dat[, .(
      compound_name = paste(sort(unique(compound_name)), collapse = "; "),
      compound_cas = paste(sort(unique(compound_cas)), collapse = "; "),
      parent_compound = paste(sort(unique(parent_compound)), collapse = "; "),
      N = .N,
      Mean = round(mean(value, na.rm = TRUE), 3),
      Median = round(median(value, na.rm = TRUE), 3),
      Min = round(min(value, na.rm = TRUE), 3),
      Max = round(max(value, na.rm = TRUE), 3),
      SD = round(sd(value, na.rm = TRUE), 3),
      Q05 = round(quantile(value, 0.05, na.rm = TRUE), 3),
      Q95 = round(quantile(value, 0.95, na.rm = TRUE), 3)
    ), by = group_var]
  })
  
  output$summary_table <- renderDT({
    stats <- summary_stats()
    datatable(stats, options = list(dom = "t"), rownames = FALSE)
  })
  
  output$density_warning_ui <- renderUI({
    dat <- tryCatch(
      filtered_data(),
      shiny.silent.error = function(e) NULL
    )
    
    if (is.null(dat)) {
      return(NULL)
    }
    
    if (nrow(dat) <= 10) {
      return(
        div(
          class = "alert alert-warning",
          paste0(
            "Concentration distribution needs more than 10 positive detections; this result has ",
            nrow(dat),
            "."
          )
        )
      )
    }
    
    NULL
  })
  
  # Density plot
  output$density_plot <- renderPlot({
    dat <- filtered_data()
    validate(
      need(
        nrow(dat) > 10,
        paste0(
          "Concentration distribution needs more than 10 positive detections; this result has ",
          nrow(dat),
          "."
        )
      )
    )
    
    ggplot(dat, aes(x = value)) +
      geom_density(
        fill   = "#2C7FB8",
        alpha  = 0.5,
        adjust = 1,
        n      = 512
      ) +
      scale_x_log10(
        breaks = scales::log_breaks(n = 6),
        labels = scales::label_number(accuracy = 0.01)
      ) +
      theme_minimal(base_size = 15) +
      labs(
        x = "Environmental concentration (µg/L, log10 scale)",
        y = "Density"
      )
  })
  
  # Data table
  output$data_table <- renderDT({
    dat <- filtered_data()
    req(nrow(dat) > 0)
    
    display_cols <- c(
      "compound_name",
      "compound_cas",
      "compound_cid",
      "parent_compound",
      "reference",
      "matrix_group",
      "value"
    )
    display_cols <- display_cols[display_cols %in% names(dat)]
    
    datatable(
      dat[, ..display_cols],
      options = list(pageLength = 10),
      rownames = FALSE
    )
  })
  
  output$dose_warning_ui <- renderUI({
    dat <- tryCatch(
      dose_filtered_data(),
      shiny.silent.error = function(e) NULL
    )
    
    if (is.null(dat)) {
      return(
        div(
          class = "alert alert-info mt-3",
          "Search for compounds in Chemical Explorer, then choose a result here."
        )
      )
    }
    
    status <- dose_input_status(dat, input$max_doses, input$spacing)
    
    if (length(status$messages) == 0) {
      return(
        div(
          class = "alert alert-success mt-3",
          paste0(nrow(dat), " positive concentration values available for dose selection.")
        )
      )
    }
    
    div(
      class = paste("alert mt-3", paste0("alert-", status$level)),
      tags$ul(lapply(status$messages, tags$li))
    )
  })
  
  # Dose selection
  dose_results <- eventReactive(input$run_doses, {
    dat <- dose_filtered_data()
    status <- dose_input_status(dat, input$max_doses, input$spacing)
    
    if (status$level == "danger") {
      return(list(
        ok = FALSE,
        message = paste(status$messages, collapse = " "),
        result = NULL,
        warnings = character()
      ))
    }
    
    loq_val <- max(0.1, min(dat$value, na.rm = TRUE) * 0.01)
    
    warnings <- character()
    
    res <- tryCatch(
      withCallingHandlers(
        select_doses(
          data             = dat,
          variable         = "value",
          central_tendency = input$central,
          min_spacing      = input$spacing,
          max_doses        = input$max_doses,
          LOQ              = loq_val
        ),
        warning = function(w) {
          warnings <<- c(warnings, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      ),
      error = function(e) e
    )
    
    if (inherits(res, "error")) {
      return(list(
        ok = FALSE,
        message = conditionMessage(res),
        result = NULL,
        warnings = warnings
      ))
    }
    
    if (length(res) == 1 && is.na(res)) {
      return(list(
        ok = FALSE,
        message = paste(c(warnings, "Dose selection did not return a dose table."), collapse = " "),
        result = NULL,
        warnings = warnings
      ))
    }
    
    res$doses$value <- round(res$doses$value, 3)
    
    list(
      ok = TRUE,
      message = NULL,
      result = res,
      warnings = warnings
    )
  })
  
  output$dose_table <- renderDT({
    res <- dose_results()
    req(res)
    validate(need(isTRUE(res$ok), res$message))
    datatable(res$result$doses, rownames = FALSE)
  })
  
  output$dose_plot <- renderPlot({
    res <- dose_results()
    req(res)
    validate(need(isTRUE(res$ok), res$message))
    print(res$result$plot)
  })
  
  # Downloads - current selected result
  output$download_filtered <- downloadHandler(
    filename = function() {
      paste0(
        "filtered_data_",
        safe_name(input$active_result),
        "_",
        Sys.Date(),
        ".csv"
      )
    },
    content = function(file) {
      dat <- filtered_data()
      req(nrow(dat) > 0)
      write_export_with_citation(dat, file)
    }
  )
  
  output$download_summary <- downloadHandler(
    filename = function() {
      paste0(
        "summary_statistics_",
        safe_name(input$active_result),
        "_",
        Sys.Date(),
        ".csv"
      )
    },
    content = function(file) {
      stats <- summary_stats()
      req(nrow(stats) > 0)
      write_export_with_citation(stats, file)
    }
  )
  
  output$download_doses <- downloadHandler(
    filename = function() {
      paste0(
        "dose_table_",
        safe_name(input$dose_active_result),
        "_",
        Sys.Date(),
        ".csv"
      )
    },
    content = function(file) {
      res <- dose_results()
      req(res)
      validate(need(isTRUE(res$ok), res$message))
      write_export_with_citation(as.data.table(res$result$doses), file)
    }
  )
  
  # Downloads - all matched results
  output$download_all_summary <- downloadHandler(
    filename = function() {
      paste0(
        "all_search_summaries_",
        safe_name(input$search_type),
        "_",
        Sys.Date(),
        ".csv"
      )
    },
    content = function(file) {
      dat <- all_summary_stats()
      req(nrow(dat) > 0)
      write_export_with_citation(dat, file)
    }
  )
  
  output$download_all_filtered <- downloadHandler(
    filename = function() {
      paste0(
        "all_filtered_search_data_",
        safe_name(input$search_type),
        "_",
        Sys.Date(),
        ".csv"
      )
    },
    content = function(file) {
      dat <- all_filtered_data()
      req(nrow(dat) > 0)
      
      export_cols <- c(
        "compound_name",
        "compound_cas",
        "compound_cid",
        "parent_compound",
        "reference",
        "matrix_group",
        "value"
      )
      export_cols <- export_cols[export_cols %in% names(dat)]
      
      write_export_with_citation(dat[, ..export_cols], file)
    }
  )
}

shinyApp(ui, server)
