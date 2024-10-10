library(shiny)
library(rhandsontable)
library(RSQLite)
library(DBI)

ui <- bootstrapPage(
  titlePanel("Sample Sheet Generating"),
  fluidRow(
    column(3, textInput("flowcell", "Flowcell:", value = "Input Flowcell")),
    column(3, dateInput("date", "Date:", value = Sys.Date()))
  ),
  fluidRow(
    column(3, selectInput("numRows", "Number of Rows:", choices = c(24, 48, 96), selected = 48)),
    column(3, selectInput("fileType", "Choose file type:", choices = c("TSV" = "tsv", "CSV" = "csv")))
  ),
  downloadButton("downloadData", "Export Samplesheet"),
  actionButton("saveToDB", "Save"),
  tags$h1(),
  rHandsontableOutput("samplesheet"),
  tags$h1()
)

server <- function(input, output, session) {
  predefined_values <- data.frame(
    I7_Index_ID = c("N701", "N702", "N703", "N704", "N705", "N706", "N707", "N708", "N709", "N710", "N711", "N712",
                    "N701", "N702", "N703", "N704", "N705", "N706", "N707", "N708", "N709", "N710", "N711", "N712"),
    index = c("TAACTCTT", "CGTAGTGA", "AGACGACG", "TCCTGAGA", "TCTCAGTC", "TAGGAGAC", "TAGGACAG", "CTCTACAC",
              "TACACTGA", "CAGTACGA", "CGTGCTAG", "GCTAGGGA", "TAACTCTT", "CGTAGTGA", "AGACGACG", "TCCTGAGA",
              "TCTCAGTC", "TAGGAGAC", "TAGGACAG", "CTCTACAC", "TACACTGA", "CAGTACGA", "CGTGCTAG", "GCTAGGGA"),
    I5_Index_ID = rep(c("S503", "S504"), each = 12),
    index2 = rep(c("TATCCTCT", "AGAGTAGA"), each = 12),
    stringsAsFactors = FALSE
  )
  # Initialize empty data frame
  generate_samplesheet <- function(numRows, flowcell, date) {
    embryo_id_placeholder <- ""  # Placeholder for Embryo_ID
    sample_id <- paste(flowcell, "-", format(date, "%Y%m%d"), "-", embryo_id_placeholder, sep = "")
    if (numRows == 24) {
      data <- predefined_values
    } else {
      data <- data.frame(
        I7_Index_ID = rep("", numRows),
        index = rep("", numRows),
        I5_Index_ID = rep("", numRows),
        index2 = rep("", numRows),
        stringsAsFactors = FALSE
      )
    }
    data.frame(
      Sample_ID = rep(sample_id, numRows),  
      Cycle_ID = rep(format(date, "%Y-%m-%d"), numRows),  # Use current date as default
      Embryo_ID = rep("", numRows),
      Cell_Type = rep("", numRows),
      DNA_ng_ul = rep("", numRows),
      Sample_Plate = rep("", numRows),
      Sample_Well = rep("", numRows),
      I7_Index_ID = data$I7_Index_ID,
      index = data$index,
      I5_Index_ID = data$I5_Index_ID,
      index2 = data$index2,
      stringsAsFactors = FALSE
    )
  }
  
  data <- reactiveVal()
  
  observe({
    req(input$flowcell, input$date, input$numRows)  # Ensure inputs are provided
    data(generate_samplesheet(input$numRows, input$flowcell, input$date))
  })
  
  # Render editable table
  output$samplesheet <- renderRHandsontable({
    rhandsontable(data(), useTypes = TRUE)
  })
  
  # Update data upon editing
  observeEvent(input$samplesheet, {
    if (!is.null(input$samplesheet)) {
      updated_data <- hot_to_r(input$samplesheet)
      
      # Update Sample_ID dynamically based on Embryo_ID
      if (!identical(updated_data$Embryo_ID, data()$Embryo_ID)) {
        updated_data$Sample_ID <- paste(input$flowcell, "-", format(input$date, "%Y%m%d"), "-", updated_data$Embryo_ID, sep = "")
      }
      
      data(updated_data)
    }
  })
  
  generate_header <- function(date, flowcell) {
    header <- data.frame(
      Key = c(
        "[Header]",
        "SampleSheetGenerating",
        "Experiment Name",
        "Date",
        "Module",
        "Workflow",
        "Library Prep Kit",
        "Index Kit",
        "Chemistry",
        "",
        "[Reads]",
        "36",
        "",
        "[Settings]",
        "FlagPCRDuplicates",
        "CustomRead2PrimerMix",
        "",
        "[Data]"
      ),
      Value = c(
        "",
        "anegin",
        flowcell,
        format(date, "%Y-%m-%d"),
        "GenerateFASTQ - 3.0.1",
        "GenerateFASTQ",
        "Custom",
        "Custom",
        "Amplicon",
        "",
        "",
        "",
        "",
        "",
        "1",
        "c3",
        "",
        ""
      ),
      stringsAsFactors = FALSE
    )
    return(header)
  }
  
  # Download data as TSV
  output$downloadData <- downloadHandler(
    filename = function() {
      paste(format(input$date, "%Y%m%d"), "_", input$flowcell, "_samplesheet.", input$fileType, sep = "")
    },
    content = function(file) {
      sample_data <- data()
      header <- generate_header(input$date, input$flowcell)
      if (input$fileType == "tsv") {
        write.table(header, file, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
        write.table(sample_data, file, sep = "\t", row.names = FALSE, quote = FALSE, append = TRUE)
      } else if (input$fileType == "csv") {
        write.table(header, file, sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)
        write.table(sample_data, file, sep = ",", row.names = FALSE, quote = FALSE, append = TRUE)
      }
    }
  )
  
  observeEvent(input$saveToDB, {
    db <- dbConnect(SQLite(), dbname = "/home/anegin97/Bioinformatics/appforclient/CNVdat/cnv_data.sqlite")
    
    table_name <- paste(format(input$date, "%Y%m%d"), "_", input$flowcell, "_samplesheet", sep = "")
    
    dbWriteTable(db, table_name, data(), row.names = FALSE, overwrite = TRUE)
    
    dbDisconnect(db)
    
    showModal(modalDialog(
      title = "Database Update",
      paste("Samplesheet saved to database with table name:", table_name),
      easyClose = TRUE,
      footer = NULL
    ))
  })
}

shinyApp(ui, server)
