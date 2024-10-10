library(readxl)
library(dplyr)
library(ggplot2)
library(DT)
library(ggtext)
library(tidyverse)
library(plotly)
library(smoother)
library(shiny)
library(shinyFiles)
library(fs)
library(shinyjs)
library(magick)

# Initialize shinyjs
shinyjs::useShinyjs()

selectizeSelectInput <- function(inputId, label, choices, selected = NULL) {
  selectizeInput(inputId, label, choices = choices, selected = selected, options = list(search = TRUE))
}

ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      body {
        background-color: #222;
        margin-top:20px; 
      }
      .dataTables_wrapper {
        background-color: #fff;
        color: #000;
        margin-top: 20px;
      }
      #table-container {
        float: right;
        margin-left: 20px;
      }
      .file-item {
        cursor: pointer;
        margin-bottom: 5px;
        padding: 5px;
        background-color: #337ab7;
        color: #fff;
        border-radius: 3px;
      }
      .file-item:hover {
        background-color: #286090;
      }
      .file-list-container {
        height: 300px; /* Adjust height as needed */
        overflow-y: auto;
      }
      .filter-box {
        margin-bottom: 10px;
      }
    "))
  ),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      style = "height: 1000px; overflow-y: auto;",
      tags$h3("CNVviso 1.0", style="color: #000; font-weight: bold;"),
      fileInput("file", "Import Data (.zip)", multiple = TRUE, accept = c(".zip")),
      textInput("fileFilter", "Sample ID"),  # Text input for filtering files
      uiOutput("fileListContainer"),  # Container for the file list
      verbatimTextOutput("directorypath"),
      tags$hr(),
      div(class = "view-box",
          tags$h5("View", style = "color: #000; font-weight: bold;"),
          checkboxInput("displayTable", "Display Data", value = FALSE),
          checkboxInput("showIdiogram", "Show Idiogram", value = FALSE)
      ),
      sliderInput("xSlider", "Chromosome size", min = 0, max = 66, value = c(0, 66))
    ),
    mainPanel(
      width = 9,
      plotlyOutput("myPlot"),
      column(7, imageOutput("idiogram"), style = "margin-top: 20px;"),
      column(5, DTOutput("table", width = "80%", height = "850px"))
    )
  )
)

server <- function(input, output, session) {
  
  observeEvent(input$file, {
    # Show loading status when a file is chosen
    shinyjs::show("loading_status")
    
    files <- input$file
    choices <- basename(files$name)
    updateSelectInput(session, "selected_file", choices = choices)
    # Dynamically generate list of files based on the filter
    output$fileListContainer <- renderUI({
      filtered_files <- if (input$fileFilter == "") {
        files$name
      } else {
        files$name[grep(input$fileFilter, files$name, ignore.case = TRUE)]
      }
      fileList <- lapply(seq_along(filtered_files), function(i) {
        tags$div(class = "file-item", onclick = paste0("Shiny.setInputValue('selected_file', '", filtered_files[i], "')"), filtered_files[i])
      })
      tags$div(class = "file-list-container", fileList)
    })
  })
  
  observeEvent(input$selected_file, {
    req(input$selected_file)
    # Extract the name of the uploaded ZIP file
    zip_file_name <- gsub(".zip", "", input$selected_file)
    # Ensure that the extraction directory exists based on the name of the ZIP file
    extraction_directory <- paste0("unzipped_files_", zip_file_name)
    if (!dir.exists(extraction_directory)) {
      dir.create(extraction_directory)
    }
    unzip(input$file$datapath[which(input$file$name == input$selected_file)], exdir = extraction_directory)
    output$directorypath <- renderPrint({
      "Process complete"
    })
    # Hide loading status when processing is complete
    shinyjs::hide("loading_status")
  })
  
  # Define reactive dataframes for CNV files
  extracted_files <- reactive({
    req(input$selected_file)
    # Extract the name of the uploaded ZIP file
    zip_file_name <- gsub(".zip", "", input$selected_file)
    list.files(paste0("unzipped_files_", zip_file_name), full.names = TRUE)
  })
  
  # Define reactive expressions to read CNR, CNS, and CALL files
  dfcnr <- reactive({
    req(extracted_files())
    file_path_cnr <- extracted_files() %>% keep(., grepl(".dedup.cnr", .)) %>% first()
    read_tsv(file_path_cnr)
  })
  
  dfcns <- reactive({
    req(extracted_files())
    file_path_cns <- extracted_files() %>% keep(., grepl(".dedup.cns", .)) %>% first()
    read_tsv(file_path_cns)
  })
  
  dfcall <- reactive({
    req(extracted_files())
    file_path_call <- extracted_files() %>% keep(., grepl(".dedup.call.cns", .)) %>% first()
    read_tsv(file_path_call)
  })
  
  idiogram <- reactive({
    req(extracted_files())
    file_path_idiogram <- extracted_files() %>% keep(.,grepl(".png",.)) %>% first()
    idiogram <- file_path_idiogram
  })
  
  output$myPlot <- renderPlotly({
    req(dfcnr(),dfcns(),dfcall())
    dfcnr<-dfcnr()
    dfcns<-dfcns()
    dfcall<-dfcall()
    dfcnr <- dfcnr %>%
      mutate(lrr = 2 * 2 ^ log2)
    for (i in 1:length(dfcnr$chromosome)) {
      dfcnr$order[i] <- i
    }
    
    dfcns <- dfcns %>%
      mutate(cncall = 2 * 2 ^ log2) %>%
      select(chromosome, start, end, cncall,weight) %>%
      filter(chromosome != "chrY")%>%
      filter(weight >= 15)
    dfcall <- dfcall %>% filter(chromosome == "chrY") %>% 
      select(chromosome,start,end,cn)
    names(dfcall)[names(dfcall) == "cn"] <- "cncall"
    dfcns<-merge(dfcall,dfcns,all=TRUE)
    df <- merge(dfcnr, dfcns, by = c("chromosome", "end"), all = TRUE) %>%
      select(-start.y)
    names(df)[names(df) == "start.x"] <- "start"
    df <- merge(df, dfcns, by.x = c("chromosome", "start"),
                by.y = c("chromosome", "start"), all = TRUE) %>%
      select(-end.y)
    names(df)[names(df) == "end.x"] <- "end"
    df[is.na(df)] <- 0
    df$cncall <- df$cncall.x + df$cncall.y
    df <- df %>% arrange(order)
    df <- df %>% select(chromosome, start, end, lrr, cncall)
    df <- df %>% mutate(lrr=ifelse(chromosome=="chrY",dfcall$cncall,lrr))
    
    CNvector <- dfcns$cncall
    CNcall <- as.numeric(df$cncall)
    
    start_index <- vector()
    end_index <- vector()
    
    for (i in 1:length(CNvector)) {
      start_index[i] <- which(CNcall == CNvector[i])[1]
    }
    start_index<-na.omit(start_index)
    for (i in 1:length(CNvector)) {
      end_index[i] <- which(CNcall == CNvector[i])[2]
    }
    end_index<-na.omit(end_index)
    for (i in 1:length(start_index)) {
      df$cncall[(start_index[i] + 1):(end_index[i] - 1)] <- df$cncall[start_index[i]]
    }
    df <- df%>%
          mutate(cncall=ifelse(chromosome!=c("chrX","chrY") & cncall==0, NA, cncall))%>%
          na.omit()
    hg38size<-data.frame(
      chromosome=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8",
                   "chr9","chr10","chr11","chr12","chr13","chr14","chr15",
                   "chr16","chr17","chr18","chr19","chr20","chr21","chr22",
                   "chrX","chrY"),
      length=c(248956422,242193529,198295559,190214555,181538259,170805979,
               159345973,145138636,138394717,135086622,133797422,133275309,
               114364328,107043718,101991189,90338345,83257441,80373285,
               58617616,64444167,46709983,50818468,156040895,57227415)
    )
    striphg38 <- hg38size %>%
      rename_all(tolower) %>%
      mutate(pos = cumsum(length / min(length))) %>%
      mutate(ymin = 0, ymax = 4,
             xmin = c(0, pos[1:23]), xmax = pos,
             fill = rep(c("xam", "trang"), length.out = nrow(hg38size))) %>%
      pivot_longer(c(ymin, ymax), values_to = "y", names_to = "yy") %>% select(-yy)
    sizepros <- hg38size %>%
      rename_all(tolower) %>%
      mutate(propos = length / min(length)) %>% select(-length) %>%
      mutate(pos = cumsum(propos)) %>%
      mutate(c = c(0, pos[1:23])) %>%
      mutate(poschr = propos / 2 + c)
    
    data <- df %>%
      inner_join(., sizepros, by = "chromosome") %>%
      group_by(chromosome) %>%
      mutate(position = cumsum(propos / length(chromosome))) %>%
      mutate(point_x = c + position) %>%
      select(chromosome, start, end, lrr, cncall, point_x)
    
    
    data <- data %>%
      mutate(lrr = ifelse(lrr > cncall + 0.1, NA, lrr)) %>% na.omit()
    data <- data %>%
      mutate(lrr = ifelse(lrr < cncall - 0.1, NA, lrr)) %>% na.omit()
    cncallcount<-table(data$cncall)
    singletons<-names(cncallcount[cncallcount==1])
    data<-data[!data$cncall %in% singletons, ]
    data <- data %>%
      group_by(chromosome,cncall) %>%
      mutate(lrr = smth.gaussian(lrr)) %>% na.omit()
    
    p<-data %>%
      ggplot(aes(x = point_x, y = lrr))+
      geom_vline(size = 0.3, xintercept = c(unique(sizepros$pos)), linetype = 2, color = "grey90") +
      geom_hline(size=0.3,yintercept=4,linetype="solid")+
      geom_vline(size=0.3,xintercept=66,linetype="solid")+
      geom_line(linewidth = 0.3, color = "green",show.legend = FALSE) +
      geom_point(aes(text=paste("Chromosome:",chromosome,"<br>Start:",start,"<br>End:",end,"<br>CN:",lrr)),size = 1.3, shape = 21, color = alpha("black",0.6), fill = "chartreuse3") +
      geom_hline(size = 0.3, yintercept = c(1, 3)) +
      geom_hline(size = 0.5, yintercept = 2, color = "green", linetype = "solid") +
      geom_vline(size = 0.3, xintercept = 61.550044) +
      scale_fill_manual(name = NULL,
                        breaks = c("xam", "trang"),
                        values = c("#F3F3F3", "#F3FAFE")) +
      scale_y_continuous(limits = c(0, 4),
                         breaks = seq(0, 4, by = 0.4),
                         expand = c(0, 0)) +
      scale_x_continuous(expand = c(0, 0),
                         breaks = c(sizepros$poschr),
                         labels = c(sizepros$chromosome),
                         guide = guide_axis(angle = 45)) +
      labs(x = "Chromosomes", y = "Copy Number") +
      theme(
        axis.text.x = element_markdown(size=12),
        axis.text.y = element_markdown(size=12),
        plot.background = element_rect(fill = "#F3FAFE"),
        panel.background = element_blank()
      ) +
      theme_classic()+
      coord_cartesian(xlim = c(input$xSlider[1], input$xSlider[2]))
    p<-ggplotly(p,tooltip = "text") %>% layout(xaxis = list(tickangle = -45))
    p
  })
  
  output$idiogram <- renderImage({
    if(input$showIdiogram){
      req(idiogram())
      # Specify the path to the local image file
      filename <- idiogram()
      
      # Read the image and resize it
      img <- image_read(filename)
      img_resized <- image_scale(img, "33%")  # Adjust the percentage as needed
      
      # Save the resized image temporarily
      temp_filename <- tempfile(fileext = ".png")
      image_write(img_resized, temp_filename)
      
      # Return a list containing the filename and the type of image
      list(src = temp_filename, contentType = "image/png")
    } else {
      # Return NULL when the checkbox is not checked
      NULL
    }
  },
  deleteFile = TRUE  # Delete the temporary file after rendering
  )
  
  output$table <- renderDT({
    if(input$displayTable){
      req(dfcnr(),dfcns(),dfcall())
      dfcnr<-dfcnr()
      dfcns<-dfcns()
      dfcall<-dfcall()
      dfcnr <- dfcnr %>%
        mutate(lrr = 2 * 2 ^ log2)
      for (i in 1:length(dfcnr$chromosome)) {
        dfcnr$order[i] <- i
      }
      
      dfcns <- dfcns %>%
        mutate(cncall = 2 * 2 ^ log2) %>%
        select(chromosome, start, end, cncall,weight) %>%
        filter(chromosome != "chrY")%>%
        filter(weight >= 15)
      dfcall <- dfcall %>% filter(chromosome == "chrY") %>% 
        select(chromosome,start,end,cn)
      names(dfcall)[names(dfcall) == "cn"] <- "cncall"
      dfcns<-merge(dfcall,dfcns,all=TRUE)
      df <- merge(dfcnr, dfcns, by = c("chromosome", "end"), all = TRUE) %>%
        select(-start.y)
      names(df)[names(df) == "start.x"] <- "start"
      df <- merge(df, dfcns, by.x = c("chromosome", "start"),
                  by.y = c("chromosome", "start"), all = TRUE) %>%
        select(-end.y)
      names(df)[names(df) == "end.x"] <- "end"
      df[is.na(df)] <- 0
      df$cncall <- df$cncall.x + df$cncall.y
      df <- df %>% arrange(order)
      df <- df %>% select(chromosome, start, end, lrr, cncall)
      df <- df %>% mutate(lrr=ifelse(chromosome=="chrY",dfcall$cncall,lrr))
      
      CNvector <- dfcns$cncall
      CNcall <- as.numeric(df$cncall)
      
      start_index <- vector()
      end_index <- vector()
      
      for (i in 1:length(CNvector)) {
        start_index[i] <- which(CNcall == CNvector[i])[1]
      }
      start_index<-na.omit(start_index)
      for (i in 1:length(CNvector)) {
        end_index[i] <- which(CNcall == CNvector[i])[2]
      }
      end_index<-na.omit(end_index)
      for (i in 1:length(start_index)) {
        df$cncall[(start_index[i] + 1):(end_index[i] - 1)] <- df$cncall[start_index[i]]
      }
      df <- df%>%
        mutate(cncall=ifelse(chromosome!=c("chrX","chrY") & cncall==0, NA, cncall))%>%
        na.omit()
      hg38size<-data.frame(
        chromosome=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8",
                     "chr9","chr10","chr11","chr12","chr13","chr14","chr15",
                     "chr16","chr17","chr18","chr19","chr20","chr21","chr22",
                     "chrX","chrY"),
        length=c(248956422,242193529,198295559,190214555,181538259,170805979,
                 159345973,145138636,138394717,135086622,133797422,133275309,
                 114364328,107043718,101991189,90338345,83257441,80373285,
                 58617616,64444167,46709983,50818468,156040895,57227415)
      )
      striphg38 <- hg38size %>%
        rename_all(tolower) %>%
        mutate(pos = cumsum(length / min(length))) %>%
        mutate(ymin = 0, ymax = 4,
               xmin = c(0, pos[1:23]), xmax = pos,
               fill = rep(c("xam", "trang"), length.out = nrow(hg38size))) %>%
        pivot_longer(c(ymin, ymax), values_to = "y", names_to = "yy") %>% select(-yy)
      sizepros <- hg38size %>%
        rename_all(tolower) %>%
        mutate(propos = length / min(length)) %>% select(-length) %>%
        mutate(pos = cumsum(propos)) %>%
        mutate(c = c(0, pos[1:23])) %>%
        mutate(poschr = propos / 2 + c)
      
      data <- df %>%
        inner_join(., sizepros, by = "chromosome") %>%
        group_by(chromosome) %>%
        mutate(position = cumsum(propos / length(chromosome))) %>%
        mutate(point_x = c + position) %>%
        select(chromosome, start, end, lrr, cncall, point_x)
      
      
      data <- data %>%
        mutate(lrr = ifelse(lrr > cncall + 0.1, NA, lrr)) %>% na.omit()
      data <- data %>%
        mutate(lrr = ifelse(lrr < cncall - 0.1, NA, lrr)) %>% na.omit()
      cncallcount<-table(data$cncall)
      singletons<-names(cncallcount[cncallcount==1])
      data<-data[!data$cncall %in% singletons, ]
      data <- data %>%
        group_by(chromosome,cncall) %>%
        mutate(lrr = smth.gaussian(lrr)) %>% na.omit()%>%ungroup()
      
      table<-data%>%select(chromosome,lrr)
      names(table)[names(table) == "lrr"] <- "cn"
      table
    }
  })
}

shinyApp(ui, server)

