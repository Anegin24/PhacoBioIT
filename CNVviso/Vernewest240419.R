library(ggplot2)
library(DT)
library(ggtext)
library(tidyverse)
library(dplyr)
library(plotly)
library(smoother)
library(shiny)
library(magick)
library(patchwork)
library(RSQLite)
library(purrr)
# Initialize shinyjs
shinyjs::useShinyjs()



selectizeSelectInput <- function(inputId, label, choices, selected = NULL) {
  selectizeInput(inputId, label, choices = choices, selected = selected, options = list(search = TRUE))
}

ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      body {
        background-color: #fff;
        margin-top:20px; 
      }
      .dataTables_wrapper {
        background-color: #fff;
        color: #000;
        margin-top: 20px;
        border-radius:3px;
      }
      #table-container {
        float: right;
        margin-left: 20px;
      }
      .main-panel {
        overflow-x: auto;
      }
      .file-item {
        cursor: pointer;
        font-size: 15px;
        margin-bottom: 5px;
        padding: 5px;
        background-color: #00a0be;
        color: #fff;
        border-radius: 3px;
        border-color: #ccc;
      }
      .file-item:hover {
        background-color: #222;
      }
      .file-list-container {
        font-site:12px;
        font-weight:bold;
        height: 300px; /* Adjust height as needed */
        overflow-y: auto;
      }
      .filter-box {
        margin-bottom: 10px;
      }
      .panel-analysis {
        background-color: #fff; /* Change background color */
        color: #fff; /* Change text color */
        font-size: 15px; /* Change font size */
        border-radius: 5px; /* Add border radius */
        padding: 10px; /* Add padding */
        margin-bottom: 20px; /* Add margin */
      }
       .panel-report {
        background-color: #fff; /* Change background color */
        color: #fff; /* Change text color */
        font-size: 15px; /* Change font size */
        border-radius: 5px; /* Add border radius */
        padding: 10px; /* Add padding */
        margin-bottom: 20px; /* Add margin */
       }
       .panel-summary {
        background-color: #fff; /* Change background color */
        color: #fff; /* Change text color */
        font-size: 15px; /* Change font size */
        border-radius: 5px; /* Add border radius */
        padding: 10px; /* Add padding */
        margin-bottom: 20px; /* Add margin */
       }
       
        .nav-tabs > li:nth-child(1) > a {
        background-color: #00a0be; 
        color: white;
        border: 1px solid #00a0be;
        border-radius: 5px;
        padding: 6px 20px;
        font-size:20px;
        font-weight:bold;
      }
    
        .nav-tabs > li:nth-child(2) > a {
        background-color: #00a0be; 
        color: white;
        border: 1px solid #00a0be; 
        border-radius: 5px;
        padding: 6px 20px;
        font-size:20px;
        font-weight:bold;
      }   
       .nav-tabs > li:nth-child(3) > a {
        background-color: #00a0be; 
        color: white;
        border: 1px solid #00a0be; 
        border-radius: 5px;
        padding: 6px 20px;
        font-size:20px;
        font-weight:bold;
       }
       .progress-bar {
       background-color: #00a0be;
       }
       .js-range-slider {
       background-color: #00a0be;
       }
    "))
  ),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      style = "height: 1155px; overflow-y: auto;",
      tags$h3("CNVviso 1.0", style="color: #000; font-weight: bold;"),
      fileInput("file", "Import Data (.zip)", multiple = TRUE, accept = c(".zip")),
      textInput("fileFilter", "Sample ID"),  # Text input for filtering files
      uiOutput("fileListContainer"),  # Container for the file list
      verbatimTextOutput("directorypath"),
      uiOutput("file_checkboxes",value=FALSE),
      actionButton("showModal", "Save Sample"),
      tags$hr(),
      div(class = "view-box",
          tags$h5("View", style = "color: #000; font-weight: bold;"),
          checkboxInput("showCircle", "Circle", value = TRUE),
          checkboxInput("displayCNV", "CNV Table", value = TRUE),
          checkboxInput("showIdiogram", "Karyotype", value = TRUE)
      ),
      sliderInput("xSlider", "Chromosome size", min = 0, max = 66, value = c(0, 66))
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        id="tabs",
        tabPanel("Analysis",
                 class="panel-analysis",
                 fluidRow(
                   column(3, plotOutput("Circle")),
                   column(9, plotlyOutput("CNVchart", height = "450px",width="100%"))
                 ),
                 tags$hr(),
                 plotlyOutput("Karyotype",height="600px",width="80%")),
        tabPanel("Report",
                 class="panel-report",
                 tags$hr(),
                 tags$h2("CNV Table",style="color: #000;"),
                 DTOutput("CNVreport",width = "80%"),
                 tags$hr(),
                 tags$h2("QC",style="color: #000;"),
                 DTOutput("QC",width = "80%")
        ),
        tabPanel("Decision Track",
                 class="panel-decision",
                 DTOutput("decision",width = "100%")
        )    
      )
    )
  )
)

server <- function(input, output, session) {
  output$file_checkboxes <- renderUI({
    req(input$file)
    file_names <- basename(input$file$name)
    checkboxGroupInput("selected_files", "Save to database",
                       choices = file_names, selected = file_names)
  })
  observeEvent(input$showModal, {
    showModal(modalDialog(
      title = "Confirm Import",
      "Are you sure you want to save the selected samples?",
      footer = tagList(
        modalButton("Cancel"),
        actionButton("ImportSample", "Import")
      )
    ))
  })
  observeEvent(input$ImportSample, {
    req(input$file)
    selected_files <- input$selected_files
    files <- input$file
    files <- files$datapath[basename(files$name) %in% selected_files]
    removeModal()
    process_file <- function(file) {
      zip_file_name <- gsub(".zip", "", basename(file))
      extraction_directory <- paste0("unzipped_files_", zip_file_name)
      unzip(file, exdir = extraction_directory)
      list.files(extraction_directory, pattern = "\\.(dedup\\.cnr|dedup\\.cns|dedup\\.call.cns|QC\\.tsv)$", full.names = TRUE)
    }
    extracted_files <- unlist(lapply(files, process_file))
    db <- dbConnect(SQLite(), dbname = "/home/anegin97/Bioinformatics/appforclient/CNVdat/cnv_data.sqlite")
    process_extracted_file <- function(file) {
      # Read the file content
      file_content <- read.table(file, header = TRUE, sep = "\t")
      
      # Extract the table name from the file name
      table_name <- basename(file)
      
      # Write the file content to the database
      dbWriteTable(db, table_name, file_content, overwrite = TRUE)
    }
    
    # Apply the process_extracted_file function to each extracted file
    lapply(extracted_files, process_extracted_file)
    dbDisconnect(db)
    lapply(files, function(file) {
      zip_file_base_name <- tools::file_path_sans_ext(basename(file))
      extraction_directory <- paste0("unzipped_files_", zip_file_base_name)
      unlink(extraction_directory, recursive = TRUE)
    })
  })
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
  
  QC <- reactive({
    req(extracted_files())
    file_path_call <- extracted_files() %>% keep(., grepl(".QC.tsv", .)) %>% first()
    read_tsv(file_path_call)%>%
      mutate(QCStatus = if_else(`Number-of-mapped-reads`/`Number-of-total-reads` > 0.65, "PASS", "FAIL"))
  })
  
  cyto <- reactive({
    read_tsv("https://raw.githubusercontent.com/Anegin24/Microbiome-analysis/main/cytoBand.txt")
  })
  annotation<-reactive({
    read_tsv("https://raw.githubusercontent.com/Anegin24/Microbiome-analysis/main/annotation.tsv")
  })
  
  data <- reactive({
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
      filter(weight >= 50)
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
    data<- data%>%
      mutate(cnline=case_when(cncall>2.3 ~ cncall,
                              cncall<1.7 ~ cncall,
                              cncall>1.7 & cncall<2.3 ~ 2))
  })
  output$CNVchart<-renderPlotly({
    req(data, input$selected_file)
    data<-data()
    name<-input$selected_file
    name <- sub("\\..*", "", name)
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
    data<-data%>%
      mutate(lrr=round(lrr,2))
    p<-data%>%
      ggplot(aes(x = point_x, y = lrr))+
      geom_vline(size = 0.3, xintercept = c(unique(sizepros$pos)), linetype = "dotted", color = "grey") +
      geom_hline(size=0.3,yintercept=4,linetype="solid")+
      geom_vline(size=0.3,xintercept=66,linetype="solid")+
      geom_hline(size=0.3,yintercept = c(1.7,2.3),linetype="dotted",color="gray")+
      geom_point(aes(text=paste("Chromosome:",chromosome,"<br>Start:",start,"<br>End:",end,"<br>CN:",lrr)),size = 1.3, shape = 21, color = alpha("black",0.6), fill = "chartreuse3") +
      geom_line(aes(x=point_x,y=cnline),linewidth = 0.9, color = "green",show.legend = FALSE) +
      geom_line(aes(x=point_x,y=lrr),linewidth = 0.1, color = "darkmagenta",show.legend=FALSE)+
      geom_hline(size = 0.3, yintercept = c(1, 3)) +
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
      labs(x = "Chromosomes", y = "Copy Number",title=name) +
      theme(
        axis.text.x = element_markdown(size=16),
        axis.text.y = element_markdown(size=16),
        axis.title = element_markdown(size=16,face="bold"),
        plot.background = element_rect(fill = "gray"),
        panel.background = element_blank()
      ) +
      theme_classic()+
      coord_cartesian(xlim = c(input$xSlider[1], input$xSlider[2]))
    p<-ggplotly(p,tooltip = "text") %>% 
      layout(xaxis = list(tickangle = -45),
             title = list(text = paste0('<span style="border: 2px solid black; padding: 5px; font-weight: bold;">', name, '</span>'), x = 0.5))
    p
  })
  output$Circle<-renderPlot({
    if(input$showCircle){
      req(data,input$selected_file)
      data<-data()
      QC<-QC()
      QCStatus<-QC$QCStatus
      name<-input$selected_file
      name <- sub("\\..*", "", name)
      data<-data[!duplicated(data$cncall), ]
      filld <- data.frame(xmin = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23),
                          xmax = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24),
                          chromosome = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6","chr7","chr8","chr9","chr10",
                                         "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21",
                                         "chr22","chrX","chrY"),
                          ymin = c(4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4),
                          ymax = c(5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5))
      filld<-merge(data,filld,by="chromosome",all=TRUE)
      filld <- filld%>%
        mutate(ymin=case_when(cncall>1.7 & cncall<2.3 ~ ymin,
                              cncall > 2.3 ~ 5,
                              cncall >0.5 & cncall < 1.7 ~ 3,
                              cncall < 0.5 ~ NA),
               ymax=case_when(cncall>1.7 & cncall <2.3 ~ ymax,
                              cncall > 2.3 ~ 6,
                              cncall> 0.5 & cncall <1.7 ~ 4,
                              cncall <0.5 ~ NA))%>%
        mutate(fill=case_when(cncall>2.3 ~ "gain",
                              cncall<1.7 ~ "loss",
                              1.7<cncall&cncall<2.3~"normal"))
      
      filld <- filld %>%
        mutate(midpoint = (xmin + xmax) / 2,
               ycenter = (ymin + ymax) / 2)
      
      ggplot(filld,aes(fill=fill)) +
        geom_segment(aes(x = xmin, xend = xmax, y = 4.5), color = "black", size = 0.5) +
        geom_segment(aes(x = xmin, xend = xmax, y = 3.5), color = "red", size = 0.5) +
        geom_segment(aes(x = xmin, xend = xmax, y = 5.5), color = "darkgreen", size = 0.5) +
        geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), color = "black") +
        scale_fill_manual(name=NULL,
                          breaks=c("gain","loss","normal"),
                          values=c("yellow4","darkred","darkblue"),
                          guide = "none")+
        coord_polar() +
        geom_text(aes(x = midpoint, y = ycenter, label = chromosome), size = 2.5, fontface = "bold",color="white") +
        scale_y_continuous(limits = c(0, 6)) +
        annotate("text", x = 0, y = 0, label = paste0(name, "\n QC: ", QCStatus), size = 3, fontface = "bold") +
        theme_void()
    }
  })
  output$Karyotype <- renderPlotly({
    if(input$showIdiogram){
      req(dfcns,dfcall,cyto)
      dfcns<-dfcns()
      dfcall<-dfcall()
      cyto<-cyto()
      print(str(dfcns))
      dfcns <- dfcns %>%
        mutate(cncall = 2 * 2 ^ log2) %>%
        select(chromosome, start, end, cncall,weight) %>%
        filter(chromosome != "chrY")%>%
        filter(weight >= 50)
      dfcall <- dfcall %>% filter(chromosome == "chrY") %>% 
        select(chromosome,start,end,cn)
      names(dfcall)[names(dfcall) == "cn"] <- "cncall"
      dfcns<-merge(dfcall,dfcns,all=TRUE)
      chrX <- dfcns %>% filter(chromosome =="chrX")
      chrY <- dfcns %>% filter(chromosome =="chrY")%>%mutate(fill=case_when(cncall>1.3~"gain"))
      if (chrY$cncall==0){
        chrX<-chrX%>%
          mutate(fill=case_when(cncall>1.7 & cncall < 2.3 ~"No change",
                                cncall<1.7 ~"loss",
                                cncall>2.3 ~"gain"))
      }
      if (chrY$cncall==1){
        chrX<-chrX%>%
          mutate(fill=case_when(cncall>1.3 ~ "gain",
                                cncall<0.7 ~ "loss",
                                cncall < 1.3 & cncall > 1.7 ~ "No change"))
      }
      sexchr<-merge(chrX,chrY,all=TRUE)
      dfcns<-dfcns%>%filter(chromosome!="chrX")%>%filter(chromosome!="chrY")
      dfcns<-dfcns%>%
        select(chromosome,start,end,cncall)%>%
        mutate(fill=case_when(cncall>2.3 ~ "gain",
                              cncall<1.7 ~ "loss",
                              1.7<cncall&cncall<2.3~"normal"))
      dfcns<-merge(dfcns,sexchr,all=TRUE)%>%
        filter(fill=="loss"|fill=="gain")
      cytoband_colors <- c(
        gpos100 = "#000000",   # Black for 100% darkly staining (gpos100)
        gpos75  = "#555555",   # Dark gray for 75% darkly staining (gpos75)
        gpos50  = "#AAAAAA",   # Medium gray for 50% darkly staining (gpos50)
        gpos25  = "#DDDDDD",   # Light gray for 25% darkly staining (gpos25)
        gneg    = "#FFFFFF",   # White for G-negative bands (gneg)
        acen    = "#FF0000",   # Red for centromeric regions (acen)
        stalk   = "#CCCCCC",   # Light gray for stalk regions (stalk) - typically lighter
        gvar    = "lightblue"    # Blue for heterochromatic regions (gvar) - or dark blue/green
      )
      x1<-data.frame(
        chromosome=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12"),
        xmin=c(1,3,5,7,9,11,13,15,17,19,21,23),
        xmax=c(2,4,6,8,10,12,14,16,18,20,22,24)
      )
      x2<-data.frame(
        chromosome=c("chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"),
        xmin=c(1,3,5,7,9,11,13,15,17,19,21,23),
        xmax=c(2,4,6,8,10,12,14,16,18,20,22,24)
      )
      x<-merge(x1,x2,all=TRUE)
      dfcns<-inner_join(dfcns,x,by="chromosome")
      cyto1<-inner_join(cyto,x1,by="chromosome")
      cyto2<-inner_join(cyto,x2,by="chromosome")
      cyto<-merge(cyto1,cyto2,all=TRUE)
      cyto1<-cyto1%>%mutate(midpoint=(xmin+xmax)/2)
      cyto2<-cyto2%>%mutate(midpoint=(xmin+xmax)/2)
      dfcns1<-dfcns%>%filter(chromosome %in% c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12"))
      dfcns1 <- dfcns1 %>%
        mutate(x = ifelse(fill == "loss", xmax + 0.2, xmin - 0.2),
               xend = x)
      cyto1plot<-cyto1%>%
        ggplot(aes(xmin=xmin,xmax=xmax,ymin=start,ymax=end,fill=dye,text=position))+
        geom_rect(color="black")+
        scale_fill_manual(name=NULL,values=cytoband_colors)+
        scale_color_manual(name=NULL,breaks=c("loss","gain"),values=c("darkred","green"))+
        scale_y_continuous(expand=c(0,0))+
        scale_x_continuous(breaks=c(cyto1$midpoint), labels=c(cyto1$chromosome))+
        theme_classic()+
        theme(axis.line=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks=element_blank(),
              axis.title=element_blank(),
              axis.text.x = element_text(size=14,face="bold"),
              legend.position="none")
      if (length(dfcns1$chromosome) > 0) {
        cyto1plot <- cyto1plot +
          geom_segment(data = dfcns1, aes(x = x, xend = xend, y = start, yend = end, color = fill), 
                       size = 1.3, inherit.aes = FALSE)
      }
      dfcns2<-dfcns%>%filter(chromosome %in% c("chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"))
      dfcns2 <- dfcns2 %>%
        mutate(x = ifelse(fill == "loss", xmax + 0.2, xmin - 0.2),
               xend = x)
      cyto2plot<-cyto2 %>%
        ggplot(aes(xmin = xmin, xmax = xmax, ymin = start, ymax = end, fill = dye,text=position)) +
        geom_rect(color = "black") +
        scale_fill_manual(name = NULL, values = cytoband_colors) +
        scale_color_manual(name=NULL,breaks=c("loss","gain"),values=c("darkred","green"))+
        scale_y_continuous(expand = c(0, 0)) +
        scale_x_continuous(breaks = cyto2$midpoint, labels = cyto2$chromosome) +
        theme_classic() +
        theme(axis.line = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks = element_blank(),
              axis.title = element_blank(),
              axis.text.x = element_text(size=14,face="bold"),
              legend.position = "none")
      if (length(dfcns2$chromosome) > 0) {
        cyto2plot <- cyto2plot +
          geom_segment(data = dfcns2, aes(x = x, xend = xend, y = start, yend = end, color = fill), 
                       size = 1.3, inherit.aes = FALSE)
      }
      p1<-ggplotly(cyto1plot,tooltip="text")
      p2<-ggplotly(cyto2plot,tooltip="text")
      subplot(p1,p2,nrows=2)
    }
  })
  CNVtable<-reactive({
    req(dfcns,dfcall)
    dfcns<-dfcns()
    dfcall<-dfcall()
    dfcns <- dfcns %>%
      mutate(cncall = 2 * 2 ^ log2) %>%
      select(chromosome, start, end, cncall,weight) %>%
      filter(chromosome != "chrY")%>%
      filter(weight >= 50)
    dfcall <- dfcall %>% filter(chromosome == "chrY") %>% 
      select(chromosome,start,end,cn)
    names(dfcall)[names(dfcall) == "cn"] <- "cncall"
    dfcns<-merge(dfcall,dfcns,all=TRUE)
    chrX <- dfcns %>% filter(chromosome =="chrX")
    chrY <- dfcns %>% filter(chromosome =="chrY")%>%mutate(Type=case_when(cncall>1.3~"Gain"))
    if (chrY$cncall==0){
      chrX<-chrX%>%
        mutate(Type=case_when(cncall>1.7 & cncall < 2.3 ~"No change",
                              cncall<1.7 ~"Loss",
                              cncall>2.3 ~"Gain"))
    }
    if (chrY$cncall==1){
      chrX<-chrX%>%
        mutate(Type=case_when(cncall>1.3 ~ "Gain",
                              cncall<0.7 ~ "Loss",
                              cncall < 1.3 & cncall > 1.7 ~ "No change"))
    }
    sexchr<-merge(chrX,chrY,all=TRUE)
    sexchr<-sexchr%>%mutate(CopyNumber=cncall)
    dfcns<-dfcns%>%filter(chromosome != "chrX")%>%filter(chromosome != "chrY")%>%
      mutate(
        CopyNumber=cncall,
        Type=case_when(CopyNumber>=2.3 ~ "Gain",
                       1.7>=CopyNumber ~ "Loss",
                       CopyNumber>1.7 & CopyNumber<2.3 ~ "No change"))
    dfcns<-merge(dfcns,sexchr,all=TRUE)
    dfcns<-dfcns%>%
      mutate(size=end-start)%>%
      select(chromosome,start,end,size,CopyNumber,Type)%>%
      filter(Type!="No change")%>%
      mutate(CopyNumber=round(CopyNumber,2))
  })
  output$CNVreport<-renderDT({
    if(input$displayCNV){
      req(CNVtable)
      CNVtable<-CNVtable()
    }
  })
  output$QC<-renderDT({
    req(QC)
    QC<-QC()
    QC
  })
  
  output$decision<-renderDT({
    req(CNVtable,annotation)
    CNVtable<-CNVtable()
    annotation<-annotation()
    annotation<-annotation%>%
      filter(Chromosome %in% CNVtable$chromosome)%>%
      select(Chromosome,Start,End,Gene,`Cyto Location`,Phenotype,Inheritance)
    # Function to filter df2 based on each row of df1
    filter_by_row <- function(df1_row, df2) {
      df2 %>%
        filter(
          Chromosome == df1_row$chromosome,
          Start >= df1_row$start,
          End <= df1_row$end
        )
    }
    annotation <- CNVtable%>%
      group_split(row_number()) %>%
      map_dfr(~ filter_by_row(.x, annotation))
  })
}

shinyApp(ui, server)
