library(ggplot2)
library(DT)
library(ggtext)
library(tidyverse)
library(plotly)
library(smoother)
library(shiny)
library(magick)
library(patchwork)
library(RSQLite)
library(purrr)


# Initialize shinyjs
shinyjs::useShinyjs()



# Function to get base names
get_base_names <- function(db) {
  tables <- dbListTables(db)
  base_names <- unique(gsub("\\..*", "", tables))
  return(base_names)
}

ui <- tagList(
  navbarPage(
  "PCG®PGS",
  tabPanel("Statistic",
           plotOutput("piechart")),
  tabPanel("Sequencing Runs",
           tags$h1("Sequencing Runs"),
           textInput("search_samples", "Batch search:", value = ""),
           uiOutput("sampleButtons"),  # Dynamic buttons for sample sheets
             tags$hr(),
             DTOutput("sampleSheet")),
  tabPanel("Sample",
  tags$h1(),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      tags$h3("CNVdatabase v1.0", style="color: #000; font-weight: bold;"),
      textInput("filter", "Filter:"),
      selectizeInput("basename", "Select Sample", choices = NULL, options = list(
        placeholder = 'Search for a basename',
        maxOptions = 1000
      )),
      actionButton("deleteTable", "Delete Sample", class = "btn-danger"),
      tags$hr(),
      tags$h4("View",style="font-weight:bold;"),
      checkboxInput("showCircle", "Circle", value = TRUE),
      checkboxInput("showIdiogram", "Karyotype", value = TRUE)
    ),
    mainPanel(
      width = 9,
      fluidRow(
        column(3, plotOutput("Circle")),
        column(9, plotlyOutput("CNVchart", height = "450px"))
      ),
      tags$hr(),
      plotlyOutput("Karyotype",height="700px",width="80%"),
      tags$hr(),
      tags$h2("CNV Table",style="color: #000;"),
      DTOutput("CNVtable",width = "80%"),
      tags$hr(),
      tags$h2("QC",style="color: #000;"),
      DTOutput("QC",width = "80%"),
      tags$h2("Decision track",style="color: #000;"),
      DTOutput("decision",width= "80%")
          )
        )
      )
    )
  )
server <- function(input, output, session) {
  # Connect to the SQLite database
  # Database connection
    db <- dbConnect(SQLite(), dbname = "/home/anegin97/Bioinformatics/appforclient/CNVdat/cnv_data.sqlite")
    sample_sheet_visibility <- reactiveVal(TRUE)
  observe({
    filter_text <- input$filter
    basenames <- get_base_names(db)
    basenames <- basenames[!grepl("_samplesheet$", basenames, ignore.case = TRUE)]
    filtered_basenames <- basenames[grepl(filter_text, basenames, ignore.case = TRUE)]
    updateSelectizeInput(session, "basename", choices = filtered_basenames, server = TRUE)
  })
  # Generate dynamic buttons for sample sheets
  output$sampleButtons <- renderUI({
    basenames <- get_base_names(db)
    sample_sheets <- basenames[grepl("samplesheet$", basenames, ignore.case = TRUE)]
    search_term <- input$search_samples
    if (search_term != "") {
      sample_sheets <- sample_sheets[grepl(search_term, sample_sheets, ignore.case = TRUE)]
    }
    buttons <- lapply(sample_sheets, function(sheet) {
      div(
        class = "col-sm-3",
        div(
          style = "display: flex; align-items: center; margin-bottom: 5px; margin-right: 0;",
          actionButton(inputId = paste0("button_", sheet), label = sheet, class = "btn-primary", style = "margin-right: 2px;"),
          actionButton(inputId = paste0("delete_button_", sheet), label = NULL, icon = icon("trash"), class = "btn-danger", style = "margin-right: 2px;"),
          downloadButton(outputId = paste0("download_tsv_", sheet), label = "TSV", icon = icon("download"), style = "margin-right: 2px;"),
          downloadButton(outputId = paste0("download_csv_", sheet), label = "CSV", icon = icon("download"))
        )
      )
    })
    fluidRow(
      do.call(tagList, buttons)
    )
  })
  observeEvent(input$search_samples, {
    basenames <- get_base_names(db)
    sample_sheets <- basenames[grepl("samplesheet$", basenames, ignore.case = TRUE)]
    search_term <- input$search_samples
    if (search_term != "") {
      sample_sheets <- sample_sheets[grepl(search_term, sample_sheets, ignore.case = TRUE)]
    }
    
    for (sheet in sample_sheets) {
      local({
        sheet_name <- sheet
        
        output[[paste0("download_tsv_", sheet_name)]] <- downloadHandler(
          filename = function() {
            paste0(sheet_name, ".tsv")
          },
          content = function(file) {
            data <- dbReadTable(db, sheet_name)
            write.table(data, file, sep = "\t", row.names = FALSE)
          }
        )
        
        output[[paste0("download_csv_", sheet_name)]] <- downloadHandler(
          filename = function() {
            paste0(sheet_name, ".csv")
          },
          content = function(file) {
            data <- dbReadTable(db, sheet_name)
            write.csv(data, file, row.names = FALSE)
          }
        )
      })
    }
  })
  
  # Observe button clicks and display the corresponding sample sheet
  observe({
    sample_sheets <- get_base_names(db)
    sample_sheets <- sample_sheets[grepl("samplesheet$", sample_sheets, ignore.case = TRUE)]
    
    lapply(sample_sheets, function(sheet) {
      observeEvent(input[[paste0("button_", sheet)]], {
        # Toggle the visibility of the sample sheet
        sample_sheet_visibility(!sample_sheet_visibility())
        
        if (sample_sheet_visibility()) {
          sample_sheet <- tryCatch({
            dbReadTable(db, sheet)
          }, error = function(e) {
            return(NULL)
          })
          
          if (is.null(sample_sheet)) {
            showNotification(paste("Sample sheet", sheet, "not found in the database"), type = "error")
          } else {
            output$sampleSheet <- renderDT({
              datatable(sample_sheet)
            })
          }
        } else {
          output$sampleSheet <- renderDT(NULL)
        }
      })
      
      observeEvent(input[[paste0("delete_button_", sheet)]], {
        showModal(modalDialog(
          title = "Confirm Delete",
          paste("Are you sure you want to delete the sample sheet:", sheet, "?"),
          easyClose = TRUE,
          footer = tagList(
            modalButton("Cancel"),
            actionButton("confirmDeleteSheet", "Delete", class = "btn-danger")
          )
        ))
        
        observeEvent(input$confirmDeleteSheet, {
          if (sheet %in% dbListTables(db)) {
            dbRemoveTable(db, sheet)
          }
          showNotification(paste("Sample sheet", sheet, "deleted."), type = "message")
          removeModal()
          updateSelectizeInput(session, "basename", choices = get_base_names(db), server = TRUE)
          output$sampleButtons <- renderUI({
            basenames <- get_base_names(db)
            sample_sheets <- basenames[grepl("samplesheet$", basenames, ignore.case = TRUE)]
            buttons <- lapply(sample_sheets, function(sheet) {
              list(
                actionButton(inputId = paste0("button_", sheet), label = sheet, class = "btn-primary"),
                actionButton(inputId = paste0("delete_button_", sheet), label = NULL, icon = icon("trash"), class = "btn-danger"),
                downloadButton(outputId = paste0("download_", sheet), label = NULL, icon = icon("download"))
              )
            })
            do.call(tagList, buttons)
          })
        })
      })
    })
  })
  observe({
    sample_sheets <- get_base_names(db)
    sample_sheets <- sample_sheets[grepl("samplesheet$", sample_sheets, ignore.case = TRUE)]
    
    lapply(sample_sheets, function(sheet) {
      output[[paste0("download_", sheet)]] <- downloadHandler(
        filename = function() {
          paste(sheet, ".tsv", sep = "")
        },
        content = function(file) {
          sample_sheet <- tryCatch({
            dbReadTable(db, sheet)
          }, error = function(e) {
            return(NULL)
          })
          
          if (!is.null(sample_sheet)) {
            write.table(sample_sheet, file, sep = "\t", row.names = FALSE, quote = FALSE)
          } else {
            showNotification(paste("Sample sheet", sheet, "not found in the database"), type = "error")
          }
        }
      )
    })
  })
  
  # Load and render tables based on selected basename
  observeEvent(input$load, {
    req(input$basename)  # Ensure a basename is selected before proceeding
    basename <- input$basename
  })
  
  dfcnr<-reactive({
    req(input$basename)
    basename<-input$basename
    dfcnrname<-paste0(basename,".dedup.cnr")
    if (dfcnrname %in% dbListTables(db)) {
      dbReadTable(db, dfcnrname)
    } else {
      NULL
    }
  })
  dfcns<-reactive({
    req(input$basename)
    basename<-input$basename
    dfcnsname<-paste0(basename,".dedup.cns")
    if (dfcnsname %in% dbListTables(db)) {
      dbReadTable(db, dfcnsname)
    } else {
      NULL
    }
  })
  dfcall<-reactive({
    req(input$basename)
    basename<-input$basename
    dfcallname<-paste0(basename,".dedup.call.cns")
    if (dfcallname %in% dbListTables(db)) {
      dbReadTable(db, dfcallname)
    } else {
      NULL
    }
  })
  QC<-reactive({
    req(input$basename)
    basename<-input$basename
    QCname<-paste0(basename,".QC.tsv")
    if (QCname %in% dbListTables(db)) {
      dbReadTable(db, QCname) %>%
        mutate(QCStatus = if_else(`Number.of.mapped.reads`/`Number.of.total.reads` > 0.75, "PASS", "FAIL"))
    } else {
      NULL
    }
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
    req(data, input$basename)
    data<-data()
    name<-input$basename
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
      req(data,input$basename,QC)
      data<-data()
      QC<-QC()
      QCStatus<-QC$QCStatus
      name<-input$basename
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
              axis.text.x = element_markdown(size=12,face="bold"),
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
              axis.text.x = element_markdown(size=12,face="bold"),
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
  
  output$DataTable <- renderDT({
    if(input$displayTable){
      req(data)
      data<-data()
      table<-data%>%select(chromosome,lrr)
      names(table)[names(table) == "lrr"] <- "cn"
      table
    }
  })
  
  CNVtable<-reactive({
      req(dfcns(),dfcall())
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
  output$CNVtable<-renderDT({
    req(CNVtable())
    CNVtable<-CNVtable()
    })
  output$QC<-renderDT({
    req(QC())
    QC<-QC()
  })
  observe({
    req(input$basename,CNVtable,QC)
    CNVtable<-CNVtable()
    QC<-QC()
    name<-input$basename
    basenames <- get_base_names(db)
    sample_sheets <- basenames[grepl("samplesheet$", basenames, ignore.case = TRUE)]
    flowcellsample <- strsplit(name, "-")[[1]][1]
    samplesheet_name <- sample_sheets[grepl(flowcellsample, sample_sheets, ignore.case = TRUE)]
    if (length(samplesheet_name) > 0) {
      samplesheet <- tryCatch({
        dbReadTable(db, samplesheet_name)
      }, error = function(e) {
        showNotification("Error reading samplesheet from the database", type = "error")
        return(NULL)
      })
      
      if (!is.null(samplesheet)) {
        # Create empty columns QCStatus and Report if they don't exist
        if (!"QCStatus" %in% colnames(samplesheet)) {
          samplesheet$QCStatus <- NA
        }
        if (!"Report" %in% colnames(samplesheet)) {
          samplesheet$Report <- NA
        }
        
        # Create the Summary data frame based on the CNVtable
        if (nrow(CNVtable) > 0) {
          Summary <- data.frame(Sample_ID = name,
                                QCStatus = QC$QCStatus,
                                Report = "Aeuploidy")
        } else {
          Summary <- data.frame(Sample_ID = name,
                                QCStatus = QC$QCStatus,
                                Report = "No Abnormalities Detected")
        }
        
        # Update existing QCStatus and Report columns in samplesheet
        samplesheet <- samplesheet %>%
          mutate(QCStatus = ifelse(Sample_ID == Summary$Sample_ID, Summary$QCStatus, QCStatus),
                 Report = ifelse(Sample_ID == Summary$Sample_ID, Summary$Report, Report))
        
        # Write the updated samplesheet back to the database
        tryCatch({
          dbWriteTable(db, samplesheet_name, samplesheet, overwrite = TRUE)
        }, error = function(e) {
          showNotification("Error writing samplesheet to the database", type = "error")
        })
      } else {
        showNotification("No matching samplesheet found in the database", type = "error")
      }
    } else {
      showNotification("No matching samplesheet found in the database", type = "error")
    }
  })
  observeEvent(input$deleteTable, {
    showModal(modalDialog(
      title = "Confirm Delete",
      "Are you sure you want to delete all tables related to this sample?",
      easyClose = TRUE,
      footer = tagList(
        modalButton("Cancel"),
        actionButton("confirmDelete", "Delete", class = "btn-danger")
      )
    ))
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

  # Confirm delete action
  observeEvent(input$confirmDelete, {
    req(input$basename)
    basename <- input$basename
    table_suffixes <- c(".dedup.cnr", ".dedup.cns", ".dedup.call", ".QC")
    tables_to_delete <- paste0(basename, table_suffixes)
    
    for (table_name in tables_to_delete) {
      if (table_name %in% dbListTables(db)) {
        dbRemoveTable(db, table_name)
      }
    }
    
    # Refresh the basenames and exclude samplesheets
    bn <- get_base_names(db)
    bn <- bn[!grepl("_samplesheet$", bn, ignore.case = TRUE)]
    
    showNotification(paste("All tables related to sample", basename, "deleted."), type = "message")
    removeModal()
    updateSelectizeInput(session, "basename", choices = bn, server = TRUE)
  })
    output$piechart <- renderPlot({
      # Function to read samplesheet tables and check for 11 columns
      read_samplesheets <- function(table_name) {
        table_data <- dbReadTable(db, table_name)
        
        required_columns <- c(
          "Sample_ID", "Cycle_ID", "Embryo_ID", "Cell_Type", 
          "DNA_ng_ul", "Sample_Plate", "Sample_Well", 
          "I7_Index_ID", "I5_Index_ID","Report","QCStatus"
        )
        
        # Check if the table has both "Report" and "QCStatus+" columns
        if (all(c("Report", "QCStatus") %in% colnames(table_data))) {
          selected_data <- table_data %>% select(all_of(required_columns))
          return(selected_data)
        } else {
          warning(paste("Warning: Table", table_name, "is missing required columns. Skipping this table."))
          return(NULL)
        }
      }
      
      
      # Get the list of tables and filter for samplesheet tables
      tables <- dbListTables(db)
      samplesheet_tables <- grep("samplesheet", tables, value = TRUE)
      
      # Read and combine samplesheet tables with 11 columns
      samplesheets_list <- lapply(samplesheet_tables, read_samplesheets)
      samplesheets_combined <- do.call(rbind, Filter(Negate(is.null), samplesheets_list))
      
      # Proceed only if data was successfully combined
      if (!is.null(samplesheets_combined) && nrow(samplesheets_combined) > 0) {
        Report <- samplesheets_combined %>%
          group_by(Report) %>%
          summarize(count = n()) %>%
          mutate(total = sum(count),
                 perc = count / total * 100) %>%
          mutate(perc = round(perc, 2),
                 labels = paste0(perc, "%"))
        
        # Plotting
        ggplot(Report, aes(x = "", y = perc, fill = Report)) +
          geom_col(color = "black") +
          geom_label(aes(label = labels), color = "white", size = 5,
                     position = position_stack(vjust = 0.5),
                     show.legend = FALSE) +
          scale_fill_manual(name = NULL,
                            breaks = c("Aeuploidy", "No Abnormalities Detected"),
                            values = c("grey", "blue")) +
          labs(title = "Embryo ratio") +
          coord_polar(theta = "y") +
          theme_void() +
          theme(legend.position = "bottom",
                plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
                legend.text = element_text(size = 16))
      } else {
        print("Debug: No valid data to plot.")
      }
    })
    
  
  # Close the database connection when the app is closed
  onStop(function() {
    DBI::dbDisconnect(db)
  })
}
shinyApp(ui, server)
