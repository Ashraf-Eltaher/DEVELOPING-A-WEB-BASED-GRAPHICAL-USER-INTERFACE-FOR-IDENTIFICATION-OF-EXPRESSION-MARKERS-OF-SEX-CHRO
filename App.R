library(shiny)
library(shinythemes)
source("Appcode.R")

genders <- c('Female','Klinefelter','Male','Turner')
ui <- fluidPage(theme = shinytheme('journal'),
                fluidRow(
                          column(2,tags$img(height = 118, width = 160, src="University-of-Bonn.png")),
                          column(2,offset = 2,tags$img(height = 90, width = 241, src="15823507_1362489713795219_6175969651089312785_n.jpg")),
                          column(2,offset = 2,tags$img(height = 100, width = 200, src="images.png"))
                        ),
                
                titlePanel(h3(tags$b("Identification of expression markers to distinguish sex chromosome abnormalities"),align = 'center')),
                br(),
                
                
                sidebarLayout(
                  sidebarPanel(
                        selectizeInput(inputId = "gender1",
                                  label = "Choose samples to compare",
                                  choices = genders,
                                  options = list(
                                    placeholder = 'Sample 1',
                                    onInitialize = I('function() { this.setValue(""); }'))),
                        conditionalPanel("input.gender1 != ''",
                            selectizeInput(inputId = "gender2",
                                      label = NULL,
                                      choices = '',
                                      options = list(
                                        placeholder = 'Sample 2',
                                        onInitialize = I('function() { this.setValue(""); }')))),
                        
                      
                    br(),
                    radioButtons('chr','Select chromosome(s)',choices = 
                                   c('all','X','autosomal')),
                    br(),
                    actionButton(inputId = "button", label = "Analyze", 
                                 style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                    
                    width = 3),
              
                  mainPanel(
                        tabsetPanel(
                          tabPanel("Exp.values",DTOutput("table")),
                          tabPanel("tSNE & PCA Plots",
                                   fluidRow(column(4,br(),br(),br(),plotOutput("tsne",width = "440px", height="365px")),
                                            column(4,offset = 3,(conditionalPanel("input.button != 0",
                                                   radioButtons('d','',choices = c('2D','3D'), selected = '3D', inline=T),offset = 2)),
                                                   conditionalPanel("input.d == '2D'",
                                                   plotOutput("plot2D",width = "400px", height="350px")),
                                                   conditionalPanel("input.d == '3D'",
                                                   rglwidgetOutput("plot3D",width = "400px", height="350px"))))),
                                   
                          tabPanel("Linear fit table",uiOutput('filter'),DTOutput("fittable"),
                                   conditionalPanel("input.button != 0",
                                    downloadButton("download", "Download",
                                                  style="color: #fff; background-color: green; border-color: Black;"))),
                          tabPanel("Volcano Plot",
                                   conditionalPanel("input.button != 0",
                                    fluidRow(column(4,numericInput('pvalue','p value threshold',value = 0.05,step = 0.01,width = 140)),
                                             column(4,numericInput('logFC','log FC threshold',value = abs(1),step = 0.1,width = 140)),
                                             column(2,style = "margin-top: 25px;",shinyWidgets::circleButton('updatePlot',icon = icon('redo'),status = 'primary', size = 'sm')))),
                              plotOutput("vplot",width = "650px", height="410px",click = 'click',brush = brushOpts('brush')),
                              conditionalPanel("input.button != 0",verbatimTextOutput('gene_name'))),
                          tabPanel("Heatmap",plotOutput("hplot",width = "550px", height="350px")),
                          tabPanel("Manhattan plot", plotOutput("mplot"),
                                   conditionalPanel("input.button != 0",
                                                    downloadButton("download2", "Download",
                                                                   style="color: #fff; background-color: green; border-color: Black;"))),
                          tabPanel("Gene ontology",tableOutput("goTable"))
                                         
                    ),width = 9
                    )))
                
                    


server <- function(input, output,session) {
  observe({
    input$gender1
    updateSelectizeInput(session, 'gender2',choices = genders[genders != input$gender1],options = list(
      placeholder = 'Sample 2',
      onInitialize = I('function() { this.setValue(""); }')))
  })
  observe({
    input$gender2
    updateSelectizeInput(session, 'gender3',choices = genders[!(genders %in% c(input$gender1,input$gender2))],options = list(
      placeholder = 'Sample 3',
      onInitialize = I('function() { this.setValue(""); }')))
  })
  
   observeEvent(input$button,{
    valuesDF <<- createValuesDF(input$gender1,input$gender2,input$chr)
    fitDF <<- fitData(valuesDF)
    manDF <- createValuesDF(input$gender1,input$gender2)
    
    output$table <- renderDT({datatable(showTable(valuesDF),style = 'bootstrap4')})
    output$tsne <- renderPlot({tsnePlot(valuesDF)})
    output$plot2D <- renderPlot({pcaPlotData2(createPCA(valuesDF))})
    output$plot3D <- renderRglwidget({
    pcaPlotData3(createPCA(valuesDF))
    rglwidget()})
    
    output$fittable <- renderDT({datatable(fitDF,caption = paste('Statistical comparison of ', input$gender1,' with ',
                                                                         input$gender2, ' expression values',sep = ''))})
    output$vplot <- renderPlot({isolate(volcanoPlotData(fitDF,input$pvalue,input$logFC))})
    output$filter <- renderUI({fluidRow(column(4,numericInput('tablePvalue','p value below',value = '',step = 0.01,width = 140)),
                                        column(4,numericInput('tableLogFC','log FC above',value = '',step = 0.1,width = 140)),
                                        column(2,style = "margin-top: 25px;",shinyWidgets::circleButton('updateTable', icon = icon('redo'),status = 'primary', size = 'sm')))})
  
                                     
    output$hplot <- renderPlot({heatmapData(valuesDF)})
    output$mplot <- renderPlot({manPlotData(manDF)})
    output$goTable  <- renderTable({goData(valuesDF)})
    })
  
  
  observeEvent(input$updatePlot,{
    output$vplot <- renderPlot({isolate(volcanoPlotData(fitDF,input$pvalue,input$logFC))})
  })

  observeEvent(input$updateTable,{
    filterDF <- filter.fitDF(valuesDF,input$tableLogFC,input$tablePvalue)
    newDF <- subset(valuesDF,rownames(valuesDF) %in% rownames(filterDF))
    output$table <- renderDT({datatable(newDF)})
    output$fittable <- renderDT({datatable(filterDF,caption = paste('Statistical comparison of ', input$gender1,' with ',
                                                                        input$gender2, ' expression values',sep = ''))})
    output$vplot <- renderPlot({volcanoPlotData(filterDF,input$pvalue,input$logFC)})
    output$pcaplot <- renderPlot({pcaPlotData(newDF,1,2)})
    output$hplot <- renderPlot({heatmapData(newDF)})
    #output$mplot <- renderPlot({manPlotData(newDF)})
    output$goTable  <- renderTable({goData(newDF)})
    #output$goPlot  <- renderPlot({goPlotData(newDF)})
  })
  
  fit.table <- reactive({filter.fitDF(valuesDF,1,0.05)})
  output$download <- downloadHandler(
    filename = function(){paste(input$gender1,input$gender2,".csv",sep = ' ')},
    content = function(file){
      write.csv(fit.table(),file)})
  
  
  observeEvent(input$click,{
  output$gene_name <- renderPrint({
    fitDF[fitDF$logFC == sort(fitDF$logFC)[MALDIquant::match.closest(input$click$x,sort(fitDF$logFC))],]
   })
  })
  
  man.table <- reactive({manTable[manTable$P.Value <0.05,]})
  output$download2 <- downloadHandler(
    filename = function(){paste(input$gender1,input$gender2, 'manhattan table',".csv",sep = ' ')},
    content = function(file){
      write.csv(man.table(),file)})
} 

shinyApp(ui = ui, server = server)