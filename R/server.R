library(shiny)
library(shinyBS)
library(RColorBrewer)
library(biomaRt)
library(Biobase)
library(reshape2)
library(ggplot2)
library(dplyr)
library(tidyr)
library(plotly)
library(shinyjs)
library(htmlwidgets)
library(DT)
library(scExtras)
library(SeuratDisk)
library(Seurat)
library(cowplot)
library(data.table)
library(NMF)
library(tibble)
library(network)
library(igraph)
library(shinyBS)
library(scExtras)
library(slingshot)
library(aws.s3)
source("functions.R")

#Specify color palette for the tSNE and UMAP plots
cpallette=c("#64B2CE", "#DA5724", "#74D944", "#CE50CA", "#C0717C", "#CBD588", "#5F7FC7",
            "#8B4484", "#D3D93E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD",
            "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#D1A33D",
            "#8A7C64", "#599861","#000099","#FFCC66","#99CC33","#CC99CC","#666666", "#695F74","#0447F9",
            "#89134F","#2CF7F0","#F72C35","#A5B617","#B05927","#B78ED8")

Sys.setenv(
  "AWS_ACCESS_KEY_ID" = "xxxxx",
  "AWS_SECRET_ACCESS_KEY" = "xxxx",
  "AWS_DEFAULT_REGION" = "us-east-1"
)
bucket="referencemapping"

server <- function(input, output,session) {
  
####### Display project list ####
   #Read the parameter file
  readexcel = reactive({
    file = read.csv("data/reflist.csv")
    return(file)
  })
  
  #Get Project list and populate drop-down
  output$reflist = renderUI({
    excel=readexcel()
    prj=as.character(excel$Reference)
    selectInput("refs","Select a reference",as.list(sort(as.character(prj))))
  })
  
  
  #display project list in Dashboard
  datasetTable <- reactive({
    user=input$username
    file=read.csv('data/reflist.csv',stringsAsFactors = F)
    return(file)
  })
  
  output$datasetTable = DT::renderDataTable({
    withProgress(session = session, message = 'Loading...',detail = 'Please Wait...',{
      DT::datatable(datasetTable(),
                    extensions = c('Buttons','Scroller'),
                    options = list(dom = 'Bfrtip',
                                   searchHighlight = TRUE,
                                   pageLength = 20,
                                   lengthMenu = list(c(30, 50, 100, 150, 200, -1), c('30', '50', '100', '150', '200', 'All')),
                                   scrollX = TRUE,
                                   buttons = c('csv')
                    ),rownames=FALSE,selection = list(mode = 'single', selected =1),escape = F)
    })
  })
  
####  load ref and query data  #######
  #Load Reference
  refload <- reactive({
        file_names <- get_bucket_df(bucket)
        inFile = paste(as.character(input$refs),'.rds',sep = '')
        #reference <- readRDS(inFile)
        reference <- s3readRDS(object = file_names[file_names$Key == inFile, "Key"], bucket = bucket)
        return(reference)
  })
  
  #Load Query
  queryload <- reactive({
    file=input$rdatafileupload
    scrna <- readRDS(file$datapath)
      return(scrna)
  })
  
  ####Get input args ####
  output$refassay = renderUI({
    withProgress(session = session, message = 'Loading..',{
    ref=refload()
    assays=names(ref@assays)
    selectInput("refassay","Select a reference assay",assays,"")
    })
  })
  
  output$queryassay = renderUI({
    withProgress(session = session, message = 'Loading..',{
    scrna=queryload()
    assays=names(scrna@assays)
    selectInput("queryassay","Select a query assay",assays,"")
    })
  })
  
  output$ndim=renderUI({
    sliderInput("ndim", "Dimension to use from reduction",min = 10, max = 100, value = 50,step=10)
  })
  
  output$ntrees=renderUI({
    sliderInput("ntrees", "Number of tress to use in nearest neighbor search",min = 5, max = 50, value = 20,step=5)
  })
  
  output$metadata = renderUI({
    withProgress(session = session, message = 'Loading..',{
    ref=refload()
    cols=colnames(ref@meta.data %>% select(starts_with(c("lineage","celltype","subtype","annotation"))))
    selectInput("metadata","Reference metadata to transfer",cols,selected=cols[1],multiple=TRUE)
    })
  })
  
  output$refred = renderUI({
    withProgress(session = session, message = 'Loading..',{
      ref=refload()
      cols=names(ref@reductions)
      selectInput("refred","Reference reduction",cols,selected=cols[1])
    })
  })
  
  output$reducmodel = renderUI({
    withProgress(session = session, message = 'Loading..',{
      ref=refload()
      cols=names(ref@reductions)
      selectInput("reducmodel","Reference reduction model",cols,selected=cols[1])
    })
  })
  ################## Map cells to reference  ###############
  mapref <- reactive({
    validate(need(input$mapbutton != 0,"Make your selections and click the Map button"))
    if(input$mapbutton == 0)
      return()
    isolate({
    withProgress(session = session, message = 'Loading Data',{
    reference=refload()
    scrna=queryload()
    })
      if(length(VariableFeatures(scrna))==0){
      withProgress(session = session, message = 'Loading Data',{
        scrna=FindVariableFeatures(scrna)
      })
    }
      withProgress(session = session, message = 'Finding Anchors',{
  anchors <- FindTransferAnchors(reference = reference,query = scrna,k.filter = NA,reference.assay = input$refassay,
                                 query.assay = input$queryassay,reference.reduction = input$refred,normalization.method = input$norm,
                                 features = intersect(rownames(x = reference), VariableFeatures(object = scrna)),dims = 1:input$ndim,
                                 n.trees = input$ntrees,mapping.score.k = 100,recompute.residuals = input$recompres)
      })
      withProgress(session = session, message = 'Mapping',{
  refdata=as.list(c(input$metadata))
  names(refdata)=input$metadata
  scrna_map<- MapQuery(
    anchorset = anchors,
    query = scrna,
    reference = reference,
    refdata = refdata,
    reference.reduction = input$refred,
    reduction.model = input$reducmodel
  )
    })
   })
  return(scrna_map)
  })

  ####### Compare Tsne plot with controls  ##########
  #generate variable list for left plot
  output$tsnea2 = renderUI({
    ref=refload()
    metadata=as.data.frame(ref@meta.data)
    var=colnames(metadata)
    selectInput("tsnea2","Select a Variable",var,"pick one")
  })
  
  #generate variable list for right plot
  output$tsneb2 = renderUI({
    scrna_map=mapref()
    metadata=as.data.frame(scrna_map@meta.data)
    var=colnames(metadata)
    selectInput("tsneb2","Select a Variable",var,"pick one")
  })
  
  output$reductiona = renderUI({
    ref=refload()
    var=names(ref@reductions)
    selectInput("reductiona","Select a dimensionality reduction",var,"pick one")
  })
  
  output$reductionb = renderUI({
    scrna_map=mapref()
    var=names(scrna_map@reductions)
    selectInput("reductionb","Select a dimensionality reduction",var,"pick one")
  })
  
  comptsne = reactive({
    ref= refload()
    scrna_map=mapref()
    plot1=DimPlot(object = ref,group.by =  input$tsnea2,label = input$checklabel1, pt.size = input$pointa2,label.size = 4,reduction = input$reductiona) + ggtitle("Reference") + theme(legend.position="bottom")
    plot2=DimPlot(object = scrna_map,group.by = input$tsneb2,label =input$checklabel2,  pt.size = input$pointa2,label.size = 7,reduction = input$reductionb) +ggtitle("Query") + theme(legend.position="bottom")
    p=plot1+plot2
    return(p)
  })
  
  #render final plot
  output$comptsne = renderPlot({
    input$mapbutton
    withProgress(session = session, message = 'Generating figure...',detail = 'Please Wait...',{
      # if(input$mapbutton == 0)
      #   return()
      # isolate({
        comptsne()
      #})
    })
  })
  
  ####### Feature plot with controls  ##########
  #generate variable list for left plot
  output$featurelist = renderUI({
    scrna_map=mapref()
    if(input$fplotoption=="gene"){
      options=sort(rownames(GetAssayData(object = scrna_map))) }
    else if(input$fplotoption=="meta"){
     options=colnames(scrna_map@meta.data %>% select(contains("predicted") & contains("score"))) }
    withProgress(session = session, message = 'Generating list...',detail = 'Please Wait...',{
      selectInput('featurelist', label='Feature Name',options,multiple=FALSE, selectize=TRUE,selected=options[1])})
  })
  
  output$featuregroupby = renderUI({
    scrna_map=mapref()
    metadata=as.data.frame(scrna_map@meta.data) %>% select(starts_with(c("var","predicted"))) %>% select(-contains("score"))
    var=colnames(metadata)
    selectInput("featuregroupby","Select metadata to group by",var,"pick one")
  })
  
  
  #Generate plot
  fplot= reactive({
    scrna_map=mapref()
    p1=DimPlot(scrna_map,group.by = input$featuregroupby ,label=F)+theme(legend.position = "bottom")
    p2=FeaturePlot(object = scrna_map,features = input$featurelist) 
    plot=p1+p2
    return(plot)
  })
  #render final plot
  output$fplot = renderPlot({
    fplot()
  })
  
  ####Download Data ####
  output$downloaddata <- downloadHandler(
    filename = function() {
      name=strsplit(input$rdatafileupload,split="[.]")[[1]][1]
      paste0(name,"_anno.RDS",sep="")
    },
    content = function(file){
      scrna_map=mapref()
      saveRDS(scrna_map,file)
    })
}#end of server