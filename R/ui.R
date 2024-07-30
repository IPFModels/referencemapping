library(shinydashboard)
library(shiny)
library(shinyBS)
library(plotly)
library(shinyjs)
library(reshape2)
library(visNetwork)
library(shinyBS)
library(bslib)
library(dashboardthemes)

options(shiny.sanitize.errors = FALSE)
options(shiny.maxRequestSize=600*1024^2) 
ui <- dashboardPage(skin = "blue",
  dashboardHeader(title = "Reference Mapping",titleWidth = 350,dropdownMenuOutput("userloggedin")),
  dashboardSidebar(width = 350,
                   div(style="overflow-y: scroll"),
            
                   tags$head(tags$style(HTML(".sidebar { height: 250vh; overflow-y: auto; }
                                             .shiny-notification{position: fixed;top: 33%;left: 45%;right: 30%;}
                                             " )
                   )),
                   sidebarMenu(
                     menuItem("Reference List", tabName = "dashboard", icon = icon("dashboard")),
                    uiOutput("reflist"),
                    fileInput('rdatafileupload', 'Upload Query RDS File'),
                    menuItem('Mapping', tabName = 'map', icon = icon('hand-o-right'),
                              menuSubItem("Preprocessing", tabName = "prep"),
                              menuSubItem("Feature Plots", tabName = "fplots")
                     )
                     
                   )#end of sidebar menu
  ),#end dashboardSidebar
  
  dashboardBody(
    # shinyDashboardThemes(
    #   theme = "blue_gradient"
    # ),
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "styles.css")
    ),
    useShinyjs(),
    tags$script("Shiny.addCustomMessageHandler('resetInputValue', function(variableName){
                Shiny.onInputChange(variableName, null);
                });
                "),
    tabItems(
      tabItem(tabName = "dashboard",
              box(
                width = 12,solidHeader = T,
                title = "Available Reference Datasets",
                DT::dataTableOutput("datasetTable")
              )
      ),
      
      ######################################################################################################################################
      
      tabItem(tabName = "prep",
              box(
                width = 12,solidHeader = TRUE,
                title = "Enter inputs for finding anchors and Map query",
                fluidRow(
                  column(4,selectInput("norm", "Normalization method",c('LogNormalize' = "LogNormalize",'SCT' = "SCT"),selected = "LogNormalize")),
                  column(4,uiOutput("refassay")),
                  column(4,uiOutput("queryassay"))
                ),
                fluidRow(
                  column(4,selectInput("recompres", "Recompute residuals",c('TRUE' = "TRUE",'FALSE' = "FALSE"),selected = "FALSE")),
                  column(4,uiOutput("ndim")),
                  column(4,uiOutput("ntrees"))
                ),
                fluidRow(
                  column(6,uiOutput("metadata")),
                  column(6,uiOutput("refred"))
                ),
                fluidRow(
                  column(6,actionButton(inputId = 'mapbutton',label = 'Map cells to reference')),
                  column(6,downloadButton('downloaddata', 'Download mapped RDS file'))
                )
              ),
              box(title = "Compare Dimension Reduction plots",solidHeader = TRUE,width=12,
                  
                  fluidRow(
                    column(6,uiOutput("tsnea2")),
                    column(6,uiOutput("tsneb2"))
                  ),
                  fluidRow(
                    column(6,uiOutput("reductiona")),
                    column(6,uiOutput("reductionb"))
                  ),
                  fluidRow(
                    column(6,checkboxInput("checklabel1", label = "Check for cell  group labelling", value = TRUE)),
                    column(6,checkboxInput("checklabel2", label = "Check for cell  group labelling", value = TRUE))
                  ),
                  sliderInput("pointa2", "Point Size:",min = 0, max = 5, value = 1,step=.25),
                  plotOutput("comptsne", height = 700)
              )
      ),
      ######################################################################################################################################
      
      tabItem(tabName = "fplots",
              box(
                width = 12,solidHeader = TRUE,title ="Feature Plots",
                radioButtons("fplotoption","Select one", c("Gene"="gene","Prediction score"="meta"),selected = "gene"),
                fluidRow(
                  column(6,uiOutput("featurelist")),
                  column(6,uiOutput("featuregroupby"))
                ),
                plotOutput("fplot", height = 700)
                ))
                
    )#end of tabitems
  )#end of dashboard body
)#end of dashboard page

