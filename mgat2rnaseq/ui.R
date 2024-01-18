#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

#library(shiny)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  
  theme = shinytheme("simplex"),
  
  # Application title
  titlePanel("Mogat2 Knockout Mice Gene Expression Under 3 Diets"),
  
  # fluidRow(
  #   column(3, uiOutput('year')),
  #   column(3, uiOutput('topic')),
  #   column(3, uiOutput('questions')),
  #   column(3, uiOutput('datatype'))
  #   
  # ),
  # 
  # 
  # 
  # fluidRow(
  #   column(4, uiOutput('useDiffButtons')),
  #   column(4, uiOutput('group1')),
  #   column(4, uiOutput('group2'))
  # ),
  # 
  
  tabsetPanel(type = 'tabs',
              
              
              tabPanel("1 Gene",
                       sidebarLayout(
                         sidebarPanel(width = 3,
                                      
                           textInput("gene", "Input Gene", value = "Mogat2", width = NULL, placeholder = NULL)
                         ),
                         mainPanel(
                           plotOutput("boxplots")
                         )
                       )
                       
              ),
              
              
              tabPanel("Table",
                       fluidRow(
                         dataTableOutput('table')
                       ),
              ),
              
              tabPanel("Heatmap",
                       fluidRow(
                         sidebarPanel(width = 3,
                                      selectInput("pathway", "select a pathway", choices = pathwaylist)
                         ),
                         mainPanel(
                           plotOutput("heatmap", width = "12in", height = "10in")
                         )
                       )
               ),
              tabPanel("Enrichment",
                       
                       
                       fluidRow(
                         column(12, style = "height: 800px;", imageOutput("gsea_chow")),
                         column(12, style = "height: 800px;", imageOutput("gsea_hf")),
                         column(12, style = "height: 800px;", imageOutput("gsea_lf")),
                       )
                       
               )
              
  ),
  
  fluidRow(
    textOutput('debug')
  )
))
