# shiny ui script for "Whither-Weather" app

library(shiny)
library(shinyalert)
library(here)
library(dygraphs)

# Define UI for application that draws a histogram
shinyUI(navbarPage("Whither Weather !!",

  tabPanel("", # user selects two cities for processing/display
    title = "Cities",
    tags$strong(verbatimTextOutput("citiesHelp")),
    tabPanel('Display length', DT::dataTableOutput('ex1'))),

  tabPanel("Temperatures",

    sidebarLayout(
      sidebarPanel(
        useShinyalert(),
        actionButton("goButton",h5(tags$b("Go")), width = 165), hr(),
        tags$strong(verbatimTextOutput("citiesSelected")), 
         dateRangeInput("dates", label = h5("Begin/End Date"),
          start = "2013-05-01", end = "2014-12-31", separator=":"),
        numericInput("searchRadius", "Search Radius (km)",
          value = 25, step = 1, min = 3, max = 50, width = NULL),
        numericInput("limit", "Max stations",
          value = 3, step = 1, min = 1, max = 10, width = NULL),
        radioButtons("metric", "Temperature Units",
                     c("Fahrenheit" = "fahr",
                       "Celsius" = "cels")),
        # widget for filter-based smoothing control
        checkboxInput("smooth", h5(tags$b("High Pass Filter"))),
        sliderInput("filtfreq", "Filter frequency:", 6, 15, 10,
          step = NULL, round = FALSE, sep = ",", pre = "",
          post = "", ticks = TRUE, animate = FALSE, width = NULL),
        actionButton("stopButton",h5(tags$b("Close")),width=165),
#        hr(),
        # htmlOutput("messages", inline = FALSE),
        # textOutput("textInfo2", inline = FALSE),
        width=3),
    
      mainPanel(
        tags$strong(verbatimTextOutput("quickInfo")),     
        dygraphOutput("dyg1"),
        br(),
        dygraphOutput("dyg2")
        )
      ) # sidebarLayout
    ), # tabPanel
  tabPanel("Help",
           htmlOutput("help")),
  tabPanel("Info",
           htmlOutput("info"))
  ) #navbarPage
) #shinyUI
