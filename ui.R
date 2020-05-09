# shiny ui script for "Whither-Weather" app

library(shiny)
library(shinyalert)
library(lubridate)
library(leaflet)
library(here)
library(shinybusy)
library(dygraphs)

endDate <- today() - 50
startDate <- as.character(endDate - 600)
endDate <- as.character(endDate)


# Define UI for application that draws a histogram
shinyUI(navbarPage("Whither Weather !!",

  tabPanel("", # user selects two cities for processing/display
    title = "Cities",
    tags$strong(verbatimTextOutput("citiesHelp")),
    tabPanel('Display length', DT::dataTableOutput('ex1'))),

  tabPanel("Temperatures",
    tags$head(
      tags$style(
        HTML(".shiny-notification {
              position:fixed;
              top: calc(40%);
              left: calc(25%);
              }
              "
            )
          )
        ),
    sidebarLayout(
      sidebarPanel(
        useShinyalert(),
        actionButton("goButton",h5(tags$b("Go")), width = 165), hr(),
        tags$strong(verbatimTextOutput("citiesSelected")), 
        dateRangeInput("dates", label = h5("Begin/End Date"),
          start = as.character(startDate),
          end = as.character(endDate), separator=":"),
        numericInput("searchRadius", "Search Radius (km)",
          value = 25, step = 1, min = 3, max = 50, width = NULL),
        numericInput("limit", "Max stations",
          value = 3, step = 1, min = 1, max = 10, width = NULL),
        radioButtons("metric", "Temperature Units",
                     c("Fahrenheit" = "fahr",
                       "Celsius" = "cels")),
        # widget for filter-based smoothing control
        checkboxInput("smooth", h5(tags$b("High Pass Smoothing Filter"))),
        sliderInput("filtfreq", "Filter Frequency Cutoff:", 6, 15, 10,
          step = NULL, round = FALSE, sep = ",", pre = "",
          post = "", ticks = TRUE, animate = FALSE, width = NULL),
        actionButton("stopButton",h5(tags$b("Close")),width=165),
        width=3),
    
      mainPanel(
        add_busy_spinner(spin = "fingerprint", color="#115546",
          position = "top-right", height='300px', width='300px'),
        tags$strong(verbatimTextOutput("quickInfo")),     
        dygraphOutput("dyg1"),
        br(),
        dygraphOutput("dyg2")
        )
      ) # sidebarLayout
    ), # tabPanel
  tabPanel("Map",
          leafletOutput("thismap"), p(),
          tags$strong(verbatimTextOutput("mapInfo")), p(),
          radioButtons("mapSelect", "Base Map",
              c("ESRI National Geographic" = "1",
                  "ESRI Topography" = "2",
                  "ESRI Mosaiced Imagery" = "3",
                  "ESRI Shaded Relief" = "4")
            )
          ),
#           actionButton("recalc", "Reload"),

  tabPanel("Help", htmlOutput("help")),
  tabPanel("Info", htmlOutput("info"))
  ) #navbarPage
) #shinyUI
