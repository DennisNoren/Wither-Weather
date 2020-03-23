# shiny server script for "Whither-Weather" app

library(shiny)
library(rnoaa)
library(dygraphs)
library(tidyverse)
library(fs)
library(lubridate)
library(here)
library(xts)
library(stringr)
library(mFilter)
library(shinyalert)

readRenviron(file.path("~", ".Renviron"))
#print(Sys.getenv("noaakey"))
# Note: need NOAA API key.
# https://www.ncdc.noaa.gov/cdo-web/webservices/v2
# Placing it in .RProfile works well for individual use, as such:
# options(noaakey = "thekeywithabout32characters")
# Placing it in .Renviron works well for shinyapps.io, as such:
# noaakey = "thekeywithabout32characters"

# 679 cities, all from CONUS
cities <- readRDS(here::here("inputData", "cities.rds"))
# 4592 stations, with TMIN and TMAX rows for each
stations <- readRDS(here::here("inputData", "stations.rds"))
# preprocessed normals for the station set, over 1 million rows,
#   which is 365 rows for each station
cumNormSS <- readRDS(here::here("inputData", "cumNormSS.rds"))

multNA <- function(x, y) { ifelse(is.na(x), NA, x * y)}

tCelsToFahr <- function(tenthCels) {tenthCels*1.8/10 + 32}
FahrToCels <- function(fahr) {(fahr-32) / 1.8}

shinyServer(function(input, output, session) {
  
  quickText <- paste0("If you have selected two cities, click 'Go' button,",
    "\n", " or you may first change any of the input controls on the left.", "\n", "There may be a 30 second or more delay to see results.", "\n",
    "Each time you change any cities or controls, click 'Go' button again.", "\n","For more detail, look at Help tab.")
  output$quickInfo <- renderText(quickText)
  output$citiesHelp <- renderText("Click to select two cities. Use page navigation and search box if desired. Then click 'Temperatures' tab.")
  
  observe({
    cit <- input$ex1_rows_selected
    cities <- cities[rev(cit), ] # reversing so latest 2
    city1 <- paste0(cities[1,1],", ", cities[1,5])
    city2 <- paste0(cities[2,1],", ", cities[2,5])
    cityText <- paste0("Cities Selected:","\n",
                       city1, "\n", city2, "\n")
    output$citiesSelected <- renderText(cityText)
  })
                                                  
    # citySelect returns a DF with the user-selected cities in Cities tab
  citySelect <-  eventReactive(input$goButton, {
    cit <- input$ex1_rows_selected
    cities <- cities[rev(cit), ] # reversing so latest 2 selected are at top
    cities
   
  })


  meanDate <- eventReactive(input$goButton, {
    as.character(mean(c(input$dates[1], input$dates[2])))
                                    })

  # getStations returns a list of DFs, each with up to k stations,
  #   where k is the maximum number to return.
  getStations <- eventReactive(input$goButton, {
    cityList <- citySelect()
    # 'top_n' grabs the bottom 2 from the list, 
    #   which are the most recently selected cities
    cityList <- top_n(cityList, 2)
#    if (is.na(cityList$id[1] || is.na(cityList$id[2]))) {
    if (nrow(cityList) < 2) {
        shinyalert(title = "Must select two cities. Return to Cities tab",
                 type = "error", closeOnClickOutside = TRUE)
    }
     
    nearStations <- meteo_nearby_stations(cityList, station_data = stations,
      var = c("TMAX", "TMIN"), lat_colname = "latitude", 
      lon_colname = "longitude", radius = input$searchRadius,
      limit = input$limit)
    if (is.na(nearStations[[1]][1]) ||
      is.na(nearStations[[2]][1])) {
      shinyalert(title = "No stations found for one or both cities; try increasing search radius", type = "error", closeOnClickOutside = TRUE)
    }

    # function meteo_nearby_stations produces list that is 
    #    alphbetically sorted by id name
    # want to keep cityList order the same, so if the two city names
    #   there are flipped alphabetically, do a swap of the list items
    #   in the nearStations list
    if (cityList$id[2] < cityList$id[1])
      nearStations <- list(nearStations[[2]], nearStations[[1]])
    numStations <- lapply(nearStations, nrow)
    nearDF <- nearStations[[1]]

    listSta <- list(nearStations, numStations)
  }) # end getStations

  getActuals <- eventReactive(input$goButton, {
    stationz <- getStations()
    dateRange <- seq.Date(from = input$dates[1],
                          to = input$dates[2], by = "days")
    # remove leap days because NOAA normals data do not include them
    dateRange <- dateRange[!(month(dateRange)==2 & day(dateRange)==29)]
    dateRange <- as.Date(dateRange)

    numDays <- length(dateRange)
    dateDF <- data.frame(dateRange = dateRange, seq <- seq(numDays))
    names(dateDF) <- c("dateRange", "seq")
    actuals <- tibble(date = as.Date(character()), station = character(),
      tmax = double(), tmid=double(), tmin=double(), city = character(),
      cityNum = integer())
    citySeq <- citySelect()
    begDate <- input$dates[1]
    endDate <- input$dates[2]
    for (citIndex in 1:2) {
      numSta <- as.integer(stationz[[2]][citIndex])
      thisCity <- as.data.frame(stationz[[1]][citIndex])
      names(thisCity) <- c("id","name","latitude",
                           "longitude","distance")
      for (station in 1:numSta) {
          staID <- str_replace(thisCity[station,1],
                  pattern = "GHCND:", replacement = "")
          actual <- meteo_tidy_ghcnd(stationid = staID,
            var=c("TMAX","TMIN"),
            date_min = begDate, date_max = endDate)

          actual$tmin <- tCelsToFahr(actual$tmin)
          actual$tmax <- tCelsToFahr(actual$tmax)
          actual$tmid <- (actual$tmax + actual$tmin) / 2
          actual$cityNum <- citIndex
          actual$city <- citySeq[citIndex, 1]
          actual <- select(.data = actual, date, station = id,
                           tmax, tmid, tmin, city, cityNum)
          actuals <- bind_rows(actuals, actual)
      } # end station loop
    } # end city loop
    return(actuals)
  }) # end getActuals reactive function
  
  # this interpolates temps among the nearby stations, for each date
  #   using inverse distance weighting, order 2
  interpolator <- reactive({ 
    actuals <- getActuals()
    gs <- getStations()
    citySeq <- citySelect()

    gs1 <- data.frame(gs[[1]][1])
    gs1 <- mutate(gs1, city = as.character(citySeq[1,1]), type = "norm")

    names(gs1) <- c("station","name","latitude","longitude",
                    "distance", "city", "type")
    gs1id <- as.vector(gs1$station)
    gs1LUT <- bind_cols(city = gs1$city, station = gs1$station,
                          invdsq = 1 / (gs1$distance^2))
    f1 <- filter(cumNormSS, station %in% gs1id)
    f1 <- mutate(f1, city = citySeq[1,1], cityNum = 1)
    
    gs2 <- as.data.frame(gs[[1]][2])
    gs2 <- mutate(gs2, city = citySeq[2,1], type = "norm")
    names(gs2) <- c("station","name","latitude","longitude",
                    "distance", "city", "type")
    gs2id <- as.vector(gs2$station)
    gs2LUT <- bind_cols(city = gs2$city, station = gs2$station,
                          invdsq = 1 / (gs2$distance^2))
    f2 <- filter(cumNormSS, station %in% gs2id)
    f2 <- mutate(f2, city = citySeq[2,1], cityNum = 2)
    
    normals <- bind_rows(f1, f2)
    gsLUT <- bind_rows(gs1LUT, gs2LUT)
    normals2 <- left_join(normals, gsLUT, by = c("city", "station"))
    normals2 <- transmute(normals2, city, date, station, invdsq,
      nMaxWt = NormalMax*invdsq, nMidWt = NormalMidrange*invdsq,
      nMinWt = NormalMin*invdsq, cityNum)
    normals2 <- group_by(normals2, city, date) %>%
      summarise(NormMax = sum(nMaxWt), NormMid = sum(nMidWt),
                NormMin = sum(nMinWt), denom = sum(invdsq),
                cityNum = first(cityNum))
    normals2 <- mutate(.data = normals2, NormalMax = NormMax/denom,
              NormalMid = NormMid/denom, NormalMin = NormMin/denom,
              type = "norm")
    actuals <- mutate(actuals, station = paste0("GHCND:", station),
          city = ifelse(cityNum == 1, citySeq[1,1], citySeq[2,1]))
    actuals2 <- left_join(actuals, gsLUT, by = c("city", "station"))
    actuals2 <- transmute(actuals2, city, date, station, invdsq,
      nMaxWt = multNA(tmax,invdsq), nMidWt = multNA(tmid,invdsq),
      nMinWt = multNA(tmin,invdsq), cityNum)
    actuals2 <- group_by(actuals2, city, date) %>%
      summarise(NormMax = sum(nMaxWt, na.rm = TRUE),
                NormMid = sum(nMidWt, na.rm = TRUE),
                NormMin = sum(nMinWt, na.rm = TRUE),
                denomMax = sum(ifelse(is.na(nMaxWt),0,invdsq)),
                denomMid = sum(ifelse(is.na(nMidWt),0,invdsq)),
                denomMin = sum(ifelse(is.na(nMinWt),0,invdsq)),
                cityNum = first(cityNum))
    actuals2 <- mutate(.data = actuals2,
      NormalMax = NormMax/denomMax,
      NormalMid = NormMid/denomMid,
      NormalMin = NormMin/denomMin,
      type = "actual")
    
    vars <- c("city", "date", "NormalMax", "NormalMid",
            "NormalMin", "cityNum", "type")
    interpolations <- bind_rows(normals2[,vars], actuals2[,vars])
  })
  
  assembler <- reactive({
    # returns a time series with normal min/max/midrange for station(s)
    #  by convention, NOAA uses dates 2010-01-01 through 2010-12-31
    # also returns time series with actual daily min/max/mid
    # it merges normals into the actuals for each year

    interpolations <- interpolator()
    normals <- filter(interpolations, type == "norm")
    normals$month <- month(normals$date)
    normals$day <- day(normals$date)
    actuals <- filter(interpolations, type == "actual")
    actuals$month <- month(actuals$date)
    actuals$day <- day(actuals$date)
    LUTvars <- c("city","NormalMax","NormalMid","NormalMin","month","day")
    normals_LUT <- normals[, LUTvars]
    actNormal <- left_join(actuals,normals_LUT, by=c("city","month","day"))
    names(actNormal) <- c("city", "date", "MaxT", "MidT", "MinT",
      "cityNum", "type", "month", "day", "MaxNorm", "MidNorm", "MinNorm") 
    
    cityList <- citySelect()
    city1 <- filter(actNormal, city == cityList[1,1])
    city2 <- filter(actNormal, city == cityList[2,1])
    list(city1, city2)
  })

 smoother <- eventReactive(input$goButton, {
  city1 <- assembler()[[1]]
  output$textInfo1 = renderText(as.character(city1[1,]))
  city2 <- assembler()[[2]]
 if (input$smooth == FALSE) # no smoothing applied
   {city1$MaxTemp <- city1$MaxT
      city1$MidTemp <- city1$MidT
      city1$MinTemp <- city1$MinT
      city2$MaxTemp <- city2$MaxT
      city2$MidTemp <- city2$MidT
      city2$MinTemp <- city2$MinT
        }
    else
    {bwf1 <- bwfilter(city1$MaxT, # Butterworth filter
                  freq = input$filtfreq, nfix = 2)
    city1$MaxTemp <- bwf1$trend
    bwf1 <- bwfilter(city1$MinT,
                  freq = input$filtfreq, nfix = 2)
    city1$MinTemp <- bwf1$trend
    bwf2 <- bwfilter(city2$MaxT,
                    freq = input$filtfreq, nfix = 2)
            city2$MaxTemp <- bwf2$trend
            bwf2 <- bwfilter(city2$MinT,
                    freq = input$filtfreq, nfix = 2)
            city2$MinTemp <- bwf2$trend
            }  
  varsID <- c("city", "date")
  varsMeas <- c("MaxTemp","MinTemp","MaxNorm","MidNorm","MinNorm")
  vars <- c(varsID, varsMeas)
  city1 <- select(.data = city1, vars)
  city2 <- select(.data = city2, vars)
  if (input$metric == "cels") {
    city1 <- city1 %>%
      mutate_at(varsMeas, FahrToCels)
    city2 <- city2 %>%
      mutate_at(varsMeas, FahrToCels)
  }
  xt1 <- xts(as.matrix(city1[,3:7]), order.by = city1$date)
  xt2 <- xts(as.matrix(city2[,3:7]), order.by = city2$date)
  xt1 <- na.approx(xt1)
  xt2 <- na.approx(xt2)
  list(xt1, xt2)

 })
  
  # build Cities tab with filter/sort/select capability
  output$ex1 <- DT::renderDataTable(DT::datatable(cities,
                options = list(pageLength = 15)))

  yLabel <- eventReactive(input$goButton, {
    ifelse (input$metric == "fahr", "Temperature Fahrenheit",
      "Temperature Celsius")
  })
  
  # temperature plot for top
  output$dyg1 <- renderDygraph({
    if (input$goButton == 0) return()
    output$quickInfo <- renderText("")
    cityList <- citySelect()
    title <- paste0(cityList[1,1],", ", cityList[1,4])
    gs <- getStations()
    city1 <- smoother()[[1]]
    msg <- paste("Number of stations = ", (gs[[2]][1]))
    dist <- paste(round((data.frame(gs[[1]][1]))[,5],1), collapse=", ")
    dygraph(city1, main = title, group= "temperatures") %>%
      dyAxis("y", yLabel()) %>%
      dyLegend(width = 500) %>%
      dyGroup(c("MaxTemp", "MinTemp"),
        color=c("red", "blue"), strokeWidth=1.5) %>%      
      dyAnnotation(x= meanDate(), text=msg,
        tooltip = paste("Distances in km =", dist),
        attachAtBottom=TRUE, width=180, height=22) %>%
      dyOptions(drawGrid = TRUE) %>%
      dyRoller(rollPeriod = 1) %>%
      dySeries(c("MaxNorm","MidNorm","MinNorm"),color="blue",strokeWidth=0)
  })

    # temperature plot for bottom
    output$dyg2 <- renderDygraph({
    if (input$goButton == 0) return()
    cityList <- citySelect()
    title <- paste0(cityList[2,1],", ", cityList[2,4])
    gs <- getStations()
    city2 <- smoother()[[2]]
    msg <- paste("Number of stations = ", (gs[[2]][2]))
    dist <- paste(round((data.frame(gs[[1]][2]))[,5],1), collapse=", ")
    dygraph(city2, main = title, group= "temperatures") %>%
      dyAxis("y", yLabel()) %>%
      dyLegend(width = 500) %>%
      dyGroup(c("MaxTemp","MinTemp"),
        color=c("red","blue"), strokeWidth=1.5) %>%      
      dyAnnotation(x= meanDate(),text=msg,
        tooltip = paste("Distances in km =", dist),
        attachAtBottom=TRUE, width=180, height=22) %>%
      dyOptions(drawGrid = TRUE) %>%
      dyRoller(rollPeriod = 1) %>%
      dySeries(c("MaxNorm","MidNorm","MinNorm"),color="blue",strokeWidth=0)
    })
    output$help <- renderText({
      path <- "inputData/WitherWeatherHelp.html"
      includeText(path)})
    output$info <- renderText({
      path <- "inputData/WitherWeatherInfo.html"
      includeText(path)})    
    session$onSessionEnded(function() { citySelect$suspend() })
    observe({ if (input$stopButton > 0) stopApp() })
  
})
