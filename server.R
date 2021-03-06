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
library(leaflet)
library(sp)
library(maps)
library(ggmap)


# Note: need NOAA API key.
# https://www.ncdc.noaa.gov/cdo-web/webservices/v2
# Placing it in .RProfile works well for individual use, as such:
# options(noaakey = "thekeywithabout32characters")
# Placing it in .Renviron works well for shinyapps.io, as such:
# noaakey = "thekeywithabout32characters"
readRenviron(file.path("~", ".Renviron"))
# print(Sys.getenv("noaakey"))

# 679 cities, all from CONUS
cities <- readRDS(here::here("inputData", "cities.rds"))
# 4592 stations, with TMIN and TMAX rows for each
stations <- readRDS(here::here("inputData", "stations.rds"))
# preprocessed normals for the station set, over 1 million rows,
#   which is 365 rows for each station
stations$label <- paste0(stations$name," (", stations$first_year, " - ",
                         stations$last_year, ")")
normals <- readRDS(here::here("inputData", "normals.rds"))

multNA <- function(x, y) { ifelse(is.na(x), NA, x * y)}

tCelsToFahr <- function(tenthCels) {tenthCels*1.8/10 + 32}
FahrToCels <- function(fahr) {(fahr-32) / 1.8}

shinyServer(function(input, output, session) {
  
  mapText <- paste0("Red markers are city centroids from dataset.","\n", 
    "Purple markers are weather stations from dataset.", "\n",
    "Use zoom control buttons, or your mouse roller.", "\n",
    "You may pan with mouse left button.", "\n",
    "Hover cursor on city to see its name.", "\n",
    "Cyan circles are 25 km and 50 km from selected cities.", "\n",
    "Yellow circles are 200 km from selected cities.")
  output$mapInfo <- renderText(mapText)
  output$citiesHelp <- renderText("Click to select two cities. Use page navigation and search box if desired. Then click 'Temperatures' tab.")
  quickText <- paste0("If you have selected two cities, click 'Go' button,","\n", 
    " or you may first change any of the input controls on the left,", "\n",
    " or look at the Map tab to see locations.", "\n",
    "There may be a 30 second or more delay to see temperature results.", "\n",
    "Results: Line charts for daily temperature min-max and daily record min-max", "\n",
    " along with the range between normal min and max as broad blue band.", "\n",
    "Each time you change any cities or controls, click 'Go' button again.", "\n",
    "For more detail, look at Help tab.")
  output$quickInfo <- renderText(quickText)
  
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
  citySelect <-  reactive( {
    cit <- input$ex1_rows_selected
    cities <- cities[rev(cit), ] # reversing so latest 2 selected are at top
    cities
  })

  # citySelectGo for getting city name in title only when temps reprocessed
  citySelectGo <-  eventReactive(input$goButton, {
    cit <- input$ex1_rows_selected
    cities <- cities[rev(cit), ] # reversing so latest 2 selected are at top
    cities
  })
  
  meanDate <- eventReactive(input$goButton, {
    as.character(mean(c(input$dates[1], input$dates[2])))})

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
  }) # end getStations eventReactive function

  output$thismap <- renderLeaflet({
    km3 <-3000; km5 <-5000; km10 <-10000; km25 <-25000; km50 <-50000; km200 <-200000
    opaq <- 1.0; trluc <- 0.7
    clr25 <- "cyan"; clr50 <- "cyan"; clr200 <- "yellow"

    us_states <- map_data("state") # map_data in package ggmap
#    https://stackoverflow.com/questions/45237646/r-leaflet-addpolygons-by-group
#    above link gave solution below
    split_data_poly = lapply(unique(us_states$group), function(x) {
      df = as.matrix(us_states[us_states$group == x, c("long", "lat")])
      polys = Polygons(list(Polygon(df)), ID = x)
      return(polys)
    })

    data_polys = SpatialPolygons(split_data_poly)

    cityPair <- citySelect()
    if (nrow(cityPair) < 2) {
      shinyalert(title = "Must select two cities. Return to Cities tab",
                 type = "error", closeOnClickOutside = TRUE)}
    baseMap <- switch(input$mapSelect,
        "1" = "Esri.NatGeoWorldMap",
        "2" = "Esri.WorldTopoMap",
        "3" = "Esri.WorldImagery",
        "4" = "Esri.WorldShadedRelief")
    leaflet(data_polys) %>%
      setView(lng = -101, lat = 40, zoom=4) %>%
      addProviderTiles(baseMap,
        providerTileOptions(noWrap = TRUE, updateWhenIdle = FALSE)) %>%
      addCircles(lng = cityPair$longitude, lat= cityPair$latitude,
        label = cityPair$id,
        labelOptions = labelOptions(noHide = T, textsize="17px")) %>%
      addCircles(lng = stations$longitude, lat = stations$latitude,
        opacity = opaq, radius=km3, color="purple", weight=2, fill=FALSE) %>%
      addCircles(lng=cities$longitude, lat=cities$latitude,
        label= cities$id, weight = 0, radius = km5, opacity = 0, color="white",
        labelOptions = labelOptions(noHide = F, textsize = "17px")) %>%
      addCircles(lng=cities$longitude, lat=cities$latitude,
        opacity = trluc, radius=km3, color="red", weight=3, fill= FALSE) %>%

      addCircles(lng = cityPair$longitude[1], lat=cityPair$latitude[1],
        weight= 3, color = clr25, radius = km25, opacity = trluc, fill = FALSE) %>%
      addCircles(lng = cityPair$longitude[1], lat=cityPair$latitude[1],
        weight= 3, color = clr50, radius = km50, opacity = opaq, fill = FALSE) %>%
      addCircles(lng = cityPair$longitude[1], lat=cityPair$latitude[1],
        weight= 3, color = clr200, radius = km200, opacity = opaq, fill = FALSE) %>%

      addCircles(lng = cityPair$longitude[2], lat=cityPair$latitude[2],
        weight= 3, color = clr25, radius = km25, opacity = trluc, fill = FALSE) %>%
      addCircles(lng = cityPair$longitude[2], lat=cityPair$latitude[2],
        weight= 3, color = clr50, radius = km50, opacity = opaq, fill = FALSE) %>%
      addCircles(lng = cityPair$longitude[2], lat=cityPair$latitude[2],
        weight= 3, color = clr200, radius = km200, opacity = opaq, fill = FALSE) %>%

      addPolygons(weight = 1, opacity = 0.8, color = "black", fill = FALSE)
    })

  getActuals <- eventReactive(input$goButton, {
    stationz <- getStations()
    if (input$dates[1] > input$dates[2]) {
      shinyalert(title = "Begin date must precede end date",
                 type = "error", closeOnClickOutside = TRUE)}
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
          if (nrow(actual) < 2) {
            shinyalert(title = "No data returned. Try changing dates, or search radius",
                       type = "error", closeOnClickOutside = TRUE)}
          
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
  }) # end getActuals eventReactive function
  
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
      distance = gs1$distance, invdsq = 1 / (gs1$distance^2))
    f1 <- filter(normals, station %in% gs1id)
    f1 <- mutate(f1, city = citySeq[1,1], cityNum = 1)
    
    gs2 <- as.data.frame(gs[[1]][2])
    gs2 <- mutate(gs2, city = citySeq[2,1], type = "norm")
    names(gs2) <- c("station","name","latitude","longitude",
                    "distance", "city", "type")
    gs2id <- as.vector(gs2$station)
    gs2LUT <- bind_cols(city = gs2$city, station = gs2$station,
      distance = gs2$distance, invdsq = 1 / (gs2$distance^2))
    f2 <- filter(normals, station %in% gs2id)
    f2 <- mutate(f2, city = citySeq[2,1], cityNum = 2)
    
    normalz <- bind_rows(f1, f2)

    gsLUT <- bind_rows(gs1LUT, gs2LUT)
    normalz2 <- left_join(normalz, gsLUT, by = c("city", "station"))
    normalz2 <- arrange(normalz2, city, distance, date)
    normalz2 <- transmute(normalz2, city, date, station, invdsq,
      nMaxWt = NormalMax*invdsq, nMidWt = NormalMidrange*invdsq,
      nMinWt = NormalMin*invdsq, cityNum, recordMax, recordMin)
    normalz2 <- group_by(normalz2, city, date) %>%
      summarise(NormMax = sum(nMaxWt), NormMid = sum(nMidWt),
                NormMin = sum(nMinWt), denom = sum(invdsq),
                cityNum = first(cityNum), recMax = first(recordMax),
                recMin = first(recordMin))
    normalz2 <- mutate(.data = normalz2, normalMax = NormMax/denom,
              normalMid = NormMid/denom, normalMin = NormMin/denom,
              type = "norm")
    actuals <- mutate(actuals, station = paste0("GHCND:", station),
          city = ifelse(cityNum == 1, citySeq[1,1], citySeq[2,1]))
    actuals2 <- left_join(actuals, gsLUT, by = c("city", "station"))
    actuals2 <- transmute(actuals2, city, date, station, invdsq,
      nMaxWt = multNA(tmax,invdsq), nMidWt = multNA(tmid,invdsq),
      nMinWt = multNA(tmin,invdsq), cityNum)
    actuals2 <- group_by(actuals2, city, date) %>%
      summarise(actMax = sum(nMaxWt, na.rm = TRUE),
                actMid = sum(nMidWt, na.rm = TRUE),
                actMin = sum(nMinWt, na.rm = TRUE),
                denomMax = sum(ifelse(is.na(nMaxWt),0,invdsq)),
                denomMid = sum(ifelse(is.na(nMidWt),0,invdsq)),
                denomMin = sum(ifelse(is.na(nMinWt),0,invdsq)),
                cityNum = first(cityNum))
    actuals2 <- mutate(.data = actuals2,
      actualMax = actMax/denomMax,
      actualMid = actMid/denomMid,
      actualMin = actMin/denomMin,
      type = "actual")
    
    varsAct <- c("city", "date", "actualMax", "actualMid",
            "actualMin", "cityNum", "type")
    varsNorm <- c("city", "date", "normalMax", "normalMid",
          "normalMin", "recMax", "recMin", "cityNum", "type")
    list(normalz2[, varsNorm], actuals2[, varsAct])
#    interpolations <- bind_rows(normalz2[,vars], actuals2[,vars])
  }) # end reactive funtion 'interpolator'
  
  assembler <- reactive({
    # returns a time series with normal min/max/midrange for station(s)
    #  by convention, NOAA uses dates 2010-01-01 through 2010-12-31
    #   for normals retrieval, without leap days, but this app's
    #   pre-processing interpolated Feb 29.
    #   for 
    # function also returns time series with actual daily min/max/mid
    # it merges normals into the actuals for each year

    interpolations <- interpolator()
#    normalz <- filter(interpolations, type == "norm")
    normalz <- interpolations[[1]]
    normalz$month <- month(normalz$date)
    normalz$day <- day(normalz$date)
#    actuals <- filter(interpolations, type == "actual")
    actualz <- interpolations[[2]]
    actualz$month <- month(actualz$date)
    actualz$day <- day(actualz$date)
    LUTvars <- c("city","normalMax","normalMid","normalMin",
              "recMax", "recMin", "month","day")
    normalz_LUT <- normalz[, LUTvars]
    actNormal <- left_join(actualz, normalz_LUT, 
                           by=c("city","month","day"))
    names(actNormal) <- c("city", "date", "MaxT", "MidT", "MinT",
      "cityNum", "type", "month", "day",
      "MaxNorm", "MidNorm", "MinNorm", "RecordMax", "RecordMin") 
    
    cityList <- citySelect()
    city1 <- filter(actNormal, city == cityList[1,1])
    city2 <- filter(actNormal, city == cityList[2,1])
    list(city1, city2)
  }) #end reactive function 'assembler'

  smoother <- eventReactive(input$goButton, {
  city1 <- assembler()[[1]]
  output$textInfo1 = renderText(as.character(city1[1,]))
  city2 <- assembler()[[2]]
 if (input$smooth == FALSE) # no smoothing applied
   {city1$DailyMax <- city1$MaxT
      city1$DailyMid <- city1$MidT
      city1$DailyMin <- city1$MinT
      city2$DailyMax <- city2$MaxT
      city2$DailyMid <- city2$MidT
      city2$DailyMin <- city2$MinT
        }
    else
    {bwf1 <- bwfilter(city1$MaxT, # Butterworth filter
                  freq = input$filtfreq, nfix = 2)
      city1$DailyMax <- bwf1$trend
      bwf1 <- bwfilter(city1$MinT, freq = input$filtfreq, nfix = 2)
      city1$DailyMin <- bwf1$trend
      bwf2 <- bwfilter(city2$MaxT, freq = input$filtfreq, nfix = 2)
      city2$DailyMax <- bwf2$trend
      bwf2 <- bwfilter(city2$MinT, freq = input$filtfreq, nfix = 2)
      city2$DailyMin <- bwf2$trend
    }  
  varsID <- c("city", "date")
  varsMeas <- c("DailyMax","DailyMin","MaxNorm","MidNorm","MinNorm",
                "RecordMax", "RecordMin")
  vars <- c(varsID, varsMeas)
  city1 <- select(.data = city1, vars)
  city2 <- select(.data = city2, vars)
  if (input$metric == "cels") {
    city1 <- city1 %>%
      mutate_at(varsMeas, FahrToCels)
    city2 <- city2 %>%
      mutate_at(varsMeas, FahrToCels)
  }
  cit1 <- data.frame(sapply(city1[,3:9], function(x) round(x, 1)))
  names(cit1) <- varsMeas
  city1 <- bind_cols(city1[,1:2], cit1)
  cit2 <- data.frame(sapply(city2[,3:9], function(x) round(x, 1)))
  names(cit2) <- varsMeas
  city2 <- bind_cols(city2[,1:2], cit2)
  xt1 <- xts(as.matrix(city1[,3:9]), order.by = city1$date)
  xt2 <- xts(as.matrix(city2[,3:9]), order.by = city2$date)
  xt1 <- na.approx(xt1)
  xt2 <- na.approx(xt2)
  list(xt1, xt2)

 }) # end eventReactive function 'smoother'
  
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
    cityList <- citySelectGo()
    title <- paste0(cityList[1,1],", ", cityList[1,4])
    gs <- getStations()
    city1 <- smoother()[[1]]
    msg <- paste("Number of stations = ", (gs[[2]][1]))
    dist <- paste(round((data.frame(gs[[1]][1]))[,5],1), collapse=", ")
    infoText <- paste("Move cursor within plot to see data values.", "\n",
      "Drag left to right to zoom dates.", "\n",
      "Drag top to bottom to zoom values.", "\n",
      "Double click to reset.")
    showNotification(infoText, duration=15, type="message")
    dygraph(city1, main = title, group= "temperatures") %>%
      dyAxis("y", yLabel()) %>%
      dyLegend(width = 700) %>%
      dyGroup(c("DailyMax", "DailyMin"),
        color=c("red", "blue"), strokeWidth=1.8) %>%
      dyGroup("RecordMax", color="orange", strokeWidth=1.5) %>%
      dyGroup("RecordMin", color="purple", strokeWidth=1.5) %>%
      dyAnnotation(x= meanDate(), text=msg,
        tooltip = paste("Distances from city centroid, in km =", dist),
        attachAtBottom=TRUE, width=200, height=22) %>%
      dyOptions(drawGrid = TRUE) %>%
      dyRoller(rollPeriod = 1) %>%
      dySeries(c("MaxNorm","MidNorm","MinNorm"),color="blue",strokeWidth=0)
  })

    # temperature plot for bottom
  output$dyg2 <- renderDygraph({
    if (input$goButton == 0) return()
    cityList <- citySelectGo()
    title <- paste0(cityList[2,1],", ", cityList[2,4])
    gs <- getStations()
    city2 <- smoother()[[2]]
    msg <- paste("Number of stations = ", (gs[[2]][2]))
    dist <- paste(round((data.frame(gs[[1]][2]))[,5],1), collapse=", ")
    dygraph(city2, main = title, group= "temperatures") %>%
      dyAxis("y", yLabel()) %>%
      dyLegend(width = 700) %>%
      dyGroup(c("DailyMax","DailyMin"),
        color=c("red","blue"), strokeWidth=1.8) %>%      
      dyGroup("RecordMax", color="orange", strokeWidth=1.5) %>%
      dyGroup("RecordMin", color="purple", strokeWidth=1.5) %>%
      dyAnnotation(x= meanDate(),text=msg,
        tooltip = paste("Distances from city centroid, in km =", dist),
        attachAtBottom=TRUE, width=200, height=22) %>%
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
