# Whither-Weather
Shiny app for processing/displaying weather data

Dennis Noren, individual non-sponsored project.

Data source is the NOAA meteorological API, accessed via rnoaa package.
Much preprocessing of weather data has already been done.
The geographic scope is that of the 48 contiguous U.S. states plus Washington DC.

The main output is dygraph plots showing "normal" temperatures for day of year
 as min, max, and midrange, with particular measured min-max values plotted
 on those. Also shown are the all-time record temperatures by day of year.
These are done relative to a set of several hundred cities. Two cities are shown top-to-bottom to allow comparison.
 
The user selects two cities from the list of cities, and a range of dates. The script calls the rnoaa function that finds the k nearest stations to both of those cities according to a search radius and the k stations (if fewer are available within radius, the script handles that).

The app also displays a map view which shows locations of all cities and weather stations in the dataset, and highlights the selected cities.

The script finds the normal temperatures for each day of the year, and does an inverse distance weighted interpolation among the stations found for that city. The set of stations and the normal temps have been preprocessed and filtered for those stations which have enough data for the normals. 

The script finds the particular highs, lows, and midranges for those stations. The same type of spatial interpolation is done among the stations associated with a city. There are however some individual days with missing values.

The script then plots the particulars and normals against each other usng the dygraphs package. Query and zoom capabilities can then be executed. A processing button is provided to allow the user to select new inputs and compute/redisplay only when requested. There are also options for displaying moving averages and a Butterworth high pass smoothing filter, as well as switching between Fahrenheit and Celsius units.

To fix/extend:

Determine why map tab view is failing on shinyapps.io, but works fine on local machine.

Change to an equal area map projection.

Fix the input data all-time record values that are obviously erroneous.

Build option for computing/displaying percentile values, using those as y axis, perhaps involving pre-processing. (long-term)

