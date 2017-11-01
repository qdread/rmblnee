# Read LAI from NEE plots.
library(XLConnect)
obsPAR <- read.csv('data/par_raw2015.csv', stringsAsFactors = FALSE)
source('calclai.r')

# Convert date and time to the correct datetime.
dt <- paste(as.character(as.Date(obsPAR$Date)),
            unlist(lapply(strsplit(as.character(obsPAR$Time), ' '), '[', 2)))
obsPAR$Datetime <- strptime(dt, '%Y-%m-%d %H:%M:%S', tz = 'America/Denver')
obsPAR$Cloud <- TRUE
obsPAR$Cloud[which(obsPAR$Sun == 'F')] <- FALSE
lats <- c(Almont = 38.715, Cinnamon = 38.992, Maxfield = 38.95, Washington_Gulch = 38.898, Almont_Low = 38.655, Judd_Falls = 38.971, Painter_Boy = 38.970, Elkton = 38.961, Almont_Obs = 38.715, Cinnamon_Obs = 38.992, Schofield = 39.034, Cinni_Top = 38.993, Brush_Creek = 38.868, Slate_River = 38.914)
obsPAR$Latitude <- lats[obsPAR$Site]
obsPAR$Below <- apply(obsPAR[,(5:9)], 1, mean)


obsPAR$LAI <- with(obsPAR, calcLAI(par_above = PAR.Above., par_below = Below, meas_day = Datetime$yday, meas_time = Datetime$hour + Datetime$min/60, latitude = Latitude, clouds = Cloud, leapyear = FALSE))