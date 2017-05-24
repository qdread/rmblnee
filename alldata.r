# Cleaned-up script for processing all the raw data from the observational gradient into a single
# combined dataframe with one row per site and all the data columns we need.

# QDR - Day 3 Mazatl (deer) in Month 1 Coatl (snake) in Year 3 Acatl (reed), aka September 3, 2015

#############################
# Edited 24 May 2017: put in file paths so that this can be run anywhere.
# Also converted xlsx files into csv files.


library(XLConnect)
library(plyr)

#########################################################################################





# Load PAR data and calculate LAI from PAR. 

obsPAR <- read.csv('data/par2015.csv')
str(obsPAR)
# Convert date and time to the correct datetime.
#dt <- paste(as.character(as.Date(obsPAR$Date, "%m/%d/%y")),
 #           unlist(lapply(strsplit(obsPAR$Time, '[:]'))))

obsPAR$Datetime <- strptime(paste(as.character(as.Date(obsPAR$Date, "%m/%d/%y")), as.POSIXct(obsPAR$Time,format="%H:%M:%S")),"%m-%d-%y %H:%M:%S", tz = 'America/Denver')

names(obsPAR)

#obsPAR$Datetime <- strptime(dt, '%Y-%m-%d %H:%M', tz = 'America/Denver')
obsPAR$Cloud <- TRUE
obsPAR$Cloud[which(obsPAR$Sun == 'F')] <- FALSE
lats <- c(Almont = 38.715, Cinnamon = 38.992, Maxfield = 38.95, Washington_Gulch = 38.898, Almont_Low = 38.655, Judd_Falls = 38.971, Painter_Boy = 38.970, Elkton = 38.961, Almont_Obs = 38.715, Cinnamon_Obs = 38.992, Schofield = 39.034, Cinni_Top = 38.993, Brush_Creek = 38.868, Slate_River = 38.914)
obsPAR$Latitude <- lats[obsPAR$Site]
obsPAR$Below <- apply(obsPAR[,(5:9)], 1, mean)
obsPAR$Par.above
#######       LAI CALCULATION FUNCTION
calcLAI <- function(par_above, par_below, meas_day, meas_time, latitude, clouds = FALSE, leapyear = FALSE) {
  
  # Calculates leaf area index given PAR above and below canopy, time and place, and sunlight conditions
  # Written by QDR on 13 June 2015
  # Sources: Decagon LP-80 manual; Decagon appendix re: beam fractionation; Wikipedia (for zenith angle calc)
  
  daysSinceEquinox <- function(yday, leapyear = FALSE) {
    yr <- if (leapyear) 366 else 365
    equinox <- yr - 285
    res <- yday - equinox + yr
    res[which(yday >= equinox)] <- yday - equinox
    #if (yday >= equinox) yday - equinox else yday - equinox + yr
  }
  
  latitude <- latitude * 2 * pi / 360 # Convert latitude to radians.
  r <- rep(0.82, length(par_above)) # Two possible values for fraction of potential PAR reaching the sensor.
  r[clouds] <- 0.2
  Fb = 1.395 + r * (-14.43 + r * (48.57 + r * (-59.024 + r * 24.835))) # Use r to calculate the beam fractionation.
  # Assume longitude is in the center of the time zone to simplify calculation of hour angle.
  hourangle <- 2 * pi * (meas_time - 12) / 24 # Hour angle of sun in radians.
  declination <- asin( sin(23.44 * 2 * pi / 360) * sin(2 * pi * (daysSinceEquinox(meas_day, leapyear) / 365.24) ) ) # Declination of sun
  theta_s <- acos( sin(latitude) * sin(declination) + cos(latitude) * cos(declination) * cos(hourangle) ) # Solar zenith angle
  A <- 0.283 + 0.785 * 0.9 + 0.159 * 0.9^2 # Constant in LAI equation
  tau <- par_below / par_above # Ratio of light passed through canopy to light above canopy
  K <- 1 / (2 * cos(theta_s) )
  LAI <- ( ( (1 - .5/K) * Fb - 1 ) * log(tau) ) / ( A * (1 - 0.47 * Fb) )
  return(LAI)
  
}



obsPAR$LAI <- with(obsPAR, calcLAI(par_above = PAR.Above., par_below = Below, meas_day = Datetime$yday, meas_time = Datetime$hour + Datetime$min/60, latitude = Latitude, clouds = Cloud, leapyear = FALSE))

LAImean <- with(obsPAR, tapply(LAI, Site, mean))
LAI95 <- with(obsPAR, tapply(LAI, Site, function(x) 1.96*sd(x)/length(x)))


#####################################################################################

# Load the NEE raw data and calculate the standardized NEE from it.

# Code to read NEE file.

read_nee <- function(filename, start = NA, end = NA, height = 1.25, length = 1.25, direc = 'data/li7500/') {
  require(ggplot2)
  needata <- read.delim(paste0(direc, filename), skip = 8, stringsAsFactors = FALSE) # Reads file.
  date_time <- read.table(paste0(direc, filename), nrow = 1, stringsAsFactors = FALSE) # Reads date and time.
  if (is.na(start)) { start <- 1 }
  if (is.na(end)) { end <- nrow(needata) }
  needata <- needata[start:end, ] # Truncates beginning and/or end if needed.
  # Plot curve.
  p <- ggplot(needata, aes(x = Relative.Time, y = CO2.umol.mol)) +
    geom_point() +
    theme_classic()
  # Calculate NEE.
  rawslope <- lm(CO2.umol.mol ~ Relative.Time, data = needata)$coeff[2]
  volume <- height * length^2
  area <- length^2
  NEE <- as.numeric( (1/0.044011) * rawslope * (volume / area) )
  cat('NEE is', NEE, 'mg C m^-2 s^-1.\n')
  # Convert date and time to R format.
  date_time <- strptime(do.call('paste', date_time[-4]), '%b %d %Y %H:%M:%S', tz = 'America/Denver')
  return(list(p=p, needata=needata, NEE=NEE, date_time=date_time))
}


NEEmetadata <- read.csv('data/NEE_metadata2015master.csv', stringsAsFactors = FALSE)
NEEs <- rep(0, nrow(NEEmetadata))


# Uses read_nee function to get slope from each reading.
for (i in 1:nrow(NEEmetadata)) {
  x <- read_nee(NEEmetadata$File.Name[i], start = NEEmetadata$start[i], end = NEEmetadata$end[i])
  NEEs[i] <- x$NEE
}

# Adds treatment information to the main dataframe.
NEEall <- transform(NEEmetadata, NEE=NEEs)
NEEall <- merge(NEEall, trt, all.x = TRUE)


NEE800 <- ddply(NEEall, .(Site, Plot), function(x) {
  dat <- data.frame(NEE = x$NEE, PAR = x$PAR.Within, Mesh = x$Mesh)
  intercept <- dat$NEE[dat$Mesh == 'dark']
  linear_fit <- lm(NEE ~ PAR + 0 + offset(rep(intercept, nrow(dat))), data = dat)
  alpha <- as.numeric(linear_fit$coeff[1])
  data.frame(NEE = alpha * 800 - intercept)
})

# Pull in elevations
elevs <- read.csv('data/elevations2015.csv', stringsAsFactors = FALSE)
NEE800 <- merge(NEE800, elevs[,1:2], all.x=T)

neeobs <- ddply(subset(NEE800, Site != 'Almont' & Site != 'Cinnamon'), .(Site), function(x) {
  NEE <- mean(x$NEE)
  NEE95 <- 1.96 * sd(x$NEE)/length(x$NEE)
  data.frame(NEE = NEE, NEE95 = NEE95, Elevation = x$Elevation[1])
})

# Combine LAI and NEE
NEEmean <- with(NEE800, tapply(NEE, Site, mean))
NEE95 <- with(NEE800, tapply(NEE, Site, function(x) 1.96*sd(x)/length(x)))
lainee <- data.frame(Site=names(LAImean),LAI=LAImean,LAI95=LAI95,NEE=NEEmean[-c(1,4)],NEE95=NEE95[-c(1,4)])

# LAI and standing biomass
sbmass <- read.csv('data/standingbiomass2015.csv', stringsAsFactors = FALSE)
sbm_mean <- ddply(sbmass, .(Site), summarize, 
                  biomass = mean(Standing.Biomass),
                  biomass95 = 1.96*sd(Standing.Biomass)/sqrt(length(Standing.Biomass)))

lainee <- merge(lainee, sbm_mean, all.x=TRUE)

###########################################################################

# Observational plant community

obs_comm <- loadWorkbook('data/obs_plant_community2015.xlsx')
obs_comm <- readWorksheet(obs_comm, sheet = getSheets(obs_comm)[1:13])
names(obs_comm) <- c('Almont','Slate_River','Washington_Gulch','Judd_Falls','Almont_Low','Brush_Creek','Painter_Boy','Schofield','Elkton','Bellview','Cinnamon','Cinni_Top','Mexican_Cut')
# Get only the "peak" measurement (late-season measurement)
obs_comm_peak <- lapply(obs_comm, function(x) x[(nrow(x)-9):nrow(x),])

# Experimental plant community
exp_comm <- loadWorkbook('./raw data/exp_plant_community2015.xlsx')
exp_comm <- readWorksheet(exp_comm, sheet = getSheets(exp_comm)[1:3])

# get rid of extraneous info in obs_comm_peak
obs_comm_peak <- 
  lapply(obs_comm_peak, function(x) {res <- x[,-c(1:5, (ncol(x)-4):ncol(x))]
                                     res[is.na(res)] <- 0
                                     res})

# Add max to the observational community (use control plots only).

max_ctrl = trt$Plot[trt$Site == 'Maxfield' & trt$Removal == 'C' & trt$Fert == 'C']
max_ctrl_comm <- exp_comm$Max[-(1:36),][max_ctrl,-(1:3)]
max_ctrl_comm[is.na(max_ctrl_comm)] <- 0
obs_comm_peak$Maxfield <- max_ctrl_comm

# Get weighted trait means for each of the observational sites.

load('./raw data/traitall.r')
trait_spmean <- ddply(trait, .(Site, Species), summarize, Height=mean(Height, na.rm=T), LMA=mean(LMA, na.rm=T), LDMC = mean(LDMC, na.rm=T), RML = mean(RML, na.rm=T))
trait_spmean$cover <- 0
for (i in 1:nrow(trait_spmean)) {
  mc <- try(mean(obs_comm_peak[[which(names(obs_comm_peak) == trait_spmean$Site[i])]][,trait_spmean$Species[i]]))
  if(!inherits(mc,'try-error')) trait_spmean$cover[i] <- mc
}

pcatraits <- trait_spmean[,3:6]
pca1 <- princomp(pcatraits[complete.cases(pcatraits), ], cor = TRUE) # This is the correct way to do the PCA. It only explains 33% of the variation in traits. PCA1 is associated weakly with shorter plants, and relatively strongly with high LMA leaves with more dry matter content, and high RML roots = "conservative" strategy!
trait_spmean$PCA1 <- NA
trait_spmean$PCA1[as.numeric(row.names(pca1$scores))] <- pca1$scores[,1]


my.wt.mean <- function(dat, wt, idx) weighted.mean(x = dat[idx], w=wt, na.rm=TRUE)

# How much do we have traits for?
totalcover <- unlist(lapply(obs_comm_peak, function(x) mean(rowSums(x))))
totaltraitcover <- ddply(trait_spmean, .(Site), summarize, cover=sum(cover))
totaltraitcover$total <- totalcover[match(totaltraitcover$Site, names(totalcover))]
totaltraitcover$prop <- with(totaltraitcover, cover/total)

traitboot_obs <- function(TABLE = trait_spmean, Trait) 
  ddply(TABLE, .(Site), function(x) {
    boot1 <- boot(data = x[,Trait], statistic = my.wt.mean, wt = x$cover, R = 49999, stype = 'i')
    #bootci1 <- boot.ci(boot1, conf=.95)
    bootsum1 <- summary(boot1)
    #ci.min <- bootci1$normal[2]
    #ci.max <- bootci1$normal[3]
    DF <- data.frame(mean = bootsum1$original, stderr = bootsum1$bootSE, median = bootsum1$bootMed)#, ci.min=ci.min, ci.max=ci.max)
    names(DF) <- paste(Trait, names(DF), sep = '.')
    return(DF)
  }, .progress = 'text')

traitboot_obs_LMA <- traitboot_obs(Trait = 'LMA')
traitboot_obs_RML <- traitboot_obs(Trait = 'RML')
traitboot_obs_LDMC <- traitboot_obs(Trait = 'LDMC')
traitboot_obs_Height <- traitboot_obs(Trait = 'Height')
traitboot_obs_PCA1 <- traitboot_obs(Trait = 'PCA1')

traitboot_obs_LMA$Site[c(1,5)] <- c('Almont_Obs','Cinnamon_Obs')
traitboot_obs_RML$Site[c(1,5)] <- c('Almont_Obs','Cinnamon_Obs')
traitboot_obs_LDMC$Site[c(1,5)] <- c('Almont_Obs','Cinnamon_Obs')
traitboot_obs_Height$Site[c(1,5)] <- c('Almont_Obs','Cinnamon_Obs')
traitboot_obs_PCA1$Site[c(1,5)] <- c('Almont_Obs','Cinnamon_Obs')

obsdat <- merge(lainee, cwmobs, all=T)
obsdat <- merge(obsdat, elevs[,1:2], all.x=T)
obsdat <- merge(obsdat, traitboot_obs_LMA)
obsdat <- merge(obsdat, traitboot_obs_RML)
obsdat <- merge(obsdat, traitboot_obs_LDMC)
obsdat <- merge(obsdat, traitboot_obs_Height)
obsdat <- merge(obsdat, traitboot_obs_PCA1)

################################################################################

# Add worldclim temp.

load('./raw data/2015dataobject.r')

latlongs <- read.csv('data/latlong.csv', stringsAsFactors = FALSE)
library(raster)
setwd('~')
wcdat <- getData('worldclim', var='bio', res=0.5, lat=latlongs$Lat, lon=latlongs$Long)
wcextract <- extract(wcdat, y=data.frame(x=latlongs$Long,y=latlongs$Lat))
latlongs$MAT <- wcextract[,1]/10
latlongs$MAP <- wcextract[,12]
latlongs$summertemp <- wcextract[,10]/10
latlongs$summerprecip <- wcextract[,18]

obsdat <- merge(obsdat, latlongs, all.x=TRUE)

#############################################################################

# Add ion exchange data, data from iButtons, and soil moisture data to obsdat.

ionexch <- read.csv('data/ionexchangeprobes.csv', stringsAsFactors = FALSE)
obsdat <- merge(obsdat, ionexch[,c(1,9:24)], all.x=TRUE)

obsmoisture <- ddply(obs_sm2015[171:230,], .(Site), function(x) data.frame(VWC = mean(x$VWC), VWC95 = 1.96*sd(x$VWC)/sqrt(length(x$VWC))))
obsdat <- merge(obsdat, obsmoisture, all.x=TRUE)

airtemps_obs2015$Site <- as.character(airtemps_obs2015$Site)
airtemps_obs2015$Site[airtemps_obs2015$Site=='Almont obs'] <- 'Almont_Obs'
airtemps_obs2015$Site[airtemps_obs2015$Site=='Cinnamon obs'] <- 'Cinnamon_Obs'
airtemps_obs2015$Site <- sub(' ', '_', airtemps_obs2015$Site)

# July temps
jultemp <- ddply(subset(airtemps_obs2015, as.Date(time) > '2015-06-30'), .(Site), summarize, ibutton_temp_July = mean(temp))
obsdat <- merge(obsdat, jultemp, all.x = TRUE)

# Soil temperatures

soiltemps_obs2015$Site <- as.character(soiltemps_obs2015$Site)
soiltemps_obs2015$Site[soiltemps_obs2015$Site=='Almont obs'] <- 'Almont_Obs'
soiltemps_obs2015$Site[soiltemps_obs2015$Site=='Cinnamon obs'] <- 'Cinnamon_Obs'
soiltemps_obs2015$Site <- sub(' ', '_', soiltemps_obs2015$Site)

# July soil temperatures
julsoiltemp <- ddply(subset(soiltemps_obs2015, as.Date(time) > '2015-06-30'), .(Site), summarize, ibutton_soiltemp_July = mean(temp))
obsdat <- merge(obsdat, julsoiltemp, all.x = TRUE)

###################################################################################

#write.table(obsdat, 'obsdat.csv', sep=',', row.names=FALSE)

# Addendum 29 Feb 2016: update NEE to new values, add FD, and add variance in trait values.

library(plyr)
NEE_bysite <- ddply(NEE800_hyperbolic2, .(Site), summarize, NEE=mean(NEE_lin_int), NEE95=qnorm(.975)*sd(NEE_lin_int)/sqrt(length(NEE_lin_int)))
# Correct NEE on obsdat
obsdat <- obsdat[, !(names(obsdat) %in% c('NEE','NEE95'))]
obsdat <- merge(obsdat, NEE_bysite, all.x=TRUE)
obsdat_old <- read.csv('obsdat.csv')
obsdatmerge <- merge(obsdat, obsdat_old[,c(1,9:23,51:55)], all=TRUE)
write.table(obsdatmerge, 'obsdat29feb.csv', sep=',', row.names=FALSE)

# Remove Brush Creek! (if desired)
obsdat_nobc <- subset(obsdat, Site != 'Brush_Creek')