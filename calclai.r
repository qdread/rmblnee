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