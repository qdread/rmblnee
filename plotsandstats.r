obsdat <- read.csv('obsdat29Feb.csv')
obsdat <- obsdat[-which(obsdat$Site=='Brush_Creek'),] # Get rid of Brush Creek

# Plot other things against LAI.

plot(LAI ~ summertemp, data=obsdat)
plot(LAI ~ summerprecip, data=obsdat)
plot(LAI ~ LMA.mean, data=obsdat)
plot(LAI ~ RML.mean, data=obsdat)
# Reciprocal transform root mass length
plot(LAI ~ I(1/RML.mean), data=obsdat)
plot(LAI ~ Height.mean, data=obsdat)

# Quadratic regression

LAItemp_quad <- lm(LAI ~ summertemp + I(summertemp^2), data=obsdat)

tplot <-
  ggplot(obsdat, aes(x=summertemp, y=LAI, ymin=LAI-LAI95, ymax=LAI+LAI95)) +
  geom_pointrange(size = 1) +
  stat_smooth(method = 'lm', formula = y ~ x + I(x^2), lwd = 2) +
  theme_classic() + 
  geom_hline(yintercept = 0, lty = 'dotted')

LAIprecip_quad <- lm(LAI ~ summerprecip + I(summerprecip^2), data=obsdat)

pplot <-
  ggplot(obsdat, aes(x=summerprecip, y=LAI, ymin=LAI-LAI95, ymax=LAI+LAI95)) +
  geom_pointrange(size = 1) +
  stat_smooth(method = 'lm', formula = y ~ x + I(x^2), lwd = 2) +
  theme_classic() + 
  geom_hline(yintercept = 0, lty = 'dotted')

# Plot plant functional proportion.

ggplot(obsdat, aes(x=tallherb, y=NEE, ymin=NEE-NEE95, ymax=NEE+NEE95)) +
  geom_pointrange(size = 1) +
  #stat_smooth(method = 'lm', formula = y ~ x + I(x^2), lwd = 2) +
  theme_classic() + 
  geom_hline(yintercept = 0, lty = 'dotted')

summary(lm(NEE~tallherb,data=obsdat))

ggplot(obsdat, aes(x=tallherb, y=LAI, ymin=LAI-LAI95, ymax=LAI+LAI95)) +
  geom_pointrange(size = 1) +
  stat_smooth(method = 'lm', lwd = 2) +
  theme_classic() + 
  geom_hline(yintercept = 0, lty = 'dotted')
#ggsave('C:/Users/Q/Dropbox/my_rmbl_stuff2015/gradient_paper//figs/LAIbytallherbs.png',height=6,width=6)
summary(lm(LAI~tallherb,data=obsdat))

ggplot(obsdat, aes(x=graminoid, y=LAI, ymin=LAI-LAI95, ymax=LAI+LAI95)) +
  geom_pointrange(size = 1) +
  stat_smooth(method = 'lm', lwd = 2) +
  theme_classic() + 
  geom_hline(yintercept = 0, lty = 'dotted')
#ggsave('C:/Users/Q/Dropbox/my_rmbl_stuff2015/gradient_paper//figs/LAIbygrasses.png',height=6,width=6)
summary(lm(LAI~graminoid,data=obsdat))


# Combine all the names from observational community and classify by plant functional type. Then add plant functional type to obsdat object.
#write.table(sort(unique(unlist(sapply(obs_comm_peak,names)))), file='~/plantfunctypes.csv',row.names=F)
# pfts <- read.csv('C:/Users/Q/Dropbox/my_rmbl_stuff2015/gradient_paper//plantfunctypes.csv')
# communities <- sapply(obs_comm_peak,colSums)
# pft_sums <- sapply(communities, function(x) tapply(x, pfts$type[match(names(x), pfts$name)], sum))
# pft_sums[is.na(pft_sums)] <- 0
# 
# pft_sums <- t(apply(pft_sums, 2, function(x) x/sum(x)))
# 
# write.table(pft_sums[1:14,][order(dimnames(pft_sums[1:14,])[[1]]),1:5], file='C:/Users/Q/Dropbox/my_rmbl_stuff2015/gradient_paper/pft_proportions.csv', sep=',', row.names=F)

# Structural equation model

obsdat <- read.csv('C:/Users/Q/Dropbox/papers/gradient_paper//obsdat29feb.csv')


model1 <- 
  'NEE ~ LAI
   LAI ~ VWC + summertemp + summerprecip + LMA.mean + RML.mean + Height.mean + tallherb'

library(lavaan)

modelfit1 <- sem(model = model1, data = obsdat)
fitMeasures(modelfit1)['bic']
inspect(modelfit1, "rsquare")
parameterEstimates(modelfit1)
modificationIndices(modelfit1, sort = TRUE)

library(semPlot)
model1plot <- semPlotModel(modelfit1)
semPaths(object = model1plot, what = 'std', residuals = FALSE, intercepts= FALSE, thresholds = FALSE, nCharNodes = 0, sizeMan = 8, edge.label.cex = 1.5, fade = FALSE, label.cex = 2.5, curvature = 2, exoVar = FALSE, exoCov=FALSE)

model2 <- 
  'NEE ~ LAI + Height.mean + LMA.mean + RML.mean + summertemp + summerprecip
LAI ~ VWC + summertemp + summerprecip + LMA.mean + RML.mean + Height.mean + tallherb'

modelfit2 <- sem(model = model2, data = obsdat)
fitMeasures(modelfit2)['bic']
inspect(modelfit2, "rsquare")
parameterEstimates(modelfit2)
modificationIndices(modelfit2, sort = TRUE)

# remove parameters.
model3 <- 
  'NEE ~ LAI + Height.mean + RML.mean + summertemp + summerprecip
LAI ~ summertemp + summerprecip + RML.mean + Height.mean'

modelfit3 <- sem(model = model3, data = obsdat)
fitMeasures(modelfit3)['bic']
inspect(modelfit3, "rsquare")
parameterEstimates(modelfit3)
modificationIndices(modelfit3, sort = TRUE)

library(semPlot)
model3plot <- semPlotModel(modelfit3)

semPaths(object = model3plot, what = 'std', residuals = FALSE, intercepts= FALSE, thresholds = FALSE, nCharNodes = 0, sizeMan = 8, edge.label.cex = 1.5, fade = FALSE, label.cex = 2.5, curvature = 2, exoVar = FALSE, exoCov=FALSE)

# Make sure this is OK, transform data.
obsdat_transform <- transform(obsdat,
                              LMA.mean = LMA.mean*100,
                              RML.mean = RML.mean*10000,
                              summerprecip = summerprecip/10
                              )

modelfit3 <- sem(model = model3, data = obsdat_transform)
fitMeasures(modelfit3)['bic']
inspect(modelfit3, "rsquare")
parameterEstimates(modelfit3)
modificationIndices(modelfit3, sort = TRUE)

model3plot <- semPlotModel(modelfit3)
#pdf('C:/Users/Q/Dropbox/papers/gradient_paper/sem01Mar.pdf', height=7, width=7)
semPaths(object = model3plot, what = 'std', residuals = FALSE, intercepts= FALSE, thresholds = FALSE, nCharNodes = 0, sizeMan = 8, edge.label.cex = 1.5, fade = FALSE, label.cex = 2.5, curvature = 2, exoVar = FALSE, exoCov=FALSE)
#dev.off()


###########################

# Further addendum: plot by functional diversity.

ggplot(obsdat, aes(x=FDlma, y=NEE, ymin=NEE-NEE95, ymax=NEE+NEE95)) +
  geom_pointrange(size = 1) +
  stat_smooth(method = 'lm', lwd = 2) +
  theme_classic() + 
  geom_hline(yintercept = 0, lty = 'dotted')
summary(lm(NEE~FDlma,data=obsdat))

ggplot(obsdat, aes(x=LMA.mean, y=NEE, ymin=NEE-NEE95, ymax=NEE+NEE95)) +
  geom_pointrange(size = 1) +
  stat_smooth(method = 'lm', lwd = 2) +
  theme_classic() + 
  geom_hline(yintercept = 0, lty = 'dotted')
summary(lm(NEE~LMA.mean,data=obsdat))

ggplot(obsdat, aes(x=FDis, y=NEE, ymin=NEE-NEE95, ymax=NEE+NEE95)) +
  geom_pointrange(size = 1) +
  stat_smooth(method = 'lm', lwd = 2) +
  theme_classic() + 
  geom_hline(yintercept = 0, lty = 'dotted')
summary(lm(NEE~FDis,data=obsdat))

model4 <- 
  'NEE ~ LAI + Height.mean + RML.mean + summertemp + summerprecip + FDis
LAI ~ summertemp + summerprecip + RML.mean + Height.mean'

modelfit4 <- sem(model = model4, data = obsdat_transform)
fitMeasures(modelfit4)['bic']
inspect(modelfit4, "rsquare")
parameterEstimates(modelfit4)
modificationIndices(modelfit4, sort = TRUE)

model4plot <- semPlotModel(modelfit4)
#pdf('C:/Users/Q/Dropbox/papers/gradient_paper/sem01Mar_withFD.pdf', height=7, width=7)
semPaths(object = model4plot, what = 'std', residuals = FALSE, intercepts= FALSE, thresholds = FALSE, nCharNodes = 0, sizeMan = 8, edge.label.cex = 1.5, fade = FALSE, label.cex = 2.5, curvature = 2, exoVar = FALSE, exoCov=FALSE)
#dev.off()

# Model with only temperature, precip, and LAI

model_abiotic <- 
  'NEE ~ LAI + summertemp + summerprecip
LAI ~ summertemp + summerprecip'

modelfit_abiotic <- sem(model = model_abiotic, data = obsdat_transform)
fitMeasures(modelfit_abiotic)['bic']
inspect(modelfit_abiotic, "rsquare")
parameterEstimates(modelfit_abiotic)
#modificationIndices(modelfit4, sort = TRUE)

model_abplot <- semPlotModel(modelfit_abiotic)
#pdf('C:/Users/Q/Dropbox/papers/gradient_paper/sem07Mar_abioticonly.pdf', height=7, width=7)
semPaths(object = model_abplot, what = 'std', residuals = FALSE, intercepts= FALSE, thresholds = FALSE, nCharNodes = 0, sizeMan = 8, edge.label.cex = 1.5, fade = FALSE, label.cex = 2.5, curvature = 2, exoVar = FALSE, exoCov=FALSE)
#dev.off()

# linear models for same data
nee.lm <- lm(NEE~LAI+summertemp+summerprecip, data=obsdat_transform)
lai.lm <- lm(LAI~summertemp+summerprecip, data=obsdat_transform)
nee.coeff <- nee.lm$coefficients
lai.coeff <- lai.lm$coefficients