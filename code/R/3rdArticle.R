# DESCRIPTION ##################################################################
# Source code used to develop the second article of my PhD research project. 
# The reference for the article is as follows:
# Samuel-Rosa, A.; Heuvelink, G.B.M.; Brus, D.; Vasques, G.M.; Anjos, L.H.C.
# Sampling for digital soil mapping in terra incognita.

# PART I - INITIAL SETTINGS ####################################################

# Workspace
rm(list = ls())
gc()

# Load packages
require(sp)
require(geoR)
require(pedometrics)
require(randomForest)
require(gstat)
require(car)

# Load auxiliary data
load("data/R/general.RData")

# Source user defined function
source("code/R/1stArticleHelper.R")
source("code/R/2ndArticleHelper.R")

# Set path to results (figures and tables)
fig_dir <- "res/fig/2ndArticle/"
tab_dir <- "res/tab/2ndArticle/"

# Initiate GRASS GIS (64) section
spgrass6::initGRASS(gisBase = "/usr/lib/grass64/", 
                    gisDbase = path.expand(dbGRASS),
                    location = "dnos-sm-rs", mapset = "predictions", 
                    pid = Sys.getpid(), override = TRUE)
system("g.region rast=dnos.raster")
spgrass6::gmeta6()

# PART II - CALIBRATION POINT SOIL DATA ########################################
# We use only the first n = 350 point soil observations available in the Santa 
# Maria dataset.
cal_data <- paste(point_data, "labData.csv", sep = "")
cal_data <- read.table(cal_data, sep = ";", header = TRUE, dec = ".", 
                       na.strings = "na", stringsAsFactors = FALSE)
id <- c("sampleid", "longitude", "latitude", "CLAY", "BUDE")
id <- match(id, colnames(cal_data))
cal_data <- cal_data[1:350, id]
coordinates(cal_data) <- ~ longitude + latitude
proj4string(cal_data) <- sirgas2000
cal_data <- spTransform(cal_data, wgs1984utm22s)
plot(cal_data, pch = 20, cex = 0.5)
rm(id)

# PART III - PARENT MODELS #####################################################
# The parent model with a linear spatial trend is defined using the best
# performing linear mixed model calibrated to the CLAY data by A. Samuel-Rosa 
# and collaborators in 2015 (Samuel-Rosa, A.; Heuvelink, G. B. M.; Vasques, G. 
# M. & Anjos, L. H. C. Do more detailed environmental covariates deliver more 
# accurate soil maps? Geoderma. v. 243-244, p. 214-227, 2015.)
# Our original idea was to use the best performing linear mixed model calibrated
# by Samuel-Rosa et al. (2015) to the three soil properties considered in their
# study: CLAY, ORCA, and ECEC. However, we have decided that working with all
# three soil properties would not add much to our study. Considering the 
# objectives of our study, it seemed to be more appropriate to use only two 
# soil properties, requiring their stochastic residuals to have different 
# quantities of spatially correlated variance. This decision would also reduce 
# the computational burden.
# We found out that all three best performing linear mixed models calibrated by
# Samuel-Rosa et al. (2015) describe processes with a large quantity of 
# spatially correlated variance in the residuals. As such, we decided to keep 
# only the linear mixed model with the largest quantity of spatially correlated
# variance in the residuals (CLAY), and select another soil variable from the
# Santa Maria dataset that contains only a small quantity of spatially 
# correlated variance in the residuals. The exploratory spatial data analysis 
# has shown that bulk density (BUDE) is that soil property.

load("data/R/2ndArticlePartI.rda")
vgmSCV(clay_lmm)
vgmSCV(orca_lmm)
vgmSCV(ecec_lmm)
rm(orca_lmm, ecec_lmm)

# We can see bellow that BUDE has an empirical distribution close to Gaussian.
# This means that is does not require transformation if we were to calibrate
# a linear mixed model (the main weakness of the BUDE data is that it is 
# available at only n = 282 observations). In fact, after we have decided to use
# only CLAY and BUDE, the idea was to calibrate linear and nonlinear spatial
# trends to both soil properties. That is where our concern with the empirical
# distribution of BUDE comes from. However, we later decided that the study did
# not require two parent models with linear and nonlinear spatial trends. It 
# seems reasonable to use only one parent model with each form of spatial trend.
# This would then be CLAY with linear and BUDE with nonlinear spatial trends.
bude_id <- which(!is.na(cal_data$BUDE))
length(bude_id)
plotHD(cal_data$BUDE[bude_id])
plotESDA(cal_data$BUDE[bude_id], cal_data$latitude[bude_id], 
         cal_data$longitude[bude_id])

# Altough we are not going to use BUDE to build a parent model with a linear
# spatial trend, we still need to calibrate a linear model to the BUDE data.
# The reason for that is that we want to use the same covariate selection 
# procedure used to calibrate the linear mixed model to the CLAY data. This will
# be of great help in case we change our minds later and decide to define two 
# parent models with each soil property (or if a reviewer requires us to do 
# that). Besides, covariate selection methods for random forest models are not 
# as straightforward as those for linear regression model.
# To select the covariates to enter the nonlinear parent model, we process the 
# calibration data using the same script used by Samuel-Rosa et al. (2015). The 
# code under sections 'COVARIATES' and 'CANDIDATE MODELS' is run in its 
# entirety. The code used by Samuel-Rosa et al. (2015) to calibrate the linear 
# mixed models and check the residuals is copied and pasted here as is.

# RUN THE CODE UNDER 'COVARIATES' AND 'CANDIDATE MODELS' FROM 
# SAMUEL-ROSA ET AL. (2015)!!!

# Update formulas
formula <- lapply(forms$main, update, BUDE ~ .)
data <- na.omit(cal_data@data)
# Build linear model series
bude_lm <- buildMS(formula, data = data, vif = TRUE, aic = TRUE)
# Get the statistics of the linear model series
bude_lm_stats <- statsMS(bude_lm, combs$num, "rmse")
# Model series plot
plotMS(bude_lm_stats, grid = c(2:6), line = "ADJ_r2", ind = 2, 
       color = c("lightyellow", "palegreen"), main = "stepwise selection")
# Get the best linear model
bude_lm <- bude_lm[tail(bude_lm_stats, 1)$id][[1]]
# Analysis of the residuals
plot(bude_lm, which = c(1:6))
plotESDA(residuals(bude_lm), lon = coordinates(cal_data)[bude_id, 1], 
         lat = coordinates(cal_data)[bude_id, 2])
# Amount of variance explained
cor(fitted(bude_lm), cal_data$BUDE[bude_id]) ^ 2





# The multiple linear regression model explains 43.49% of the variance. The 
# residuals empirical distribution is close to Gaussian. Their spatial 
# dependence is low. The lack of short-distance point-pairs limits the 
# estimation of the nugget effect. We fit a linear mixed model to the BUDE data
# setting fix.kappa = FALSE. The REML estimate of kappa is equal to 2.756, 
# which is equivalent to a Gaussian model. The amount of spatially correlated
# variance i 28%. We re-calibrate the linear mixed model using a Gaussian
# covariance model. The amount of spatially correlated variance falls to 
# 26.18%. The linear mixed model explains 61.18% of the variance of the 
# calibration data.

# # Calculate and check the empirical variogram
# lags <- vgmLags(cal_data@coords, n = 7, cutoff = 1)
# v <- fitVariog(y = "BUDE", model = bude_lm, lambda = 1,
#                data = as.data.frame(cal_data), breaks = lags)
# plot(v, pch = 20)
# # Estimate model parameters using REML
# ini.cov.pars <- as.matrix(expand.grid(c(50, 75, 100), c(400, 500, 600)))
# nugget <- c(275, 300, 325)
# bude_lmm <- fitREML(y = "BUDE", model = bude_lm, nugget = nugget,
#                     data = as.data.frame(cal_data), ini.cov.pars = ini.cov.pars, 
#                     cov.model = "gaussian")
# summary(bude_lmm)
# lines(bude_lmm)
# # Spatially correlated variance
# vgmSCV(bude_lmm)
# # Amount of variance explained
# cor(fitted(bude_lmm), cal_data$BUDE[bude_id]) ^ 2

# PART IV - REGRESSION-KRIGING MODELS
# We calibrate RANDOM FOREST models using the existing point soil data to
# model the deterministic spatial trend. The random forest models are calibrated
# using the same covariates included in the best performing linear mixed models.
# Covariate values are retrieved from the trend matrix of the linear mixed
# models so that we do not have to sample the raster data in the GRASS GIS
# database. The random forest algorithm is set to grow 500 trees. At each split,
# it uses 03 covariates for CLAY, and 04 for BUDE. These values are close to the
# default values (1/3), and changing them has only a minor effect on the amount
# of variance explained by the random forest model. For CLAY, the random forest
# model explains 58.44% of the variance of the calibration data. For BUDE, the
# amount of variance explained is only 37.56%.
# CLAY
# set.seed(2001)
# clay_rkm <- list()
# clay_rkm$clay_rfm <- 
#   randomForest::randomForest(x = trend.matrix(clay_lmm), mtry = 3,
#                              y = cal_data$CLAY, ntree = 500)
# cor(clay_rkm$clay_rfm$predicted, cal_data$CLAY) ^ 2

# We now use the set of covariates selected to enter de linear regression model
# to calibrate a random forest to the BUDE data. The random forest model
# explains 37.56% of the variance of the BUDE data. We also produce a figure
# showing the bias of the random forest model: low values are over predicted,
# while high values are under predicted.

set.seed(2001)
bude_rkm <- list()
bude_rkm$bude_rfm <- 
  randomForest::randomForest(x = bude_lm$model[, -1], mtry = 4,
                             y = bude_lm$model[, 1], ntree = 500)
cor(bude_rkm$bude_rfm$predicted, cal_data$BUDE[bude_id]) ^ 2

bias <- data.frame(obs = cal_data$BUDE[bude_id], 
                   est = bude_rkm$bude_rfm$predicted,
                   res = bude_rkm$bude_rfm$predicted - cal_data$BUDE[bude_id])
bias <- bias[order(bias$obs), ]
dev.off()
png("tmp/random-forest-bias.png", width = 150, height = 150, units = "mm",
    res = 300)
plot(y = bias$obs, x = bias$est, ylab = "Observed", xlab = "Fitted",
     xlim = c(60, 200), ylim = c(60, 200),
     main = "Random forest model calibrated to bulk density data",
     sub = "Black: 1:1 line, Red: linear model fit")
abline(a = 0, b = 1)
abline(lm(bias$obs ~ bias$est), col = "red")
dev.off()




# ORCA
# orca_rkm <- list()
# orca_rkm$orca_rfm <- 
#   randomForest::randomForest(x = trend.matrix(orca_lmm), mtry = 4,
#                              y = cal_data$ORCA, ntree = 500)
# ECEC
# ecec_rkm <- list()
# ecec_rkm$ecec_rfm <- 
#   randomForest::randomForest(x = trend.matrix(ecec_lmm), mtry = 4,
#                              y = cal_data$ECEC, ntree = 500)





# Spatial exploratory data analysis of the random forest residuals
# The residuals of the random forest models are defined as the difference 
# between the predicted and observed values. We check the residuals regarding
# their empirical distribution and find that the BUDE residuals are Gaussian
# distributed, with a slight negative mean.

# BUDE
cal_data$BUDEres <- NA
cal_data$BUDEres[bude_id] <- bude_rkm$bude_rfm$predicted - cal_data$BUDE[bude_id]
plotHD(na.omit(cal_data$BUDEres), BoxCox = FALSE)



# The CLAY residuals are slightly skewed to the left due to a few 
# under-predictions. However, this skewness should not be a problem for 
# calibrating a Gaussian model. The mean CLAY residual is positive.
# # CLAY
# cal_data$CLAYres <- clay_rkm$clay_rfm$predicted - cal_data$CLAY
# plotHD(cal_data$CLAYres, BoxCox = FALSE)
# 
# ORCA
# cal_data$ORCAres <- orca_rkm$orca_rfm$predicted - cal_data$ORCA
# plotHD(cal_data$ORCAres, BoxCox = FALSE)
# ECEC
# cal_data$ECECres <- ecec_rkm$ecec_rfm$predicted - cal_data$ECEC
# plotHD(cal_data$ECECres, BoxCox = FALSE)
# 
# Check range
# sapply(cal_data@data[, c("CLAYres", "ORCAres", "ECECres")], range)
# 







# Next, we check the random forest residuals regarding their spatial correlation
# structure using bubble plots, empirical isotropic and anysotropic variograms,
# and variogram maps. We use exponentially spaced lag-distance classes. The main
# results are as follows:
# BUDE: the residuals present spatial autocorrelation, being difficult to
#    identify the type of correlation structure (exponential or Gaussian); the
#    nugget variance is much larger than zero; there also is a slightly 
#    departure from stationarity.

lags <- vgmLags(cal_data@coords, n = 8, cutoff = 1)[-1]

# BUDE
plotESDA(na.omit(cal_data$BUDEres), lat = coordinates(cal_data)[bude_id, 2],
         lon = coordinates(cal_data)[bude_id, 1], lags = lags)
tmp <- na.omit(as.data.frame(cal_data))
coordinates(tmp) <- colnames(coordinates(cal_data))
v <- gstat::variogram(BUDE ~ 1, data = tmp, boundaries = lags)
plot(v, pch = 20, col = "black")
rm(tmp)








# CLAY: the residuals present spatial autocorrelation with a sharp increase 
#    of the semi-variance at short separation distances, suggesting an 
#    exponential covariance structure; the nugget variance appears to be only
#    slightly larger than zero; there is a slightly departure from stationarity,
#    but it is hard to tell if this departure is due to non-stationarity in the
#    covariance structure or in the spatial mean.
# # CLAY
# plotESDA(cal_data$CLAYres, lat = coordinates(cal_data)[, 2],
#          lon = coordinates(cal_data)[, 1], lags = lags)
# v <- gstat::variogram(CLAYres ~ 1, data = cal_data, boundaries = lags)
# plot(v, pch = 20, col = "black")
# 
# ORCA
# plotESDA(cal_data$ORCAres, lat = cal_data$latitude, 
#          lon = cal_data$longitude, lags = lags)
# v <- gstat::variogram(ORCAres ~ 1, data = cal_data, alpha = alpha, 
#                       boundaries = lags)
# plot(v, pch = 20, col = "black")
# ECEC
# plotESDA(cal_data$ECECres, lat = cal_data$latitude, 
#          lon = cal_data$longitude, lags = lags)
# v <- gstat::variogram(ECECres ~ 1, data = cal_data, alpha = alpha, 
#                       boundaries = lags)
# plot(v, pch = 20, col = "black")
# rm(v, lags, alpha)






# In the traditional regression-kriging setting, we would model the stochastic
# residuals using a geostatistical model assuming a constant spatial mean 
# beta = 0. However, our approach consists of calibrating a linear mixed model
# to the BUDE data using the random forest predictions as a covariate. This step
# makes the prediction step simpler, and enables us to correct some of the 
# existing bias of the random forest model. In their study, Samuel-Rosa et al.
# (2015) used an exponential covariance model for CLAY. Here, we use a Gaussian
# covariance model for BUDE.
# Despite the small evidence for a departure from stationarity, we use an
# isotropic model to reduce the computational burden in next processing steps. 

require(geoR)
lags <- vgmLags(cal_data@coords)

# BUDE
v <- geoR::variog(coords = coordinates(cal_data)[bude_id, ], 
                  data = na.omit(cal_data$BUDEres), breaks = lags)
plot(v, pch = 20)
ini.cov.pars <- as.matrix(expand.grid(c(50, 75, 100), c(750, 1000, 1250)))
geodata <- as.geodata(cal_data, data.col = "BUDE",
                      covar.col = colnames(bude_lm$model)[-1])
bude_rkm$bude_vgm <- likfit(geodata = geodata, 
                            lik.method = "REML", nugget = c(350, 375, 400),
                            ini.cov.pars = ini.cov.pars, cov.model = "gaussian")
rm(ini.cov.pars)
lines(bude_rkm$bude_vgm, col = "blue")
summary(bude_rkm$bude_vgm)
vgmSCV(bude_rkm$bude_vgm)
cor(fitted(bude_rkm$bude_vgm), na.omit(cal_data$BUDEres)) ^ 2
# save image
dev.off()
png("res/fig/2ndArticle/bude-res-vario.png")
plot(v, pch = 20, main = "BUDE residuals")
lines(bude_rkm$bude_vgm)
dev.off()




# The nugget variance of the CLAY residuals is regarded as fixed and determined
# visually because geoR::likfit always estimates nugget = 0. For the purposes 
# of our study, it is important to have some nugget variance. Besides, our
# calibration data is sub-optimal for estimating the nugget due to the lack of
# short distance points. We set the nugget of the CLAY model as to produce an
# amount of spatially correlated variance of 85%. The geostatistical model
# explains 99.83% and 29.91% of the calibration data, respectively.
# 
# # CLAY
# v <- geoR::variog(coords = coordinates(cal_data), data = cal_data$CLAYres,
#                   breaks = lags)
# plot(v, pch = 20)
# ini.cov.pars <- as.matrix(expand.grid(c(3000, 3500, 4000), c(400, 500, 600)))
# clay_rkm$clay_vgm <- likfit(as.geodata(cal_data, data.col = "CLAYres"), 
#                             lik.method = "REML", nugget = c(500, 1000, 1500),
#                             ini.cov.pars = ini.cov.pars)
# nugget <- clay_rkm$clay_vgm$sigmasq * 0.15
# clay_rkm$clay_vgm <- likfit(as.geodata(cal_data, data.col = "CLAYres"), 
#                             lik.method = "REML", nugget = nugget,
#                             ini.cov.pars = ini.cov.pars, fix.nugget = TRUE)
# lines(clay_rkm$clay_vgm)
# summary(clay_rkm$clay_vgm)
# vgmSCV(clay_rkm$clay_vgm)
# cor(fitted(clay_rkm$clay_vgm), cal_data$CLAYres) ^ 2
# rm(ini.cov.pars, nugget)
# # save image
# dev.off()
# png("res/fig/2ndArticle/clay-res-vario.png")
# plot(v, pch = 20, main = "CLAY residuals")
# lines(clay_rkm$clay_vgm)
# dev.off()
# 
# ORCA
# v <- geoR::variog(coords = coordinates(cal_data), data = cal_data$ORCAres,
#                   breaks = lags)
# plot(v, pch = 20)
# orca_rkm$orca_vgm <- fit_reml(obj = cal_data, z = "ORCAres")
# lines(orca_rkm$orca_vgm)
# summary(orca_rkm$orca_vgm)
# save image
# dev.off()
# png("res/fig/2ndArticle/orca-res-vario.png")
# plot(v, pch = 20, main = "ORCA residuals")
# lines(orca_rkm$orca_vgm)
# dev.off()
# 
# ECEC
# v <- geoR::variog(coords = coordinates(cal_data), data = cal_data$ECECres,
#                   breaks = lags)
# plot(v, pch = 20)
# ecec_rkm$ecec_vgm <- fit_reml(obj = cal_data, z = "ECECres")
# lines(ecec_rkm$ecec_vgm)
# summary(ecec_rkm$ecec_vgm)
# save image
# dev.off()
# png("res/fig/2ndArticle/ecec-res-vario.png")
# plot(v, pch = 20, main = "ECEC residuals")
# lines(ecec_rkm$ecec_vgm)
# dev.off()

# Check residuals beta and range
clay_rkm$clay_vgm$beta / diff(range(cal_data$CLAYres)) * 100
bude_rkm$bude_vgm$beta / diff(range(na.omit(cal_data$BUDEres))) * 100
# orca_rkm$orca_vgm$beta / diff(range(cal_data$ORCAres))
# ecec_rkm$ecec_vgm$beta / diff(range(cal_data$ECECres))

# Ratio of spatial dependence
rsd(clay_rkm$clay_vgm)
rsd(bude_rkm$bude_vgm)
# rsd(orca_rkm$orca_vgm)
# rsd(ecec_rkm$ecec_vgm)

# PART VI - SPATIAL PREDICTIONS ################################################

# CLAY

# BUDE

# PART V - SAVE/LOAD DATA ######################################################
ls()
save(bude_id, bude_lm, bude_lmm, bude_rkm, clay_lmm, clay_rkm, cal_data,
     file = "data/R/2ndArticlePartII.rda")
load(file = "data/R/2ndArticlePartII.rda")

# PART VI - OPTIMIZATION OF SAMPLE CONFIGURATIONS ##############################
# We optimize sample configurations regarding five (05) single objective 
# functions (DIST, CORR, PPL, MSSD, MKV), and three (03) multi-objective 
# functions (ACDC and SPAN). We use three (03) sample sizes: 100, 200, and 400.
# For objective functions that depend on the covariates, such as DIST, CORR,
# ACDC, and SPAN, different sample configurations are optimized for CLAY, ORCA,
# and ECEC.

require(spsann)
require(parallel)

# initial settings
pts <- c(100, 200, 400)
iter <- 100000
x_minmax <- c(5, 5230)
y_minmax <- c(5, 5970)
models <- c("clay_lmm", "bude_lmm")
# prepare candi
candi <- spgrass6::readRAST6(vname = "dnos.raster")
gridded(candi) <- FALSE
candi <- coordinates(candi)
colnames(candi) <- c("x", "y")
# prepare covars
covars <- unique(c(trend.terms(clay_lmm), trend.terms(bude_lmm)))
covars <- na.omit(spgrass6::readRAST6(vname = covars)@data)
id <- which(sapply(covars, is.integer))
covars[, id] <- lapply(covars[, id], as.factor)

# DIST
# Optimize the sample configuration
osc_dist <- list(CLAY = list(), BUDE = list())
for (j in 1:length(models)) {
  for (i in 1:length(pts)) {
    set.seed(2001)
    osc_dist[[j]][[i]] <- 
      optimDIST(points = pts[i], candi = candi, track = TRUE,
                covars = covars[, trend.terms(eval(parse(text = models[j])))],
                iterations = iter, x.max = x_minmax[2], y.max = y_minmax[2], 
                x.min = x_minmax[1], y.min = y_minmax[1])
  }
}
# Check and save osc
plot.OSC(osc_dist$CLAY[[1]])
save(osc_dist, file = "data/R/2ndArticle_osc_dist.rda")
rm(osc_dist)

# CORR
# Optimize the sample configuration
# There is a bug in optimCORR causing the following error:
# Error in if (new_energy <= old_energy) { :
#   missing value where TRUE/FALSE needed
# This bug affects optimACDC and optimSPAN. It was found with j = 2 and i = 1,
# when 4% of the optimization was completed.
osc_corr <- list(CLAY = list(), BUDE = list())
for (j in 1:length(models)) {
  for (i in 1:length(pts)) {
    set.seed(2001)
    osc_corr[[j]][[i]] <- 
      optimCORR(points = pts[i], candi = candi, track = TRUE,
                covars = covars[, trend.terms(eval(parse(text = models[j])))],
                iterations = iter, x.max = x_minmax[2], y.max = y_minmax[2], 
                x.min = x_minmax[1], y.min = y_minmax[1])
  }
}
# Check and save osc
plot.OSC(osc_corr$CLAY[[1]])
save(osc_corr, file = "data/R/2ndArticle_osc_corr.rda")
rm(osc_corr)

# ACDC
# Calculate nadir and utopia points
load("data/R/2ndArticle_osc_corr.rda")
load("data/R/2ndArticle_osc_dist.rda")
nadir_acdc <- list(CLAY = list(), BUDE = list())
for (j in 1:length(models)) {
  for (i in 1:length(pts)) {
    nadir_acdc[[j]][[i]] <- 
      nadirACDC(osc = list(DIST = osc_dist[[j]][[i]], 
                           CORR = osc_corr[[j]][[i]]), candi = candi,
                covars = covars[, trend.terms(eval(parse(text = models[j])))])
  }
}
rm(osc_corr, osc_dist)
# Optimize the sample configuration
osc_acdc <- list(CLAY = list(), BUDE = list())
for (j in 1:length(models)) {
  for (i in 1:length(pts)) {
    nadir <- list(user = list(DIST = nadir_acdc[[j]][[i]]["CORR", "DIST"],
                              CORR = nadir_acdc[[j]][[i]]["DIST", "CORR"]))
    utopia <- list(user = list(DIST = nadir_acdc[[j]][[i]]["DIST", "DIST"],
                               CORR = nadir_acdc[[j]][[i]]["CORR", "CORR"]))
    set.seed(2001)
    osc_acdc[[j]][[i]] <- 
      optimACDC(points = pts[i], candi = candi, track = TRUE,
                covars = covars[, trend.terms(eval(parse(text = models[j])))],
                iterations = iter, x.max = x_minmax[2], y.max = y_minmax[2], 
                x.min = x_minmax[1], y.min = y_minmax[1], 
                nadir = nadir, utopia = utopia)
  }
}
# Check and save osc
plot.OSC(osc_acdc$CLAY[[1]])
save(osc_acdc, nadir_acdc, file = "data/R/2ndArticle_osc_acdc.rda")
rm(osc_acdc, nadir_acdc)

# PPL
# Set-up clusters
nclus <- 3
clus <- rep("localhost", nclus)
cl <- makeCluster(clus, type = "SOCK")
clusterEvalQ(cl, require(spsann))
clusterExport(cl, list("pts", "iter", "x_minmax", "y_minmax", "candi"))
# Optimize the sample configuration
osc_ppl <- parLapply(cl, pts, function (pts) {
  set.seed(2001)
  optimPPL(points = pts, candi = candi, track = TRUE,
           iterations = iter, x.max = x_minmax[2], y.max = y_minmax[2], 
           x.min = x_minmax[1], y.min = y_minmax[1])
})
# Stop cluster
stopCluster(cl)
# Check and save osc
plot.OSC(osc_ppl[[1]])
save(osc_ppl, file = "data/R/2ndArticle_osc_ppl.rda")
rm(osc_ppl)

# MSSD
# Set-up clusters
nclus <- 3
clus <- rep("localhost", nclus)
cl <- makeCluster(clus, type = "SOCK")
clusterEvalQ(cl, require(spsann))
clusterExport(cl, list("pts", "iter", "x_minmax", "y_minmax", "candi"))
# Optimize the sample configuration
osc_mssd <- parLapply(cl, pts, function (pts) {
  set.seed(2001)
  optimMSSD(points = pts, candi = candi, track = TRUE, 
            iterations = iter, x.max = x_minmax[2], y.max = y_minmax[2], 
            x.min = x_minmax[1], y.min = y_minmax[1])
})
# Stop cluster
stopCluster(cl)
# Check and save osc
plot.OSC(osc_mssd[[1]])
save(osc_mssd, file = "data/R/2ndArticle_osc_mssd.rda")
rm(osc_mssd)

# SPAN
# Calculate nadir and utopia points
nadir_span <- list(CLAY = list(), BUDE = list())
for (j in 1:length(models)) {
  for (i in 1:length(pts)) {
    nadir_span[[j]][[i]] <- 
      nadirSPAN(osc = list(DIST = osc_dist[[j]][[i]], CORR = osc_corr[[j]][[i]],
                           PPL = osc_ppl[[i]], MSSD = osc_mssd[[i]]),
                candi = candi, x.max = x_minmax[2], y.max = y_minmax[2], 
                covars = covars[, trend.terms(eval(parse(text = models[j])))])
  }
}
# Optimize the sample configuration
osc_span <- list(CLAY = list(), BUDE = list())
for (j in 1:length(models)) {
  for (i in 1:length(pts)) {
    nadir <- list(user = list(DIST = max(nadir_span[[j]][[i]][, "DIST"]),
                              CORR = max(nadir_span[[j]][[i]][, "CORR"]),
                              PPL = max(nadir_span[[j]][[i]][, "PPL"]),
                              MSSD = max(nadir_span[[j]][[i]][, "MSSD"])))
    utopia <- list(user = list(DIST = min(nadir_span[[j]][[i]][, "DIST"]),
                               CORR = min(nadir_span[[j]][[i]][, "CORR"]),
                               PPL = min(nadir_span[[j]][[i]][, "PPL"]),
                               MSSD = min(nadir_span[[j]][[i]][, "MSSD"])))
    set.seed(2001)
    osc_span[[j]][[i]] <- 
      optimSPAN(points = pts[i], candi = candi, track = TRUE, plotit = TRUE,
                covars = covars[, trend.terms(eval(parse(text = models[j])))],
                iterations = iter, x.max = x_minmax[2], y.max = y_minmax[2], 
                x.min = x_minmax[1], y.min = y_minmax[1],
                nadir = nadir, utopia = utopia)
  }
}
# plot.OSC(osc_span$CLAY[[3]])
save(osc_span, file = "data/R/2ndArticle_osc_span.rda")
rm(osc_span)

# MKV
# Optimize the sample configuration
# osc_mkv <- list(CLAY = list(), BUDE = list())
# for (j in 1:length(models)) {
#   for (i in 1:length(pts)) {
#     set.seed(2001)
#     osc_mkv[[j]][[i]] <- 
#       optimMKV(points = pts[i], candi = candi, nmax = 20, track = TRUE,
#                covars = covars[, trend.terms(eval(parse(text = models[j])))],
#                eqn = ..., vgm = ...,
#                iterations = iter, x.max = x_minmax[2], y.max = y_minmax[2], 
#                x.min = x_minmax[1], y.min = y_minmax[1])
#   }
# }
# plot.OSC(osc_mkv$CLAY[[3]])
# save(osc_mkv, file = "data/R/2ndArticle_osc_mkv.rda")
# rm(osc_mkv)

rm(candi, covars)
gc()

# PART VII - SEQUENTIAL GAUSSIAN SIMULATIONS ###################################

# Prepare prediction grid
system("r.mask -o input=buffer_BASIN_10")
dnos_raster <- spgrass6::readRAST6(vname = "dnos.raster")
gridded(dnos_raster) <- FALSE
dnos_raster <- coordinates(dnos_raster)
colnames(dnos_raster) <- c("x", "y")
# data("meuse.grid")
# dnos_raster <- meuse.grid[, 1:2]
# rm(meuse.grid)

# Prepare covariates
covars <- unique(c(trend.terms(clay_lmm), trend.terms(bude_lmm)))
covars <- na.omit(spgrass6::readRAST6(vname = covars)@data)
# set.seed(2001)
# covars <- covars[sample(1:nrow(covars), nrow(dnos_raster)), ]
# id <- which(sapply(covars, is.integer))
# covars <- covars[order(covars[, -id][, 1]), ]

# Settings
seed <- c(2001, 1984, 1964, 1500, 1822)
nsim <- 1000 / length(seed)
nmax <- 20

# CLAY - linear spatial trend
g_dummy <- gstat(formula = update(clay_lmm$trend, z ~ .),
                 dummy = TRUE, beta = clay_lmm$beta, nmax = nmax,
                 maxdist = clay_lmm$practicalRange,
                 model = as.vgm.variomodel(clay_lmm))
for (i in 1:length(seed)) {
  set.seed(seed[i])
  clay_lmm_sim <- 
    predict(g_dummy, nsim = nsim, 
            newdata = newdata4sim(dnos_raster, covars, clay_lmm), 
            debug.level = -1)@data
  file <- paste("data/R/2ndArticle_clay_lmm_sim_", seed[i], ".rda", sep = "")
  save(clay_lmm_sim, file = file)
  rm(clay_lmm_sim)
  gc()
}
rm(g_dummy)

# BUDE - linear spatial trend
g_dummy <- gstat(formula = update(bude_lmm$trend, z ~ .),
                 dummy = TRUE, beta = bude_lmm$beta, nmax = nmax,
                 maxdist = bude_lmm$practicalRange,
                 model = as.vgm.variomodel(bude_lmm))
for (i in 1:length(seed)) {
  set.seed(seed[i])
  bude_lmm_sim <- 
    predict(g_dummy, nsim = nsim, 
            newdata = newdata4sim(dnos_raster, covars, bude_lmm), 
            debug.level = -1)@data
  file <- paste("data/R/2ndArticle_bude_lmm_sim_", seed[i], ".rda", sep = "")
  save(bude_lmm_sim, file = file)
  rm(bude_lmm_sim)
  gc()
}
rm(g_dummy)

# CLAY - non-linear spatial trend
g_dummy <- gstat(formula = z ~ 1, dummy = TRUE, beta = 0, nmax = nmax,
                 maxdist = clay_rkm$clay_vgm$practicalRange,
                 model = as.vgm.variomodel(clay_rkm$clay_vgm))
for (i in 1:length(seed)) {
  set.seed(seed[i])
  clay_rkm_sim <- 
    predict(g_dummy, nsim = nsim, 
            newdata = newdata4sim(dnos_raster, trend = "RF"), 
            debug.level = -1)@data + 
    randomForest::predict(clay_rkm$clay_rfm, 
                          newdata = covars[, trend.terms(clay_lmm)])
  file <- paste("data/R/2ndArticle_clay_rkm_sim_", seed[i], ".rda", sep = "")
  save(clay_rkm_sim, file = file)
  rm(clay_rkm_sim)
  gc()
}
rm(g_dummy)

# BUDE - non-linear spatial trend
g_dummy <- gstat(formula = z ~ 1, dummy = TRUE, beta = 0, nmax = nmax,
                 maxdist = bude_rkm$bude_vgm$practicalRange,
                 model = as.vgm.variomodel(bude_rkm$bude_vgm))
for (i in 1:length(seed)) {
  set.seed(seed[i])
  bude_rkm_sim <- 
    predict(g_dummy, nsim = nsim, 
            newdata = newdata4sim(dnos_raster, trend = "RF"), 
            debug.level = -1)@data +
    randomForest::predict(bude_rkm$bude_rfm, 
                          newdata = covars[, trend.terms(bude_lmm)])
  file <- paste("data/R/2ndArticle_bude_rkm_sim_", seed[i], ".rda", sep = "")
  save(bude_rkm_sim, file = file)
  rm(bude_rkm_sim)
  gc()
}
rm(g_dummy)


# Check simulations
file <- paste("data/R/2ndArticle_bude_lmm_sim_", seed[i], ".rda", sep = "")
load(file)
tmp <- cbind(dnos_raster, bude_lmm_sim)
gridded(tmp) <- ~ x + y
spplot(tmp)

# PART VIII - SAMPLING #########################################################

# Settings
seed <- c(2001, 1984, 1964, 1500, 1822)
nsim <- 1000 / length(seed)
pts <- c(100, 200, 400)
vars <- c("clay", "bude")
model <- c("lmm", "rkm")
load("data/R/2ndArticle_osc_dist.rda")
load("data/R/2ndArticle_osc_ppl.rda")

# Prepare covariates and validation data
system("r.mask -o input=buffer_BASIN_10")
dnos_raster <- spgrass6::readRAST6(vname = "dnos.raster")
gridded(dnos_raster) <- FALSE
dnos_raster <- coordinates(dnos_raster)
colnames(dnos_raster) <- c("x", "y")
covars <- unique(c(trend.terms(clay_lmm), trend.terms(bude_lmm)))
covars <- na.omit(spgrass6::readRAST6(vname = covars)@data)
set.seed(2001)
id <- sample(1:nrow(dnos_raster), 100)
val_data <- cbind(id, dnos_raster[id, ], covars[id, ])
rm(id, dnos_raster)

i <- j <- k <- 1

tmp <- matrix(nrow = 200, ncol = 2)
l <- 1

# CLAY - nonlinear spatial trend (DIST)
for (i in 1:length(seed)) {
  file <- paste("data/R/2ndArticle_clay_rkm_sim_", seed[i], ".rda", sep = "")
  load(file)
  for (j in 1:length(vars)) {
    for (k in 1:length(pts)) {
      for (l in 1:ncol(clay_rkm_sim)) {
        
        print(l)
        
        # calibration data
        y <- clay_rkm_sim[osc_dist[[j]][[k]][, "id"], l]
        # y <- clay_rkm_sim[osc_ppl[[k]][, "id"], l]
        if (anyNegative(y)) y[whichNegative(y)] <- 0
        x <- covars[osc_dist[[j]][[k]][, "id"], trend.terms(clay_lmm)]
        coords <- osc_dist[[j]][[k]][, 2:3]
        # x <- covars[osc_ppl[[k]][, "id"], trend.terms(clay_lmm)]
        # coords <- osc_ppl[[k]][, 2:3]
        
        # regression-kriging model
        set.seed(2001)
        rkm <- randomForest(x = x, y = y, mtry = 3, ntree = 500)
        res <- rkm$predicted - y
        rkm <- predict(rkm, newdata = val_data[, trend.terms(clay_lmm)])
        
        # We set the initial values for the covariance parameters using the
        # calibration data and the knowledge of the soil data generating 
        # process. The initial values of tausq and sigmasq are set using the
        # total variance of the random forest residuals. According to the soil
        # data generating process, tausq and sigmasq comprise 15% and 85% of 
        # the total variance of the random forest residuals, respectively. The
        # initial value of range (phi) is set computing a variogram of the 
        # random forest residuals up to the 

        # v <- variog(geodata = list(coords = coords, data = res), 
                    # breaks = vgmLags(coords), messages = FALSE)
        # plot(v, ylim = c(0, 8000))

        icp <- vgmICP(z = res, coords = coords, scv = vgmSCV(clay_rkm$clay_vgm),
                      by = list(tausq = 300, sigmasq = 200, phi = 200))
        vgm <- likfit(geodata = list(coords = coords, data = res), 
                      lik.method = "REML", nugget = icp$tausq,
                      ini.cov.pars = icp$cov.pars, messages = FALSE)
        
        # print(vgmSCV(vgm))
        # lines(vgm)
        # abline(v = vgm$phi, col = "blue", lwd = 2)
        # abline(v = unique(icp$cov.pars[, "phi"]), col = "red")
        # abline(h = icp$tausq, col = "red")
        # abline(h = unique(icp$cov.pars[, "sigmasq"]) + icp$tausq[2], 
               # col = "green")
        # l <- l + 1
        # Sys.sleep(2)
        
        
        # If the calibration point data is insuficient to estimate the true
        # nugget variance (tausq) and the true spatially correlated variance
        # (sigmasq), geoR::likfit estimates tausq = sigmasq, and a range (phi)
        # equal to zero. This is equivalent to a pure nugget effect model.
        if (vgm$phi == 0) vgm$cov.model <- "pure.nugget"
        
        locations <- data.frame(coords, z = res, x)
        coordinates(locations) <- ~ x + y
        newdata <- val_data[, c("x", "y", trend.terms(clay_lmm))]
        coordinates(newdata) <- ~ x + y
        rkm <- rkm + 
          krige(formula = z ~ 1, locations = locations, newdata = newdata, 
                nmax = 20, model = as.vgm.variomodel(vgm), beta = 0,
                debug.level = 0)$var1.pred
        
        # linear mixed model
        y <- y + 1
        lambda <- powerTransform(y, lower = 0)$lambda
        locations$z <- bcPower(y, lambda)
        res <- lm(update(clay_lmm$trend, z ~ .), data = locations@data)$resid
        
        # v <- variog(geodata = list(coords = coords, data = res), 
                    # breaks = vgmLags(coords), messages = FALSE)
        # plot(v)
        
        icp <- vgmICP(z = res, coords = coords, 
                      scv = vgmSCV(clay_rkm$clay_vgm),
                      by = list(tausq = round(var(res) / 10),
                                sigmasq = round(var(res) / 10), phi = 200))
        lmm <- likfit(geodata = list(coords = coords, data = locations$z,
                                     covariates = x), 
                      lik.method = "REML", nugget = icp$tausq,
                      ini.cov.pars = icp$cov.pars, 
                      trend = clay_lmm$trend, messages = FALSE)
        
        # print(vgmSCV(lmm))
        # lines(lmm)
        # abline(v = lmm$phi, col = "blue", lwd = 2)
        # abline(v = unique(icp$cov.pars[, "phi"]), col = "red")
        # abline(h = icp$tausq, col = "red")
        # abline(h = unique(icp$cov.pars[, "sigmasq"]) + icp$tausq[2], 
               # col = "green")
        # l <- l + 1
        # Sys.sleep(2)
        
        if (lmm$phi == 0) lmm$cov.model <- "pure.nugget"
        lmm <- krige(formula = update(clay_lmm$trend, z ~ .), 
                     locations = locations, newdata = newdata, nmax = 20, 
                     model = as.vgm.variomodel(lmm), beta = lmm$beta,
                     debug.level = 0)@data
        for (i in 1:nrow(lmm)) {
          lmm$var1.pred[i] <- 
            mean(rboxcox(n = 20000, lambda = lambda, lambda2 = 1, verbose = F,
                         mean = lmm$var1.pred[i], sd = lmm$var1.var[i] ^ 0.5))
        }
        if (anyNegative(lmm$var1.pred)) {
          lmm$var1.pred[whichNegative(lmm$var1.pred)] <- 0
        }

#         obs <- clay_rkm_sim[val_data$id, l]
#         if (anyNegative(obs)) {
#           obs[whichNegative(obs)] <- 0
#         }
#         
#         tmp[l, ] <- c(round(mean(rkm - obs)), round(mean(lmm$var1.pred - obs)))
#         
      }
    }
  }
}

require(car)


#  We fix kappa because it is often poorly identified, and because the maximum
#  likelihood estimators for phi (range) and kappa tend to be strongly 
#  correlated (Diggle & Ribeiro Jr., 2007).
# , divided by a constant c = 100. The
# division by a constant is performed to scale the data values to help in the
# convergence of numerical optimization algorithms (this is a recommendation of
# the authors of the geoR package).
#  Despite the
# empirical distribution is close to the normal, we transform the residuals 
# using the Box-Cox family of power transformations to guarantee that we use 
# the same approach used when fitting the linear mixed models. Because the 
# Box-Cox transformation is not defined for negative data, the residuals are
# linearly transformed to a positive-valued variable using 
# x* = x + abs(min(x)) + 1.
