# Federal Rural University of Rio de Janeiro
# Graduate School in Agronomy - Soil Science
# Project: Contribution to the Construction of Models to Predict Soil Properties
# Author: Alessandro Samuel-Rosa
# Wageningen, August of 2014
#
# Supervisory committee:
# Lúcia Helena Cunha dos Anjos (UFRRJ)
# Gustavo de Mattos Vasques (Embrapa)
# Gerard B M Heuvelink (ISRIC)
#
# Working script: sample size and sampling design
#
# e-mail: alessandrosamuel@yahoo.com.br
# homepage: soil-scientist.net
#
# SETTINGS #####################################################################
rm(list = ls())
gc()
options(device = x11, stringsAsFactors = FALSE)
require(sp)
require(rgdal)
require(maptools)
require(raster)
require(spgrass6)
require(lattice)
require(latticeExtra)
require(grid)
require(spcosa)
require(clhs)
require(ggplot2)
require(geoR)
require(RandomFields)
require(pedometrics)
source("/home/alessandro/PROJECTS/pedometrics/pedometrics/cooking/readRAST.R")
source("/home/alessandro/PROJECTS/pedometrics/pedometrics/cooking/stratify.R")
source("/home/alessandro/PROJECTS/pedometrics/pedometrics/cooking/form2vect.R")
source("/home/alessandro/PROJECTS/pedometrics/pedometrics/cooking/spSample.R")
load("sm-dnos-general.RData")
load("sm-dnos-phd-chap2.RData")
# load("sm-dnos-phd-chap1.RData")
# load("sm-dnos-phd-chap1-final-models.RData")
# sim_predictors <- preds$main[[32]]
# sim_predictors <- strsplit(sim_predictors, " + ", fixed = TRUE)[[1]]
# sim_models <- list()
# sim_models$clay <- clay_best_lmm
# sim_models$soc <- carbon_best_lmm
# sim_models$ecec <- ecec_best_lmm
initGRASS(gisBase = "/usr/lib/grass64/", gisDbase = GRASSgisDbase,
          location = "dnos-sm-rs", mapset = "predictions",
          pid = Sys.getpid(), override = TRUE)
system("g.region rast=dnos.raster")
system("g.remove MASK")
#
# SAMPLE #######################################################################
sim_sample               <- list()
sim_sample$size          <- c(25, 50, 100, 200, 400, 800)
sim_sample$spcosa_strata <- list()
sim_sample$spcosa        <- list()
sim_sample$cLHS          <- list()
sim_sample$spcLHS        <- list()
sim_sample$spcosa_vario  <- list()
sim_sample$cLHS_vario    <- list()
sim_sample$spcLHS_vario  <- list()
sim_sample$validation    <- list()
#
# SAMPLE - spcosa --------------------------------------------------------------
region <- readVECT6("buffer_BASIN_10")
set.seed(2014)
for (i in 1:length(sim_sample$size)) {
  sim_sample$spcosa_strata[i] <- stratify(region, nTry = 10,
                                          nStrata = sim_sample$size[i],
                                          maxIterations = 10000)
}
sim_sample$spcosa <- lapply(sim_sample$spcosa_strata, spSample)
rm(region)
gc()
# 
# SAMPLE - cLHS ----------------------------------------------------------------
system("r.mask buffer_BASIN_10")
maps <- readRAST(sim_predictors[15:20], out = "RasterStack")
set.seed(2014)
for (i in length(sim_sample$size)) {
  sim_sample$cLHS[i] <- spSample(maps, n = sim_sample$size[i], type = "cLHS")
}
rm(maps)
gc()
system("g.remove MASK")
# 
# SAMPLE - spcLHS --------------------------------------------------------------
#
# SAMPLE - spcosa_vario --------------------------------------------------------
#
# SAMPLE - cLHS_vario ----------------------------------------------------------
#
# SAMPLE - spcLHS_vario --------------------------------------------------------
#
# SAMPLE - validation ----------------------------------------------------------
region <- readVECT6("buffer_BASIN_10")
set.seed(2014)
sim_sample$validation <- spSample(region, n = 1000, type = "random")
plot(sim_sample$validation, pch = 20, cex = 0.3, col = "red")
plot(region, add = TRUE)
rm(region)
gc()
# 
# SIMULATE REALITIES ###########################################################
#
# SIMULATE REALITIES - clay ----------------------------------------------------
clay_sim_sample              <- list()
clay_sim_sample$spcosa       <- list()
clay_sim_sample$cLHS         <- list()
clay_sim_sample$spcLHS       <- list()
clay_sim_sample$spcosa_vario <- list()
clay_sim_sample$cLHS_vario   <- list()
clay_sim_sample$spcLHS_vario <- list()
clay_sim_sample$validation   <- list()
#
trend <- form2vect(sim_models$clay$trend)
beta <- sim_models$clay$beta
names(beta) <- c("Intercept", trend)
trend <- readRAST(trend)
trend <- over(x = sim_sample$validation, y = trend)
trend <- as.matrix(cbind(rep(1, length(sim_sample$validation)), trend))
trend <- trend %*% beta
#
set.seed(2014)
clay_sim_sample$validation <- grf(n = length(sim_sample$validation), 
                                  nsim = 1, mean = trend,
                                  grid = coordinates(sim_sample$validation),
                                  cov.model = sim_models$clay$cov.model, 
                                  cov.pars = c(sim_models$clay$sigmasq, 
                                               sim_models$clay$phi),
                                  nugget = sim_models$clay$nugget, 
                                  lambda = sim_models$clay$lambda)
summary(clay_sim_sample$validation)


#
# SIMULATE REALITIES - soc -----------------------------------------------------
soc_sim_sample              <- list()
soc_sim_sample$spcosa       <- list()
soc_sim_sample$cLHS         <- list()
soc_sim_sample$spcLHS       <- list()
soc_sim_sample$spcosa_vario <- list()
soc_sim_sample$cLHS_vario   <- list()
soc_sim_sample$spcLHS_vario <- list()
soc_sim_sample$validation   <- list()
#
trend <- form2vect(sim_models$soc$trend)
beta <- sim_models$soc$beta
names(beta) <- c("Intercept", trend)
trend <- readRAST(trend)
trend <- over(x = sim_sample$validation, y = trend)
trend <- as.matrix(cbind(rep(1, length(sim_sample$validation)), trend))
trend <- trend %*% beta
#
# SIMULATE REALITIES - ecec ----------------------------------------------------
ecec_sim_sample              <- list()
ecec_sim_sample$spcosa       <- list()
ecec_sim_sample$cLHS         <- list()
ecec_sim_sample$spcLHS       <- list()
ecec_sim_sample$spcosa_vario <- list()
ecec_sim_sample$cLHS_vario   <- list()
ecec_sim_sample$spcLHS_vario <- list()
ecec_sim_sample$validation   <- list()
#
trend <- form2vect(sim_models$ecec$trend)
beta <- sim_models$ecec$beta
names(beta) <- c("Intercept", trend)
trend <- readRAST(trend)
trend <- over(x = sim_sample$validation, y = trend)
trend <- as.matrix(cbind(rep(1, length(sim_sample$validation)), trend))
trend <- trend %*% beta
#
# TABLES #######################################################################
#
# SAVE #########################################################################
save(sim_predictors, sim_models, sim_sample, ecec_sim_sample, soc_sim_sample,
     clay_sim_sample,
     file = "sm-dnos-phd-chap2.RData")
# End!