# Weird low CLAY values #######################################################################################
# 26 March 2016
# Check weird low CLAY values of soils derived from igneous rocks.
soil_data <- cbind(
  read.table("~/projects/dnos-sm-rs/dnos-sm-rs-general/data/labData.csv", sep = ";", header = TRUE, 
             na.strings = "na"),
  read.table("~/projects/dnos-sm-rs/dnos-sm-rs-general/data/fieldData.csv", sep = ";", header = TRUE,
             na.strings = "na"))
id <- which(soil_data$CLAY <= 100 & soil_data$PARENT == "igneous")
lattice::bwplot(
  CLAY ~ LAND | PARENT, data = soil_data[-id, ],
  ylab = expression(paste('Argila (g ', kg^-1, ')', sep = '')),
  strip = lattice::strip.custom(bg = "lightgray"), scales = list(x = list(rot = 45)),
  par.settings = list(
    fontsize = list(text = 14, points = 8), box.rectangle = list(col = "black"), 
    box.umbrella = list(col = "black"), plot.symbol = list(col = "black", cex = 0.7),
    box.dot = list(cex = 0.7),
    layout.widths = list(left.padding = 0, right.padding = 0),
    layout.heights = list(top.padding = 0, bottom.padding = 0)))
write.csv(soil_data[id, ], path.expand("~/tmp/soil_data.csv"))

# Open data in QGIS. Use satellite image to check the location.
# - PabloP4: Observation from Profile 4 of Miguel (2010). The description of this profile has many other erros 
#   that had to be corrected before such as the sum of particle size distribution. This is just another error.



# - sm-dnos-302: SAND is very high. The observation locations looks like a place of accumulation of sediments.
#   I guess that we have made an incorrect identification of the parent material!
# - 28/02/11: Six observations made in the same date, specifically, in the last field campaing (sm-dnos-325, 
#   sm-dnos-327, sm-dnos-328, sm-dnos-330, sm-dnos-331). Some of these might be derived from sediments of 
#   igneous rocks. All have high ORCA: CLAY determination was preceded by H2O2 treatment. SAND is relatively 
#   low. Have we underestimated CLAY due to H2O2? Have we made an incorrect identification of the parent 
#   material?
# - sm-dnos-167: SAND is relatively high. Have we made an incorrect identification of the parent material? The
#   observation locations looks like a place of accumulation of sediments.

