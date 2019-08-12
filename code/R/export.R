rgrass7::initGRASS(
  gisBase = "/usr/lib/grass76/", gisDbase = path.expand("~/dbGRASS"), location = "dnos-sm-rs", 
  mapset = "predictions", pid = Sys.getpid(), override = TRUE)

library(sp)
library(sf)
library(magrittr)

basin <- 
  rgrass7::readVECT(vname = "buffer_BASIN_10") %>% 
  sf::st_as_sf()

# basin10plus30m.geojson
rgrass7::readVECT(vname = "buffer_BASIN_10") %>% 
  sf::st_as_sf() %>% 
  sf::st_transform(crs = 4674) %>% 
  write_sf(
    dsn = "/home/alessandro/oCloud/dnos-sm-rs/vector/basin10plus30m.geojson", 
    layer_options = c("WRITE_BBOX=YES", "DESCRIPTION=Boundary of the DNOS basin computed with r.watershed using a 5-m resolution DEM derived from 10-m equidistant contour lines, later added a buffer of 30 m to compensate for positional errors of the source elevation data. Coordinate Reference System: EPSG:4674 (SIRGAS 2000)"))

# geology25k
rgrass7::readVECT(vname = "GEO_25") %>% 
  sf::st_as_sf() %>% 
  sf::st_intersection(basin, .) %>% 
  sf::st_transform(crs = 4674) %>% 
  dplyr::select(-cat, -cat.1, -id) %>% 
  write_sf(
    dsn = "/home/alessandro/oCloud/dnos-sm-rs/vector/geology25k.geojson", 
    layer_options = c("WRITE_BBOX=YES", "DESCRIPTION=Geological map published at a cartographic scale of 1:25.000. Coordinate Reference System: EPSG:4674 (SIRGAS 2000)"))

# geology50k
rgrass7::readVECT(vname = "GEO_50") %>% 
  sf::st_as_sf() %>% 
  sf::st_intersection(basin, .) %>% 
  sf::st_transform(crs = 4674) %>% 
  dplyr::select(-cat, -cat.1, -id) %>% 
  write_sf(
    dsn = "/home/alessandro/oCloud/dnos-sm-rs/vector/geology50k.geojson", 
    layer_options = c("WRITE_BBOX=YES", "DESCRIPTION=Geological map published at a cartographic scale of 1:50.000. Coordinate Reference System: EPSG:4674 (SIRGAS 2000)"))

# deposits25k
rgrass7::readVECT(vname = "DEP_25") %>% 
  sf::st_as_sf() %>% 
  sf::st_intersection(basin, .) %>% 
  sf::st_transform(crs = 4674) %>% 
  dplyr::select(-cat, -cat.1, -id) %>% 
  write_sf(
    dsn = "/home/alessandro/oCloud/dnos-sm-rs/vector/deposits25k.geojson", 
    layer_options = c("WRITE_BBOX=YES", "DESCRIPTION=Map of quaternary deposits published at a cartographic scale of 1:25.000. Coordinate Reference System: EPSG:4674 (SIRGAS 2000)"))

# landuse1980
tmp <- rgrass7::readVECT(vname = "LU1980")
tmp %>% 
  slot("polygons") %>% 
  lapply(maptools::checkPolygonsHoles) %>% 
  sp::SpatialPolygons(proj4string = sp::CRS("+proj=utm +zone=22 +south +datum=WGS84 +units=m +no_defs")) %>% 
  sp::SpatialPolygonsDataFrame(data = tmp@data) %>% 
  sf::st_as_sf() %>% 
  sf::st_intersection(basin, .) %>% 
  sf::st_transform(crs = 4674) %>% 
  dplyr::select(-cat, -cat.1, -id) %>% 
  write_sf(
    dsn = "/home/alessandro/oCloud/dnos-sm-rs/vector/landuse1980.geojson", 
    layer_options = c("WRITE_BBOX=YES", "DESCRIPTION=Land use map of the year of 1980. Coordinate Reference System: EPSG:4674 (SIRGAS 2000)"))

# landuse2009
tmp <- rgrass7::readVECT(vname = "LU2009", with_c = TRUE)
tmp %>% 
  slot("polygons") %>% 
  lapply(maptools::checkPolygonsHoles) %>% 
  sp::SpatialPolygons(proj4string = sp::CRS("+proj=utm +zone=22 +south +datum=WGS84 +units=m +no_defs")) %>% 
  sp::SpatialPolygonsDataFrame(data = tmp@data) %>% 
  sf::st_as_sf() %>% 
  sf::st_intersection(basin, .) %>% 
  sf::st_transform(crs = 4674) %>% 
  dplyr::select(-cat, -cat.1, -Id) %>% 
  write_sf(
    dsn = "/home/alessandro/oCloud/dnos-sm-rs/vector/landuse2009.geojson", 
    layer_options = c("WRITE_BBOX=YES", "DESCRIPTION=Land use map of the year 2009. Coordinate Reference System: EPSG:4674 (SIRGAS 2000)"))

# pedology100k
rgrass7::readVECT(vname = "SOIL_100") %>% 
  sf::st_as_sf() %>% 
  sf::st_intersection(basin, .) %>% 
  sf::st_transform(crs = 4674) %>% 
  dplyr::select(-cat, -cat.1, -id) %>% 
  write_sf(
    dsn = "/home/alessandro/oCloud/dnos-sm-rs/vector/pedology100k.geojson", 
    layer_options = c("WRITE_BBOX=YES", "DESCRIPTION=Pedological map published at a cartographic scale of 1:100.000. Coordinate Reference System: EPSG:4674 (SIRGAS 2000)"))

# pedology25k
rgrass7::readVECT(vname = "SOIL_25") %>% 
  sf::st_as_sf() %>% 
  sf::st_intersection(basin, .) %>% 
  sf::st_transform(crs = 4674) %>% 
  dplyr::select(-cat, -cat.1, -Id) %>% 
  write_sf(
    dsn = "/home/alessandro/oCloud/dnos-sm-rs/vector/pedology25k.geojson", 
    layer_options = c("WRITE_BBOX=YES", "DESCRIPTION=Pedological map published at a cartographic scale of 1:100.000. Coordinate Reference System: EPSG:4674 (SIRGAS 2000)"))

# topodata
topodata(
  sheet = "29S54_", layer = "all", destfolder = path.expand("~/oCloud/dnos-sm-rs/raster"))
topodata_files <- list.files(path.expand("~/oCloud/dnos-sm-rs/raster"), full.names = TRUE)
for (i in topodata_files) {
  gdalUtils::gdalwarp(
    srcfile = i,
    dstfile = gsub(pattern = ".tif", replacement = "_CUT.tif", i),
    s_srs = "EPSG:4674", t_srs = "EPSG:4674",
    crop_to_cutline = TRUE, co = "COMPRESS=DEFLATE",
    wo = c("CUTLINE_ALL_TOUCHED=TRUE", "NUM_THREADS=ALL_CPUS", "OPTIMIZE_SIZE=TRUE"),
    cutline = path.expand("~/oCloud/dnos-sm-rs/vector/basin10plus30m.geojson"))
  file.remove(i)
  file.rename(
    from = gsub(pattern = ".tif", replacement = "_CUT.tif", i),
    to = gsub("_CUT.tif", ".tif", gsub(pattern = ".tif", replacement = "_CUT.tif", i)))
}
# gdalUtils::gdalwarp(srcfile = "/home/alessandro/oCloud/dnos-sm-rs/raster/29S54_DD.tif",
#                     dstfile = "/home/alessandro/oCloud/dnos-sm-rs/raster/29S54_DD_CUT.tif", 
#                     s_srs = "EPSG:4674", t_srs = "EPSG:4674", dstalpha = TRUE,
#                     crop_to_cutline = TRUE, co = "COMPRESS=DEFLATE",
#                     wo = c("CUTLINE_ALL_TOUCHED=TRUE", "NUM_THREADS=ALL_CPUS", "OPTIMIZE_SIZE=TRUE"),
#                     cutline = "/home/alessandro/oCloud/dnos-sm-rs/vector/basin10plus30m.geojson")
# system("g.list vector")
