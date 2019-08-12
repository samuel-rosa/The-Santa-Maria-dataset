library(sf)
library(magrittr)
library(dplyr)

# Point soil observations: soildata and soilsamples
soildata <- febr::febr(
  dataset = "ctb0003", merge = TRUE, variable = "all", 
  standardization = list(repetition =  "combine", crs = "EPSG:4674"))
febr::febr2xlsx(soildata, file = path.expand("~/oCloud/dnos-sm-rs/vector/soildata.xlsx"))
write.table(
  soildata[[2]], file = path.expand("~/oCloud/dnos-sm-rs/vector/soilsamples.csv"), sep = "\tab", 
  row.names = FALSE, dec = ",")

# Initiate GRASS GIS
rgrass7::initGRASS(
  gisBase = "/usr/lib/grass76/", gisDbase = path.expand("~/dbGRASS"), location = "dnos-sm-rs", 
  mapset = "predictions", pid = Sys.getpid(), override = TRUE)

# Base data
basin <- 
  rgrass7::readVECT(vname = "buffer_BASIN_10") %>% 
  sf::st_as_sf()
points <- 
  febr::observation("ctb0003") %>% 
  febr::febr2spdf() %>% 
  sf::st_as_sf() %>% 
  sf::st_transform(crs = 32722)
full_hull <- 
  points %>% 
  sf::st_union() %>% 
  sf::st_convex_hull() %>% 
  sf::st_union(basin)

# fullhull
full_hull %>% 
  sf::st_transform(crs = 4674) %>% 
  sf::st_write(
    dsn = "/home/alessandro/oCloud/dnos-sm-rs/vector/fullhull.geojson", delete_dsn = TRUE,
    layer_options = c("WRITE_BBOX=YES", "DESCRIPTION=Convex hull of the DNOS basin and soil observation points. Coordinate Reference System: EPSG:4674 (SIRGAS 2000)"))

# basin10plus30m.geojson
rgrass7::readVECT(vname = "buffer_BASIN_10") %>% 
  sf::st_as_sf() %>% 
  sf::st_transform(crs = 4674) %>% 
  write_sf(
    dsn = "/home/alessandro/oCloud/dnos-sm-rs/vector/basin10plus30m.geojson", delete_dsn = TRUE,
    layer_options = c("WRITE_BBOX=YES", "DESCRIPTION=Boundary of the DNOS basin computed with r.watershed using a 5-m resolution DEM derived from 10-m equidistant contour lines, later added a buffer of 30 m to compensate for positional errors of the source elevation data. Coordinate Reference System: EPSG:4674 (SIRGAS 2000)"))

# geology25k
rgrass7::readVECT(vname = "GEO_25") %>% 
  sf::st_as_sf() %>% 
  sf::st_intersection(., full_hull) %>% 
  sf::st_transform(crs = 4674) %>% 
  dplyr::select(-cat, -id) %>% 
  sf::write_sf(
    dsn = "/home/alessandro/oCloud/dnos-sm-rs/vector/geology25k.geojson", delete_dsn = TRUE,
    layer_options = c("WRITE_BBOX=YES", "DESCRIPTION=Geological map published at a cartographic scale of 1:25.000. Coordinate Reference System: EPSG:4674 (SIRGAS 2000)"))

# geology50k
rgrass7::readVECT(vname = "GEO_50") %>% 
  sf::st_as_sf() %>% 
  sf::st_intersection(., full_hull) %>% 
  sf::st_transform(crs = 4674) %>% 
  dplyr::select(-cat, -id) %>% 
  sf::write_sf(
    dsn = "/home/alessandro/oCloud/dnos-sm-rs/vector/geology50k.geojson", delete_dsn = TRUE,
    layer_options = c("WRITE_BBOX=YES", "DESCRIPTION=Geological map published at a cartographic scale of 1:50.000. Coordinate Reference System: EPSG:4674 (SIRGAS 2000)"))

# deposits25k
rgrass7::readVECT(vname = "DEP_25") %>% 
  sf::st_as_sf() %>% 
  sf::st_intersection(., full_hull) %>% 
  sf::st_transform(crs = 4674) %>% 
  dplyr::select(-cat, -id) %>% 
  write_sf(
    dsn = "/home/alessandro/oCloud/dnos-sm-rs/vector/deposits25k.geojson", delete_dsn = TRUE,
    layer_options = c("WRITE_BBOX=YES", "DESCRIPTION=Map of quaternary deposits published at a cartographic scale of 1:25.000. Coordinate Reference System: EPSG:4674 (SIRGAS 2000)"))

# landuse1980
tmp <- rgrass7::readVECT(vname = "LU1980")
tmp %>% 
  slot("polygons") %>% 
  lapply(maptools::checkPolygonsHoles) %>% 
  sp::SpatialPolygons(proj4string = sp::CRS("+proj=utm +zone=22 +south +datum=WGS84 +units=m +no_defs")) %>% 
  sp::SpatialPolygonsDataFrame(data = tmp@data) %>% 
  sf::st_as_sf() %>% 
  sf::st_intersection(., full_hull) %>% 
  sf::st_transform(crs = 4674) %>% 
  dplyr::select(-cat, -id) %>% 
  write_sf(
    dsn = "/home/alessandro/oCloud/dnos-sm-rs/vector/landuse1980.geojson", delete_dsn = TRUE,
    layer_options = c("WRITE_BBOX=YES", "DESCRIPTION=Land use map of the year of 1980. Coordinate Reference System: EPSG:4674 (SIRGAS 2000)"))

# landuse2009
tmp <- rgrass7::readVECT(vname = "LU2009", with_c = TRUE)
tmp %>%
  slot("polygons") %>%
  lapply(maptools::checkPolygonsHoles) %>%
  sp::SpatialPolygons(proj4string = sp::CRS("+proj=utm +zone=22 +south +datum=WGS84 +units=m +no_defs")) %>%
  sp::SpatialPolygonsDataFrame(data = tmp@data) %>%
  sf::st_as_sf() %>% 
  dplyr::filter(!is.na(land_use)) %>% 
  sf::st_intersection(., full_hull) %>% 
  sf::st_transform(crs = 4674) %>%
  dplyr::select(-cat, -Id) %>% 
  # dplyr::filter(sf::st_geometry_type(.) == "GEOMETRYCOLLECTION") %>%
  write_sf(
    dsn = "/home/alessandro/oCloud/dnos-sm-rs/vector/landuse2009.geojson", delete_dsn = TRUE,
    layer_options = c("WRITE_BBOX=YES", "DESCRIPTION=Land use map of the year 2009. Coordinate Reference System: EPSG:4674 (SIRGAS 2000)"))

# pedology100k
rgrass7::readVECT(vname = "SOIL_100") %>% 
  sf::st_as_sf() %>% 
  sf::st_intersection(., full_hull) %>% 
  sf::st_transform(crs = 4674) %>% 
  dplyr::select(-cat, -id) %>% 
  write_sf(
    dsn = "/home/alessandro/oCloud/dnos-sm-rs/vector/pedology100k.geojson", delete_dsn = TRUE, 
    layer_options = c("WRITE_BBOX=YES", "DESCRIPTION=Pedological map published at a cartographic scale of 1:100.000. Coordinate Reference System: EPSG:4674 (SIRGAS 2000)"))

# pedology25k
rgrass7::readVECT(vname = "SOIL_25") %>% 
  sf::st_as_sf() %>% 
  sf::st_intersection(., full_hull) %>% 
  sf::st_transform(crs = 4674) %>% 
  dplyr::select(-cat, -Id) %>% 
  write_sf(
    dsn = "/home/alessandro/oCloud/dnos-sm-rs/vector/pedology25k.geojson", delete_dsn = TRUE,
    layer_options = c("WRITE_BBOX=YES", "DESCRIPTION=Pedological map published at a cartographic scale of 1:100.000. Coordinate Reference System: EPSG:4674 (SIRGAS 2000)"))

# topodata
topodata(
  sheet = "29S54_", layer = "all", destfolder = path.expand("~/oCloud/dnos-sm-rs/raster"))
topodata_files <- list.files(path.expand("~/oCloud/dnos-sm-rs/raster"), full.names = TRUE)
projwin <- path.expand("~/oCloud/dnos-sm-rs/vector/fullhull.geojson") %>% sf::st_read() %>% sf::st_bbox()
for (i in topodata_files) {
  gdalUtils::gdal_translate(
    src_dataset = i, 
    dst_dataset = gsub(pattern = ".tif", replacement = "_CUT.tif", i), 
    projwin = projwin[c(1, 4, 3, 2)], co = "COMPRESS=DEFLATE")
  file.remove(i)
  file.rename(
    from = gsub(pattern = ".tif", replacement = "_CUT.tif", i),
    to = gsub("_CUT.tif", ".tif", gsub(pattern = ".tif", replacement = "_CUT.tif", i)))
}
