# Fix GRASS 64 projection information
wgs1984utm22s <- sp::CRS("+init=epsg:32722")

# Start with PERMANENT mapset
spgrass6::initGRASS(
  gisBase = "/usr/lib/grass64/", gisDbase = path.expand("~/dbGRASS"),location = "dnos-sm-rs", 
  mapset = "PERMANENT", pid = Sys.getpid(), override = TRUE)

# Overwrite projection
system("g.proj --quiet -c epsg=32722")
spgrass6::gmeta()

# Set region parameters
system("g.region rast=dnos.raster@predictions")
spgrass6::gmeta()

# Now with 'predictions' mapset
spgrass6::initGRASS(
  gisBase = "/usr/lib/grass64/", gisDbase = path.expand("~/dbGRASS"),location = "dnos-sm-rs", 
  mapset = "predictions", pid = Sys.getpid(), override = TRUE)
system("g.region rast=dnos.raster")
spgrass6::gmeta()

# Try load file
tmp <- spgrass6::readRAST("ACC_10")

