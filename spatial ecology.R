set.seed(11)

suppressPackageStartupMessages({
  library(terra)
  library(sf)
})

# Data folder + fixed filenames
data_dir <- normalizePath("D:/spatial ecology", winslash="/", mustWork=TRUE)
lcm_file  <- file.path(data_dir, "LCMUK.tif")
dem_file  <- file.path(data_dir, "demScotland.tif")
scot_file <- file.path(data_dir, "scotSamp.shp")
occ_file  <- file.path(data_dir, "Melesmeles.csv")
stopifnot(file.exists(lcm_file), file.exists(dem_file), file.exists(scot_file), file.exists(occ_file))

# Output folder
out_dir_w1 <- file.path(data_dir, "outputs_week1")
out_dir_w3 <- file.path(data_dir, "outputs_ppm")
dir.create(out_dir_w1, showWarnings=FALSE, recursive=TRUE)
dir.create(out_dir_w3, showWarnings=FALSE, recursive=TRUE)

# Terra temp (avoid RAM spikes)
dir.create(file.path(data_dir, "terra_tmp"), showWarnings=FALSE)
terra::terraOptions(tempdir=file.path(data_dir, "terra_tmp"), memfrac=0.6)
# Load LCM, project points, crop raster to points (+5 km)
LCM <- terra::rast(lcm_file)
if (terra::nlyr(LCM) > 1) LCM <- LCM[[1]]

meles_proj <- sf::st_transform(meles_ll, terra::crs(LCM))

xy <- sf::st_coordinates(meles_proj)
stopifnot(nrow(xy) > 0)

ext_new <- terra::ext(min(xy[,1])-5000, max(xy[,1])+5000,
                      min(xy[,2])-5000, max(xy[,2])+5000)

LCM_crop <- terra::crop(LCM, ext_new)
# Aggregate LCM to 100 m first (modal)
# LCM is 25 m, fact=4 gives 100 m
LCM100 <- terra::aggregate(LCM_crop, fact = 4, fun = "modal")
LCM_fac <- terra::as.factor(LCM100)
# Reclassify to broadleaf binary raster 
broadleaf_code <- 2  # TODO: change after checking levels(LCM_fac) if needed

lv <- levels(LCM_fac)[[1]]
codes <- lv[,1]
reclass <- ifelse(codes == broadleaf_code, 1, 0)
RC <- apply(cbind(codes, reclass), 2, as.numeric)

broadleaf <- terra::classify(LCM_fac, RC)