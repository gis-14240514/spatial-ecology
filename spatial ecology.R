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