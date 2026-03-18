set.seed(11)

# load packages
pkgs <- c("terra", "sf")
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install) > 0) install.packages(to_install, dependencies = TRUE)
invisible(lapply(pkgs, library, character.only = TRUE))

# file paths
data_dir <- normalizePath("E:/spatial ecology", winslash = "/", mustWork = TRUE)
out_dir  <- file.path(data_dir, "outputs_week1")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

dir.create(file.path(data_dir, "terra_tmp"), showWarnings = FALSE)
terra::terraOptions(tempdir = file.path(data_dir, "terra_tmp"), memfrac = 0.6)

lcm_file <- file.path(data_dir, "LCMUK.tif")
occ_file <- file.path(data_dir, "Melesmeles.csv")

stopifnot(file.exists(lcm_file), file.exists(occ_file))

# read badger records
meles <- read.csv(occ_file, stringsAsFactors = FALSE)

# clean coordinates
meles$Longitude <- as.numeric(trimws(as.character(meles$Longitude)))
meles$Latitude  <- as.numeric(trimws(as.character(meles$Latitude)))
meles <- meles[!is.na(meles$Longitude) & !is.na(meles$Latitude), ]

# drop very uncertain records
if ("Coordinate" %in% names(meles)) {
  meles$Coordinate <- as.numeric(trimws(as.character(meles$Coordinate)))
  meles <- meles[is.na(meles$Coordinate) | meles$Coordinate < 1000, ]
}

if ("Coordinate.uncertainty.in.metres" %in% names(meles)) {
  meles$Coordinate.uncertainty.in.metres <-
    as.numeric(trimws(as.character(meles$Coordinate.uncertainty.in.metres)))
  meles <- meles[
    is.na(meles$Coordinate.uncertainty.in.metres) |
      meles$Coordinate.uncertainty.in.metres < 1000, ]
}

# drop unconfirmed records
if ("Identification" %in% names(meles)) {
  meles <- meles[is.na(meles$Identification) | meles$Identification != "Unconfirmed", ]
}

if ("Identification.verification.status" %in% names(meles)) {
  meles <- meles[
    is.na(meles$Identification.verification.status) |
      tolower(trimws(meles$Identification.verification.status)) != "unconfirmed", ]
}

if (nrow(meles) == 0) stop("No badger records left after cleaning.")

# make sf points
meles_ll <- sf::st_as_sf(
  meles,
  coords = c("Longitude", "Latitude"),
  crs = 4326
)

# use the same practical extent
studyExtent <- terra::ext(-4.2, -2.7, 56.5, 57.5)
meles_vect_ll <- terra::vect(meles_ll)
meles_crop_ll <- terra::crop(meles_vect_ll, studyExtent)

if (nrow(meles_crop_ll) == 0) {
  stop("No badger points fall inside the study extent.")
}

# read LCM and project points
LCM <- terra::rast(lcm_file)
if (terra::nlyr(LCM) > 1) LCM <- LCM[[1]]

meles_proj <- terra::project(meles_crop_ll, LCM)

# crop LCM around the points
xy <- terra::crds(meles_proj)
stopifnot(nrow(xy) > 0)

extent_new <- terra::ext(
  min(xy[, 1]) - 5000,
  max(xy[, 1]) + 5000,
  min(xy[, 2]) - 5000,
  max(xy[, 2]) + 5000
)

LCM_crop <- terra::crop(LCM, extent_new)

# broadleaf = class 1 in this version
LCM_crop <- terra::as.factor(LCM_crop)
lv <- levels(LCM_crop)[[1]]
codes <- lv[, 1]

RCmatrix <- cbind(codes, ifelse(codes == 1, 1, 0))
RCmatrix <- apply(RCmatrix, 2, as.numeric)

broadleaf <- terra::classify(LCM_crop, RCmatrix)

terra::writeRaster(
  broadleaf,
  file.path(out_dir, "01_broadleaf_binary.tif"),
  overwrite = TRUE
)

# create background points
bg_v <- terra::spatSample(
  broadleaf,
  size = 1000,
  as.points = TRUE,
  method = "random",
  na.rm = TRUE
)

Pres <- data.frame(terra::crds(meles_proj), Pres = 1)
Abs  <- data.frame(terra::crds(bg_v), Pres = 0)
melesData <- rbind(Pres, Abs)

melesSF <- sf::st_as_sf(
  melesData,
  coords = c("x", "y"),
  crs = terra::crs(broadleaf)
)

# function to get broadleaf % in each buffer
cell_area <- prod(terra::res(broadleaf))

landBuffer <- function(points_sf, radius_m, broadleaf_rast) {
  buf <- sf::st_buffer(points_sf, dist = radius_m)
  bl_sum <- terra::extract(broadleaf_rast, terra::vect(buf), fun = "sum", na.rm = TRUE)
  bl_area <- bl_sum[, 2] * cell_area
  buf_area <- as.numeric(sf::st_area(buf))
  (bl_area / buf_area) * 100
}

# test 100 to 2000 m
radii <- seq(100, 2000, by = 100)

pct_list <- vector("list", length(radii))
for (i in seq_along(radii)) {
  r <- radii[i]
  pct_list[[i]] <- landBuffer(melesSF, r, broadleaf)
}

pct_mat <- do.call(cbind, pct_list)
colnames(pct_mat) <- paste0("radius_", radii)

glmData <- as.data.frame(pct_mat)
glmData$Pres <- melesData$Pres
glmData <- glmData[complete.cases(glmData), ]

# fit one GLM for each radius
glmRes <- data.frame(radius = radii, logLik = NA_real_)

for (i in seq_along(radii)) {
  r <- radii[i]
  xname <- paste0("radius_", r)
  
  m <- glm(
    as.formula(paste0("Pres ~ ", xname)),
    family = "binomial",
    data = glmData
  )
  
  glmRes$logLik[i] <- as.numeric(logLik(m))
}

opt_row <- glmRes[which.max(glmRes$logLik), , drop = FALSE]
opt_radius <- opt_row$radius

# save outputs
write.csv(glmRes, file.path(out_dir, "02_logLik_by_radius.csv"), row.names = FALSE)
write.csv(opt_row, file.path(out_dir, "03_optimum_radius.csv"), row.names = FALSE)

png(file.path(out_dir, "04_logLik_curve.png"), width = 900, height = 600, res = 150)
plot(
  glmRes$radius, glmRes$logLik,
  type = "b", pch = 19,
  xlab = "Buffer radius (m)",
  ylab = "log-likelihood",
  main = "Characteristic scale of broadleaf woodland for Meles meles"
)
abline(v = opt_radius, lty = 2, col = "red")
text(opt_radius, max(glmRes$logLik),
     labels = paste0("opt = ", opt_radius, " m"),
     pos = 4)
dev.off()

message("Optimum radius = ", opt_radius, " m")
message("Outputs written to: ", out_dir)