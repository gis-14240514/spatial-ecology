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

set.seed(11)

# load packages
pkgs <- c("terra", "sf", "pROC")
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install) > 0) install.packages(to_install, dependencies = TRUE)

suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(pROC)
})

# file paths
data_dir <- normalizePath("E:/spatial ecology", winslash = "/", mustWork = TRUE)
out_dir  <- file.path(data_dir, "outputs_week2_800m")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

dir.create(file.path(data_dir, "terra_tmp"), showWarnings = FALSE)
terra::terraOptions(tempdir = file.path(data_dir, "terra_tmp"), memfrac = 0.6)

lcm_file  <- file.path(data_dir, "LCMUK.tif")
dem_file  <- file.path(data_dir, "demScotland.tif")
scot_file <- file.path(data_dir, "scotSamp.shp")
occ_file  <- file.path(data_dir, "Melesmeles.csv")

stopifnot(file.exists(lcm_file), file.exists(dem_file),
          file.exists(scot_file), file.exists(occ_file))

# fixed radii
radius_broadleaf <- 800
radius_urban <- 2300

# read study area
scot <- sf::st_read(scot_file, quiet = TRUE)
scot_buf <- sf::st_buffer(scot, dist = 1000)

# read rasters
LCM <- terra::rast(lcm_file)
if (terra::nlyr(LCM) > 1) LCM <- LCM[[1]]
DEM <- terra::rast(dem_file)

LCM <- terra::mask(terra::crop(LCM, terra::vect(scot_buf)), terra::vect(scot))
DEM <- terra::mask(terra::crop(DEM, terra::vect(scot_buf)), terra::vect(scot))

# keep 100 m version
LCM100 <- terra::aggregate(LCM, fact = 4, fun = "modal")
LCM_fac <- terra::as.factor(LCM100)

# this version used class 2 for broadleaf
broadleaf_code <- 2

lv <- levels(LCM_fac)[[1]]
codes <- lv[, 1]

broadleaf <- terra::classify(
  LCM_fac,
  apply(cbind(codes, ifelse(codes == broadleaf_code, 1, 0)), 2, as.numeric)
)

urb <- rep(0, length(codes))
if (length(codes) >= 2) urb[(length(codes) - 1):length(codes)] <- 1

urban <- terra::classify(
  LCM_fac,
  apply(cbind(codes, urb), 2, as.numeric)
)

# make circular weights
makeWeights <- function(radius_m, template_rast) {
  cell <- terra::res(template_rast)[1]
  nPix <- round(radius_m / cell)
  nPix <- (nPix * 2) + 1
  wm <- matrix(1:nPix^2, nrow = nPix, ncol = nPix)
  
  cx <- ceiling(ncol(wm) / 2)
  cy <- ceiling(nrow(wm) / 2)
  focalCell <- wm[cx, cy]
  indFocal <- which(wm == focalCell, arr.ind = TRUE)
  
  d <- vector("list", length = nPix^2)
  for (i in 1:(nPix^2)) {
    ind <- which(wm == i, arr.ind = TRUE)
    dx <- abs(ind[1, 1] - indFocal[1, 1]) * cell
    dy <- abs(ind[1, 2] - indFocal[1, 2]) * cell
    d[[i]] <- sqrt(dx^2 + dy^2)
  }
  
  wm[] <- unlist(d)
  wm[wm > radius_m] <- NA
  
  wm_norm <- wm
  wm_norm[!is.na(wm_norm)] <- 1 / length(wm_norm[!is.na(wm_norm)])
  wm_norm
}

w_bl  <- makeWeights(radius_broadleaf, broadleaf)
w_urb <- makeWeights(radius_urban, broadleaf)

# convert to % cover
bl_pct  <- 100 * terra::focal(broadleaf, w = w_bl,  fun = "sum", na.rm = TRUE)
urb_pct <- 100 * terra::focal(urban,     w = w_urb, fun = "sum", na.rm = TRUE)
elevR   <- terra::resample(DEM, bl_pct)

covs <- c(bl_pct, urb_pct, elevR)
names(covs) <- c("bl_pct", "urb_pct", "elev")
covs <- terra::mask(terra::crop(covs, terra::vect(scot)), terra::vect(scot))

terra::writeRaster(covs, file.path(out_dir, "P1_covariates_week2_800m.tif"), overwrite = TRUE)

# read badger records
meles <- read.csv(occ_file, stringsAsFactors = FALSE)

meles$Longitude <- as.numeric(trimws(as.character(meles$Longitude)))
meles$Latitude  <- as.numeric(trimws(as.character(meles$Latitude)))
meles <- meles[!is.na(meles$Longitude) & !is.na(meles$Latitude), ]

if ("Coordinate" %in% names(meles)) {
  meles$Coordinate <- as.numeric(trimws(as.character(meles$Coordinate)))
  meles <- meles[is.na(meles$Coordinate) | meles$Coordinate < 1000, ]
}

if ("Identification" %in% names(meles)) {
  meles <- meles[is.na(meles$Identification) | meles$Identification != "Unconfirmed", ]
}

pres_ll <- sf::st_as_sf(meles, coords = c("Longitude", "Latitude"), crs = 4326)
pres_sf <- sf::st_transform(pres_ll, terra::crs(covs))
pres_sf <- pres_sf[scot, ]

stopifnot(nrow(pres_sf) >= 5)

# make background points
bg_n <- max(2000, 5 * nrow(pres_sf))
bg_v <- terra::spatSample(covs[[1]], size = bg_n, as.points = TRUE,
                          na.rm = TRUE, method = "random")

pres_v <- terra::vect(pres_sf)

# extract values at presence and background points
pres_dat <- terra::extract(covs, pres_v)
bg_dat   <- terra::extract(covs, bg_v)

pres_dat$Pres <- 1
bg_dat$Pres   <- 0

# keep coordinates for spatial CV
pres_xy <- as.data.frame(terra::crds(pres_v))
bg_xy   <- as.data.frame(terra::crds(bg_v))
names(pres_xy) <- c("x", "y")
names(bg_xy)   <- c("x", "y")

pres_dat <- cbind(pres_dat, pres_xy)
bg_dat   <- cbind(bg_dat, bg_xy)

dat <- rbind(pres_dat, bg_dat)
dat <- dat[, c("bl_pct", "urb_pct", "elev", "Pres", "x", "y")]
dat <- dat[complete.cases(dat), ]

write.csv(dat, file.path(out_dir, "P2_training_table.csv"), row.names = FALSE)

# full GLM
glm_mcl <- glm(Pres ~ bl_pct + urb_pct + elev, family = "binomial", data = dat)
capture.output(summary(glm_mcl), file = file.path(out_dir, "P3_GLM_summary.txt"))

# coefficients as OR
coefs <- summary(glm_mcl)$coefficients
ci <- suppressMessages(confint(glm_mcl))

coef_tbl <- data.frame(
  term = rownames(coefs),
  estimate = coefs[, "Estimate"],
  se = coefs[, "Std. Error"],
  z = coefs[, "z value"],
  p = coefs[, "Pr(>|z|)"],
  ci_low = ci[, 1],
  ci_high = ci[, 2]
)

coef_tbl$OR <- exp(coef_tbl$estimate)
coef_tbl$OR_low <- exp(coef_tbl$ci_low)
coef_tbl$OR_high <- exp(coef_tbl$ci_high)

write.csv(coef_tbl, file.path(out_dir, "P3b_coefficients_OR_CI.csv"), row.names = FALSE)

pt <- coef_tbl[coef_tbl$term != "(Intercept)", , drop = FALSE]

png(file.path(out_dir, "P3c_coef_forest.png"), width = 1000, height = 650, res = 150)
par(mar = c(5, 10, 4, 2))
y <- seq_len(nrow(pt))
plot(pt$OR, y, log = "x", pch = 19, yaxt = "n",
     xlab = "Odds Ratio (log scale)", ylab = "",
     main = "GLM coefficients (Odds Ratios with 95% CI)")
axis(2, at = y, labels = pt$term, las = 1)
segments(pt$OR_low, y, pt$OR_high, y)
abline(v = 1, lty = 2)
dev.off()

# full ROC
dat$pred <- predict(glm_mcl, type = "response")

roc_obj <- pROC::roc(dat$Pres, dat$pred, quiet = TRUE)
auc_val <- as.numeric(pROC::auc(roc_obj))
write.csv(data.frame(AUC = auc_val), file.path(out_dir, "P4_AUC.csv"), row.names = FALSE)

png(file.path(out_dir, "P5_ROC.png"), width = 900, height = 650, res = 150)
plot(roc_obj, main = paste0("ROC (GLM, 800 m) — AUC = ", round(auc_val, 3)))
dev.off()

# best threshold by Youden
best <- pROC::coords(
  roc_obj,
  x = "best",
  best.method = "youden",
  ret = c("threshold", "sensitivity", "specificity"),
  transpose = FALSE
)
thr_best <- as.numeric(best["threshold"])

write.csv(data.frame(best), file.path(out_dir, "P5b_best_threshold_youden.csv"), row.names = FALSE)

# confusion matrices
confusion_at <- function(thr, obs, pred) {
  pred_class <- ifelse(pred >= thr, 1, 0)
  TP <- sum(pred_class == 1 & obs == 1)
  TN <- sum(pred_class == 0 & obs == 0)
  FP <- sum(pred_class == 1 & obs == 0)
  FN <- sum(pred_class == 0 & obs == 1)
  acc <- (TP + TN) / (TP + TN + FP + FN)
  data.frame(threshold = thr, TP = TP, TN = TN, FP = FP, FN = FN, Accuracy = acc)
}

write.csv(confusion_at(0.5, dat$Pres, dat$pred),
          file.path(out_dir, "P6_confusion_threshold_0.5.csv"),
          row.names = FALSE)

write.csv(confusion_at(thr_best, dat$Pres, dat$pred),
          file.path(out_dir, "P6b_confusion_threshold_youden.csv"),
          row.names = FALSE)

# random 5-fold CV
k <- 5
set.seed(11)
fold_id_random <- sample(rep(1:k, length.out = nrow(dat)))

cv_random <- data.frame(
  fold = 1:k,
  AUC = NA_real_,
  threshold = NA_real_,
  sensitivity = NA_real_,
  specificity = NA_real_
)

for (i in 1:k) {
  train_dat <- dat[fold_id_random != i, ]
  test_dat  <- dat[fold_id_random == i, ]
  
  m_i <- glm(Pres ~ bl_pct + urb_pct + elev, family = "binomial", data = train_dat)
  pred_i <- predict(m_i, newdata = test_dat, type = "response")
  
  roc_i <- pROC::roc(test_dat$Pres, pred_i, quiet = TRUE)
  auc_i <- as.numeric(pROC::auc(roc_i))
  
  best_i <- pROC::coords(
    roc_i,
    x = "best",
    best.method = "youden",
    ret = c("threshold", "sensitivity", "specificity"),
    transpose = FALSE
  )
  
  cv_random$AUC[i] <- auc_i
  cv_random$threshold[i] <- as.numeric(best_i["threshold"])
  cv_random$sensitivity[i] <- as.numeric(best_i["sensitivity"])
  cv_random$specificity[i] <- as.numeric(best_i["specificity"])
}

write.csv(cv_random, file.path(out_dir, "P7_random_cv_results.csv"), row.names = FALSE)
write.csv(
  data.frame(
    mean_AUC = mean(cv_random$AUC, na.rm = TRUE),
    sd_AUC = sd(cv_random$AUC, na.rm = TRUE)
  ),
  file.path(out_dir, "P7b_random_cv_summary.csv"),
  row.names = FALSE
)

# spatial 5-fold CV
set.seed(11)
km <- kmeans(scale(dat[, c("x", "y")]), centers = 5, nstart = 50)
dat$sp_fold <- km$cluster

cv_spatial <- data.frame(
  fold = 1:k,
  AUC = NA_real_,
  threshold = NA_real_,
  sensitivity = NA_real_,
  specificity = NA_real_
)

for (i in 1:k) {
  train_dat <- dat[dat$sp_fold != i, ]
  test_dat  <- dat[dat$sp_fold == i, ]
  
  if (length(unique(test_dat$Pres)) < 2) next
  
  m_i <- glm(Pres ~ bl_pct + urb_pct + elev, family = "binomial", data = train_dat)
  pred_i <- predict(m_i, newdata = test_dat, type = "response")
  
  roc_i <- pROC::roc(test_dat$Pres, pred_i, quiet = TRUE)
  auc_i <- as.numeric(pROC::auc(roc_i))
  
  best_i <- pROC::coords(
    roc_i,
    x = "best",
    best.method = "youden",
    ret = c("threshold", "sensitivity", "specificity"),
    transpose = FALSE
  )
  
  cv_spatial$AUC[i] <- auc_i
  cv_spatial$threshold[i] <- as.numeric(best_i["threshold"])
  cv_spatial$sensitivity[i] <- as.numeric(best_i["sensitivity"])
  cv_spatial$specificity[i] <- as.numeric(best_i["specificity"])
}

write.csv(cv_spatial, file.path(out_dir, "P8_spatial_cv_results.csv"), row.names = FALSE)
write.csv(
  data.frame(
    mean_AUC = mean(cv_spatial$AUC, na.rm = TRUE),
    sd_AUC = sd(cv_spatial$AUC, na.rm = TRUE)
  ),
  file.path(out_dir, "P8b_spatial_cv_summary.csv"),
  row.names = FALSE
)

fold_sf <- sf::st_as_sf(dat, coords = c("x", "y"), crs = terra::crs(covs))

png(file.path(out_dir, "P8c_spatial_folds.png"), width = 900, height = 700, res = 150)
plot(sf::st_geometry(scot), col = "grey95", border = "grey50",
     main = "Spatial cross-validation folds (800 m GLM)")
cols <- c("red", "blue", "green4", "orange", "purple")
for (i in 1:k) {
  plot(sf::st_geometry(fold_sf[dat$sp_fold == i, ]),
       add = TRUE, col = cols[i], pch = 16, cex = 0.35)
}
legend("topright", legend = paste("Fold", 1:k), col = cols, pch = 16, bty = "n")
dev.off()

# partial response curves
base <- data.frame(
  bl_pct = median(dat$bl_pct, na.rm = TRUE),
  urb_pct = median(dat$urb_pct, na.rm = TRUE),
  elev = median(dat$elev, na.rm = TRUE)
)

make_response <- function(varname, n = 200) {
  rng <- range(dat[[varname]], na.rm = TRUE)
  xseq <- seq(rng[1], rng[2], length.out = n)
  newd <- base[rep(1, n), , drop = FALSE]
  newd[[varname]] <- xseq
  p <- predict(glm_mcl, newdata = newd, type = "response")
  data.frame(x = xseq, p = p, var = varname)
}

resp_bl  <- make_response("bl_pct")
resp_urb <- make_response("urb_pct")
resp_e   <- make_response("elev")

resp_all <- rbind(resp_bl, resp_urb, resp_e)
write.csv(resp_all, file.path(out_dir, "P9_partial_response_data.csv"), row.names = FALSE)

png(file.path(out_dir, "P9b_response_bl_pct.png"), width = 900, height = 650, res = 150)
plot(resp_bl$x, resp_bl$p, type = "l",
     xlab = "Broadleaf (%)", ylab = "Predicted probability",
     main = "Partial response: Broadleaf (%)")
dev.off()

png(file.path(out_dir, "P9c_response_urb_pct.png"), width = 900, height = 650, res = 150)
plot(resp_urb$x, resp_urb$p, type = "l",
     xlab = "Urban (%)", ylab = "Predicted probability",
     main = "Partial response: Urban (%)")
dev.off()

png(file.path(out_dir, "P9d_response_elev.png"), width = 900, height = 650, res = 150)
plot(resp_e$x, resp_e$p, type = "l",
     xlab = "Elevation", ylab = "Predicted probability",
     main = "Partial response: Elevation")
dev.off()

# probability map
pred_prob <- terra::predict(covs, glm_mcl, type = "response", na.rm = TRUE)
terra::writeRaster(pred_prob, file.path(out_dir, "P10_pred_prob.tif"), overwrite = TRUE)

png(file.path(out_dir, "P11_pred_prob.png"), width = 900, height = 650, res = 150)
plot(pred_prob, main = "Predicted probability (GLM, 800 m)")
dev.off()

capture.output(sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))
message("Week2 GLM with CV (800 m) completed. Outputs in: ", out_dir)

set.seed(11)

# make sure spatstat loads from user library
user_lib <- file.path(Sys.getenv("USERPROFILE"), "Documents", "R", "win-library",
                      paste0(R.version$major, ".", R.version$minor))
dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(user_lib, setdiff(.libPaths(), user_lib)))

# load packages
pkgs <- c("terra", "sf", "raster",
          "spatstat.geom", "spatstat.model", "spatstat.explore")
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install) > 0) install.packages(to_install, dependencies = TRUE)

suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(raster)
  library(spatstat.geom)
  library(spatstat.model)
  library(spatstat.explore)
})

stopifnot(requireNamespace("spatstat.utils", quietly = TRUE))

# file paths
data_dir <- normalizePath("E:/spatial ecology", winslash = "/", mustWork = TRUE)
out_dir  <- file.path(data_dir, "outputs_ppm_800m")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

dir.create(file.path(data_dir, "terra_tmp"), showWarnings = FALSE)
terra::terraOptions(tempdir = file.path(data_dir, "terra_tmp"), memfrac = 0.6)

lcm_file  <- file.path(data_dir, "LCMUK.tif")
dem_file  <- file.path(data_dir, "demScotland.tif")
occ_file  <- file.path(data_dir, "Melesmeles.csv")
scot_file <- file.path(data_dir, "scotSamp.shp")

stopifnot(file.exists(lcm_file), file.exists(dem_file),
          file.exists(occ_file), file.exists(scot_file))

# fixed radii
radius_broadleaf <- 800
radius_urban <- 2300

# read study area
scot <- sf::st_read(scot_file, quiet = TRUE)
if (is.na(sf::st_crs(scot))) stop("scotSamp.shp has no CRS.")
scot_ll <- sf::st_transform(scot, 4326)

# read badger records
meles <- read.csv(occ_file, stringsAsFactors = FALSE)

meles$Longitude <- as.numeric(trimws(as.character(meles$Longitude)))
meles$Latitude  <- as.numeric(trimws(as.character(meles$Latitude)))
meles <- meles[!is.na(meles$Longitude) & !is.na(meles$Latitude), ]

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

if ("Identification" %in% names(meles)) {
  meles <- meles[is.na(meles$Identification) | meles$Identification != "Unconfirmed", ]
}

if ("Identification.verification.status" %in% names(meles)) {
  meles <- meles[
    is.na(meles$Identification.verification.status) |
      tolower(trimws(meles$Identification.verification.status)) != "unconfirmed", ]
}

if (nrow(meles) == 0) stop("No badger records left after cleaning.")

meles_ll <- sf::st_as_sf(meles, coords = c("Longitude", "Latitude"), crs = 4326)
meles_ll <- suppressWarnings(sf::st_intersection(meles_ll, scot_ll))

if (nrow(meles_ll) < 5) stop("Too few badger records inside scotSamp.shp.")

# read rasters
LCM <- terra::rast(lcm_file)
if (terra::nlyr(LCM) > 1) LCM <- LCM[[1]]
DEM <- terra::rast(dem_file)

scot_proj  <- sf::st_transform(scot_ll, terra::crs(LCM))
meles_proj <- sf::st_transform(meles_ll, terra::crs(LCM))

LCM_crop <- terra::mask(
  terra::crop(LCM, terra::vect(scot_proj)),
  terra::vect(scot_proj)
)

DEM_crop <- terra::mask(
  terra::crop(DEM, terra::vect(scot_proj)),
  terra::vect(scot_proj)
)

# this version used class 1 for broadleaf
LCM_fac <- terra::as.factor(LCM_crop)
lv <- levels(LCM_fac)[[1]]
codes <- lv[, 1]

broadleaf <- terra::classify(
  LCM_fac,
  apply(cbind(codes, ifelse(codes == 1, 1, 0)), 2, as.numeric)
)

urb <- rep(0, length(codes))
if (length(codes) >= 2) urb[(length(codes) - 1):length(codes)] <- 1

urban <- terra::classify(
  LCM_fac,
  apply(cbind(codes, urb), 2, as.numeric)
)

# make circular weights
makeWeights <- function(radius_m, template_rast) {
  cell <- terra::res(template_rast)[1]
  nPix <- round(radius_m / cell)
  nPix <- (nPix * 2) + 1
  wm <- matrix(1:nPix^2, nrow = nPix, ncol = nPix)
  
  cx <- ceiling(ncol(wm) / 2)
  cy <- ceiling(nrow(wm) / 2)
  focalCell <- wm[cx, cy]
  indFocal <- which(wm == focalCell, arr.ind = TRUE)
  
  d <- vector("list", length = nPix^2)
  for (i in 1:(nPix^2)) {
    ind <- which(wm == i, arr.ind = TRUE)
    dx <- abs(ind[1, 1] - indFocal[1, 1]) * cell
    dy <- abs(ind[1, 2] - indFocal[1, 2]) * cell
    d[[i]] <- sqrt(dx^2 + dy^2)
  }
  
  wm[] <- unlist(d)
  wm[wm > radius_m] <- NA
  
  wm_norm <- wm
  wm_norm[!is.na(wm_norm)] <- 1 / length(wm_norm[!is.na(wm_norm)])
  wm_norm
}

broadleaf_opt <- terra::focal(
  broadleaf,
  w = makeWeights(radius_broadleaf, broadleaf),
  fun = "sum",
  na.rm = TRUE
)

urban_2300 <- terra::focal(
  urban,
  w = makeWeights(radius_urban, broadleaf),
  fun = "sum",
  na.rm = TRUE
)

DEM_res <- terra::resample(DEM_crop, broadleaf_opt)

env_stack <- c(broadleaf_opt, urban_2300, DEM_res)
names(env_stack) <- c("broadleaf_opt", "urban_2300", "elev")

env_stack <- terra::mask(
  terra::crop(env_stack, terra::vect(scot_proj)),
  terra::vect(scot_proj)
)

terra::writeRaster(
  env_stack,
  file.path(out_dir, "P1_covariates_stack.tif"),
  overwrite = TRUE
)

# convert raster to spatstat image
raster.as.im <- function(r) {
  stopifnot(inherits(r, "Raster"))
  ex <- raster::extent(r)
  ncolr <- raster::ncol(r)
  nrowr <- raster::nrow(r)
  rx <- raster::xres(r)
  ry <- raster::yres(r)
  
  xx <- seq(from = raster::xmin(ex) + rx / 2,
            to   = raster::xmax(ex) - rx / 2,
            length.out = ncolr)
  yy <- seq(from = raster::ymin(ex) + ry / 2,
            to   = raster::ymax(ex) - ry / 2,
            length.out = nrowr)
  
  v <- raster::getValues(r)
  mat <- matrix(v, nrow = nrowr, ncol = ncolr, byrow = TRUE)
  mat <- mat[nrowr:1, , drop = FALSE]
  
  spatstat.geom::im(mat, xcol = xx, yrow = yy)
}

broadleafIm <- raster.as.im(raster::raster(env_stack[["broadleaf_opt"]]))
urbanIm     <- raster.as.im(raster::raster(env_stack[["urban_2300"]]))
elevIm      <- raster.as.im(raster::raster(env_stack[["elev"]]))

# build analysis window
valid_mask_r <- !is.na(raster::raster(env_stack[[1]]))
valid_mask_im <- raster.as.im(valid_mask_r)
W <- spatstat.geom::as.owin(valid_mask_im)

# fill edge NA values
broadleafIm$v[!is.finite(broadleafIm$v)] <- 0
urbanIm$v[!is.finite(urbanIm$v)]         <- 0
elev_fill <- median(elevIm$v[is.finite(elevIm$v)], na.rm = TRUE)
elevIm$v[!is.finite(elevIm$v)] <- elev_fill

# make ppp object
xy <- sf::st_coordinates(meles_proj)
pppMeles <- spatstat.geom::ppp(xy[, 1], xy[, 2], window = W)
pppMeles <- spatstat.geom::as.ppp(pppMeles)
pppMeles <- pppMeles[W]

png(file.path(out_dir, "P2_ppp_qc.png"), width = 900, height = 650, res = 150)
plot(pppMeles, main = "Badger ppp — QC")
dev.off()

# rescale to km
pppMeles    <- spatstat.geom::rescale(pppMeles, 1000)
broadleafIm <- spatstat.geom::rescale(broadleafIm, 1000)
urbanIm     <- spatstat.geom::rescale(urbanIm, 1000)
elevIm      <- spatstat.geom::rescale(elevIm, 1000)

# tune quadrature density
ndTry <- seq(100, 1000, by = 100)
aic_tbl <- data.frame(nd = ndTry, AIC = NA_real_)

for (i in seq_along(ndTry)) {
  Q_i <- spatstat.geom::quadscheme(pppMeles, method = "grid", nd = ndTry[i])
  fit_i <- spatstat.model::ppm(Q_i ~ broadleafIm + elevIm + urbanIm)
  aic_tbl$AIC[i] <- AIC(fit_i)
}

write.csv(aic_tbl, file.path(out_dir, "P3_quadrature_AIC.csv"), row.names = FALSE)

best_nd <- aic_tbl$nd[which.min(aic_tbl$AIC)]
Q <- spatstat.geom::quadscheme(pppMeles, method = "grid", nd = best_nd)

png(file.path(out_dir, "P4_quadrature_AIC_plot.png"), width = 900, height = 600, res = 150)
plot(aic_tbl$nd, aic_tbl$AIC, type = "b", pch = 19,
     xlab = "nd (grid quadrature)",
     ylab = "AIC",
     main = "Quadrature tuning by AIC")
abline(v = best_nd, lty = 2)
dev.off()

# fit ppm
ppm_fit <- spatstat.model::ppm(
  Q ~ poly(broadleafIm, 3) + poly(elevIm, 2) + poly(urbanIm, 2) + x + y
)

capture.output(summary(ppm_fit), file = file.path(out_dir, "P5_ppm_summary.txt"))

envK <- spatstat.explore::envelope(
  ppm_fit, spatstat.explore::Kest,
  nsim = 39, VARIANCE = TRUE, nSD = 1, global = TRUE
)

png(file.path(out_dir, "P6_ppm_Kest_envelope.png"), width = 900, height = 650, res = 150)
plot(envK, main = "Kest envelope — ppm")
dev.off()

# fit Thomas kppm
kppm_thomas <- spatstat.model::kppm(
  Q ~ poly(broadleafIm, 3) + poly(elevIm, 2) + poly(urbanIm, 2) + x + y,
  clusters = "Thomas"
)

capture.output(summary(kppm_thomas), file = file.path(out_dir, "P7_kppm_thomas_summary.txt"))

envKinh <- spatstat.explore::envelope(
  kppm_thomas, spatstat.explore::Kinhom,
  nsim = 39, VARIANCE = TRUE, nSD = 1, global = TRUE
)

png(file.path(out_dir, "P8_kppm_Kinhom_envelope.png"), width = 900, height = 650, res = 150)
plot(envKinh, main = "Kinhom envelope — Thomas kppm")
dev.off()

# ROC and intensity map
roc_obj <- spatstat.explore::roc(kppm_thomas)
auc_val <- spatstat.model::auc.kppm(kppm_thomas)

write.csv(
  data.frame(
    AUC_obs = auc_val["obs"],
    AUC_theo = auc_val["theo"]
  ),
  file.path(out_dir, "P9_auc_kppm.csv"),
  row.names = FALSE
)

png(file.path(out_dir, "P10_ROC_curve.png"), width = 900, height = 650, res = 150)
plot(roc_obj, main = "ROC — Thomas kppm")
dev.off()

pred_im <- predict(kppm_thomas)
pred_rast <- terra::rast(pred_im)

terra::writeRaster(
  pred_rast,
  file.path(out_dir, "P11_intensity_thomas.tif"),
  overwrite = TRUE
)

png(file.path(out_dir, "P12_intensity_map.png"), width = 900, height = 650, res = 150)
plot(pred_rast, main = "Predicted intensity — Thomas kppm")
dev.off()

capture.output(sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))
message("DONE. Outputs in: ", out_dir)