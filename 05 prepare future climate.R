
# Load libraries ----------------------------------------------------------
require(pacman)
p_load(terra, sf, fs, tidyverse, spatialEco, climateStability, parallelDist, glue, outliers, dismo, scales, glue, rnaturalearthdata, rnaturalearth, openxlsx)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# Load data ---------------------------------------------------------------

root <- '//catalogue/workspace-cluster9/2024/WCR_Sep/p1/tif/tc_2c'
fles <- dir_ls(root, regexp = '.nc$')
fles <- as.character(fles)
fles <- grep(paste0(c('ppt', 'tmin', 'tmax'), collapse = '|'), fles, value = T)

mskr <- terra::rast('./tif/pca/pca-bsl_run1.tif')
extn <- ext(mskr)

## Vector data 
wrld <- ne_countries(returnclass = 'sf')
wrld <- vect(wrld)
wrld <- terra::crop(wrld, extn)

# Function ----------------------------------------------------------------
calc.avrg <- function(month){
  
  ## To start
  cat('To process: ', month, '\n')
  
  ## 
  lyer <- map(.x = 1:length(fles), .f = function(i){
    rst <- terra::rast(fles[i], lyr = month)
    return(rst)
  }) %>% 
    reduce(., c)
  
  ## 
  prec <- lyer[[grep('ppt', names(lyer), value = F)]]
  prec <- mean(prec)
  names(prec) <- glue('prec_{month}')
  
  tmin <- lyer[[grep('tmin', names(lyer), value = F)]]
  tmin <- mean(tmin)  
  names(tmin) <- glue('tmin_{month}')    
  
  tmax <- lyer[[grep('tmax', names(lyer), value = F)]]
  tmax <- mean(tmax)  
  names(tmax) <- glue('tmax_{month}')    
  
  ## 
  prec <- terra::crop(prec, wrld)
  prec <- terra::mask(prec, wrld)
  
  tmin <- terra::crop(tmin, wrld)
  tmin <- terra::mask(tmin, wrld)    
  
  tmax <- terra::crop(tmax, wrld)
  tmax <- terra::mask(tmax, wrld)
  
  ##
  stck <- c(prec, tmin, tmax)
  return(stck)
  
}

# To apply the function ---------------------------------------------------
rstr <- map(1:12, calc.avrg)
rstr <- reduce(rstr, c)

## Filtering the variables in different object
prec <- rstr[[grep('prec', names(rstr), value = F)]]
tmin <- rstr[[grep('tmin', names(rstr), value = F)]]
tmax <- rstr[[grep('tmax', names(rstr), value = F)]]

## To write the raster
dir.create('./tif/tc_2c')
terra::writeRaster(x = prec, filename = './tif/tc_2c/prec.tif', overwrite = TRUE)
terra::writeRaster(x = tmin, filename = './tif/tc_2c/tmin.tif', overwrite = TRUE)
terra::writeRaster(x = tmax, filename = './tif/tc_2c/tmax.tif', overwrite = TRUE)

# To make the bioclimatic variables ---------------------------------------

## Raster to matrix
prec.mtrx <- terra::as.matrix(prec)
tmin.mtrx <- terra::as.matrix(tmin)
tmax.mtrx <- terra::as.matrix(tmax)

## Coordinates
crds <- terra::crds(prec, na.rm = F)

## To calculate the bioclimatic variables
bioc <- dismo::biovars(prec = prec.mtrx, tmin = tmin.mtrx, tmax = tmax.mtrx)
bioc <- cbind(crds, bioc)
bioc <- terra::rast(bioc, type = 'xyz')

## To write the raster
dir.create('./tif/tc_2c')
terra::writeRaster(x = bioc, filename = './tif/tc_2c/bioc.tif', overwrite = TRUE)

# To normalice the bios ---------------------------------------------------
bioc.mtrx <- terra::as.data.frame(bioc, xy = T)
crds <- bioc.mtrx[,1:2]
bioc.mtrx <- bioc.mtrx[,3:ncol(bioc.mtrx)]
nrml.bioc <- outliers::scores(x = bioc.mtrx, type = 'z')

# Table to raster ---------------------------------------------------------
nrml.rstr <- terra::rast(cbind(crds, nrml.bioc), type = 'xyz', crs = 'EPSG:4326')
terra::writeRaster(x = nrml.rstr, filename = './tif/tc_2c/ext_1/bioc_nrml.tif', overwrite = T)

# To do the PCA for the future --------------------------------------------

## PCA model
pca_model <- readRDS(file = './rds/pca_model-run1.rrds')

## Project
pca_projection <- predict(pca_model, newdata = nrml.bioc[,1:ncol(nrml.bioc)])
crds.pca <- as_tibble(cbind(crds[,1:2], pca_projection$coord))
rstr.pca <- terra::rast(crds.pca, type = 'xyz', crs = 'EPSG:4326')

full.crd <- as.matrix(crds(rstr.pca, na.rm = F))

## Raster to matrix 
pcas.mtrx <- as.matrix(rstr.pca, na.rm = T)

## To calculate the euclidean distance

## Points
pnts <- read_csv('./tbl/points_crds.csv', show_col_types = FALSE)
gids <- pull(pnts, id)
pnts <- as_tibble(cbind(pnts, terra::extract(rstr.pca, pnts[,2:3])))
pnts <- dplyr::select(pnts, -ID)

### Function to use
euc_clm <- function(gid){
  
  ## Filtering 
  cat('To process: ', gid, '\n')
  pnt <- filter(pnts, id == gid)
  pnt <- dplyr::select(pnt, id, Longitude, Latitude, Dim.1, Dim.2)
  bse <- dplyr::select(pnt, starts_with('Dim'))
  pca <- rstr.pca
  pca <- pca[[1:2]]
  
  ## Apply the euclidean distance 
  dst <- dista(xnew = as.matrix(bse), as.matrix(pcas.mtrx[,1:2]), type = 'euclidean')
  dst <- as.vector(dst)
  
  ## Table to raster
  fnl <- cbind(full.crd, euc = dst)
  fnl <- terra::rast(fnl, type = 'xyz')
  names(fnl) <- gid
  
  ## Table to raster 
  terra::writeRaster(x = fnl, filename = glue('./tif/euclidean/t2c/run-1_euc-t2c_{gid}.tif'), overwrite = TRUE)
  return(fnl)
  
}

## To apply the function 
eucl <- map(gids, euc_clm)
eucl <- reduce(eucl, c)
terra::writeRaster(x = eucl, filename = './tif/euclidean/t2c/run-1_ecu-t2c_allPoints_rescale.tif')


# To normalice the raster results -----------------------------------------
eucl <- terra::rast('./tif/euclidean/t2c/run-1_ecu-t2c_allPoints_rescale.tif')
rstr <- eucl

scle <- map(.x = 1:nlyr(rstr), .f = function(i){
  
  ## To select the raster
  cat('To start', i, '\n')
  rst <- rstr[[i]]
  nme <- names(rst)
  
  ## To normalice the raster
  scl <- rescale0to1(rst)
  scl <- raster.invert(scl)
  
  ## Finish 
  cat('Done!\n')
  return(scl)
  
})

scle <- reduce(scle, c)
terra::writeRaster(scle, './tif/euclidean/t2c/run-1_ecu-t2c_allPoints_rescale_2.tif')
