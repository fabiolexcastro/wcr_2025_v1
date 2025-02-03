
# Load libraries ----------------------------------------------------------
require(pacman)
p_load(terra, sf, fs, tidyverse, parallelDist, Rfast, glue, outliers, spatialEco, climateStability, dismo, scales, glue, rnaturalearthdata, rnaturalearth, openxlsx)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# Load data ---------------------------------------------------------------

## Points
pnts <- read_csv('./tbl/points_crds.csv', show_col_types = FALSE)
gids <- pull(pnts, id)

## Raster 
bioc <- terra::rast('./tif/tc_baseline/ext_1/bioc_normal.tif')
pcas <- terra::rast('./tif/pca/pca-bsl_run1.tif')

# PCA matrix --------------------------------------------------------------
pcas.mtrx <- terra::as.data.frame(pcas, xy = T, na.rm = T)

# Extract the values for the points ---------------------------------------
pnts <- as_tibble(cbind(pnts, terra::extract(pcas, pnts[,2:3])))
pnts <- dplyr::select(pnts, -ID)
  
# Function ----------------------------------------------------------------
euc_clm <- function(gid){
  
  ## Filtering 
  cat('To process: ', gid, '\n')
  pnt <- filter(pnts, id == gid)
  pnt <- dplyr::select(pnt, id, Longitude, Latitude, Dim.1, Dim.2)
  bse <- dplyr::select(pnt, starts_with('Dim'))
  pca <- pcas
  pca <- pca[[1:2]]
  
  ## Apply the euclidean distance 
  dst <- dista(xnew = as.matrix(bse), as.matrix(pcas.mtrx[,3:4]), type = 'euclidean')
  dst <- as.vector(dst)
  
  ## Table to raster
  fnl <- cbind(crds(pcas), dst)
  fnl <- terra::rast(fnl, type = 'xyz')
  
  ## Table to raster 
  terra::writeRaster(x = fnl, filename = glue('./tif/euclidean/run-1_euc-bsl_{gid}.tif'), overwrite = TRUE)
  return(fnl)
  
}

##
rstr <- map(pull(pnts, id), euc_clm)
rstr <- reduce(rstr, c)
names(rstr) <- pull(pnts, id)

## 
terra::writeRaster(x = rstr, filename = './tif/euclidean/run-1_euc-bsl_allPoints.tif', overwrite = TRUE)




