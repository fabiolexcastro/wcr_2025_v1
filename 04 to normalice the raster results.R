
# Load libraries ----------------------------------------------------------
require(pacman)
p_load(terra, sf, fs, tidyverse, spatialEco, climateStability, parallelDist, glue, outliers, dismo, scales, glue, rnaturalearthdata, rnaturalearth, openxlsx)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# Load data ---------------------------------------------------------------
rstr <- terra::rast('./tif/euclidean/run-1_euc-bsl_allPoints.tif')

# To normalice ------------------------------------------------------------
scle <- map(.x = 1:nlyr(rstr), .f = function(i){
  
  ## To select the raster
  cat('To start!\n')
  rst <- rstr[[i]]
  nme <- names(rst)
  
  ## To normalice the raster
  scl <- rescale0to1(rst)
  scl <- raster.invert(scl)

  ## Finish 
  cat('Done!\n')
  return(scl)
  
})

##
scle <- reduce(scle, c)
terra::writeRaster(x = scle, filename = './tif/euclidean/run-1_ecu-bsl_allPoints_rescale.tif', overwrite = TRUE)


