

# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(terra, sf, fs, tidyverse, geodata, glue, rnaturalearthdata, rnaturalearth)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn =1)

# Load data ---------------------------------------------------------------

## Extent
extn <- c(-156, 156, -30, 30)

## Vector data
wrld <- ne_countries(returnclass = 'sf', scale = 50)
wrld <- vect(wrld)
wrld <- terra::crop(wrld, extn)

# Raster data -------------------------------------------------------------

## Download baseline climate 
bioc.bsln <- worldclim_global(var = 'bioc', res = 2.5, path = './tmpr')
bioc.bsln <- terra::crop(bioc.bsln, extn)
names(bioc.bsln) <- glue('bioc_{1:19}')
dir_create('./tif/wc_25')
terra::writeRaster(bioc.bsln, filename = glue('./tif/wc_25/bioc.tif'), overwrite = TRUE)

## Download future data
gcms <- c("ACCESS-CM2", "ACCESS-ESM1-5", "AWI-CM-1-1-MR", "BCC-CSM2-MR", "CanESM5", "CanESM5-CanOE", "CMCC-ESM2", "CNRM-CM6-1", "CNRM-CM6-1-HR", "CNRM-ESM2-1", "EC-Earth3-Veg", "EC-Earth3-Veg-LR", "FIO-ESM-2-0", "GFDL-ESM4", "GISS-E2-1-G", "GISS-E2-1-H", "HadGEM3-GC31-LL", "INM-CM4-8", "INM-CM5-0", "IPSL-CM6A-LR", "MIROC-ES2L", "MIROC6", "MPI-ESM1-2-HR", "MPI-ESM1-2-LR", "MRI-ESM2-0", "UKESM1-0-LL")
ssps <- c('370')
prds <- c('2021-2040', '2041-2060')

### Loop to download
map(.x = 1:length(gcms), .f = function(i){
  
  try(
    expr = {
      
      ## To download
      cat('To process: ', gcms[i], '\n')
      prec <- cmip6_world(model = gcms[i], ssp = ssps[1], time = '2041-2060', var = 'prec', path = './tmpr', res = 2.5)
      prec <- terra::crop(prec, ext(-156, 156, -30, 30))
      
      tmin <- cmip6_world(model = gcms[i], ssp = ssps[1], time = '2041-2060', var = 'tmin', path = './tmpr', res = 2.5)
      tmin <- terra::crop(tmin, ext(-156, 156, -30, 30))
      
      tmax <- cmip6_world(model = gcms[i], ssp = ssps[1], time = '2041-2060', var = 'tmax', path = './tmpr', res = 2.5)
      tmax <- terra::crop(tmax, ext(-156, 156, -30, 30))
      
      ## To write the raster 
      dout <- glue('./tif/c6_30s/ssp370-2041-2060/{gcms[i]}')
      dir_create(dout)
      terra::writeRaster(x = prec, filename = glue('{dout}/prec.tif'), overwrite = TRUE)
      terra::writeRaster(x = tmin, filename = glue('{dout}/tmin.tif'), overwrite = TRUE)
      terra::writeRaster(x = tmax, filename = glue('{dout}/tmax.tif'), overwrite = TRUE)
      rm(prec, tmin, tmax)
      gc(reset = T)
      cat('Done!\n')
      
    }
    
  )

})

# Check the downloaded files ----------------------------------------------

dirs <- dir_ls('./tif/c6_30s/ssp370-2041-2060') %>% as.character()
fles.dirs <- map(dirs, dir_ls)
fles.dirs

## Loop to create bioclimatic variables 
map(.x = 2:length(fles.dirs), .f = function(i){
  
  ## Files
  fls <- as.character(dir_ls(dirs[i]))
  ppt <- rast(grep('prec', fls, value = T))
  tmn <- rast(grep('tmin', fls, value = T))
  tmx <- rast(grep('tmax', fls, value = T))
  
  ## To matrix 
  ppt.mtx <- as.matrix(as.data.frame(ppt, xy = F, na.rm = T))
  tmn.mtx <- as.matrix(as.data.frame(tmn, xy = F, na.rm = T))
  tmx.mtx <- as.matrix(as.data.frame(tmx, xy = F, na.rm = T))
  
  ## Coordinates 
  crd <- as.data.frame(ppt[[1]], xy = T, na.rm = T)
  
  ## To create the bioclimatic variables
  bio.mtx <- dismo::biovars(prec = ppt.mtx, tmin = tmn.mtx, tmax = tmx.mtx)
  bio.mtx <- as.data.frame(bio.mtx)
  
  ## Join bios with the coordinates
  bio.mtx <- cbind(crd[,1:2], bio.mtx)
  
  ## Matrix to raster 
  bio.rst <- terra::rast(bio.mtx, type = 'xyz', crs = 'EPSG:4326')
  
  ## To write the raster
  terra::writeRaster(x = bio.rst, filename = glue('{dirs[i]}/bioc.tif'), overwrite = TRUE)
  cat('Done!\n')
  rm(ppt, tmn, tmx, ppt.mtx, tmn.mtx, tmx.mtx, crd, bio.mtx, bio.rst)
  gc(reset = T)
  
  
})
