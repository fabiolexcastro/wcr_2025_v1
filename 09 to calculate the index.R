

## Fabio Alexander Castro - Llanos 
## Alliance Bioversity - CIAT 
## July 30 th - 2024

## Get the index worldwide

# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(terra, fs, sf, openxlsx, readxl, geodata, ggspatial, ggrepel, rnaturalearthdata, rnaturalearth, climateStability, cptcity, spatialEco, parallelDist, tidyr, tidyverse, gridExtra, glue, gtools, readxl, dtw, dtwclust)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# Functions to use --------------------------------------------------------
calc.ndd.cdd <- function(yr){
  
  ## To list the files
  cat('>>> Process: ', yr, '\n')
  fls <- dir_ls(root.prec) %>% 
    as.character() %>% 
    grep(paste0('.', yr, '.'), ., value = T)
  
  ## To apply by each month
  rst <- map(.x = 1:12, .f = function(m){
    
    ## Filtering the month and reclassified by negative values
    cat('To processing: ', month.abb[m], '\n')
    mnt <- ifelse(m < 10, paste0('0', m), m)
    trr <- rast(grep(paste0('.', yr, '.', mnt, '.'), fls, value = T))
    trr[trr < 0] <- 0
    trr <- terra::crop(trr, ext(mskr))
    trr <- terra::resample(trr, mskr, method = 'bilinear')
    trr <- terra::mask(trr, wrld)
    
    ## Number of dry days
    ndd <- terra::app(x = trr, fun = function(x){ ndd <- sum(x < 1, na.rm = T); return(ndd)})
    names(ndd) <- glue('ndd_{yr}-{mnt}')
    ndd <- terra::mask(ndd, wrld)
    
    ## Number of consecutive dry days
    calc_cdd <- function(PREC, p_thresh = 1){
      runs <- rle(PREC < p_thresh)
      cons_days <- max(runs$lengths[runs$values==1], na.rm=TRUE)
      return(cons_days)
    }
    cdd <- terra::app(x = trr, fun = function(x){y = calc_cdd(PREC = x, p_thresh = 1); return(y)})
    cdd <- terra::mask(x = cdd, mask = wrld)
    names(cdd) <- glue('cdd_{yr}-{mnt}')
    
    ## To write the rasters 
    terra::writeRaster(x = ndd, filename = glue('./tif/index_world_v1/baseline/ndd_{yr}-{mnt}.tif'), overwrite = TRUE)
    terra::writeRaster(x = cdd, filename = glue('./tif/index_world_v1/baseline/cdd_{yr}-{mnt}.tif'), overwrite = TRUE)
    
    ## Return the two final rasters
    rm(trr, ndd, cdd); gc(reset = TRUE)
    cat('Done :)\n')
    
  })
  
}
calc.ntx <- function(yr){
  
  ## To list the files
  cat('>>> Process: ', yr, '\n')
  fls <- dir_ls(root.tmax) %>% as.character() %>% grep('temperature', ., value = T) %>% grep('maximum', ., value = T)
  fls <- grep(paste0('_', yr), fls, value = T)
  print(length(fls))
  
  ## To apply the function for each month
  map(.x = 1:12, .f = function(m){
    
    ## Filtering the month
    cat('To process: ', month.abb[m], '\n')
    mnt <- ifelse(m < 10, paste0('0', m), as.character(m))
    fle <- grep(paste0('_', yr, mnt), fls, value = T)
    
    ## To read the raster and convert the units
    rst <- rast(fle)
    rst <- terra::crop(rst, ext(mskr))
    rst <- terra::resample(rst, mskr, method = 'bilinear')
    rst <- rst - 273.15
    
    ## Threshold
    thr <- 30
    
    ## To calculate the index and write the final 
    ntx <- terra::app(x = rst, fun = function(x){ntxv = sum(x >= thr, na.rm = T); return(ntxv)})
    ntx <- terra::mask(ntx, wrld)
    names(ntx) <- glue('ntx30_{yr}-{mnt}')
    
    ## To write the final raster and finish
    terra::writeRaster(x = ntx, filename = glue('./tif/index_world_v1/baseline/ntx30_{yr}-{mnt}.tif'), overwrite = TRUE)
    rm(rst, ntx); gc(reset = TRUE)
    cat('Done!\n')
    
  })
  
  ## Finish
  cat('Done!\n')
  
}

# Load data ---------------------------------------------------------------

## Tabular data 
pnts <- as_tibble(read.xlsx('../V1/r/tbl/points/IMLVT Site info.xlsx'))
crds <- dplyr::select(pnts, X1, id, Longitude, Latitude)
data <- read_csv('./data/tbl/crds_index.csv', show_col_types = F)

## Directories (raster data)
root.tmax <- '//catalogue/WFP_ClimateRiskPr1/1.Data/AgERA5/2m_temperature-24_hour_maximum/'
root.tmin <- '//catalogue/WFP_ClimateRiskPr1/1.Data/AgERA5/2m_temperature-24_hour_minimum/'
root.dewp <- ' '
root.prec <- '//CATALOGUE/WFP_ClimateRiskPr1/1.Data/Chirps'

## Vector data 
wrld <- vect(ne_countries(returnclass = 'sf', scale = 50))
idn0 <- gadm(country = 'IDN', level = 0, path = './tmpr')

# To create a mask --------------------------------------------------------

## Precipitation
mask.prec <- dir_ls(root.prec) %>% .[1] %>% rast() %>% ifel(. < 0, 0, .)
mask.prec <- mask.prec * 0 + 1; names(mask.prec) <- 'mask'

## Temperature 
mask.tasm <- dir_ls(root.tmax) %>% .[1] %>% rast()

### To check the resolution 
res(mask.prec)
res(mask.tasm)

## Study area --------------------------------------------------------------
mask.tasm <- terra::crop(mask.tasm, ext(mask.prec))
mask.tasm <- terra::mask(mask.tasm, wrld)
mask.tasm <- mask.tasm * 0 + 1
mask.prec <- terra::crop(mask.prec, ext(mask.prec))
mask.prec <- terra::mask(mask.prec, wrld)

## To make a resample  -----------------------------------------------------
mask.prec <- terra::resample(mask.prec, mask.tasm, method = 'bilinear')
mask.prec <- terra::mask(mask.prec, wrld)

## Sum the two rasters -----------------------------------------------------
mskr <- mask.tasm + mask.prec
dout <- glue('./tif/base')
dir_create(dout)
terra::writeRaster(x = mskr, filename = paste0(dout, '/', 'mask_v1.tif'), overwrite = TRUE)

mskr <- terra::rast('./tif/base/mask_v1.tif')

# To calculate NDD / CDD --------------------------------------------------
map(.x = 2016:2023, .f = calc.ndd.cdd)

# To calculate NTX --------------------------------------------------------
map(2016:2023, calc.ntx)
