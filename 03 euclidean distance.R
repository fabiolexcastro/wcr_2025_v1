

# Load libraries ----------------------------------------------------------
require(pacman)
p_load(terra, sf, fs, tidyverse, parallelDist, Rfast, glue, outliers, spatialEco, climateStability, dismo, scales, glue, rnaturalearthdata, rnaturalearth, openxlsx)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# Load data ---------------------------------------------------------------

## Points 
pnts <- as_tibble(read.xlsx('./tbl/IMLVT Site info.xlsx'))
pnts <- mutate(pnts, id = gsub('Bulegeni\t\t_UGA', 'Bulegeni_UGA', id))
gids <- unique(pnts$id)

## Raster pca
fle.bsln <- './tif/output/pca/pca-bsl_run1.tif'
fls.ftre <- dir_ls('./tif/c6_25/ssp370-2041-2060') %>% map(., dir_ls) %>% unlist() %>% as.character() %>% grep('rpca', ., value = T) 
gcms <- basename(dirname(fls.ftre))

# Function to make the euclidean ------------------------------------------
euc.clm <- function(gid, fle, prd){
  
  # gid <- gids[1]
  # fle <- fle.bsln
  # prd <- 'bsl'
  
  ## To read as a raster
  cat('To processs: ', gid, ' ', prd, '\n')
  rst <- terra::rast(fle)
  pnt <- filter(pnts, id == gid)
  
  ## Extract the baseline value 
  bse <- cbind(pnt[,c('id', 'Longitude', 'Latitude')], terra::extract(rast(fle.bsln), pnt[,c('Longitude', 'Latitude')]))
  bse <- dplyr::select(bse, -ID)
  
  ## Raster to matrix 
  mtx <- terra::as.matrix(as.data.frame(rst, na.rm = T, xy = F))
  
  ## To apply the euclidean distance 
  dst <- dista(as.matrix(bse[,4:ncol(bse)]), as.matrix(mtx), type = 'euclidean')
  dst <- as.vector(dst)
  
  ## Vector to raster 
  fnl <- cbind(crds(rst), dst)
  fnl <- terra::rast(fnl, type = 'xyz', crs = 'EPSG:4326')
  names(fnl) <- glue('euc_{prd}_{gid}')
  
  ## To write the raster 
  terra::writeRaster(x = fnl, filename = glue('./tif/output/euc/euc-run1_{gid}_{prd}.tif'), overwrite = TRUE)
  cat('Done!\n')
  
  ## To clean the object for optimize the RAM
  rm(rst, pnt, bse, mtx, dst, fnl)
  gc(reset = T)
  
}

# To apply the function ---------------------------------------------------

## Baseline
euc.bsl <- map(.x = gids, .f = function(g){euc.clm(gid = g, fle = fle.bsln, prd = 'bsl')})
euc.bsl <- dir_ls('./tif/output/euc') %>% as.character() %>% grep('bsl', ., value = T) %>% rast()

## Future 
euc.ftr <- map(.x = gcms, .f = function(gc){
  
  ## GCMe
  gcme <- grep(paste0('/', gc, '/'), fls.ftre, value = T)
  
  ## Apply
  rstr <- map(.x = gids, .f = function(g){
    euc.clm(gid = g, fle = gcme, prd = basename(dirname(gcme)))
  })
  
  ## Finish
  cat('Finish: ', gcme, '\n')
  
})

fls.ftre[4] %>% rast()
