
# Load libraries ----------------------------------------------------------
require(pacman)
p_load(terra, sf, fs, tidyverse, tidyterra, rmapshaper, parallelDist, Rfast, glue, outliers, spatialEco, climateStability, dismo, scales, glue, rnaturalearthdata, rnaturalearth, openxlsx)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# Load data ---------------------------------------------------------------

## Points 
pnts <- as_tibble(read.xlsx('./tbl/IMLVT Site info.xlsx'))
pnts <- mutate(pnts, id = gsub('Bulegeni\t\t_UGA', 'Bulegeni_UGA', id))
gids <- unique(pnts$id)

## List the files 
fles <- as.character(dir_ls('./tif/output/euc'))

## Vector data 
wrld <- ne_countries(returnclass = 'sf', scale = 50)
cntn <- ms_dissolve(wrld, field = 'continent')

# To draw the maps --------------------------------------------------------
# gid <- gids[1]
# prc <- 0.99

# To make the binary analysis ---------------------------------------------

## Function
make.bins <- function(gid, prc){
  
  ##
  cat('To process: ', gid, '\n')
  fls <- grep(gid, fles, value = T)
  
  ## Current
  bsl <- rast(grep('bsl', fls, value = T))
  
  ## Future 
  ftr <- rast(fls[-grep('bsl', fls, value = F)])
  
  ## Point
  pnt <- filter(pnts, id == gid)
  
  ## To normalice and intert the raster
  
  ### Baseline
  bsl.nrm <- rescale0to1(bsl)
  bsl.nrm <- raster.invert(bsl.nrm)
  
  ### Future 
  ftr.nrm <- map(.x = 1:nlyr(ftr), .f = function(i){print(i); rescale0to1(ftr[[i]])})
  ftr.nrm <- reduce(ftr.nrm, c)
  
  ## To extract the percentile 
  thr <- as.numeric(terra::global(x = bsl.nrm, fun = stats::quantile, probs = prc, na.rm = T))
  
  ## To binarize 
  mtx <- matrix(c(0, thr, 0, thr, 1, 1), ncol = 3, byrow = TRUE)
  bsl.bin <- terra::classify(bsl.nrm, mtx, include.lowest = T)
  ftr.bin <- terra::classify(ftr.nrm, mtx, include.lowest = T)
  
  ## To make a stack 
  stk <- c(bsl.bin, ftr.bin)
  
  ## To write the raster normaliced and binarized (current and future)
  terra::writeRaster(x = stk, filename = glue('./tif/output/euc-nrm-bin/euc-run1_{gid}_{prc}.tif'), overwrite = TRUE)
  
  ## To compile the rasters binnary
  ftr.avg <- mean(ftr.nrm)
  ftr.avg.bin <- terra::classify(ftr.avg, mtx, include.lowest = T)
  ftr.gcm.bin <- terra::classify(ftr.nrm, mtx, include.lowest = T)
  
  ## To write the rasters
  terra::writeRaster(x = ftr.gcm.bin, filename = glue('./tif/output/euc-nrm-bin-gcm/euc-run1_{gid}_{prc}_gcms.tif'), overwrite = TRUE)
  rm(ftr.avg, ftr.avg.bin, ftr.gcm.bin, stk, bsl.bin, ftr.bin)
  gc(reset = T)
  cat('Done!\n')
  
}

## To apply the function
map(.x = gids, .f = function(gd){make.bins(gid = gd, prc = 0.90)})

# To read the results -----------------------------------------------------
rstr.bsln <- as.character(dir_ls('./tif/output/euc-nrm-bin')) %>% grep('0.9.tif', ., value = T) %>% map(rast, lyr = 1) %>% reduce(., c)
fles.ftre <- as.character(dir_ls('./tif/output/euc-nrm-bin-gcm'))
fles.ftre <- fles.ftre %>% grep('0.9_', ., value = T) 

# To draw the maps -------------

## Function
make.map <- function(gid){
  
  # gid <- gids[1]
  
  ##
  cat('Point: ', gid, '\n')
  rst.ftr <- rast(grep(gid, fles.ftre, value = T))
  rst.bsl <- rstr.bsln[[grep(gid, names(rstr.bsln))]]
  
  ## Future modal 
  rst.ftr.mdl <- modal(rst.ftr)
  
  ## 
  pnt <- filter(pnts, id == gid)
  
  ## 
  rst.bsl <- as.factor(rst.bsl)
  levels(rst.bsl) <- data.frame(id = c(0, 1), class = c('Not similarity', 'Similarity'))
  
  ##
  rst.ftr.mdl <- as.factor(rst.ftr.mdl)
  levels(rst.ftr.mdl) <- data.frame(id = c(0, 1), class = 'Not similarity',  'Similarity')
  
  ## Extract GCM names
  extract_gcm_names <- function(r){names(r) %>% strsplit("_") %>% sapply(function(x) x[2]) %>% unique()}
  gcm_names <- extract_gcm_names(rst.ftr)
  
  ## Baseline map con geom_spatraster
  g.bsl <- ggplot() + 
    geom_spatraster(data = rst.bsl, aes(fill = class)) +
    scale_fill_manual(values = c('Not Similarity' = 'grey80', 
                                 'Similarity' = 'forestgreen'),
                      na.value = NA) +
    geom_sf(data = wrld, fill = NA, col = 'grey30') +
    labs(x = '', y = '', fill = '') +
    ggtitle(label = gid) +
    coord_sf() +
    theme_minimal() +
    theme(
      plot.title = element_text(face = 'bold', size = 14),
      strip.text = element_text(face = 'bold'),
      legend.position = 'bottom'
    )
  
  ## Make one stack 
  stk <- c(rst.bsl, rst.ftr.mdl)
  names(stk) <- c('Baseline', 'Future')
  stk <- as.factor(stk)
  
  ## Raster tot able
  stk_df <- as.data.frame(stk, xy = TRUE, na.rm = TRUE) %>% tidyr::pivot_longer(cols = c(Baseline, Future), names_to = "period", values_to = "value")
  stk_df$value <- factor(stk_df$value, levels = c(0, 1), labels = c('Not similarity', 'Similarity'))
  
  ## To draw the map
  g.bsl.ftr <- ggplot() + 
    geom_raster(data = stk_df, aes(x = x, y = y, fill = value)) +
    scale_fill_manual(values = c('Not similarity' = 'grey80', 'Similarity' = 'forestgreen'), na.value = NA) +
    geom_sf(data = wrld, fill = NA, col = 'grey30') +
    geom_point(dat = pnt, aes(x = Longitude, y = Latitude), col = 'brown') +
    labs(x = '', y = '', fill = '') +
    ggtitle(label = paste(gid)) +
    facet_wrap(~period, ncol = 2, nrow = 1) +
    coord_sf() +
    theme_minimal() +
    theme(
      plot.title = element_text(face = 'bold', size = 14, hjust = 0.5),
      strip.text = element_text(face = 'bold'),
      legend.position = 'bottom'
    )
  
  # To save the maps 
  ggsave(plot = g.bsl.ftr, filename = glue('./png/maps/prc90/bsl-ftr_{gid}.jpg'), units = 'in', width = 10, height = 5, dpi = 300)
  rm(rst.bsl, rst.ftr, tbl.bsl, tbl.ftr, g.bsl, g.ftr)
  gc(reset = T)
  
}

## To apply the function
map(gids, make.map)
