

# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(terra, fs, sf, tidyverse, gtools, glue, ggspatial, RColorBrewer, rnaturalearthdata, rnaturalearth)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# Load data ---------------------------------------------------------------

## Vector data
wrld <- ne_countries(returnclass = 'sf', scale = 50)

## Baseline data
fles.bsln <- dir_ls('./tif/tc_baseline/ext_3', regexp = '.tif$') %>% as.character()
prec.bsln <- grep('prec', fles.bsln, value = T) %>% rast()
tmin.bsln <- grep('tmin', fles.bsln, value = T) %>% rast()
tmax.bsln <- grep('tmax', fles.bsln, value = T) %>% rast()

## Future data
fles.ftre <- as.character(dir_ls('./tif/tc_2c/ext_3'))
prec.ftre <- grep('prec', fles.ftre, value = T) %>% rast()
tmin.ftre <- grep('tmin', fles.ftre, value = T) %>% rast()
tmax.ftre <- grep('tmax', fles.ftre, value = T) %>% rast()

## Points 
pnts <- read_csv('./tbl/points_crds.csv', show_col_types = FALSE)

# To extract the values ---------------------------------------------------

## Function to use
extr.vles <- function(stk){
  
  cat('To start the analysis\n')
  vls <- cbind(pnts, terra::extract(stk, pnts[,c('Longitude', 'Latitude')]))
  vls <- as_tibble(vls)
  return(vls)
  
}

## Baseline
vles.bsln <- list(prec.bsln, tmin.bsln, tmax.bsln) %>%
  map(., extr.vles) %>% 
  reduce(., inner_join) %>% 
  gather(vari, value, -c(id, Longitude, Latitude, ID)) %>% 
  separate(data = ., col = vari, into = c('Variable', 'Month'), sep = '_') %>% 
  spread(Variable, value) %>% 
  mutate(period = 'Baseline')

## Future
vles.ftre <- list(prec.ftre, tmin.ftre, tmax.ftre) %>%
  map(., extr.vles) %>% 
  reduce(., inner_join) %>% 
  gather(vari, value, -c(id, Longitude, Latitude, ID)) %>% 
  separate(data = ., col = vari, into = c('Variable', 'Month'), sep = '_') %>% 
  spread(Variable, value) %>% 
  mutate(period = 'Future')

vles.full <- rbind(vles.bsln, vles.ftre)  

## To write the table 
write.csv(vles.full, './tbl/values/values_monthly_bsl-ftr.csv', row.names = FALSE)


  