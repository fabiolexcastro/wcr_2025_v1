

# Load libraries ----------------------------------------------------------
require(pacman)
p_load(terra, sf, fs, tidyverse, spatialEco, RColorBrewer, OpenStreetMap, climateStability, parallelDist, glue, outliers, dismo, scales, glue, rnaturalearthdata, rnaturalearth, openxlsx)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# Load data ---------------------------------------------------------------

## Raster data 
bsln <- terra::rast('./tif/euclidean/run-1_ecu-bsl_allPoints_rescale.tif') 
ftre <- terra::rast('./tif/euclidean/t2c/run-1_ecu-t2c_allPoints_rescale.tif')
gids <- names(bsln)

## Tabular data 
pnts <- read_csv('./tbl/points_crds.csv', show_col_types = FALSE)

## Vector data 
wrld <- ne_countries(returnclass = 'sf', scale = 50)  

# To make the map  --------------------------------------------------------

## Function
make.map <- function(gid){
  
  ## To make a stack from the two rasters
  cat('To process: ', gid, '\n')
  bsl <- bsln[[gid]]
  ftr <- ftre[[gid]]  
  stk <- c(bsl, ftr)
  names(stk) <- c('Baseline', 'Future 2°C')
  
  ## Raster to table
  tbl <- terra::as.data.frame(stk, xy = T) %>% 
    as_tibble() %>% 
    gather(var, value, -c(x, y)) %>% 
    mutate(var = factor(var, levels = c('Baseline', 'Future 2°C')))
  
  ## Point tibble
  pnt <- filter(pnts, id == gid)
  
  ## Open Street Map
  lat1 <- -30 ; lat2 <- 24
  lon1 <- -110 ; lon2 <- 155
  map <- openmap(c(lat2,lon1), c(lat1,lon2), zoom = 2, type = c("osm", "stamen-toner", "stamen-terrain","stamen-watercolor", "esri","esri-topo", 'esri-physical', 'esri-shaded')[1], mergeTiles = TRUE)
  map <- openproj(map, projection = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
  
  ## To build the map 
  gmp <- autoplot(map) + 
    geom_tile(data = tbl, aes(x = x, y = y, fill = value)) +
    facet_wrap(~var, ncol = 1, nrow = 2) +
    scale_fill_gradientn(colors = brewer.pal(n = 9, name = 'RdYlGn')) +
    geom_sf(data = wrld, fill = NA, col = 'grey30', inherit.aes = FALSE) + 
    geom_point(data = pnt, aes(x = Longitude, y = Latitude), col = 'brown', size = 1) +
    coord_sf(xlim = ext(stk)[1:2], ylim = ext(stk)[3:4]) + 
    labs(x = 'Lon', y = 'Lat', fill = 'Suitability') +
    ggtitle(label = glue('{gid}')) +
    theme_minimal() +
    theme(
      strip.text = element_text(face = 'bold', hjust = 0.5), 
      plot.title = element_text(face = 'bold', hjust = 0.5, size = 22),
      axis.text.y = element_text(angle = 90, hjust = 0.5, size = 6),
      axis.text.x = element_text(size = 6),
      axis.title = element_text(size = 8), 
      legend.position = 'bottom'
    ) +
    guides(fill = guide_legend( 
      direction = 'horizontal',
      keyheight = unit(1.15, units = "mm"),
      keywidth = unit(15, units = "mm"),
      title.position = 'top',
      title.hjust = 0.5,
      label.hjust = .5,
      nrow = 1,
      byrow = T,
      reverse = F,
      label.position = "bottom"
    )) 
  
  ## To save the map
  ggsave(plot = gmp, filename = glue('./png/maps/euclidean/{gid}_prds_raw.jpg'), units = 'in', width = 7, height = 8, dpi = 300)
  
  ## To reclassify the raster
  
  ### Threshold
  thr <- as.numeric(quantile(x = pull(tbl, value), probs = 0.9, na.rm = T))
  tbl <- tbl %>% mutate(value_rcl = ifelse(value < thr, NA, value))
  tbl <- tbl %>% mutate(class = ifelse(is.na(value_rcl), 'Not suitable', 'Suitable'))
  tbl <- tbl %>% mutate(class = factor(class, levels = c('Not suitable', 'Suitable')))
  
  ## To make the classify map
  grc <- autoplot(map) + 
    geom_tile(data = tbl, aes(x = x, y = y, fill = class)) +
    scale_fill_manual(values = c('grey70', 'darkgreen')) +
    facet_wrap(~var, ncol = 1, nrow = 2) +
    geom_sf(data = wrld, fill = NA, col = 'grey50', inherit.aes = FALSE) + 
    geom_point(data = pnt, aes(x = Longitude, y = Latitude), col = 'brown', size = 1) +
    coord_sf(xlim = ext(stk)[1:2], ylim = ext(stk)[3:4]) + 
    labs(x = 'Lon', y = 'Lat', fill = 'Suitability') +
    ggtitle(label = glue('{gid}')) +
    theme_minimal() +
    theme(
      strip.text = element_text(face = 'bold', hjust = 0.5), 
      plot.title = element_text(face = 'bold', hjust = 0.5, size = 22),
      axis.text.y = element_text(angle = 90, hjust = 0.5, size = 6),
      axis.text.x = element_text(size = 6),
      axis.title = element_text(size = 8), 
      legend.position = 'bottom'
    ) +
    guides(fill = guide_legend( 
      direction = 'horizontal',
      keyheight = unit(1.15, units = "mm"),
      keywidth = unit(15, units = "mm"),
      title.position = 'top',
      title.hjust = 0.5,
      label.hjust = .5,
      nrow = 1,
      byrow = T,
      reverse = F,
      label.position = "bottom"
    )) 
  
  ## To save the map
  ggsave(plot = grc, filename = glue('./png/maps/euclidean/{gid}_prds_rcl.jpg'), units = 'in', width = 10, height = 6.5, dpi = 300)
  cat('Done!\n')
  
  
}

## To apply the function for making the maps
map(gids[30:length(gids)], make.map)
gids
