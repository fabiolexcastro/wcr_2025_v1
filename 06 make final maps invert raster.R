

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

# To invert the raster ----------------------------------------------------
bsln <- spatialEco::raster.invert(bsln)
ftre <- spatialEco::raster.invert(ftre)

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

# Climatic monthly analysis -----------------------------------------------

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

# To get the threshold ----------------------------------------------------

## Function to use
thrs.get <- function(gid){
  
  ## To make a stack from the two rasters
  cat('To process: ', gid, '\n')
  bsl <- bsln[[gid]]
  ftr <- ftre[[gid]]  
  stk <- c(bsl, ftr)
  names(stk) <- c('Baseline', 'Future 2°C')

  ## Raster to table 
  tbl <- terra::as.data.frame(stk, xy = TRUE)
  tbl <- gather(tbl, var, value, -c(x, y))
  tbl <- as_tibble(tbl)
  tbl <- filter(tbl, var == 'Baseline')
  
  ### Threshold
  thr <- quantile(x = pull(tbl, value), probs = c(0.1, 0.05, 0.01), na.rm = T)
  thr <- rownames_to_column(as.data.frame(thr))
  thr <- mutate(thr, percentil = parse_number(rowname), threshold = thr)
  thr <- mutate(thr, id = gid)
  thr <- dplyr::select(thr, id, percentil, threshold)
  
  thr.10 <- thr %>% filter(percentil == 10) %>% pull(3)
  thr.05 <- thr %>% filter(percentil ==  5) %>% pull(3)
  thr.01 <- thr %>% filter(percentil ==  1) %>% pull(3)
  
  ## Tolerance
  tolerance <- 0.000001
  
  ## Get the point 
  pnt <- pnts %>% filter(id == gid)
  
  ## Get the coordinate
  pnt.10 <- bsl %>% setNames('Baseline') %>% as.data.frame(xy = TRUE) %>% filter(abs(Baseline - thr.10) <= tolerance)
  pnt.05 <- bsl %>% setNames('Baseline') %>% as.data.frame(xy = TRUE) %>% filter(abs(Baseline - thr.05) <= tolerance)
  pnt.01 <- bsl %>% setNames('Baseline') %>% as.data.frame(xy = TRUE) %>% filter(abs(Baseline - thr.01) <= tolerance)
  
  ## Sampling with n > 1
  pnt.10 <- pnt.10 %>% sample_n(size = 1, replace = FALSE)
  pnt.05 <- pnt.05 %>% sample_n(size = 1, replace = FALSE)
  pnt.01 <- pnt.01 %>% sample_n(size = 1, replace = FALSE)
  
  ## Tidy the tables (thresholds + points) / Coordinates
  pnt.th <- rbind(pnt.10, pnt.05, pnt.01) %>% mutate(percentil = c(10, 5, 1)) %>% mutate(id = gid, .before = x)
  pnt <- pnt %>% setNames(c('id', 'x', 'y'))
  pnt <- pnt %>% mutate(Baseline = pull(terra::extract(bsl, pnt[,c(2,3)]), 2))
  pnt <- pnt %>% mutate(percentil = 100)
  pnt <- rbind(pnt.th, pnt)
  
  ## Extract the values for the points
  pnt.vls <- list(
    cbind(pnt, terra::extract(prec.bsln, pnt[,c('x', 'y')])),
    cbind(pnt, terra::extract(tmin.bsln, pnt[,c('x', 'y')])),
    cbind(pnt, terra::extract(tmax.bsln, pnt[,c('x', 'y')]))
  ) %>% 
    reduce(., inner_join) %>% 
    dplyr::select(-ID) %>% 
    gather(var, value, -c(id, x, y, Baseline, percentil)) %>% 
    as_tibble() %>% 
    separate(data = ., col = 'var', into = c('Variable', 'Month'), sep = '_') %>% 
    inner_join(., tibble(Month = as.character(1:12), month_abb = month.abb), by = 'Month') %>% 
    mutate(month_abb = factor(month_abb, levels = month.abb)) %>% 
    dplyr::select(-Month) %>% 
    spread(Variable, value) %>% 
    mutate(type = ifelse(percentil == 100, 'Trial', percentil), 
           type = factor(type, levels = c('Trial', '10', '5', '1')))
  
  ## To make the graph 
  rlc <- 10
  
  ggp <- ggplot(data = pnt.vls, aes(x = month_abb)) +
    geom_col(aes(y = prec), fill = 'grey50') +
    geom_line(aes(y = tmin * rlc, group = 1), col = 'grey30') +
    geom_line(aes(y = tmax * rlc, group = 1), col = 'grey30') +
    facet_wrap(.~type) +
    scale_y_continuous(sec.axis = sec_axis(~./rlc, name = 'Temperature ºC')) +
    ggtitle(label = gid) +
    labs(x = 'Month', y = 'Prec. (mm)') +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = 'bold', size = 16),
      axis.text.y = element_text(angle = 90, hjust = 0.5),
      strip.text = element_text(size = 14, face = 'bold')
    )
  
  ## To save the graph 
  ggsave(plot = ggp, filename = glue('./png/graphs/climatogram_percentiles_invert/climatogram_{gid}.jpg'), units = 'in', width = 9, height = 7, dpi = 300, create.dir = TRUE)
  
  ## Finish 
  cat('Done!\n')
  return(thr)
  
}

## To apply the function
thr.tbl <- map(gids, function(i){
  try(
    expr = {
      thrs.get(gid = i)
    }
  )
})
thr.tbl <- as_tibble(bind_rows(thr.tbl[map_lgl(thr.tbl, is.data.frame)]))
write.csv(thr.tbl, './tbl/values/thresholds_v1.csv', row.names = FALSE)

# To make the map for different percentiles -------------------------------

make.map.rcl <- function(gid){
  
  # gid <- gid[1]

  ## Filtering 
  cat('To process: ', gid, '\n')
  bsl <- bsln[[gid]]
  ftr <- ftre[[gid]]
  stk <- c(bsl, ftr)
  names(stk) <- c('Baseline', 'Future')
  
  ## Point
  pnt <- filter(pnts, id == gid)
  
  ## Threshold
  thr <- thr.tbl %>% filter(id == gid)
  
  ## To reclassify 
  tbl <- map_dfr(.x = 1:nrow(thr), .f = function(i){
    
    ## Filtering threshold
    th <- thr[i,]
    pr <- pull(th, percentil)
    th <- pull(th, threshold)
    
    ## To classify
    rs <- stk
    rc <- terra::ifel(rs < th, 0, 1)
    
    ## Raster to table 
    tb <- terra::as.data.frame(rc, xy = T)
    tb <- as_tibble(tb)
    tb <- gather(tb, var, value, -c(x, y))
    tb <- mutate(tb, var = factor(var, levels = c('Baseline', 'Future')))
    tb <- mutate(tb, percentil = pr)
    
    ## Finish
    cat('Done!\n')
    return(tb)
      
  }) %>% 
    mutate(percentil = factor(percentil, levels = c('90', '95', '99')),
           class = ifelse(value == 0, 'Unsuitable', 'Suitable'),
           class = factor(class, levels = c('Unsuitable', 'Suitable')))

  ## To make the map 
  
  ## Open Street Map
  lat1 <- -30 ; lat2 <- 24
  lon1 <- -110 ; lon2 <- 155
  map <- openmap(c(lat2,lon1), c(lat1,lon2), zoom = 2, type = c("osm", "stamen-toner", "stamen-terrain","stamen-watercolor", "esri","esri-topo", 'esri-physical', 'esri-shaded')[1], mergeTiles = TRUE)
  map <- openproj(map, projection = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
  
  grc <- autoplot(map) + 
    geom_tile(data = tbl, aes(x = x, y = y, fill = class)) +
    scale_fill_manual(values = c('grey70', 'darkgreen'), na.translate = FALSE) +
    facet_wrap(~var + percentil, ncol = 2, nrow = 3) +
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
  ggsave(plot = grc, filename = glue('./png/maps/euclidean_percentiles/{gid}_prds_percentiles.jpg'), units = 'in', width = 15, height = 8, dpi = 300, create.dir = T)
  cat('Done!\n')
    
}

map(.x = 2:length(gids), .f = function(z){
  try(
    expr = {
      make.map.rcl(gid = gids[z])
    }
  )
})
  












