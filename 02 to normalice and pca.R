

# Load libraries ----------------------------------------------------------
require(pacman)
p_load(terra, sf, fs, tidyverse, geodata, parallelDist, factoextra, FactoMineR, corrr, ggcorrplot, outliers, dismo, scales, glue, rnaturalearthdata, rnaturalearth, openxlsx)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# Load data ---------------------------------------------------------------

## Points 
pnts <- as_tibble(read.xlsx('./tbl/IMLVT Site info.xlsx'))
vctr <- terra::vect(pnts, c('Longitude', 'Latitude'), crs = 'EPSG:4326')

## Raster
bsln <- terra::rast('./tif/wc_25/bioc.tif')
ftre.fles <- dir_ls('./tif/c6_25/ssp370-2041-2060') %>% 
  map(., dir_ls) %>% 
  unlist() %>% 
  as.character() %>% 
  grep('bioc.tif', ., value = T) 


# To normalice ------------------------------------------------------------

## Bios baseline [statistics]
bsln.stts <- map_dfr(.x = 1:nlyr(bsln), .f = function(i){
  r <- bsln[[i]] 
  a <- terra::global(x = r, fun = mean, na.rm = T)
  s <- terra::global(x = r, fun = sd, na.rm = T)
  f <- tibble(variable = names(r), mean = as.numeric(a), sdt = as.numeric(s))
  return(f)
})
bsln.stts

## To normalice baseline
bsln.nrml <- map(.x = 1:nlyr(bsln), .f = function(i){
  cat('To process: ', names(bsln[[i]]), '\n')  
  r <- bsln[[i]]
  n <- names(r)
  d <- filter(bsln.stts, variable == n)
  f <- (r - pull(d, mean)) / pull(d, sdt)
  return(f)
})
bsln.nrml <- reduce(bsln.nrml, c)

terra::writeRaster(x = bsln.nrml, glue('./tif/wc_25/bioc_nrml.tif'), overwrite = TRUE)
bsln.nrml <- terra::rast('./tif/wc_25/bioc_nrml.tif')

## To normalice future 
map(.x = 1:length(ftre.fles), .f = function(i){
  
  cat('To process: ', basename(dirname(ftre.fles[i])), '\n')
  rstr <- rast(ftre.fles[i])
  names(rstr) <- glue('bioc_{1:19}')
  
  ## Loop by each layers
  rstr.nrml <- map(.x = 1:19, .f = function(j){
    
    cat('Layer: ', j, '\t')
    rst <- rstr[[j]]
    nme <- names(rst)
    dfm <- filter(bsln.stts, variable == nme)
    fnl <- (r - pull(d, mean)) / pull(d, sdt)
    rm(rst, nme, dfm)
    gc(reset = T)
    return(fnl)
    
  }) %>% 
    reduce(., c)
  
  ## To write the raster
  dire <- dirname(ftre.fles[i])
  terra::writeRaster(x = rstr, filename = glue('{dire}/bioc_nrml.tif'), overwrite = TRUE)
  rm(rstr)
  gc(reset = T)
  
})

## Read the normalice future datasets
ftre.fles <- dir_ls('./tif/c6_25/ssp370-2041-2060') %>% 
  map(., dir_ls) %>% 
  unlist() %>% 
  as.character() %>% 
  grep('bioc_nrml', ., value = T) 

# To make the PCA ---------------------------------------------------------

## Seed 
set.seed(123)

## Sample size 
land_mask <- !is.na(bsln.nrml[[1]])
valid_cells <- which(values(land_mask) == TRUE) 
sample_cells <- sample(valid_cells, min(50000, length(valid_cells)))
sampled_data <- as.data.frame(terra::extract(bsln.nrml, sample_cells))

## Get the coordinates from sample cells 
coords <- terra::xyFromCell(bsln.nrml, sample_cells)
sampled_df <- data.frame(coords, sampled_data)

## PCA Model 
pca_model <- PCA(sampled_df[, -c(1,2)], scale = FALSE, graph = FALSE)
saveRDS(object = pca_model, file = './rds/pca_model-run1.rds')

# Visualizar la varianza explicada
cntr.dims <- pca_model$eig %>% 
  as.data.frame() %>% 
  as_tibble() %>%
  mutate(porc = round(`percentage of variance`, 1))
cntr.dims <- mutate(cntr.dims[1:10,], dimension = 1:10)
cntr.dims <- mutate(cntr.dims, porc_cum = round(`cumulative percentage of variance`, 0))

gg.cntr <- fviz_eig(pca_model) + 
  ggtitle(label = 'Scree plot') + 
  geom_text(data = cntr.dims, aes(x = dimension, y = porc, label = porc_cum), vjust = -1.0, col = 'grey60') + 
  theme_minimal(base_family = 'Segoe UI') +
  theme(plot.title = element_text(face = 'bold', hjust = 0.5, size = 14), 
        axis.text.y = element_text(angle = 90, hjust = 0.5)) 

gg.cntr
ggsave(plot = gg.cntr, filename = './png/graphs/screeplot_run1.jpg', units = 'in', width = 6, height = 5, dpi = 300, create.dir = T)
dir_create('./tbl/pca')
write.csv(cntr.dims, './tbl/pca/run-1_dims.csv', row.names = FALSE)

## Project baseline
bsln.nrml <- terra::as.data.frame(bsln.nrml, xy = T, na.rm = T)
pca_projection <- predict(pca_model, newdata = bsln.nrml[,3:ncol(bsln.nrml)])
crds.pca.bsln <- as_tibble(cbind(bsln.nrml[,1:2], pca_projection$coord))
rstr.pca.bsln <- terra::rast(crds.pca.bsln, type = 'xyz', crs = 'EPSG:4326')

## Project future 
rstr.pca.ftre <- map(.x = 1:length(ftre.fles), .f = function(i){
  cat('To process: ', i, '\n')
  tble <- ftre.fles[i] %>% rast() %>% as.data.frame(., xy = T, na.rm = T)
  gcme <- basename(dirname(ftre.fles[i]))
  pcap <- predict(pca_model, newdata = tble[,3:ncol(tble)])
  crds <- as_tibble(cbind(tble[,1:2], pcap$coord))
  rstr <- terra::rast(crds, type = 'xyz', crs = 'EPSG:4326')
  names(rstr) <- glue('pca-{gcme}_{1:5}')
  dire <- dirname(ftre.fles[i])
  terra::writeRaster(x = rstr, filename = glue('{dire}/rpca.tif'), overwrite = TRUE)
  rm(tble, gcme, pcap, crds, rstr, dire)
  gc(reset = T) 
})

# To save the results -----------------------------------------------------

##
dir_create('./tif/output/pca')
terra::writeRaster(x = rstr.pca.bsln, filename = './tif/output/pca/pca-bsl_run1.tif', overwrite = TRUE)

##
ftre.fles <- dir_ls('./tif/c6_25/ssp370-2041-2060') %>% map(., dir_ls) %>% unlist() %>% as.character() %>% grep('rpca', ., value = T) 
