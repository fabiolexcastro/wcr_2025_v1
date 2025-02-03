

## Fabio Castro - Llanos 
## Alliance Bioversity - CIAT
## Normalice bioclimatic data and make the PCA analysis

# Load libraries ----------------------------------------------------------
require(pacman)
p_load(terra, sf, fs, tidyverse, parallelDist, factoextra, FactoMineR, corrr, ggcorrplot, outliers, dismo, scales, glue, rnaturalearthdata, rnaturalearth, openxlsx)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# Load data ---------------------------------------------------------------

## Spatial data 
bioc <- terra::rast('./tif/tc_baseline/ext_1/bioc.tif')
wrld <- vect(rnaturalearth::ne_countries(scale = 50, type = 'countries'))

plot(bioc[[1]])

## Tabular data 
pnts <- as_tibble(openxlsx::read.xlsx(xlsxFile = './tbl/IMLVT Site info.xlsx'))
pnts <- dplyr::select(pnts, id, Longitude, Latitude)
pnts <- mutate(pnts, id = gsub('Bulegeni\t\t_UGA', 'Bulegeni', id))
write.csv(pnts, './tbl/points_crds.csv', row.names = FALSE)

# To normalice the datasets -----------------------------------------------

## Raster to table 
bioc.tble <- terra::as.data.frame(bioc, xy = T) %>% as_tibble()

## To scale 
nrml.tble <- outliers::scores(x = bioc.tble[,3:ncol(bioc.tble)], type = 'z')
nrml.tble <- cbind(bioc.tble[,1:2], nrml.tble)
nrml.tble <- as_tibble(nrml.tble)

# Normaliced - Table to raster --------------------------------------------
bioc.nrml <- terra::rast(nrml.tble, type = 'xyz')
terra::writeRaster(x = bioc.nrml, filename = './tif/tc_baseline/ext_1/bioc_normal.tif', overwrite = TRUE)
bioc.nrml <- terra::rast('./tif/tc_baseline/ext_1/bioc_normal.tif')

## Raster to table
nrml.bioc <- terra::as.data.frame(bioc.nrml, xy = T)

# PCA analysis ------------------------------------------------------------

## Seed 
set.seed(123)

## Sample size 
land_mask <- !is.na(bioc.nrml[[1]])
valid_cells <- which(values(land_mask) == TRUE)
sample_cells <- sample(valid_cells, min(5000, length(valid_cells)))
sampled_data <- as.data.frame(terra::extract(bioc.nrml, sample_cells))

## Get the coordinates from sample cells
coords <- terra::xyFromCell(bioc.nrml, sample_cells)
sampled_df <- data.frame(coords, sampled_data)

## PCA Model 
pca_model <- PCA(sampled_df[, -c(1,2)], graph = FALSE)

pca_model <- readRDS(file = './rds/pca_model-run1.rrds')

# Visualizar la varianza explicada
cntr.dims <- pca_model$eig %>% as.data.frame() %>% as_tibble() %>% mutate(porc = round(`percentage of variance`, 1))
cntr.dims <- mutate(cntr.dims[1:10,], dimension = 1:10)
cntr.dims <- mutate(cntr.dims, porc_cum = round(`cumulative percentage of variance`, 0))

gg.cntr <- fviz_eig(pca_model) + 
  ggtitle(label = 'Scree plot') + 
  geom_text(data = cntr.dims, aes(x = dimension, y = porc, label = porc_cum), vjust = -1.0, col = 'grey60') + 
  theme_minimal(base_family = 'Segoe UI') +
  theme(plot.title = element_text(face = 'bold', hjust = 0.5, size = 14), 
        axis.text.y = element_text(angle = 90, hjust = 0.5)) 

ggsave(plot = gg.cntr, filename = './png/graphs/screeplot_run1.jpg', units = 'in', width = 6, height = 5, dpi = 300, create.dir = T)

## Project
pca_projection <- predict(pca_model, newdata = nrml.bioc[,3:ncol(nrml.tble)])
crds.pca <- as_tibble(cbind(nrml.tble[,1:2], pca_projection$coord ))
rstr.pca <- terra::rast(crds.pca, type = 'xyz', crs = 'EPSG:4326')

# To save the results -----------------------------------------------------

## Raster
dir_create('./tif/pca')
terra::writeRaster(x = rstr.pca, filename = './tif/pca/pca-bsl_run1.tif', overwrite = TRUE)

## Model
dout <- glue('./rds')
saveRDS(object = pca_model, file = './rds/pca_model-run1.rrds')
