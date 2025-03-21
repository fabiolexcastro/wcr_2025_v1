


## Fabio Alexander Castro - Llanos 
## Alliance Bioversity - CIAT 
## February 14 th - 2024

## Get the index worldwide

## id	Longitude	Latitude	bioc_1	bioc_10	bioc_11	bioc_12	bioc_13	bioc_14	bioc_15	bioc_16	bioc_17	bioc_18	bioc_19	bioc_2	bioc_21	bioc_22	
## bioc_23	bioc_24	bioc_25	bioc_26	bioc_27	bioc_28	bioc_29	bioc_3	bioc_4	bioc_5	bioc_6	bioc_7	bioc_8	bioc_9
## ntx30	rnge	tmax	tmin	VPD	cdd	ndd	ndd5	prec	srad	t10	ttrm


# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(terra, fs, sf, openxlsx, readxl, geodata, ggspatial, dismo, ggrepel, rnaturalearthdata, rnaturalearth, climateStability, cptcity, spatialEco, parallelDist, tidyr, tidyverse, gridExtra, glue, gtools, readxl, dtw, dtwclust)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# ETP variables -----------------------------------------------------------
etpvars <- function(x){
  
  p <- matrix(nrow = 1, ncol = 9)
  colnames(p) = paste('bio', 25:33, sep = '')
  
  tavg <- x[25:36] ### Temp
  prec <- x[13:24] ### PREC
  pet <- x[1:12]  ### PET
  
  ### if the values are NA the bios are NA
  if(all(is.na(x))) { 
    return(p)
  } else {
    
    window <- function(x)  { 
      lng <- length(x)
      x <- c(x,  x[1:3])
      m <- matrix(ncol = 3, nrow = lng)
      for (i in 1:3) { m[,i] <- x[i:(lng+i-1)] }
      apply(m, MARGIN = 1, FUN = sum)
    }
    
    ### BIO_23: Annual PET
    p[,1] <- sum(pet)
    ### BIO_24: PET seasonality (Coefficient of Variation)
    p[,2] <- cv(pet)
    ### BIO_25: MAX PET
    p[,3] <- max(pet)
    ### BIO_26: Min PET
    p[,4] <- min(pet)
    ### BIO_27: Range of PET (PETmax-PETmin)
    p[,5] <- p[,3]-p[,4]
    
    wet <- window(prec)
    hot <- window(tavg)/3
    pet2 <- c(pet,pet[1:2])
    
    ### BIO_28: PET of wettest quarter
    p[,6] <- sum(pet2[c(which.max(wet):(which.max(wet)+2))])
    ### BIO_29:	PET of driest quarter
    p[,7] <- sum(pet2[c(which.min(wet):(which.min(wet)+2))])
    ### BIO_30:	PET of warmest quarter
    p[,8] <- sum(pet2[c(which.max(hot):(which.max(hot)+2))])
    ### BIO_31:	PET of coldest quarter
    p[,9] <- sum(pet2[c(which.min(hot):(which.min(hot)+2))])
    
  }
  
  round(p, digits = 2)
  return(p)
  
} 

# Load data ---------------------------------------------------------------

## Tabular data 
pnts <- read_csv('./tbl/points_crds.csv', show_col_types = FALSE)
crds <- dplyr::select(pnts, id, Longitude, Latitude)

## Directories (raster data)
root.tmax <- '//catalogue/WFP_ClimateRiskPr1/1.Data/AgERA5/2m_temperature-24_hour_maximum/'
root.tmin <- '//catalogue/WFP_ClimateRiskPr1/1.Data/AgERA5/2m_temperature-24_hour_minimum/'
root.dewp <- '//catalogue/WFP_ClimateRiskPr1/1.Data/AgERA5/2m_dewpoint_temperature/'
root.prec <- '//catalogue/WFP_ClimateRiskPr1/1.Data/AgERA5/Precipitation_flux'
root.srad <- '//CATALOGUE/WFP_ClimateRiskPr1/1.Data/AgERA5/solar_radiation_flux'

## Vector data 
wrld <- vect(ne_countries(returnclass = 'sf', scale = 50))
idn0 <- gadm(country = 'IDN', level = 0, path = './tmpr')

## Masking
mskr <- terra::rast('./tif/base/mask_v1.tif')

# Functions  --------------------------------------------------------------
make.clmt <- function(yr){
  
  ## Filtering
  # yr <- 2015
  
  cat('To process: ', yr, '\n')
  fls.ppt <- dir_ls(root.prec) %>% as.character() %>% grep(paste0('_', yr), ., value = T )
  fls.tmx <- dir_ls(root.tmax) %>% as.character() %>% grep(paste0('_', yr), ., value = T )
  fls.tmn <- dir_ls(root.tmin) %>% as.character() %>% grep(paste0('_', yr), ., value = T )
  fls.dew <- dir_ls(root.dewp) %>% as.character() %>% grep(paste0('_', yr), ., value = T)
  fls.srd <- dir_ls(root.srad) %>% as.character() %>% grep(paste0('_', yr), ., value = T)
  
  ## By each month 
  map(.x = 1:12, .f = function(m){
    
    ## Month
    cat('To process: ', month.abb[m], '\n')
    mnt <- ifelse(m < 10, paste0('0', m), as.character(m))
    
    ## Filtering
    ppt <- grep(paste0('_', yr, mnt), fls.ppt, value = T)
    tmx <- grep(paste0('_', yr, mnt), fls.tmx, value = T)
    tmn <- grep(paste0('_', yr, mnt), fls.tmn, value = T)
    dew <- grep(paste0('_', yr, mnt), fls.dew, value = T)
    srd <- grep(paste0('_', yr, mnt), fls.srd, value = T)
    
    ## To read the raster and convert the units
    ppt <- rast(ppt)
    tmx <- rast(tmx)
    tmn <- rast(tmn)
    dew <- rast(dew)
    srd <- rast(srd)
    
    # ## To make the resample
    ppt.rsm <- terra::crop(ppt, ext(mskr))
    ppt.rsm <- terra::resample(ppt.rsm, mskr, method = 'bilinear')
    ppt.rsm <- terra::mask(ppt.rsm, wrld)
    
    dew.rsm <- terra::crop(dew, ext(mskr))
    dew.rsm <- terra::mask(dew.rsm, ext(mskr))
    dew.rsm <- terra::mask(dew.rsm, wrld)

    tmn.rsm <- terra::crop(tmn, ext(mskr))
    tmn.rsm <- terra::resample(tmn.rsm, mskr, method = 'bilinear')
    tmn.rsm <- terra::mask(tmn.rsm, wrld)

    tmx.rsm <- terra::crop(tmx, ext(mskr))
    tmx.rsm <- terra::resample(tmx.rsm, mskr, method = 'bilinear')
    tmx.rsm <- terra::mask(tmx.rsm, wrld)
    
    srd.rsm <- terra::crop(srd, ext(mskr))
    srd.rsm <- terra::resample(srd.rsm, mskr, method = 'bilinear')
    srd.rsm <- terra::mask(srd.rsm, wrld)
    
    ## To write the raster
    out <- glue('./tif/climate/baseline/daily')
    dir_create(out)
    terra::writeRaster(x = ppt.rsm, filename = glue('./tif/climate/baseline/daily/prec_{yr}-{mnt}.tif'), overwrite = TRUE)
    terra::writeRaster(x = tmn.rsm, filename = glue('./tif/climate/baseline/daily/tmin_{yr}-{mnt}.tif'), overwrite = TRUE)
    terra::writeRaster(x = tmx.rsm, filename = glue('./tif/climate/baseline/daily/tmax_{yr}-{mnt}.tif'), overwrite = TRUE)
    terra::writeRaster(x = dew.rsm, filename = glue('./tif/climate/baseline/daily/dewp_{yr}-{mnt}.tif'), overwrite = TRUE)
    terra::writeRaster(x = srd.rsm, filename = glue('./tif/climate/baseline/daily/srad_{yr}-{mnt}.tif'), overwrite = TRUE)
    cat('Done!\n')
    
    ## To remove 
    rm(ppt.rsm, tmn.rsm, tmx.rsm, ppt, tmx, tmn, mnt)
    rm(dew.rsm, dew)
    gc(reset = TRUE)
    
  })
  
  ## Finish 
  cat('Done: ', yr, '\n')
  
}

# To apply the function ---------------------------------------------------
map(2015:2022, make.clmt)

# Daily to monthly --------------------------------------------------------

##
clma <- dir_ls('./tif/climate/baseline/daily', regexp = '.tif$')

##
yr2mn <- function(yr){
  
  ## Greeping climate
  cat('To process: ', yr, '\n')
  fls <- grep(paste0('_', yr, '-'), clma, value = T)
  ppt <- grep('prec_', fls, value = T)
  tmn <- grep('tmin_', fls, value = T)
  tmx <- grep('tmax_', fls, value = T)
  
  ## By each month 
  rsl <- map(.x = 1:12, .f = function(m){
    
    ## Filtering and read as a raster
    mnt <- ifelse(m < 10, paste0('0', m), as.character(m))
    pr <- rast(grep(paste0(yr, '-', mnt), ppt, value = T))
    tx <- rast(grep(paste0(yr, '-', mnt), tmx, value = T))
    tn <- rast(grep(paste0(yr, '-', mnt), tmn, value = T))
    
    ## To summarise 
    pr <- ifel(pr < 0, 0, pr)
    pr <- sum(pr)
    tx <- mean(tx)
    tn <- mean(tn)
    
    ## To convert 
    tx <- tx - 273.15
    tn <- tn - 273.15
    
    ## To change the names 
    names(pr) <- glue('prec_{yr}-{mnt}')
    names(tx) <- glue('tmax_{yr}-{mnt}')
    names(tn) <- glue('tmin_{yr}-{mnt}')
    
    ## Finish 
    return(list(pr, tn, tx))
    
  })
  
  ppt <- map(rsl, 1) %>% reduce(., c)
  tmn <- map(rsl, 2) %>% reduce(., c)
  tmx <- map(rsl, 3) %>% reduce(., c)
  
  ## To write the raster 
  terra::writeRaster(x = ppt, filename = glue('./tif/climate/baseline/monthly/prec_{yr}.tif'), overwrite = TRUE)
  terra::writeRaster(x = tmn, filename = glue('./tif/climate/baseline/monthly/tmin_{yr}.tif'), overwrite = TRUE)
  terra::writeRaster(x = tmx, filename = glue('./tif/climate/baseline/monthly/tmax_{yr}.tif'), overwrite = TRUE)
  
  ## Finish
  rm(rsl, ppt, tmn, tmx)
  gc(reset = T)
  cat('Done!\n')
  
}

##
map(2015:2022, yr2mn)

# To calculate thermal amplitude -------------------------------------------
root <- './tif/climate/baseline/monthly'
fles <- dir_ls(root)
tmin <- grep('tmin', fles, value = T) %>% as.character() %>% rast()
tmax <- grep('tmax', fles, value = T) %>% as.character() %>% rast()
prec <- grep('prec', fles, value = T) %>% as.character() %>% rast()

## To calculate the differnece 
dfrn <- (tmax - tmin)
names(dfrn) <- gsub('tmax', 'rnge', names(dfrn))
terra::writeRaster(x = dfrn, filename = './tif/climate/baseline/monthly/range-tasm.tif', overwrite = TRUE)
dfrn <- terra::rast('./tif/climate/baseline/monthly/range-tasm.tif')

## Annual range 
names(dfrn)
year <- 2015:2022

dfrn.year <- reduce(map(.x = year, .f = function(i){mean(dfrn[[grep(i, names(dfrn))]])}), c)
names(dfrn.year) <- glue('range-tasm_{2015:2022}')
writeRaster(x = dfrn.year, filename = './tif/climate/baseline/range-tasm_year.tif', overwrite = TRUE)

idn0 <- gadm(country = 'IDN', level = 0, path = './tmpr')
plot(idn0)
plot(dfrn.year[[1]], add = TRUE)

# To exract the values ----------------------------------------------------
rnge.vles <- as_tibble(cbind(crds, terra::extract(dfrn.year, crds[,c('Longitude', 'Latitude')])))
write.csv(rnge.vles, './tbl/points/points_range-tasm_annual.csv', row.names = FALSE)

# To calculate srad average -----------------------------------------------
fles.srad <- dir_ls('./tif/climate/baseline/daily', regexp = '.tif$') %>% as.character() %>% grep('srad', ., value = T)
nmes <- basename(fles.srad)
srad <- map(fles.srad, rast)
srad <- map(.x = 1:length(srad), .f = function(i){
  print(i)
  sr <- srad[[i]] %>% sum()
  return(sr)
})
srad <- reduce(srad, c)
names(srad) <- gsub('.tif$', '', nmes)
terra::writeRaster(x = srad, filename = './tif/climate/baseline/monthly/srad.tif', overwrite = TRUE)

# To calculate t mean  ----------------------------------------------------
tavg <- (tmax + tmin) / 2
names(tavg) <- gsub('tmax_', 'tavg_', names(tavg))
tavg.vles <- as_tibble(cbind(crds, terra::extract(tavg, crds[,c('Longitude', 'Latitude')])))
tavg.vles <- tavg.vles %>% 
  gather(var, value, -c(id, ID, Longitude, Latitude)) %>% 
  mutate(year = str_sub(var, 6, 9)) %>% 
  group_by(id, Longitude, Latitude, ID, year) %>% 
  dplyr::summarise(value = mean(value, na.rm = T)) %>% 
  ungroup() %>% 
  spread(year, value) %>% 
  setNames(c('id', 'Longitude', 'Latitude', 'ID', glue('tavg_{2015:2022}')))

write.csv(tavg.vles, './tbl/points/points_tavg-anual.csv', row.names = FALSE)

# T10 ---------------------------------------------------------------------
calc.t10 <- function(yr){
  
  # yr <- 2015
  
  ## To list the files
  cat('>>> Process: ', yr, '\n')
  root <- './tif/climate/baseline/daily'
  fles <- dir_ls(root, regexp = '.tif$') %>% as.character() %>% grep('/tmin', ., value = T) %>% grep(yr, ., value = T)
  
  map(.x = 1:12, .f = function(m){
    
    ## Filtering the month
    cat('To process: ', month.abb[m], '\n')
    mnt <- ifelse(m < 10, paste0('0', m), as.character(m))
    fle.tmn <- grep(paste0('_', yr, '-', mnt), fles, value = T) 
    rst <- rast(fle.tmn)
    
    ## To calculate the index
    rst <- rst - 273.15
    t10 <- terra::app(x = rst, fun = function(x){ t10 <- sum(x < 10, na.rm = T); return(t10)})
    names(t10) <- glue('t10_{yr}-{mnt}')
    t10 <- terra::mask(t10, wrld)
    
    ## To write the final raster and finish
    terra::writeRaster(x = t10, filename = glue('./tif/index_world/baseline/t10_{yr}-{mnt}.tif'), overwrite = TRUE)
    rm(rst, t10); gc(reset = TRUE)
    cat('Done!\n')
    
  })
  
}

t10 <- map(2015:2022, calc.t10)
t10 <- dir_ls('./tif/index_world/baseline', regexp = '.tif$')
t10 <- as.character(t10)
t10 <- grep('t10_', t10, value = T)
t10 <- rast(t10)

# CDD ---------------------------------------------------------------------
prec.day <- dir_ls('./tif/climate/baseline/daily', regexp = '.tif$')
prec.day <- grep('prec', prec.day, value = T)
prec.day <- as.character(prec.day)
years <- 2015:2022

map(.x = 1:length(years), .f = function(i){
  
  ## To list the files
  cat('>>> Process: ', i, '\n')
  yr <- years[i]
  
  ## To apply by each month
  rst <- map(.x = 1:12, .f = function(m){
    
    ## Filtering the month and reclassified by negative values
    cat('To processing: ', month.abb[m], '\n')
    mnt <- ifelse(m < 10, paste0('0', m), m)
    trr <- rast(grep(paste0('_', yr, '-', mnt, '.'), prec.day, value = T))
    trr[trr < 0] <- 0
    
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
    terra::writeRaster(x = ndd, filename = glue('./tif/index_world/baseline/ndd_{yr}-{mnt}.tif'), overwrite = TRUE)
    terra::writeRaster(x = cdd, filename = glue('./tif/index_world/baseline/cdd_{yr}-{mnt}.tif'), overwrite = TRUE)
    
    ## Return the two final rasters
    rm(trr, ndd, cdd); gc(reset = TRUE)
    cat('Done :)\n')
    
  })
  
})

# NDD 5 -------------------------------------------------------------------
calc.nd5 <- function(yr){
  
  ## To list the files
  cat('>>> Process: ', yr, '\n')
  fls <- dir_ls(root) %>% 
    as.character() %>% 
    grep(paste0('_', yr, '-'), ., value = T) %>% 
    grep('prec', ., value = T) %>% 
    grep('.tif$', ., value = T)
  
  ## To apply by each month
  rst <- map(.x = 1:12, .f = function(m){
    
    ## Filtering the month and reclassified by negative values
    cat('To processing: ', month.abb[m], '\n')
    mnt <- ifelse(m < 10, paste0('0', m), m)
    trr <- rast(grep(paste0('_', yr, '-', mnt, '.'), fls, value = T))
    trr[trr < 0] <- 0
    
    ## Number of dry days < 5 mm
    ndd <- terra::app(x = trr, fun = function(x){ ndd <- sum(x < 5, na.rm = T); return(ndd)})
    names(ndd) <- glue('ndd5_{yr}-{mnt}')
    ndd <- terra::mask(ndd, wrld)
    
    ## To write the rasters 
    terra::writeRaster(x = ndd, filename = glue('./tif/index_world/baseline/ndd5_{yr}-{mnt}.tif'), overwrite = TRUE)
    
    ## Return the two final rasters
    rm(trr, ndd); gc(reset = TRUE)
    cat('Done :)\n')
    
  })
  
}
map(2016:2022, calc.nd5)

cdd <- as.character(dir_ls('./tif/index_world/baseline', regexp = '.tif$'))
ndd <- grep('ndd_', cdd, value = T)
nd5 <- grep('ndd5', cdd,  value = T)
cdd <- grep('cdd', cdd, value = T)

ndd <- rast(ndd)
nd5 <- rast(nd5)
cdd <- rast(cdd)

# VPD ---------------------------------------------------------------------
calc.vpd <- function(yr){
  
  # yr <- 2015
  
  ## To list the files
  cat('>>> Process: ', yr, '\n')
  root <- './tif/climate/baseline/daily'
  fls.tas <- dir_ls(root) %>% grep('.tif$', ., value = T) %>% as.character() %>%  grep(paste0('_', yr), ., value = T) %>% as.character() %>% grep('tm', ., value = T)
  fls.dwp <- dir_ls(root) %>% grep('.tif$', ., value = T) %>% as.character() %>%  grep(paste0('_', yr), ., value = T) %>% as.character() %>% grep('de', ., value = T)
  
  map(.x = 1:12, .f = function(m){
    
    ## Filtering the month
    cat('To process: ', month.abb[m], '\n')
    mnt <- ifelse(m < 10, paste0('0', m), as.character(m))
    fle.tmx <- grep(paste0('_', yr, '-', mnt), fls.tas, value = T) %>% grep('tmax_', ., value = T)
    fle.tmn <- grep(paste0('_', yr, '-', mnt), fls.tas, value = T) %>% grep('tmin_', ., value = T)
    fle.dew <- grep(paste0('_', yr, '-', mnt), fls.dwp, value = T) 
    
    ## To read the raster and convert the units
    rst.tmx <- rast(fle.tmx)
    rst.tmn <- rast(fle.tmn)
    rst.tav <- (rst.tmx + rst.tmn) / 2
    rst.dew <- rast(fle.dew)
    
    ## Es
    es <- 0.61078*exp(17.2694*((rst.tav-273.16)/(rst.tav-35.86))) #es.2 <- app(rst.t, function(x){0.61078*exp(17.2694*((x-273.16)/(x-35.86)))})
    
    ## Ea
    ea <- 0.61078*exp(17.2694*((rst.dew-273.16)/(rst.dew-35.86)))
    
    ## To calculate VPD 
    VPD <- es - ea
    VPD <- mean(VPD)
    names(VPD) <- glue('VPD_{yr}-{mnt}')
    
    VPD %>% crop(., idn0) %>% mask(., idn0) %>% plot()
    
    ## To write the final raster and finish
    terra::writeRaster(x = VPD, filename = glue('./tif/index_world/baseline/vpd_{yr}-{mnt}.tif'), overwrite = TRUE)
    rm(es, ea, VPD); gc(reset = TRUE)
    cat('Done!\n')
    
  })
  
}
map(2015:2022, calc.vpd)
vpd <- as.character(dir_ls('./tif/index_world/baseline', regexp = '.tif$'))
vpd <- grep('vpd', vpd, value = T)
vpd <- rast(vpd)

# TTRM --------------------------------------------------------------------
fles <- dir_ls(root, regexp = '.tif$')
fles <- as.character(fles)
fles.tmax <- fles[grep('tmax_', fles)]
fles.tmin <- fles[grep('tmin_', fles)]
fles.tmax <- rast(fles.tmax)
fles.tmin <- rast(fles.tmin)  

##
make.ttrm <- function(yr){
  
  ## To start the analysis
  cat('To process: ', yr, '\n')
  
  ## Filter the year
  tmax <- fles.tmax[[grep(yr, time(fles.tmax))]]
  tmin <- fles.tmin[[grep(yr, time(fles.tmin))]]
  
  ## To convert the units
  tmax <- tmax - 273.15
  tmin <- tmin - 273.15
  
  ## To calculate the tavg
  tavg <- (tmax + tmin) / 2

  ## To calculate the time thermic variable
  ttr.mnt <- map(.x = 1:12, .f = function(m){
    
    cat('To process the month: ', month.abb[m], '\n')
    mnt <- ifelse(m < 10, paste0('0', m), as.character(m))
    tav <- tavg[[grep(paste0('-', mnt, '-'), time(tavg))]]
    ttr <- tavg - 10
    ttr <- sum(ttr)
    names(ttr) <- glue('ttrm_{yr}-{mnt}')
    return(ttr)
    
  }) %>% 
    reduce(., c)
  
  ttr.sum <- sum(ttr.mnt)
  names(ttr.sum) <- glue('ttr_{yr}')
  
  ## To write the raster 
  terra::writeRaster(x = ttr.sum, filename = glue('./tif/index_world/baseline/ttrm_{yr}.tif'), overwrite = TRUE)
  
  rm(tmax, tmin, tavg, ttr.mnt)
  gc(reset = TRUE)
  return(ttr.sum)

}

##
ttrm <- map(2015:2022, make.ttrm)
ttrm <- rast(grep('ttrm', dir_ls('./tif/index_world/baseline', regexp = '.tif$'), value = T))
ttrm <- mean(ttrm)
names(ttrm) <- 'ttrm'

# To calculate bioclimatic variables --------------------------------------

fles.mnth <- dir_ls('./tif/climate/baseline/monthly', regexp = '.tif$')
fles.mnth <- as.character(fles.mnth)

prec.rstr <- grep('prec', fles.mnth, value = T)
tmin.rstr <- grep('tmin', fles.mnth, value = T)
tmax.rstr <- grep('tmax', fles.mnth, value = T)

## Srad 
srad <- '//catalogue/workspace-cluster9/DATA/ET_SolRad'
srad <- dir_ls(srad, regexp = 'et_solrad')
srad <- mixedsort(srad)
srad <- rast(as.character(srad))

## Function to calculate bioclimatic variables
make.bioc <- function(yr){
  
  # yr <- 2017
  
  ## To read as a raster
  cat('To start the analysis\n')
  prec <- rast(grep(paste0('_', yr, '.tif'), prec.rstr, value = T))
  tmin <- rast(grep(paste0('_', yr, '.tif'), tmin.rstr, value = T))
  tmax <- rast(grep(paste0('_', yr, '.tif'), tmax.rstr, value = T))
  tavg <- (tmax + tmin) / 2
  
  ## To calculate the bioclimatic variables 
  prec.mt <- as.matrix(prec)
  tmin.mt <- as.matrix(tmin)
  tmax.mt <- as.matrix(tmax)
  
  ## Make bioclimatic variables 
  bioc.mt <- dismo::biovars(prec = prec.mt, tmin = tmin.mt, tmax = tmax.mt)
  crds.mt <- terra::crds(prec, na.rm = FALSE)
  bioc.mt <- cbind(crds.mt, bioc.mt)
  bioc.rs <- rast(bioc.mt, type = 'xyz')
  
  ## To estimate the ETP 
  srad <- '//catalogue/workspace-cluster9/DATA/ET_SolRad'
  srad <- dir_ls(srad, regexp = 'et_solrad')
  srad <- mixedsort(srad)
  srad <- rast(as.character(srad))
  
  srad <- terra::resample(srad, prec, method = 'bilinear')
  srad <- terra::crop(srad, prec)
  srad <- terra::mask(srad, prec)
  
  ## To calculate the ETP
  etps <- 0.0013 * 0.408 * srad * (tavg + 17) * (tmax - tmin - 0.0123 * prec) ^ 0.76
  names(etps) <- glue('etps_{1:12}')
  etps <- etps * c(31,29,31,30,31,30,31,31,30,31,30,31)
  
  for(i in 1:12){
    etps[[i]][is.na(etps[[i]])] <- 0
  }
  
  etps <- terra::crop(etps, mskr)
  etps <- terra::mask(etps, mskr)
  
  ## To calculate bioclimatic variables for ETPs
  ETPAndPrec <- cbind(as.matrix(etps), as.matrix(prec), as.matrix(tavg))
  etpbios    <- t(apply(ETPAndPrec, 1, etpvars))
  nms <- paste0('bioc', 21:29)
  bio.etp <- map(.x = 1:ncol(etpbios), .f = function(i){
    print(i)
    rsl <- mskr
    terra::values(rsl) <- etpbios[,i]
    return(rsl)
  })
  bio.etp <- reduce(bio.etp, c)
  names(bio.etp) <- nms
  
  ## Add to just one stack 
  bios <- c(bioc.rs, bio.etp)
  terra::writeRaster(x = bios, filename = glue('tif/climate/baseline/bios_{yr}.tif'), overwrite = TRUE)
  rm(bio.etp, etps, srad, prec.mt, tmin.mt, tmax.mt, bioc.mt, crds.mt, bios)
  cat('Done!\n')
  
}

map(2015:2022, make.bioc)

## To calculate the average 
bios <- dir_ls('./tif/climate/baseline', regexp = '.tif$')
bios <- as.character(grep('bios_', bios, value = T))
bios <- rast(bios)
bios.nmes <- unique(names(bios))

bios.avrg <- map(.x = 1:length(bios.nmes), .f = function(i){
  
  cat('To process: ', bios.nmes[i], '\n')
  bio <- bios[[grep(paste0(paste0(bios.nmes[i], '$'), collapse = '|'), names(bios))]]
  bio <- mean(bio)
  names(bio) <- bios.nmes[i]
  return(bio)  

})

bios.avrg <- bios.avrg %>% reduce(., c)
terra::writeRaster(bios.avrg, './tif/climate/baseline/bios_years-avrg.tif', overwrite = TRUE)
bios.avrg <- setNames(bios.avrg, c(paste0('bioc_', 1:19), paste0('bioc_', 21:29)))

# To process temperature tables -------------------------------------------


##
bios.vles <- gather(as_tibble(cbind(crds, terra::extract(bios.avrg, crds[,c('Longitude',  'Latitude')]))), var, value, -c(id, Longitude, Latitude, ID))
bios.vles <- bios.vles %>% spread(var, value) %>% dplyr::select(-ID)

##
tmin.vles <- gather(as_tibble(cbind(crds, terra::extract(tmin, crds[,c('Longitude',  'Latitude')]))), var, value, -c(id, Longitude, Latitude, ID))
tmin.vles <- separate(data = tmin.vles, col = 'var', sep = '_', into = c('variable', 'date'))
tmin.vles <- mutate(tmin.vles, year = str_sub(date, 1, 4), month = str_sub(date, 6, 7))
tmin.vles <- tmin.vles %>% group_by(id, Longitude, Latitude, variable, month) %>% reframe(value = mean(value, na.rm = T))
tmin.vles <- tmin.vles %>% group_by(id, Longitude, Latitude, variable) %>% reframe(value = mean(value, na.rm = T))

##
tmax.vles <- gather(as_tibble(cbind(crds, terra::extract(tmax, crds[,c('Longitude',  'Latitude')]))), var, value, -c(id, Longitude, Latitude, ID))
tmax.vles <- separate(data = tmax.vles, col = 'var', sep = '_', into = c('variable', 'date'))
tmax.vles <- mutate(tmax.vles, year = str_sub(date, 1, 4), month = str_sub(date, 6, 7))
tmax.vles <- tmax.vles %>% group_by(id, Longitude, Latitude, variable, month) %>% reframe(value = mean(value, na.rm = T))
tmax.vles <- tmax.vles %>% group_by(id, Longitude, Latitude, variable) %>% reframe(value = mean(value, na.rm = T))

##
tavg.vles <- gather(as_tibble(cbind(crds, terra::extract(tavg, crds[,c('Longitude',  'Latitude')]))), var, value, -c(id, Longitude, Latitude, ID))
tavg.vles <- separate(data = tavg.vles, col = 'var', sep = '_', into = c('variable', 'date'))
tavg.vles <- mutate(tavg.vles, year = str_sub(date, 1, 4), month = str_sub(date, 6, 7))
tavg.vles <- tavg.vles %>% group_by(id, Longitude, Latitude, variable, month) %>% reframe(value = mean(value, na.rm = T))
tavg.vles <- tavg.vles %>% group_by(id, Longitude, Latitude, variable) %>% reframe(value = mean(value, na.rm = T))

##
prec.vles <- gather(as_tibble(cbind(crds, terra::extract(prec, crds[,c('Longitude',  'Latitude')]))), var, value, -c(id, Longitude, Latitude, ID))
prec.vles <- separate(data = prec.vles, col = 'var', sep = '_', into = c('variable', 'date'))
prec.vles <- mutate(prec.vles, year = str_sub(date, 1, 4), month = str_sub(date, 6, 7))
prec.vles <- prec.vles %>% group_by(id, Longitude, Latitude, variable, month) %>% reframe(value = mean(value, na.rm = T))
prec.vles <- prec.vles %>% group_by(id, Longitude, Latitude, variable) %>% reframe(value = mean(value, na.rm = T))

##
srad.vles <- gather(as_tibble(cbind(crds, terra::extract(srad, crds[,c('Longitude',  'Latitude')]))), var, value, -c(id, Longitude, Latitude, ID))
srad.vles <- separate(data = srad.vles, col = 'var', sep = '_', into = c('variable', 'date'))
srad.vles <- mutate(srad.vles, year = str_sub(date, 1, 4), month = str_sub(date, 6, 7))
srad.vles <- srad.vles %>% group_by(id, Longitude, Latitude, variable, month) %>% reframe(value = mean(value, na.rm = T))
srad.vles <- srad.vles %>% group_by(id, Longitude, Latitude, variable) %>% reframe(value = mean(value, na.rm = T))

##
t10.vles <- gather(as_tibble(cbind(crds, terra::extract(t10, crds[,c('Longitude',  'Latitude')]))), var, value, -c(id, Longitude, Latitude, ID))
t10.vles <- separate(data = t10.vles, col = 'var', sep = '_', into = c('variable', 'date'))
t10.vles <- mutate(t10.vles, year = str_sub(date, 1, 4), month = str_sub(date, 6, 7))
t10.vles <- t10.vles %>% group_by(id, Longitude, Latitude, variable, month) %>% reframe(value = mean(value, na.rm = T))
t10.vles <- t10.vles %>% group_by(id, Longitude, Latitude, variable) %>% reframe(value = mean(value, na.rm = T))

##
ndd.vles <- gather(as_tibble(cbind(crds, terra::extract(ndd, crds[,c('Longitude',  'Latitude')]))), var, value, -c(id, Longitude, Latitude, ID))
ndd.vles <- separate(data = ndd.vles, col = 'var', sep = '_', into = c('variable', 'date'))
ndd.vles <- mutate(ndd.vles, year = str_sub(date, 1, 4), month = str_sub(date, 6, 7))
ndd.vles <- ndd.vles %>% group_by(id, Longitude, Latitude, variable, year) %>% reframe(value = sum(value, na.rm = T))
ndd.vles <- ndd.vles %>% group_by(id, Longitude, Latitude, variable) %>% reframe(value = mean(value, na.rm = T))

##
cdd.vles <- gather(as_tibble(cbind(crds, terra::extract(cdd, crds[,c('Longitude',  'Latitude')]))), var, value, -c(id, Longitude, Latitude, ID))
cdd.vles <- mutate(cdd.vles, value = ifelse(is.infinite(value), 0, value))
cdd.vles <- separate(data = cdd.vles, col = 'var', sep = '_', into = c('variable', 'date'))
cdd.vles <- mutate(cdd.vles, year = str_sub(date, 1, 4), month = str_sub(date, 6, 7))
cdd.vles <- cdd.vles %>% group_by(id, Longitude, Latitude, variable, year) %>% reframe(value = sum(value, na.rm = T))
cdd.vles <- cdd.vles %>% group_by(id, Longitude, Latitude, variable) %>% reframe(value = mean(value, na.rm = T))

##
cdd.vles <- gather(as_tibble(cbind(crds, terra::extract(cdd, crds[,c('Longitude',  'Latitude')]))), var, value, -c(id, Longitude, Latitude, ID))
cdd.vles <- mutate(cdd.vles, value = ifelse(is.infinite(value), 0, value))
cdd.vles <- separate(data = cdd.vles, col = 'var', sep = '_', into = c('variable', 'date'))
cdd.vles <- mutate(cdd.vles, year = str_sub(date, 1, 4), month = str_sub(date, 6, 7))
cdd.vles <- cdd.vles %>% group_by(id, Longitude, Latitude, variable, year) %>% reframe(value = sum(value, na.rm = T))
cdd.vles <- cdd.vles %>% group_by(id, Longitude, Latitude, variable) %>% reframe(value = mean(value, na.rm = T))

##
nd5.vles <- gather(as_tibble(cbind(crds, terra::extract(nd5, crds[,c('Longitude',  'Latitude')]))), var, value, -c(id, Longitude, Latitude, ID))
nd5.vles <- separate(data = nd5.vles, col = 'var', sep = '_', into = c('variable', 'date'))
nd5.vles <- mutate(nd5.vles, year = str_sub(date, 1, 4), month = str_sub(date, 6, 7))
nd5.vles <- nd5.vles %>% group_by(id, Longitude, Latitude, variable, year) %>% reframe(value = sum(value, na.rm = T))
nd5.vles <- nd5.vles %>% group_by(id, Longitude, Latitude, variable) %>% reframe(value = mean(value, na.rm = T))

##
vpd.vles <- gather(as_tibble(cbind(crds, terra::extract(vpd, crds[,c('Longitude',  'Latitude')]))), var, value, -c(id, Longitude, Latitude, ID))
vpd.vles <- separate(data = vpd.vles, col = 'var', sep = '_', into = c('variable', 'date'))
vpd.vles <- mutate(vpd.vles, year = str_sub(date, 1, 4), month = str_sub(date, 6, 7))
vpd.vles <- vpd.vles %>% group_by(id, Longitude, Latitude, variable, month) %>% reframe(value = mean(value, na.rm = T))
vpd.vles <- vpd.vles %>% group_by(id, Longitude, Latitude, variable) %>% reframe(value = mean(value, na.rm = T))

##
ttr.vles <- gather(as_tibble(cbind(crds, terra::extract(ttrm, crds[,c('Longitude',  'Latitude')]))), variable, value, -c(id, Longitude, Latitude, ID))
ttr.vles <- ttr.vles %>% dplyr::select(-ID)

##
rng.vles <- rnge.vles %>% dplyr::select(-ID) %>% gather(var, value, -c(id, Longitude, Latitude))
rng.vles <- rng.vles %>% separate(col = 'var', into = c('variable', 'year'), sep = '_')
rng.vles <- rng.vles %>% group_by(id, Longitude, Latitude) %>% reframe(value = mean(value, narm = T))
rng.vles <- rng.vles %>% rename(rnge = value)

##
tmin.vles <- tmin.vles %>% spread(variable, value)
tmax.vles <- tmax.vles %>% spread(variable, value)
tavg.vles <- tavg.vles %>% spread(variable, value)
prec.vles <- prec.vles %>% spread(variable, value)
srad.vles <- srad.vles %>% spread(variable, value)
t10.vles  <- t10.vles %>% spread(variable, value)
vpd.vles  <- vpd.vles %>% spread(variable, value)
ndd.vles  <- ndd.vles %>% spread(variable, value)
cdd.vles  <- cdd.vles %>% spread(variable, value)
ttr.vles  <- ttr.vles %>% spread(variable, value)
nd5.vles  <- nd5.vles %>% spread(variable, value)


## Join these tables into only one
clma.vles <- list(bios.vles, tmin.vles, tmax.vles, tavg.vles, prec.vles, srad.vles, t10.vles, rng.vles, vpd.vles, ndd.vles, ttr.vles, nd5.vles, cdd.vles) %>% 
  reduce(., inner_join, by = c('id', 'Longitude', 'Latitude'))

colnames(clma.vles)

write.csv(clma.vles, './tbl/values/values_baseline.csv', row.names = FALSE)


# id	Longitude	Latitude	bioc_1	bioc_10	bioc_11	bioc_12	bioc_13	bioc_14	bioc_15	bioc_16	bioc_17	
# bioc_18	bioc_19	bioc_2	bioc_21	bioc_22	bioc_23	bioc_24	bioc_25	bioc_26	bioc_27	bioc_28	bioc_29	
# bioc_3	bioc_4	bioc_5	bioc_6	bioc_7	bioc_8	bioc_9	ntx30	rnge	tmax	tmin	VPD	cdd	ndd	
# ndd5	prec	srad	t10	ttrm

##
vars <- colnames(clma.vles)[4:ncol(clma.vles)]

# To make the graph  ------------------------------------------------------

make.graph <- function(vrb){
  
  # vrb <- vars[1]
  
  ##
  cat('To process: ', vrb, '\n')
  tbl <- dplyr::select(clma.vles, id, Longitude, Latitude, vrb)
  colnames(tbl) <- c('id', 'Longitude', 'Latitude', 'value')
  
  ##
  library(ggrepel)
  
  ggp <- ggplot(data = tbl, aes(x = id, y = value)) +
    geom_point() +
    labs(x = '', y = vrb) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 0.5),
      axis.text.y = element_text(angle = 90, hjust = 0.5)
    )
  
  ggsave(plot = ggp, filename = glue('./png/graphs/points_values/{vrb}.jpg'), units = 'in', width = 7, height = 5, dpi = 300)
  return(ggp)
  
  
  
}

ggps <- map(vars, make.graph)

# Save all plots to a single PDF
pdf("./png/graphs/points_values/all_plots.pdf", width = 7, height = 5)
walk(ggps, print)
dev.off()

