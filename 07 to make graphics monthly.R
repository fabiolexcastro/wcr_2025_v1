

# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(terra, fs, sf, tidyverse, gtools, glue, ggspatial, RColorBrewer, rnaturalearthdata, rnaturalearth)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# Load data ---------------------------------------------------------------
tble <- read_csv('./tbl/values/values_monthly_bsl-ftr.csv', show_col_types = FALSE)
gids <- pull(tble, id) %>% unique()

# To make the graph -------------------------------------------------------
make.graph <- function(gid){

  ## To make the filter
  cat('To process: ', gid, '\n')
  tbl <- filter(tble, id == gid)
  tbl <- rename(tbl, month = Month)
  
  ## Add the abbreviate month column
  lbl <- tibble(month = 1:12, month_abb = month.abb)
  tbl <- inner_join(tbl, lbl, by = 'month')
  tbl <- mutate(tbl, month_abb = factor(month_abb, levels = month.abb))
  
  ## Factor the period 
  tbl <- mutate(tbl,  period = factor(period, levels = c('Baseline', 'Future')))
  rlc <- 35
  
  ## Add tavg
  tbl <- mutate(tbl, tavg = (tmax + tmin) / 2)
  tbl <- mutate(tbl, period = ifelse(period == 'Baseline', 'Current', 'Future'))
  tbl <- mutate(tbl, period = factor(period, levels = c('Current', 'Future')))
  
  ## To make the graph 
  ggp <- ggplot() +
    geom_col(data = tbl, aes(x = month_abb, y = prec, fill = period), position = 'dodge') +
    geom_line(data = tbl, aes(x = month_abb, y = tmin * rlc, col = period, group = period, linetype = 'A')) +
    geom_line(data = tbl, aes(x = month_abb, y = tavg * rlc, col = period, group = period, linetype = 'B')) +
    geom_line(data = tbl, aes(x = month_abb, y = tmax * rlc, col = period, group = period, linetype = 'C')) +
    scale_fill_manual(name = 'Temp. / Prec.',
                      values = c('Current' = '#63BE5C', 'Future' = '#009933'),
                      labels = c('Current' = 'Current', 'Future' = 'Future'),
                      breaks = c('Current', 'Future')) +
    scale_color_manual(name = 'Temperature',
                       values = c('Current' = '#63BE5C', 'Future' = '#009933'),
                       labels = c('Current' = 'Current', 'Future' = 'Future'),
                       breaks = c('Current', 'Future')) +
    scale_linetype_manual(name = ' ', 
                          values = c("A" = 1, 'B' = 2, 'C' = 3), 
                          labels = c("A" = "Min", 'B' = 'Avg', 'C' = 'Max')) +
    scale_y_continuous(sec.axis = sec_axis(~./rlc, name = 'Temperature ºC'),
                       limits = c(0, 1250)) +
    ggtitle(label = gid) +
    labs(x = '', y = 'Precipitation (mm)', caption = 'Source: Terraclimate') +
    theme_minimal() + 
    theme(
      legend.position = 'bottom',
      plot.title = element_text(face = 'bold', hjust = 0.5, size = 18)
    ) +
    guides(linetype = guide_legend(nrow = 2, keywidth = 3, order = 4, title.position = 'top', size = 15),
           color = guide_legend(nrow = 2, keywidth = 3, order = 3, title.position = 'top', size = 15),
           fill = guide_legend(order = 1, title.position = 'top', size = 15),
           size = guide_legend(order = 2, nrow = 2, title.position = 'top', size = 15)) 
  
  ## To save the graph
  ggsave(plot = ggp, filename = glue('./png/graphs/climatogram/{gid}.jpg'), units = 'in', width = 7, height = 6, dpi = 300, create.dir = T)
  cat('Done!\n')
  
}

## 
map(gids, make.graph)
map(gids[16:length(gids)], make.graph)
