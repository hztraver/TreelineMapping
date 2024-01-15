########## LVIS VALIDATION #########
## Canopy height error metrics
## Boxplot showing distribution of LVIS height within canopy present vs absent classes

library(data.table); library(dplyr); library(terra); library(ggplot2); library(tidyr)

# 30m gridded LVIS data joined with modeled canopy presence, height & standard error from Random Forest
dt = fread("data/LVIS.csv")

# Binary canopy presence/absence using 50% threshold in Random Forest classification probability
dt$Pres = ifelse(dt$can.prob > 0.5, 1, 0)
# Binary disturbance from Composite 2 Change - pixels can have multiple disturbance events over the time series (i.e. values > 1)
dt$dist2 = ifelse(dt$dist == 0, 0, 1)

#### Canopy height ####

# Filter to pixels with a canopy present and no disturbance history
dt_samp = dt %>% filter(Pres == 1) %>% filter(dist2 == 0)
dt_samp$cluster = kmeans(dt_samp$LVIS98, 3)$cluster # cluster by LVIS 98th percentile height

set.seed(1234) # stratified sample by cluster
dt_samp = dt_samp %>% group_by(., cluster) %>% slice_sample(., n = 10000)

# Pixel error 
dt_samp$bias = dt_samp$height - dt_samp$LVIS98
dt_samp$mae = abs(dt_samp$height - dt_samp$LVIS98)
dt_samp$mape = abs((dt_samp$LVIS98 - dt_samp$height)/dt_samp$LVIS98 )*100
dt_samp$rmse = (dt_samp$height - dt_samp$LVIS98)^2

# Height bins
dt_samp = dt_samp %>% mutate(Canopy_class = case_when(LVIS98 <= 2 ~ 'a.<2m',
                                LVIS98 > 2 & LVIS98 <= 5 ~ 'b.2-5m',
                                LVIS98 > 5 & LVIS98 <= 10 ~ 'c.5-10m',
                                LVIS98 > 10 ~ 'e.>10m'))

# Error metrics by height bins
dt_samp %>%
  group_by(Canopy_class) %>% 
  summarise(mae = mean(mae),
            mape = mean(mape),
            rmse = sqrt(mean(rmse)),
            bias = mean(bias),
            prop.n = n()/nrow(dt_samp))

# Regression plot
ggplot() +
  geom_hex(aes(y = height, x = LVIS98), data = dt_samp) +
  scale_fill_continuous(high = "#0C1CAE", low = '#DFE1F7') +
  xlim(c(0.5,25)) + 
  ylim(c(0.5,25)) +
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed', color = 'red')+
  ylab('Predicted Canopy Height (m)')+
  xlab('LVIS 98. Canopy Height (m)')+
  theme_classic()+
  theme(legend.position = c(0.2,0.8), axis.title = element_text(size = 18),
        axis.text = element_text(size = 16), legend.text = element_text(size = 12),
        legend.title = element_blank())

#### Canopy presence/absence ####

## Boxplot
dt_long = dt %>% 
  pivot_longer(., cols = LVIS95:LVIS100,
               names_to = "metric", values_to = "LVIS_height")

dt_long$Pres = factor(dt_long$Pres)
levels(dt_long$Pres) = c('Canopy Absent', 'Canopy Present')

ggplot(aes(x = factor(metric, levels = c('LVIS100', 'LVIS98', 'LVIS95')), y = LVIS_height, 
           fill = Pres), data = dt_long) +
  geom_boxplot(outlier.alpha = 0) +
  xlab('LVIS Height Metric')+
  ylab('LVIS Height (m)')+
  ylim(c(0,15)) +
  theme_bw()+
  theme(legend.position = 'bottom', axis.title = element_text(size = 18),
        axis.text = element_text(size = 16), legend.text = element_text(size = 12),
        legend.title = element_blank())
