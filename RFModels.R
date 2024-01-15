library(terra) ; library(dplyr) ; library(data.table) ; library(ggplot2); library(ranger)

UTM_list = paste0('UTM_', 8:15, 'N')

#read in canopy height ~ BAP + DEM data across all UTM Zones
csv_paths = paste0('E:/UTM_Zones/', UTM_list, '/ICESamples_V2.csv')
dt = lapply(csv_paths, function(x) fread(file = x)) %>%
  rbindlist()

#points overlapping across UTM Zones
#remove duplicate samples
dt = distinct(dt, SampleID, .keep_all = TRUE)

#### CANOPY CLASSIFICATION #####

dt_samp = dt %>%
  filter(., Dist == 0) %>%
  filter(if_all(ends_with("median_buff.cv"), ~ . < 0.1)) %>%
  select(-contains(c('.cv', 'Dist')))

set.seed(1234)
dt_samp = dt_samp %>%
  group_by(Canopy) %>%
  slice_sample(n = 50000) %>%
  ungroup()

# RF variabales
variables = c('Canopy', names(dt_samp[,13:37]))
rf_data = dt_samp[, variables]
rf_data$Canopy = as.factor(rf_data$Canopy)

#get sample IDs for points in model
class.ID = dt_samp$SampleID

#Ranger
canopy_class = ranger(Canopy ~.,
                      data = rf_data,
                      num.trees = 1001,
                      importance = 'permutation',
                      scale.permutation.importance = T,
                      probability = T,
                      mtry = 5)


d = data.frame(pred = canopy_class$predictions, obs = rf_data$Canopy)
d$pred = ifelse(d$pred.1 > 0.7, 1,0)

errorTab = function(tab) { 
  overall = 1 - (tab[1,1] + tab[2,2])/ (tab[1,1] + tab[2,2] + tab[1,2] + tab[2,1])
  class1 = 1 - (tab[2,2] / (tab[2,2] + tab[2,1]))
  class0 = 1 - (tab[1,1] / (tab[1,1] + tab[1,2]))
  
  return(list = c(OverallError = overall, Class1Error = class1, Class0Error = class0))
  
  }

tab = table(d$pred, d$obs)
errorTab(tab)

#### CANOPY HEIGHT REGRESSION #####

#filter by disturbance, undisturbed only
#filter by coefficient of variation, <5% change in spatial buffer
#filter canopy height < 20m 
dt.sub = dt %>%
  filter(., Dist == 0) %>%
  filter(if_all(ends_with("median_buff.cv"), ~ . < 0.1)) %>%
  filter(h_canopy_20m < 20) %>%
  select(-contains(c('.cv', 'Dist'))) %>%
  filter(., Canopy == 1)

#select sample of ICESat points 
#k-means clustering of canopy height 
  set.seed(1234)
  dt.sub$cluster = kmeans(dt.sub$h_canopy_20m, centers = 3)$cluster
  
  #sample by cluster
  set.seed(1234)
  dt_samp = dt.sub %>%
    group_by(cluster) %>%
    slice_sample(n = 10000)
  
  #RF input variables
  variables = c('h_canopy_20m', names(dt.sub[,13:37]))
  rf_data = dt_samp[, variables]
  
#get SampleID of points in regression
reg.ID = dt_samp$SampleID
  
  #regression
  mod = ranger(h_canopy_20m ~.,
                      data = rf_data, 
                      num.trees = 1000,
                      importance = 'permutation',
                      mtry = 10,
                      min.node.size = 1, 
                      keep.inbag = TRUE)

#function to get error stats and obs vs predicted plot from RF model
# rfMod = random forest model
# data = input RF dataset
errorRF = function(rfMod, data) {
  
  df = data.frame(preds = rfMod$predictions,
                   obs = data$h_canopy_20m)
  mod = lm(obs ~ preds, data = df)
  
  x = c(mae = mean(abs(df$preds - df$obs)),
         mape = mean(abs((df$obs - df$preds)/df$obs))*100, 
         rmse = sqrt(mean((df$obs - df$preds)^2)),
         bias = mean(df$preds- df$obs),
         r2 = summary(mod)$r.squared)
  
  plot = ggplot(aes(x = preds, y = obs), data = df)+
    geom_point(size = 2, alpha = 0.3)+
    geom_abline(intercept = 0, slope = 1, linetype = 'dashed', color = 'red')+
    #scale_x_continuous(limits = c(0,14), breaks = seq(0,14, by=2))+
    #scale_y_continuous(limits = c(0,14), breaks = seq(0,14, by=2))+
    xlab('Predicted Max. Canopy Height (m)')+
    ylab('Observed Max. Canopy Height (m)')+
    theme_bw()+
    theme(legend.position = c(0.8,0.2), axis.title = element_text(size = 14),
          axis.text = element_text(size = 12), legend.text = element_text(size = 10),
          legend.title = element_blank())
  
  return(list(Error = x, Plot = plot))
}

###### TRAIN GLOBAL MODEL ########
errorRF(mod, rf_data)$Error
errorRF(mod, rf_data)$Plot
#varImpPlot(mod[1]$model)

save(list = c('mod', 'canopy_class'), file = "E:/RandomForest/RFModels_V2.RData")
load(file = "E:/MS1/RandomForest/RFModels_V2.RData")

##### APPLY MODEL TO ALL ICESAT POINTS #########

dt.pred = dt %>%
  #filter(., Canopy == 1) %>%
  filter(., Dist == 0) %>%
  select(-contains(c('.cv', 'Dist')))

#assign k-means groups 
dt.pred$cluster = kmeans(dt.pred$h_canopy_20m, centers = 3)$cluster 

#canopy presence/absence
variables = c('Canopy', names(dt.sub[,13:37]))

#predict canopy presence on all undisturbed points
can.p = predict(canopy_class, dt.pred[, ..variables])$predictions
can.p.class = ifelse(can.p[,2] > 0.5, 1,0)

#predict canopy height on all undisturbed points
variables = c('h_canopy_20m', names(dt.sub[,13:37]))
can.h = predict(mod, dt.pred[, ..variables])

#merge predictions
merge = as.data.table(cbind(dt.pred, can.p = can.p.class, can.h = can.h$predictions))

#take balanced samples of canopy presence/absence and canopy height
set.seed(1234)
merge.p = merge %>% group_by(Canopy) %>%
  slice_sample(n = 100000)

set.seed(1234)
merge.h = merge %>%
  filter(Canopy == 1) %>%
  group_by(cluster) %>% 
  slice_sample(n=40000)

#Canopy presence error 11%
#01 error: 12%
#00 error: 9%
errorTab(table(merge.p$Canopy, merge.p$can.p))

#Canopy height error stats
#mae       mape       rmse       bias         r2 
#1.5007122 60.4847184  2.0925585  0.4020529  0.5391998 

errorSum = function(preds, obs) {
  
  df = data.frame(preds = preds,
                   obs = obs)
  mod = lm(obs ~ preds, data = df)
  
  x = c(mae = mean(abs(df$preds - df$obs)),
         mape = mean(abs((df$obs - df$preds)/df$obs))*100, 
         rmse = sqrt(mean((df$preds - df$obs)^2)),
         bias = mean(df$preds- df$obs),
         r2 = summary(mod)$r.squared)
  
  return(x)
}

errorSum(merge.h$can.h, merge.h$h_canopy_20m)

#error stats by 
merge.h = merge.h %>% mutate(Canopy_class = case_when(h_canopy_20m <= 2 ~ 'a.<2m',
                                                  h_canopy_20m > 2 & h_canopy_20m <= 5 ~ 'b.2-5m',
                                                  h_canopy_20m > 5 & h_canopy_20m <= 10 ~ 'c.5-10m',
                                                  h_canopy_20m > 10 ~ 'e.>10m'))

merge.h$bias = merge.h$can.h - merge.h$h_canopy_20m
merge.h$mae = abs(merge.h$can.h - merge.h$h_canopy_20m)
merge.h$mape = abs((merge.h$h_canopy_20m - merge.h$can.h)/merge.h$h_canopy_20m )*100
merge.h$rmse = (merge.h$can.h - merge.h$h_canopy_20m)^2

summary_error = merge.h %>%
  filter(Canopy == 1) %>%
  group_by(Canopy_class) %>% 
  summarise(mae = mean(mae),
            mape = mean(mape),
            rmse = sqrt(mean(rmse)),
            bias = mean(bias),
            n = n(),
            prop = n()/nrow(merge.h))

summary_error

### DENSITY PLOT
ggplot() +
  geom_hex(aes(x = can.h, y = h_canopy_20m), data = merge.h) +
  scale_fill_continuous(trans = "log10", high = "#0C1CAE", low = '#DFE1F7') +
  xlim(c(0,25)) + 
  ylim(c(0,25)) +
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed', color = 'red')+
  xlab('Predicted Canopy Height (m)')+
  ylab('Observed Canopy Height (m)')+
  theme_bw()+
  theme(legend.position = c(0.8,0.2), axis.title = element_text(size = 14),
        axis.text = element_text(size = 12), legend.text = element_text(size = 10),
        legend.title = element_blank())
