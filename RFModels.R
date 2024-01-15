########## RANDOM FOREST MODEL TRAINING ###########
## out of bag error metrics for classifcation and regression 

library(terra) ; library(dplyr) ; library(data.table) ; library(ggplot2); library(ranger)

## dataset of ICESat ATL08 20m Spatial Support Units - joined to RF predictors
dt = fread("data/RFpredictors.csv") %>% 
  filter(., Dist == 0) %>% # filter disturbed
  filter(if_all(ends_with("median_buff.cv"), ~ . < 0.1)) %>% # spectrally 'homogeneous' SSUs < 10% variation in buffer
  select(-contains(c('.cv', 'Dist'))) # remove the disturbance and cv attributes 

#### CANOPY PRESENCE/ABSENCE CLASSIFICATION #####

# balanced sample for training the classification model
set.seed(1234)
dt_samp = dt %>%
  group_by(Canopy) %>%
  slice_sample(n = 50000) %>%
  ungroup()

# RF spectral and topographic predictors
variables = c('Canopy', names(dt_samp[,13:37]))
dt_samp = dt_samp[, variables]
dt_samp$Canopy = as.factor(dt_samp$Canopy)

# Train RF classification using ranger package
# this will take several minutes
canopy_class = ranger::ranger(Canopy ~.,
                      data = dt_samp,
                      num.trees = 1001,
                      importance = 'permutation',
                      scale.permutation.importance = T,
                      probability = T,
                      mtry = 5)

# out of bag classification error
# make confusion matrix of predicted and observed
# 1 = canopy present, 0 = canopy absent
d = data.frame(pred = canopy_class$predictions, obs = dt_samp$Canopy)
d$pred = ifelse(d$pred.1 > 0.5, 1,0) # if RF probability of class 1 > 0.5 assign canopy present class

# function to print overall and class error
errorTab = function(tab) { 
  overall = 1 - (tab[1,1] + tab[2,2])/ (tab[1,1] + tab[2,2] + tab[1,2] + tab[2,1])
  class1 = 1 - (tab[2,2] / (tab[2,2] + tab[2,1]))
  class0 = 1 - (tab[1,1] / (tab[1,1] + tab[1,2]))
  
  return(list = c(OverallError = overall, Class1Error = class1, Class0Error = class0))
  
  }

errorTab(table(d$pred, d$obs))

#### CANOPY HEIGHT REGRESSION #####

dt_samp = dt %>%
  filter(h_canopy_20m < 20) %>% # preferentially train on lower canopies
  filter(., Canopy == 1) # only SSUs with a canopy present

#k-means clustering to stratify sample by height
  set.seed(1234)
  dt_samp$cluster = kmeans(dt_samp$h_canopy_20m, centers = 3)$cluster
  
  # sample by cluster
  set.seed(1234)
  dt_samp = dt_samp %>%
    group_by(cluster) %>%
    slice_sample(n = 30000)
  
  #RF input variables
  variables = c('h_canopy_20m', names(dt_samp[,13:37]))
  dt_samp = dt_samp[, variables]
  
  # Train RF regression using ranger package
  canopy_reg = ranger::ranger(h_canopy_20m ~.,
                      data = dt_samp, 
                      num.trees = 1000,
                      importance = 'permutation',
                      mtry = 10,
                      min.node.size = 1, 
                      keep.inbag = TRUE)

# function to get error stats (MAE, MAPE, RMSE, BIAS, R2) and obs vs predicted plot from RF model
# inputs: rfMod = ranger object, data = input RF dataset
errorMetricsRF = function(rfMod, data) {
  
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

errorMetricsRF(canopy_reg, dt_samp)$Error
errorMetricsRF(canopy_reg, dt_samp)$Plot # the plot will take some time to render
