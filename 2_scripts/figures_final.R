## Code for Multiplex Canada Clinical risk prediction model outputs 

## By Fio Vialard adapted from Cindy Leung Soo's code

## Created on February 20th 


# 1. Set-up ---------------------------------------------------------------

# packages

require(tidyverse)
require(bayesplot)
# install.packages("gganimate")
library(gganimate)
theme_set(bayesplot::theme_default(base_family = "sans"))

df_perf <- read.csv("3_intermediate/projpred/submodels_performance.csv")
df_perf <- df_perf %>%
  mutate(data=case_when(
    data %in% "training" ~ "Development",
    T ~ "Validation"
  ))

colnames <- c("Past drug injection", "Type of past STI test", "Work status", "Education status", "Gender",  "Alcohol consumption", "STBBI testing history", 
              "Monthly income status","Past Syphilis infection", "Reference model")

# 2. Ribbon plot -------------------------------------

## function

ribbon_plot <- function(data,metric, colnames, ylab){
  df<- data %>% filter(performance==metric) %>% #create the dataframe with the needed performance metric
    mutate(model=as.factor(model), 
           value = rep(1:length(colnames), 2)) # assign a numbered value for each model
  plot <- df %>% 
    ggplot(aes(x=value, y=Q50 , fill = data))+ # call the models on the x axis and the performance metric on the y axis
    geom_ribbon(aes(ymin = Q5.5, ymax = Q94.5), alpha = 0.4, colour = NA)+ # add the 89%CrI
    #geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5), alpha = 0.2, colour = NA)+ # add the 95%CrI
    geom_line(aes(colour = data, group = data)) +  # Use colour and fill aesthetics
    scale_fill_manual(values = c("#ACCBF9", "#8170a2")) + # fill line and ribbon with chosen colour
    scale_colour_manual(values = c("#ACCBF9", "#8170a2")) +  
    scale_x_continuous(breaks =1:length(colnames), # assign x-label to predictors in each model
                       labels = str_wrap(colnames,10))+ # wrap the names so they don't take too much space
    # geom_vline(xintercept = which.max(data$Q50), colour = "#36325c", alpha = 0.8)+ # add a line to show the most performant model
    labs(y = paste(ylab),
      # x = "STBBI predictor", # add legend
         fill = "",
         colour = "")+
    theme(legend.position = "none",
          axis.text = element_text(size=24))
  # +
  #   theme(legend.position = "top") # put the legend at the top
  # 
  return(plot)
}


# 3. Point range ----------------------------------------------------------



## function




point.range <-  function(data,metric, colnames, ylab){
  df<- data %>% filter(performance==metric) %>% #create the dataframe with the needed performance metric
    mutate(model=as.factor(model), 
           value = rep(1:length(colnames), 2)) # assign a numbered value for each model
plot  <- df %>% 
    ggplot(aes(x=value, y=Q50, color = data))+ 
    # call the models on the x axis and the performance metric on the y axis
    geom_pointrange(aes(
      ymin= Q5.5,
      ymax = Q94.5),
      alpha=0.4,
      size=1,
      stroke = 1,
      linewidth=0.75,
      position=position_dodge(width = 0.5))+
    scale_color_manual(values=c("#ACCBF9", "#8170a2"))+
    scale_x_continuous(breaks =1:length(colnames), # assign x-label to predictors in each model
                       labels = str_wrap(colnames,10))+ 
    labs(y = paste(ylab),
         # x = "STBBI predictor", # add legend
         fill = "",
         colour = "")+
    theme(legend.position = "none",
          axis.text = element_text(size=24))
  return(plot)
}






## AUC
auc.plot<- point.range(df_perf, "auc", colnames, ylab="AUC")+
  theme(legend.position = "top", legend.text = element_text(size=28))+ # add legend 
  #xlab(NULL)+ # remove the x label
  theme(axis.title = element_text(size=28))

auc.plot

## Sensitivity 
sens.plot <- point.range(df_perf, metric="sensitivity", colnames, ylab = "Sensitivity") +
  ylim(0.75,0.95)+ # add a limit so the plots are on the same y-axis
  xlab(NULL)+
  theme(axis.title = element_text(size=28))
sens.plot

## Specificity
spec.plot <- point.range(df_perf, metric="specificity", colnames, ylab = "Specificity")+
  # ylim(0.1,0.5) +
  xlab("STBBI predictor")+
  theme(axis.title = element_text(size=28))

spec.plot

perf.proj<- bayesplot_grid(auc.plot, sens.plot, spec.plot)

perf.proj

## All 3 plots together 

ggsave("4_outputs/performance_projections_poster.png",perf.proj, device = "png", width = 16, height= 14)



# 3. presentation plots  --------------------------------------------------

auc.plot <- auc.plot + 
  theme(axis.text = element_text(size=18),
        axis.title = element_text(size=20))+
  xlab("Submodels")

auc.plot

ggsave("4_outputs/auc.png", auc.plot, device = "png", width=14, height = 10)

sens.plot <- sens.plot +
  theme(legend.position = "Top",
    axis.text = element_text(size=18),
        axis.title = element_text(size=20))
sens.plot


spec.plot <- spec.plot+
  theme(axis.text = element_text(size=18),
        axis.title = element_text(size=20))

spec.plot

sp.sn.proj<- bayesplot_grid(sens.plot, spec.plot)

sp.sn.proj


ggsave("4_outputs/sp.sn.png", sp.sn.proj, device = "png", width=14, height = 10)
