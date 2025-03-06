
# Functions for Bayesian clinical risk prediction model development and validation 

# By Fio Vialard

# includes functions for:
# 1) Individual Bayesian sn, sp and AUC
# 2) Multiple Bayesian sn, sp, AUC (using for loops of functions 1)
# 3) Plots for histograms (single and grouped), counts 
# 4) Measuring sn, sp, and AUC from training and testing datasets (using loops of functions 2) 
# 5) Missing data histograms with p-value of test (for imputation)
# 6) Prdiction matrix visualization of imputation 
# 7) making quick 2 by 2 tables



## package 
#install.packages("caret")
require(caret)
require(pROC)
require(arsenal)

## 1.Functions for sn, sp, ER, AUC  ------------------------------------------------------------------

## by Cindy Leung Soo

# Function for computing Error rate
error.rate <- function(predictor, response){
  ER <- as.numeric(count(predictor != response)/length(response))
  ER
}

# Function for computing bayesian Error rate
bayes.ER <- 
  function(response, predictor, iter = 100){
    v.ER <- numeric(iter)
    for(i in 1:iter){
      v.ER[i] <- error.rate(response, predictor[i,])
    }
    v.ER
  }

# Function for computing AUC
bayes.auc <- function(response, predictor, iter = 100){
    v.auc <- numeric(iter)
    for(i in 1:iter){
      v.auc[i] <- auc(response, predictor[i,]) %>% suppressMessages() %>% as.numeric()
    }
    v.auc
  }

# Function for computing Sensitivity
bayes.sens <- 
  function(response, predictor, iter = 100){
    v.sensitivity <- numeric(iter)
    for(i in 1:iter){
      v.sensitivity[i] <- sensitivity(as.factor(predictor[i,]), as.factor(response))
    }
    v.sensitivity
  }

# Function for computing specificity
bayes.spec <- 
  function(response, predictor, iter = 100){
    v.specificity <- numeric(iter)
    for(i in 1:iter){
      v.specificity[i] <- specificity(as.factor(predictor[i,]), as.factor(response))
    }
    v.specificity
  }


# 2.Functions for multiple  -------------------------------------------------------
## by Cindy Leung Soo


# Sensitivity
mult.sens<- function(outcome, l.pred){
  l.isens <- list()
for (i in 1:length(l.pred)) {
  l.isens[[i]] <- bayes.sens(outcome, l.pred[[i]]$pred, iter = 2000)
}
# Storing medians and intervals in vector
df.isens <- data.frame(matrix(ncol = 5, nrow = length(l.pred),
                              dimnames = list(colnames)))
colnames(df.isens) <- c("Q50", "Q5.5", "Q94.5", "Q2.5", "Q97.5")
for(i in 1:length(l.pred)){
  df.isens[i,] <- quantile(l.isens[[i]], probs = c(0.5, 0.055, 0.945, 0.025, 0.975)) %>% as.vector()
}
df.isens <- df.isens %>% rownames_to_column()
return(df.isens)
}



# Specificity
mult.spec<- function(outcome, l.pred){

l.ispec <- list()
# Computing specificities
for (i in 1:length(l.pred)) {
  l.ispec[[i]] <- bayes.spec(outcome, l.pred[[i]]$pred, iter = 2000)
}
# Storing medians and intervals in vector
df.ispec <- data.frame(matrix(ncol = 5, nrow = length(l.pred),
                              dimnames = list(colnames)))

colnames(df.ispec) <- c("Q50", "Q5.5", "Q94.5", "Q2.5", "Q97.5")
for(i in 1:length(l.pred)){
  df.ispec[i,] <- quantile(l.ispec[[i]], probs = c(0.5, 0.055, 0.945, 0.025, 0.975)) %>% as.vector()
}

df.ispec <- df.ispec %>% rownames_to_column()
return(df.ispec)
}


# AUC
mult.auc<- function(outcome, l.pred){
  

l.iauc <- list()
# Computing AUC
for (i in 1:length(l.pred)) {
  l.iauc[[i]] <- bayes.auc(outcome, l.pred[[i]]$fitted, iter = 2000)
}
# Storing in dataframe
df.iauc <- data.frame(matrix(ncol = 5, nrow = length(l.pred),
                             dimnames = list(colnames)))
colnames(df.iauc) <- c("Q50", "Q5.5", "Q94.5", "Q2.5", "Q97.5")

for(i in 1:length(l.pred)){
  df.iauc[i,] <- quantile(l.iauc[[i]], probs = c(0.5, 0.055, 0.945, 0.025, 0.975)) %>% as.vector()
}

df.iauc <- df.iauc %>% rownames_to_column()
return(df.iauc)
}

## will need to debug
mult.ER<- function(outcome, l.pred){
  
# Error rate
l.ier <- list()
# Computing error rates
for (i in 1:length(l.pred)) {
  l.ier[[i]] <- bayes.ER(outcome, l.pred[[i]]$pred, iter = 2000)
}
return(l.ier)
}
 
# 3.Plots -------------------------------------------------------------------


histo.plot <- function(data, xvar, label){
data.gather <- data %>% 
  gather(key = xvar, value = "count") %>%
  group_by(xvar) %>%
  summarise(count = sum(count))

plot <- data.gather %>% 
  mutate(xvar = fct_reorder(xvar, desc(count))) %>%
  ggplot(aes(x=xvar, y = count)) +
  geom_bar(stat= "identity", fill = "lightblue") +
  geom_text(label=data.gather$count, nudge_y = -0.5) +
  xlab(label) +
  ylab("Number of CRPMs")+
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=15))

return(plot)
}

histo.plot.switch <- function(data, xvar, label){
  data.gather <- data %>% 
    gather(key = xvar, value = "count") %>%
    group_by(xvar) %>%
    summarise(count = sum(count))
  
  plot <- data.gather %>% 
    mutate(xvar = fct_reorder(xvar, count)) %>%
    ggplot(aes(x=xvar, y = count)) +
    geom_bar(stat= "identity", fill = "#FFFFCCFF") +
    geom_text(label=data.gather$count, nudge_y = -0.5) +
    xlab(label) +
    ylab("Number of CRPMs")+
    scale_y_continuous(breaks= seq(0,12, by=2), # add breaks at every 2 lines
                       expand = expansion(mult = c(0, 0.05)),
                       limits = c(0, 12))+ # remove floating bar
    theme(axis.text = element_text(size=11),
          axis.title = element_text(size=12))+
    coord_flip()
  
  return(plot)
}



count.plot <- function(data, variable, label){
  data.grouped <- data %>%
    select({{ variable }}) %>%
    group_by({{ variable }}) %>%
    summarise(count = n()) 
  
  plot <- data.grouped %>%
    mutate(xvar = fct_reorder({{ variable }}, desc(count))) %>%
    ggplot(aes(x=xvar, y = count)) +
    geom_bar(stat= "identity", fill = "lightblue") +
    geom_text(aes(label=count), nudge_y = -0.5) +
    xlab(label) +
    ylab("Number of CRPMs")+
    theme(axis.text = element_text(size=12),
          axis.title = element_text(size=15))
  
  return(plot)
}


grouped_plot <- function(data, variable, group, label, legend) {
  # Group data and calculate counts
  data.grouped <- data %>%
    select({{ variable }}, {{ group }}) %>% # the {{ }} code is to ensure that the function treats the called variable as a column and not an object
    group_by({{ variable }}, {{ group }}) %>%
    summarise(count = n(), .groups = "drop") 
  
  # Calculate total stacked counts for each variable to reorder bars
  data.grouped <- data.grouped %>%
    group_by({{ variable }}) %>%
    mutate(total_count = sum(count)) %>%
    ungroup() %>%
    mutate(xvar = fct_reorder({{ variable }}, total_count)) # Reorder by total count
  
  # Plot the data
  plot <- data.grouped %>%
    ggplot(aes(x = xvar, y = count, fill = {{ group }})) +
    geom_bar(stat = "identity", position = "stack") +
    geom_text(aes(label = count), 
              position = position_stack(vjust = 0.75)) +
    xlab(label) +
    ylab("Number of CRPMs") +
    scale_y_continuous(breaks = seq(0, 12, by = 2), # set number of breaks
                       expand = expansion(mult = c(0, 0.05)), # remove the floating bar
                       limits = c(0, 12))+ # set limits of axis
    scale_fill_manual(values = c("#007FFFFF", "#4CC3FFFF", "#99EDFFFF"),
                      name = legend) +
    theme(axis.text = element_text(size = 11),
          axis.title = element_text(size = 12),
          legend.title=element_text(size=14, face = "bold"), # change the appearance of the legend
          legend.text = element_text(size=12))+
    coord_flip()

  return(plot)
}
# 4. sub-model eval -------------------------------------------------------

# create a function that re-fits the initial formula
# we need the initial reference model, and the predictor ranking
# start with one model, works 
# add a for loop to do so for all submodels


submodels_train <- function(ref.model, predictor_ranking, new_data){
  submodels_list <- list()  # Create an empty list to store submodels
  
  for(i in 1:length(predictor_ranking)){
    submodel <- NULL
    pred <- NULL
    fitted <- NULL
    predictors <- head(predictor_ranking, i) # for each predictor, call it and the previous ones
    formula <- reformulate(setdiff(predictors, "stbbi.pos"), response = "stbbi.pos") # re-set the formula
    
    submodel <- update(ref.model, formula, verbose = F) # update the reference model with new formula
    pred <- posterior_predict(submodel, newdata = new_data) #create prediction object with newdata
    fitted <- posterior_epred(submodel, newdata = new_data)# create fitted values object with newdata
    new_sub <- submodel
    new_sub$pred <- pred # add new values to the submodel
    new_sub$fitted <- fitted
    submodels_list[[i]] <- new_sub # Add submodel to the list
  }
  
  ref.model$pred <- posterior_predict(ref.model, newdata = new_data)
  ref.model$fitted <- posterior_epred(ref.model, newdata = new_data)
  
  submodels_list[[length(predictor_ranking)+1]] <- ref.model
  
    return(submodels_list)
}

submodels_train2 <- function(ref.model, predictor_ranking, new_data, outcome="stbbi.pos"){
  submodels_list <- list()  # Create an empty list to store submodels
  
  for(i in 1:length(predictor_ranking)){
    submodel <- NULL
    pred <- NULL
    fitted <- NULL
    predictors <- head(predictor_ranking, i) # for each predictor, call it and the previous ones
    formula <- reformulate(setdiff(predictors, outcome), response = outcome) # re-set the formula
    
    submodel <- update(ref.model, formula, verbose = F) # update the reference model with new formula
    pred <- posterior_predict(submodel, newdata = new_data) #create prediction object with newdata
    fitted <- posterior_epred(submodel, newdata = new_data)# create fitted values object with newdata
    new_sub <- submodel
    new_sub$pred <- pred # add new values to the submodel
    new_sub$fitted <- fitted
    submodels_list[[i]] <- new_sub # Add submodel to the list
  }
  
  ref.model$pred <- posterior_predict(ref.model, newdata = new_data)
  ref.model$fitted <- posterior_epred(ref.model, newdata = new_data)
  
  submodels_list[[length(predictor_ranking)+1]] <- ref.model
  
  return(submodels_list)
}



submodels_test <- function(submodels_list, new_data){
  submodels_list_new <- list()
  for(i in 1:length(submodels_list)){
    new_sub <- NULL
    pred <- NULL
    fitted <- NULL
    pred <- posterior_predict(submodels_list[[i]], newdata = new_data) #create prediction object with newdata
    fitted <- posterior_epred(submodels_list[[i]], newdata = new_data)# create fitted values object with newdata
    new_sub <- submodels_list[[i]]
    new_sub$pred <- pred # add new values to the submodel
    new_sub$fitted <- fitted
    submodels_list_new[[i]] <- new_sub
  }
  return(submodels_list_new)
}


# 5. Missing data corr graphs  --------------------------------------------


get.correlation.site <- function(dataset = NA,
                                 variable = NA,
                                 title = "Title") {
  plotlist <-list()
  for(i in 1:length(variable)){
  na <- dataset %>% # determine missingness
    mutate(missingness = as.factor(ifelse(is.na(!!sym(variable[[i]])),
                                          "Missing", "Not missing")))
  
  # t-test
  X1 <- na %>%
    filter(site == "REZO" & missingness == "Missing") %>%
    nrow() # number of event for prop 1
  X2 <- na %>%
    filter(site == "RECAP" & missingness == "Missing") %>%
    nrow()
  n <- nrow(na)
  t.test <- prop.test(x = c(X1, X2), 
                      n = c(n, n), 
                      alternative = "two.sided")
  
  #create label for p-value
  p.value <- ifelse(t.test$p.value<0.05, "p<0.05", "p>0.05")
  # color 
  color <- ifelse(t.test$p.value<0.05, "gold", "white")
  
  
  corr.plot <- na %>%
    ggplot(aes(x = site)) +
    geom_bar(fill = "lightblue") +
    labs(title = title[[i]],
         x = "Site", y = "Frequency") +
    theme_minimal() +
    facet_wrap(~missingness, ncol = 2)+
    geom_label(label = p.value, y = 200, x=2, fill = color) # create plot
  
  plotlist[[i]] <- corr.plot
  }
  
  return(plotlist)
}

get.correlation.inf <- function(dataset = NA,
                                variable = NA,
                                title = "Title") {
  na <- dataset %>% # determine missingness
    mutate(missingness = as.factor(ifelse(is.na(variable),
                                          "Missing", "Not missing")))
  # t-test
  X1 <- na %>%
    filter(stbbi.pos == "Infected" & missingness == "Missing") %>%
    nrow() # number of event for prop 1
  X2 <- na %>%
    filter(stbbi.pos == "Non-infected" & missingness == "Missing") %>%
    nrow()
  n <- nrow(na)
  t.test <- prop.test(x = c(X1, X2), 
                      n = c(n, n), 
                      alternative = "two.sided")
  
  #create label for p-value
  p.value <- ifelse(t.test$p.value<0.05, "p<0.05", "p>0.05")
  # color 
  color <- ifelse(t.test$p.value<0.05, "gold", "white")
  
  
  corr.plot <- na %>%
    ggplot(aes(x = stbbi.pos)) +
    geom_bar(fill = "lightblue") +
    labs(title = title,
         x = "Infection", y = "Frequency") +
    theme_minimal() +
    facet_wrap(~missingness, ncol = 2)+
    geom_label(label = p.value, y = 200, x=2, fill = color) # create plot
  
  return(corr.plot)
}


# 6. Prediction matrix for imputation -------------------------------------

prediction.matrix <- function(imp.pred
){
  names <- colnames(imp.pred)
  
  long.format <- imp.pred %>% 
    as_tibble() %>%
    rowid_to_column(var="Predictors") %>% #assign names to rows
    mutate(Predictors = factor(names, ordered = T, levels = names)) %>%
    gather (key = "Predicted", value = "Z", -1) %>%  # change to long format
    mutate(Predicted = factor(Predicted, ordered = T, levels = names))
  # levels to keep it in the specific order 
  
  pred.plot <- long.format %>% 
    ggplot(aes(Predictors, Predicted, fill = as.factor(Z))) +
    geom_tile() +
    theme(axis.text.x = element_text( 
      angle = 90)) +
    labs(fill= "Did prediction occur?") + 
    scale_fill_brewer()
  
  return(pred.plot)
}



# 7. 2 by 2 tables --------------------------------------------------------

# using arsenal package


two_by_two <-  function(xvar, 
                    yvar,
                    data = df){
  formula <- reformulate(xvar, response = yvar)
  
  two_by_two_tab<- tableby(formula, data)
  
  return(two_by_two_tab)
}
