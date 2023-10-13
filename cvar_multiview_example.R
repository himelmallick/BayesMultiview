###################
# Clear workspace #
###################

rm(list = ls())

##################
# Load libraries #
##################

library(multiview)

####################
# Load source code #
####################

############################
# Source codes from GitHub #
############################

devtools::source_url("https://github.com/himelmallick/BayesMultiview/blob/master/cvar_multiview.R?raw=TRUE")

##########################
# Toy Example (Gaussian) #
##########################

# Generate data based on a factor model
set.seed(1)
x = matrix(rnorm(100*20), 100, 20)
z = matrix(rnorm(100*20), 100, 20)
U = matrix(rnorm(100*5), 100, 5)
for (m in seq(5)){
  u = rnorm(100)
  x[, m] = x[, m] + u
  z[, m] = z[, m] + u
  U[, m] = U[, m] + u}
x = scale(x, center = TRUE, scale = FALSE)
z = scale(z, center = TRUE, scale = FALSE)
beta_U = c(rep(0.1, 5))
y = U %*% beta_U + 0.1 * rnorm(100)
x_list<-list(x=x,z=z)
DD1<-cvar.multiview(x_list, y = y)
fit1<-DD1$opt_cv
fit1
plot(fit1)

# Extract coefficients
coef(fit1, s=DD1$opt_lambda, alpha = DD1$opt_alpha, rho = DD1$opt_rho)

# Extract ordered coefficients
coef_ordered(fit1,  s=DD1$opt_lambda, alpha = DD1$opt_alpha, rho = DD1$opt_rho)

# Make predictions
predict(fit1, newx = list(x[1:5, ],z[1:5,]), s=DD1$opt_lambda, alpha = DD1$opt_alpha, rho = DD1$opt_rho)

###################
# Clear workspace #
###################

rm(DD1); rm(fit1)

##########################
# Toy Example (Binomial) #
##########################

by = 1 * (y > median(y)) 
DD2<-cvar.multiview(x_list, by, 
                   family = binomial(),  
                   alpha = c(0, 0.5, 1), 
                   rho = c(0, 0.1, 0.25, 0.5, 1)) 
fit2<-DD2$opt_cv
fit2
plot(fit2)

# Extract coefficients
coef(fit2, s=DD2$opt_lambda, alpha = DD2$opt_alpha, rho = DD2$opt_rho)

# Extract ordered coefficients
coef_ordered(fit2, s=DD2$opt_lambda, alpha = DD2$opt_alpha, rho = DD2$opt_rho)

# Make predictions
predict(fit2, newx = list(x[1:5, ],z[1:5,]), s=DD2$opt_lambda, alpha = DD2$opt_alpha, rho = DD2$opt_rho, type = "response")

###################
# Clear workspace #
###################

rm(DD2); rm(fit2)

#####################################
# Toy Example (Binomial + AUC Loss) #
#####################################

by = 1 * (y > median(y)) 
DD3<-cvar.multiview(x_list, by, 
                    family = binomial(),  
                    alpha = c(0, 0.5, 1), 
                    rho = c(0, 0.1, 0.25, 0.5, 1),
                    type.measure = 'auc') 
fit3<-DD3$opt_cv
fit3
plot(fit3)

# Extract coefficients
coef(fit3, s=DD3$opt_lambda, alpha = DD3$opt_alpha, rho = DD3$opt_rho)

# Extract ordered coefficients
coef_ordered(fit3, s=DD3$opt_lambda, alpha = DD3$opt_alpha, rho = DD3$opt_rho)

# Make predictions
predict(fit3, newx = list(x[1:5, ],z[1:5,]), s=DD3$opt_lambda, alpha = DD3$opt_alpha, rho = DD3$opt_rho, type = "response")

###################################
# Toy Example (IntegratedLearner) #
###################################

library(IntegratedLearner)
library(SuperLearner)
library(tidyverse)

################################
# Prepare dataset in IL format #
################################

# Sample_metadata (Must have variables Y, subjectID, and rownames)
sample_metadata<-as.data.frame(by)
colnames(sample_metadata)<-'Y'
sample_metadata$subjectID<-paste('Subject', 1:nrow(sample_metadata), sep = '_')
rownames(sample_metadata)<-sample_metadata$subjectID

# Feature table (features in rows and samples in columns)
feature_table_t<-as.data.frame(Reduce(cbind, x_list))
feature_table<-as.data.frame(t(feature_table_t))
colnames(feature_table)<-paste('Subject', 1:ncol(feature_table), sep = '_')
rownames(feature_table)<-c(paste('View_1_Feature_', 1:20, sep = ''), paste('View_2_Feature_', 1:20, sep = ''))

# Feature_metadata table
rowID<-c(rep('View_1', 20), rep('View_2', 20))
feature_metadata<-cbind.data.frame(featureID = rownames(feature_table), featureType = rowID)
rownames(feature_metadata)<-rownames(feature_table)

# Final sanity check
all(rownames(feature_metadata)==rownames(feature_table))
all(rownames(sample_metadata)==colnames(feature_table))

########################################
# Run IntegratedLearner (Default BART) #
########################################

IL_BART<-IntegratedLearner(feature_table,
                  sample_metadata, 
                  feature_metadata,
                  family=binomial())
                                    
IL_BART$AUC.train

##############################
# Run IntegratedLearner (RF) #
##############################

IL_RF<-IntegratedLearner(feature_table,
                      sample_metadata, 
                      feature_metadata,
                      base_learner = 'SL.randomForest',
                      family=binomial())

IL_RF$AUC.train

