##################
# Load libraries #
##################

library(multiview)

###################################################################
# Multiview CV with Simultaneous Tuning of Alpha, Lambda, and Rho #
###################################################################

cvar.multiview<-function(x_list, 
                         y, 
                         nfolds = 10, 
                         alpha = c(0, seq(0.05, 0.95, 0.05), 1), 
                         rho = c(0, 0.1, 0.25, 0.5, 1, 5, 10), 
                         family = gaussian(), 
                         seed = 1234,
                         type.measure = 'deviance',
                         lambda.choice = 'lambda.min', 
                         verbose = TRUE)
{
  
  ################################
  # Set seed for reproducibility #
  ################################
  
  set.seed(seed)
  foldid <- sample(rep(seq(nfolds), length.out = length(y)))
    
  ####################################
  # Loop over alpha-rho combinations #
  ####################################
  
  alpha_rho<-apply(expand.grid(alpha, rho), 1, paste, collapse="_")
  a1<-lapply(alpha_rho, .cvfunc, xmat = x_list, ymat = y, foldid = foldid, family = family, type.measure = type.measure, verbose = verbose)
  DD <- plyr::ldply(a1, extractMultiviewInfo)
  DD$alpha<-as.numeric(sapply(strsplit(alpha_rho, '_'), "[[" , 1))
  DD$rho<-as.numeric(sapply(strsplit(alpha_rho, '_'), "[[" , 2))
  
  ###################
  # Extract results #
  ###################
  
  if (lambda.choice == 'lambda.1se') {
    min_index<-which(DD$cv.1se==min(DD$cv.1se))
  }
  else {
    min_index<-which(DD$cv.1se==min(DD$cv.1se))
  }
  opt_alpha<-DD$alpha[min_index]
  opt_rho<-DD$rho[min_index]
  opt_cv<-a1[[min_index]]
  if (lambda.choice == 'lambda.1se') {
    opt_lambda<-opt_cv$lambda.1se
  }
  else {
    opt_lambda<-opt_cv$lambda.min
  }
  
  ##########
  # Return #
  ##########
  
  return(list(opt_cv = opt_cv, opt_lambda = opt_lambda, opt_alpha = opt_alpha, opt_rho = opt_rho))
}  

##########################################################################
# Helper Function to Perform Multiview CV for Each Pair of Alpha and Rho #
##########################################################################

.cvfunc <- function(a_r, xmat, ymat, foldid, family, type.measure, verbose)
{
  a<-as.numeric(sapply(strsplit(a_r, '_'), "[[" , 1))
  r<-as.numeric(sapply(strsplit(a_r, '_'), "[[" , 2))
  multiviewFit<-multiview::cv.multiview(x_list = xmat, y = ymat, alpha = a, rho = r, family = family, foldid = foldid, type.measure = type.measure)
  if(verbose) cat('Multiview CV completed for alpha_rho:', a_r, "\n")
  return(multiviewFit)
}

############################################################################################
# Helper Function to Extract Optimal Alpha and Lambda Corresponding to the Lowest CV Error #
############################################################################################

extractMultiviewInfo <- function(object)
{
  # Find Lambdas
  lambda.min <- object$lambda.min
  lambda.1se <- object$lambda.1se
  
  # Determine Optimal Lambda
  which.min <- which(object$lambda == lambda.min)
  which.1se <- which(object$lambda == lambda.1se)
  
  # Summarize Optimal Lambda and Alpha
  data.frame(lambda.min = lambda.min, cv.min = object$cvm[which.min],
             lambda.1se = lambda.1se, cv.1se = object$cvm[which.1se])
}


