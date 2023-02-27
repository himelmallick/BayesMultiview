########################################################
# Generate simulated multi-omics data for benchmarking #
########################################################

# Simulates three inter-connected omics layers using InterSIM
# Always generates 658 features: gene expression (131), methylation (367), and protein (160)
# Generates realistic effect sizes from log-normal distribution using Splatter
# snr is used to control the signal to noise ratio

library(InterSIM)
library(tidyverse)
library(splatter) # Bioconductor package if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager"); BiocManager::install("splatter")

gen_simmba<-function(nsample, # Sample size
                     snr = 1, # Signal to noise ratio
                     p.train = 0.7, # Train-test split
                     de.prob = 0.1, # DE probability across all modalities
                     de.downProb = 0.5, # Down-regulation probability
                     de.facLoc = 1, # DE factor location 
                     de.facScale = 0.4, # DE factor scale
                     nrep = 100, 
                     seed = 1234){
                         
  
  # Set random seed
  set.seed(seed)
  
  # Initialize list of datasets
  trainDat<-testDat<-vector("list", nrep)
  names(trainDat)<- names(testDat)<-paste('Rep', 1: nrep, sep = '_')

  # Loop over nrep
  for (k in 1:nrep){
    
    # Generate X 
    pcl<-trigger_InterSIM(n = nsample)
    X<-pcl$feature_table %>% t() %>% as.matrix()
    
    # Initialize coefficients
    nfeature<-ncol(X)
    
    # Generate realistic effect sizes using Splatter
    de.facs <- splatter:::getLNormFactors(n.facs = nfeature, 
                                          sel.prob = de.prob, 
                                          neg.prob = de.downProb,
                                          fac.loc = de.facLoc, 
                                          fac.scale = de.facScale)
    
    # Convert to LFCs
    beta0<-log2(de.facs) # Assuming LFC is log-normally distributed
    
    # Back-calculate residual sigma2 based on beta0 and snr
    mu = as.matrix(X)%*%beta0
    sigma2 = as.vector(var(mu)/snr)
    
    # Generate Y
    Y = X%*%beta0 + rnorm(nsample)*sigma2 
  
    # Insert Y into the simulated datasets
    pcl$sample_metadata$Y<-as.vector(Y)
    
    # Training / test set splitting
    train<-test<-pcl
    tr.row <- sample(1L:nsample, round(nsample * p.train), replace = FALSE)
    train$sample_metadata<-pcl$sample_metadata[tr.row, , drop = FALSE]
    test$sample_metadata<-pcl$sample_metadata[-tr.row, , drop = FALSE]
    train$feature_table<-pcl$feature_table[, tr.row, drop = FALSE]
    test$feature_table<-pcl$feature_table[, -tr.row, drop = FALSE]
    
    # Store train and test data
    trainDat[[k]]<-train; rm(train)
    testDat[[k]]<-test; rm(test)
    
  }

  # Return synthetic data and simulation parameters
  return(list(trainDat = trainDat,
              testDat = testDat,
              snr = snr,
              p.train = p.train,
              de.prob = de.prob,
              de.downProb = de.downProb, 
              de.facLoc = de.facLoc, 
              de.facScale = de.facScale,
              nrep = nrep,
              seed = seed))
}

##########################
# A wrapper for InterSIM #
##########################

# Currently only accepts the desired number of samples

trigger_InterSIM<-function(n){
  
  # Generate null multi-omics dataset
  sim.data <- InterSIM(n.sample = n, 
                       cluster.sample.prop = c(0.51,0.49), # c(0.5,0.5) does not work!!!
                       delta.methyl = 0, # Must be zero for null data simulation
                       delta.expr = 0, # Must be zero for null data simulation
                       delta.protein = 0) # Must be zero for null data simulation
  
  ####################
  # Process features #
  ####################
  
  names<-c('methyl', 'expr', 'protein')
  list_features<-vector("list", length(names))
  names(list_features)<-names
  list_features[[1]]<-t(sim.data$dat.methyl)
  list_features[[2]]<-t(sim.data$dat.expr)
  list_features[[3]]<-t(sim.data$dat.protein)
  
  #############################
  # Format generated datasets #
  #############################
  
  feature_table<-as.data.frame(Reduce(rbind, list_features))
  colnames(feature_table)<-str_to_title(colnames(feature_table)) 
  
  sample_metadata<-sim.data$clustering.assignment
  colnames(sample_metadata)<-c('subjectID', 'Y')
  sample_metadata$subjectID<-str_to_title(sample_metadata$subjectID)
  rownames(sample_metadata)<-sample_metadata$subjectID
  
  rowID<-rep(names, sapply(list_features, nrow))
  feature_metadata<-cbind.data.frame(featureID = rownames(feature_table), featureType = rowID)
  rownames(feature_metadata)<-feature_metadata$featureID
  
  #####################
  # Save them as list #
  #####################
  
  pcl<-list(feature_table= feature_table, 
            sample_metadata = sample_metadata, 
            feature_metadata = feature_metadata)

  ##########
  # Return #
  ##########
  
  return(pcl)
  
}

# Inspiration for snr calculation
# https://github.com/dingdaisy/cooperative-learning/
# https://github.com/nanxstats/msaenet


# Test Run
DD<-gen_simmba(nsample = 200, nrep = 2)
all(colnames(DD$trainDat$Rep_1$feature_table)==rownames(DD$trainDat$Rep_1$sample_metadata))
all(colnames(DD$trainDat$Rep_2$feature_table)==rownames(DD$trainDat$Rep_2$sample_metadata))
all(rownames(DD$trainDat$Rep_1$feature_table)==rownames(DD$trainDat$Rep_2$feature_tablea))
nsample
