#########################################################
# Script for Analysis of transformations for automated, #
#       high-throughput clustering of CyTOF data        #
#########################################################

## Required packages
library(flowCore)
library(MASS)
library(forecast)
library(AID)
library(scales)
library(FlowSOM)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(pheatmap)

## Required scripts
source("run_FlowSOM.R")
source("helper_match_evaluate_multiple.R")  


###########################
# LOAD BENCHMARK DATASETS #
###########################

## Function for loading the .fcs files given a path
load_data <- function(path){
  fcs_data <- read.FCS(path, transformation = FALSE, truncate_max_range = FALSE)
  data <- exprs(fcs_data)
  dim(data)
  return(data)
}

## Paths of .fcs files
files <- list(
  Levine32 = "data/Levine_32dim_notransform.fcs",
  Levine13 = "data/Levine_13dim_notransform.fcs",
  Samusik39 = "data/Samusik_all_notransform.fcs"
)

## Load data
data <- lapply(files, load_data)

# lapply(data, dim)

## Select marker columns only
marker_cols <- list(
  Levine32 = 5:36, 
  Levine13 = 1:13, 
  Samusik39 = 9:47
)

## "data.s" is a list of 5 datasets (only marker columns)
data.s <- data
for (i in 1:length(data)){
  data.s[[i]] <- data[[i]][ , marker_cols[[i]] ]
}

## Function to retrieve (a list of) all clusters of a given "dataset"
##  (Clusters ordered decreasingly by size)
##  If "useNaN = TRUE", non-assigned cells (NaN) are included as clusters in the list
##  If "transform = TRUE", returned clusters are TRANSFORMED with ASINH (cofactor = 5)

get_clusters <- function(index, transform = TRUE, useNaN = FALSE){
  dataset <- data.s[[index]]
  labels <- data[[index]][, "label"]
  
  if (useNaN == TRUE){
    # Get decreasing order (including "NaN" and excluding "NA")
    indexes <- names(table(labels, exclude = NA))[order( table(labels, exclude = NA), decreasing = TRUE )]
    
  } else{
    # Get decreasing order (only assigned cells)
    indexes <- order( table(labels), decreasing = TRUE)
  }
  
  clusters <- lapply( indexes, function(u) { dataset[which(labels == u),]} )
  names(clusters) <- indexes
  
  cof <- 5    # CyTOF data
  
  if (transform == TRUE)
    clusters <- lapply( clusters, function(u) asinh(u/cof) )
  
  return(clusters)
}

## Function to retrieve the clusters in a given dataset, given the labels
## The order of events in dataset and their labels must match
return_clusters <- function(dataset, labels, orderit = FALSE){
  
  # Only for flow cytometry (rare) datasets (just one cluster each):
  # Change label of unassigned cells (label = 0) to "NaN"
  if ( names(table(labels)[1]) == "0" ){
    labels[labels == "0"] <- NaN
  } 
  if (orderit == FALSE){
    indexes <- names(table(labels))
  } else{
    indexes <- order( table(labels), decreasing = TRUE)
  }
  
  clusters <- lapply( indexes, function(u) { dataset[which(labels == u),]} )
  names(clusters) <- indexes
  
  return(clusters) 
}


#####################
# ESTIMATION METHOD #
#####################

## This function is used because, for now, goal is to only estimate lambda
make_positive <- function(dataset){
  p.dataset <- dataset
  eps <- 0
  if (sum(dataset < 0)){
    min <- abs(min(dataset))
    eps <- min + 1e-10  # data have to be strictly positive so add a small increment 
                        # to avoid min(dataset) = zero
    p.dataset <- dataset + min
  }
  p.dataset[p.dataset == 0] <- 1e-10  # data have to be strictly positive so add a small increment 
                                      # to avoid min(dataset) = zero
  return(list(data = p.dataset, eps = eps))
}

## Function to estimate the parameters (lambda, epsilon) of the 
## BoxCox transformation necessary to achieve normality of "dataset"
## using a specific "method" for the estimation

## Possible methods:
##  - "loglik": maximize profile log-likelihood
##  - "sw": maximize Shapiro-Wilk test statistic
##          to use this method, size of "dataset" has to be between 3 and 5000

estimate_trans <- function(dataset, method){
  
  lambda <- NA
  r <- make_positive(dataset) # make data positive for the BoxCox estimation
  p.dataset <- r$data # if data are already positive, p.dataset = dataset
  eps <- r$eps        # and eps = 0
  
  if (method == "loglik"){
    out.fore <- BoxCox.lambda(p.dataset, method = "loglik", 
                              lower = 0, upper = 3) # we only consider positive lambdas (for now?)
    lambda <- out.fore
  } else if (method == "sw"){
    if ( length(p.dataset) > 5000 ){        # size of "dataset" has to be between 3 and 5000
      p.dataset <- p.dataset[sample(4500)]  # so draw a sample of size = 4500
    }
    out.aid <- boxcoxnc(p.dataset, 
                        lam = seq(-4,4,0.1)  # such a large range to avoid error "Enlarge the range
                                             # of the lambda in a positive/negative direction"
                        )#plotit = FALSE, verbose = FALSE)
    lambda <- out.aid$lambda.hat
  }
  return( c(lambda = lambda, eps = eps) )
}

### MAPPING BETWEEN BOXCOX AND ASINH

## Calculate sum of absolute deviation between transformations 
## for a specific cofactor (in asinh) and lambda (in x^lambda)

sumsq <- function(cofactor, lambda, x) {
  # sum( abs(asinh(x/cofactor) - BoxCox(x, lambda = lambda) ) )
  # sum( abs(asinh(x/cofactor) - (x^lambda)) )
  if (lambda == 0 ){
    sum( abs(asinh(x/cofactor) - log(x)) )
  } else{
    sum( abs(asinh(x/cofactor) - (x^lambda)) )
  }
}


## Given a lambda, map it to the "closest" cofactor in asinh(x/cofactor)
## Mapping through min of the absolute deviation between both transformations
map_cofactor <- function(lambda, cofactors=1:200, x, plot=FALSE){
  res <- sapply(cofactors, sumsq, lambda = lambda, x = x)
  the_cof <- cofactors[which.min(res)]
  
  if (plot == TRUE){
    plot(cofactors, res, xlab = paste0("Cofactor (Min=",the_cof,")"), 
         ylab = "Abs. Deviation", main = paste0("Cofactor2 from est. Lambda = ",lambda))
    abline(v=the_cof,col="blue")
  }
  return( c(cofactor = the_cof, min = min(res)) )
}

#############################
# SIMULATION TO TEST METHOD #
#############################

## Applies the inverse of the BoxCox transformation (with parameters lambda and epsilon [lambda2 or shift])
## To apply the "original" BoxCox transformation use eps = 0
## This function is taken from and almost identical to InvBoxCox() {forecast}

InvBoxCoxEps <- function(x, lambda, eps){
  if (lambda < 0) x[x > -1/lambda] <- NA # Not sure how to deal with this??
  if (lambda == 0){
    out <- exp(x) - eps
  } else if (lambda > 0) {
    xx <- x * lambda + 1
    out <- ( sign(xx) * abs(xx)^(1/lambda) ) - eps
  }
  return(out)
}

## Outputs a dataset of size "n" following a normal distribution with 
## a specified "mean" and "sd"

get_sample <- function(n, mean, sd){
  dataset <- rnorm(n, mean = mean, sd = sd)
  return(dataset)
}

## SIMULATIONS

# Possible methods ( see description below in estimate_trans() ):
  #   "loglik",
  #   "sw"

# Choose values

method <- "loglik" 
lambda <- 0.25 # testing for one lambda
eps <- 0      # for now, estimate only lambda
nrep <- 1000

set.seed(123)

### APPLY ESTIMATION
df <- data.frame(lambda=rep(NA,nrep),eps=NA,n=NA,mean=NA,sd=NA,nneg=NA,nzer=NA,
                 nneg.t=NA,lower.t=NA,upper.t=NA,lower=NA,upper=NA)

for (i in 1:nrep){
  
  if( i %% 100 == 0 )
    cat("i:",i,"\n")
  
  ## Get dataset
  n <- 10000
  mean <- sample( seq(0,6, length = 5), 1 )
  sd <- sample( seq(0.1,2, length = 5), 1 )
  
  true <- get_sample(n, mean, sd)
  lower.t <- range(true)[1] # min of sample dataset
  upper.t <- range(true)[2] # max of sample dataset
  nneg.t <- sum(true < 0)   # number of negatives of sample dataset
  
  ## Inverse-transform dataset
  dataset <- InvBoxCoxEps(true, lambda, 0)
  
  nneg <- sum(dataset < 0)
  nzer <- sum(dataset == 0)
  lower <- range(dataset)[1]
  upper <- range(dataset)[2]
  
  ## Estimate
  r <- estimate_trans(dataset, method)
  
  df$lambda[i] <- r[1]
  df$eps[i] <- r[2]
  df$mean[i] <- mean
  df$sd[i] <- sd
  df$n[i] <- n
  df$nneg[i] <- nneg
  df$nneg.t[i] <- nneg.t
  df$nzer[i] <- nzer
  df$lower[i] <- lower
  df$upper[i] <- upper
  df$lower.t[i] <- lower.t
  df$upper.t[i] <- upper.t
  
}

df$mean <- as.factor(round(df$mean,3))
df$sd <- as.factor(round(df$sd,3))
df$nneg <- as.factor(df$nneg)

dff <- df
dff$group <- with(dff, interaction(df$mean, df$sd, df$n))
dff$meanneg <- ave(dff$nneg.t, dff$group, mean)
dff$meanneg <- as.factor(round(dff$meanneg,0))

p1 <- ggplot(dff, aes(x = sd, y = lambda)) + geom_violin() + facet_grid(~mean) + coord_cartesian(ylim = c(0, 1))
plot(p1)



#######################################
# RUN AND EVALUATE FLOWSOM CLUSTERING #
#######################################

clus_and_eval <- function(x, labels_true) {
  F1 <- matrix( NA, nrow = length(table(labels_true)), ncol = nruns )
  pre <- matrix( NA, nrow = length(table(labels_true)), ncol = nruns )
  rec <- matrix( NA, nrow = length(table(labels_true)), ncol = nruns )
  labmat <- matrix( NA, nrow = length(labels_true), ncol = nruns )
  for (i in 1:nruns) {
    cat('Run',i,'\n')
    labels <- run_FlowSOM(x, k)
    scores <- helper_match_evaluate_multiple(clus_algorithm = labels, clus_truth = labels_true)
    labmat[,i] <- labels
    F1[,i] <- scores$F1
    pre[,i] <- scores$pr
    rec[,i] <- scores$re
    
  }	
  return(list(labels = labmat,
              F1 = F1,
              pre = pre, 
              rec = rec))
}


###################
# APPLYING METHOD #
###################

## ESTIMATE TRANSFORMATIONS

df.list <- list()
for (j in 1:3){
  clusters <- get_clusters(j, transform = TRUE, useNaN = FALSE)     # transformed
  o.clusters <- get_clusters(j, transform = FALSE, useNaN = FALSE)  # raw
  
  cof <- 5
  
  ## Get marker-wise means for each cluster
  
  ## 'means' is a 'clusters x markers' matrix (after transposing) containing the mean value
  ## of each marker for each cluster
  means <- sapply(clusters, function(u){
    apply(u, 2, mean)
  })
  means <- t(means)
  
  ## Re-order columns of 'means' dataframe (markers) by decreasing maximum mean
  new <- means[, order(apply(means,2,max), decreasing = TRUE)]
  
  hmarkers <- colnames(new)     # Examine all markers
  
  df <- data.frame(lambda=rep(NA,length(hmarkers)),eps=NA,# lambda.abs=NA,
                   cofactor=NA,minsq=NA,cluster=NA,marker=NA,mean=NA,n=NA)
  
  for (i in 1:length(hmarkers)){
    marker <- hmarkers[i]                   # For each marker
    maxmean <- max(means[,marker])          # Select the maximum mean
    pp <- names(which.max(means[,marker]))  # Pick the cluster correspoding
                                            # to maxmean (positive population)
    
    ## Positive population for that marker
    cluster <- clusters[[ pp ]][ ,marker] # Select "cluster" by name of cluster ("pp")
                                          # and name of marker ("marker")
    n <- length(cluster)
    o.cluster <- sinh(cluster)*cof  # Inverse of asinh(x / cof); i.e. original datavalues
    
    r <- estimate_trans(o.cluster, "loglik") # Estimate transf introducing a shift
    m <- map_cofactor(lambda=r[1],x=abs(o.cluster),plot=FALSE)
    
    df$lambda[i] <- r[1]
    df$eps[i] <- r[2]
    df$cofactor[i] <- m[1]
    df$minsq[i] <- m[2]
    df$mean[i] <- maxmean
    df$cluster[i] <- pp
    df$marker[i] <- marker
    df$n[i] <- n
  }
  df.list[[j]] <- df
}


## APPLY TRANSFORMATIONS

dt.list <- list()
for (j in 1:3){
  
  cof <- 5
  
  dataset <- data.s[[j]]          # raw dataset without labels
  df <- df.list[[j]]
  labels <- data[[j]][,"label"]
  
  data_transf <- list()
  
  ## Raw data
  data_transf[[1]] <- dataset  
  
  ## Standard transform (asinh/cof)
  data_transf[[2]] <- asinh(dataset/cof)  
  
  ## Get range of asinh (original) dataset to normalize the scales of the other transformations
  z <- data_transf[[2]]
  r.dat <- apply(z,2,function(u){
    round(range(u),1)
  })
  rownames(r.dat) <- c("min","max")
  
  ## TRANSFORM + NO SCALING
  
  ## Shift + Power transform (est. lambda)
  data_transf[[3]] <- sapply(1:nrow(df),function(u){
    dataset <- data.s[[j]][,u]
    lambda <- df$lambda[ df$marker == colnames(data.s[[j]])[u] ]
    p.dataset <- make_positive(dataset)$data # if data are already positive, p.dataset = dataset
    t.dataset <- p.dataset^lambda
    if (lambda == 0)
      t.dataset <- log(p.dataset)
    return(t.dataset)
  })
  
  ## Asinh with mapped cofactor from est. lambda
  data_transf[[4]] <- sapply(1:nrow(df),function(u){
    dataset <- data.s[[j]][,u]
    marker <- colnames(data.s[[j]])[u]
    cofactor <- df$cofactor[ df$marker == marker ]
    t.dataset <- asinh(dataset/cofactor)
    return(t.dataset)
  })
  
  ## Shift + BoxCox transform
  data_transf[[5]] <- sapply(1:nrow(df),function(u){
    dataset <- data.s[[j]][,u]
    marker <- colnames(data.s[[j]])[u]
    lambda <- df$lambda[ df$marker == marker ]
    p.dataset <- make_positive(dataset)$data # if data are already positive, p.dataset = dataset
    t.dataset <- BoxCox(p.dataset, lambda = lambda)
    return(t.dataset)
  })
  
  ## TRANSFORM + SCALING WITH ASINH 5 RANGE 
  
  ## Shift + Power transform (est. lambda)
  data_transf[[6]] <- sapply(1:nrow(df),function(u){
    marker <- colnames(data.s[[j]])[u]
    dataset <- data_transf[[3]][,u]

    ## for scaling with range
    ra <- as.numeric(r.dat[,colnames(r.dat) == marker])
    t.dataset <- rescale( dataset, to = ra )
    if (sum(ra != round(range(t.dataset),1)) > 0){
      # cat("Dataset",j,"Power","Asinh5:",ra,marker,":",range(t.dataset),
      #     "\n",file = "ranges.txt", append = TRUE)  
    }
    return(t.dataset)
  })
  
  ## Asinh with mapped cofactor from est. lambda
  data_transf[[7]] <- sapply(1:nrow(df),function(u){
    marker <- colnames(data.s[[j]])[u]
    dataset <- data_transf[[4]][,u]
    
    ## for scaling with range
    ra <- as.numeric(r.dat[,colnames(r.dat) == marker])
    t.dataset <- rescale( dataset, to = ra )
    if (sum(ra != round(range(t.dataset),1)) > 0){
      # cat("Dataset",j,"AsinhCof","Asinh5:",ra,marker,":",range(t.dataset),
      #     "\n",file = "ranges.txt", append = TRUE)  
    }
    return(t.dataset)
  })
  
  ## Shift + BoxCox transform
  data_transf[[8]] <- sapply(1:nrow(df),function(u){
    marker <- colnames(data.s[[j]])[u]
    dataset <- data_transf[[5]][,u]
    
    ## for scaling with range
    ra <- as.numeric(r.dat[,colnames(r.dat) == marker])
    t.dataset <- rescale( dataset, to = ra )
    if (sum(ra != round(range(t.dataset),1)) > 0){
      # cat("Dataset",j,"BoxCox","Asinh5:",ra,marker,":",range(t.dataset),
      #     "\n",file = "ranges.txt", append = TRUE)  
    }
    return(t.dataset)
  })
  
  for (u in 2:length(data_transf))
    dimnames(data_transf[[u]]) <- dimnames(data_transf[[1]])
  
  names(data_transf) <- c("Raw", "Asinh5", "Power", "AsinhCof", "BoxCox", 
                          "PowerSc", "AsinhCofSc", "BoxCoxSc")
  
  
  dt.list[[j]] <- data_transf
}


## EVALUATE ON FLOWSOM CLUSTERING

# number of meta-clusters k
k <- 40

# number of repetitions of the analysis
nruns <- 10

output <- list()
for (j in 1:3){
  
  data_transf <- dt.list[[j]]
  labels <- data[[j]][, "label"]  # true labels of dataset

  set.seed(123)
  res <- list()
  cat( system.time( res <- lapply(data_transf, function(x) { clus_and_eval(x, labels) }) ), '\n' )
  output[[j]] <- res
}

names(output) <- c("Levine32", "Levine13", "Samusik39")

save(output, file = paste0("output_",nruns,"run_",k,"k_scaled.RData"))

load("output_10run_40k_scaled.RData")


##########################################
# COMPARE RESULTS AND PLOT DISTRIBUTIONS #
##########################################

names <- c("Levine32", "Levine13", "Samusik39")

## Barplots of results (lambdas and cofactors) from estimation method
pdf("plots/est-barplots.pdf")
par(mfrow = c(3,1))
for (i in 1:3){
  df <- df.list[[i]]
  barplot(table(df$lambda), xlab = "Lambda", ylab = "Frequency", main = names[i])
}
for(i in 1:3){
  df <- df.list[[i]]
  barplot(table(df$cofactor), xlab = "Cofactor", ylab = "Frequency", main = names[i])
}
dev.off()

## Plot comparing lambdas with positive population (log) size and mean
pdf("plots/lambda-mean-size.pdf")
par(mfrow = c(2,3))
for (i in 1:3){
  df <- df.list[[i]]
  plot(df$mean, df$lambda, pch = 20, ylab = "Lambda", xlab = "Pop. Mean", main = names[i])
}
for (i in 1:3){
  df <- df.list[[i]]
  plot(log(df$n), df$lambda, pch = 20, ylab = "Lambda", xlab = "Log(Pop. Size)", main = names[i])
}
dev.off()

## Plots comparing effect of all transformations on clusters for each marker

filenames <- list(
  Levine32 = "plots/Levine32-transf-marker-sc.pdf", 
  Levine13 = "plots/Levine13-transf-marker-sc.pdf", 
  Samusik39 = "plots/Samusik39-transf-marker-sc.pdf"
)


for (j in 1:3){
  plots.l <- list()
  df <- df.list[[j]]
  data_transf <- dt.list[[j]] # get all transformed forms for the given dataset
  labels <- data[[j]][, "label"] # get true gates for the given dataset
  # events in data_transf must be in same 
  # order as in labels
  
  ## For each transformed set, get all clusters
  cdata_transf <- lapply(data_transf, function(u){
    return_clusters(u, labels, orderit = FALSE) 
  })
  
  m.data <- melt(cdata_transf)
  colnames(m.data) <- c("Cell", "Marker", "Value", "Cluster", "Transf")
  n <- length(cdata_transf[[1]])
  m.data$Cluster <- factor(m.data$Cluster, levels = 1:n)
  m.data$Transf <- factor(m.data$Transf)
  # m.data$Marker <- factor(m.data$Marker)
  
  ## reorder to better display comparisons between transformations
  m.data$Transf <- factor(m.data$Transf, levels(m.data$Transf)[c(8,1,2,4,6,3,5,7)])
  
  
  col.vec <- c("#F8766D", "#E38900", "#C49A00", "#99A800", "#53B400", "#00BC56",
               "#00C094", "#00BFC4", "#00B6EB", "#06A4FF", "#A58AFF", "#DF70F8", 
               "#FB61D7", "#FF66A8")
  if (j > 1){}
    col.vec <- c("#F8766D", "#ED813E", "#DE8C00", "#CD9600", "#B79F00", "#9DA700", 
                  "#7CAE00", "#49B500", "#00BA38", "#00BE67", "#00C08B", "#00C1A9",
                  "#00BFC4", "#00BBDC", "#00B4F0", "#00A9FF", "#619CFF", "#9F8CFF",
                  "#C77CFF", "#E36EF6", "#F564E3", "#FF61CC", "#FF64B0", "#FF6C91")
  
  ## Plot one marker per page (all transformations)
  plots.l <- lapply(colnames(data_transf[[1]]), function(u) {
    lambda <- df$lambda[ df$marker == u ]
    cofactor <- df$cofactor[ df$marker == u ]
    m <- m.data[m.data$Marker == u,]
    maxc <- as.numeric(df$cluster[ df$marker == u ])
    col.vec[maxc] <- "black"
    p <- ggplot(m, aes(x = Value, color = Cluster, group = Cluster)) + 
      geom_density() + facet_wrap(~ Transf, scales = "free") +
      scale_color_manual(values=col.vec) +
      theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      ggtitle(paste0(u, " Lambda=", lambda, " Cofactor=", cofactor))
  })
  
  plots <- marrangeGrob(plots.l, nrow=1, ncol=1)
  ggsave(file = filenames[[j]], plots)
}


## Plots comparing effect of all transformations on cluster used for estimation 
## for each marker

filenames2 <- list(
  Levine32 = "plots/Levine32-transf-marker-cluster.pdf", 
  Levine13 = "plots/Levine13-transf-marker-cluster.pdf", 
  Samusik39 = "plots/Samusik39-transf-marker-cluster.pdf"
)


for (j in 1:3){
  plots.l <- list()
  df <- df.list[[j]]
  data_transf <- dt.list[[j]] # get all transformed forms for the given dataset
  labels <- data[[j]][, "label"] # get true gates for the given dataset
  # events in data_transf must be in same 
  # order as in labels
  
  ## For each transformed set, get all clusters
  cdata_transf <- lapply(data_transf, function(u){
    return_clusters(u, labels, orderit = FALSE) 
  })
  
  m.data <- melt(cdata_transf)
  colnames(m.data) <- c("Cell", "Marker", "Value", "Cluster", "Transf")
  n <- length(cdata_transf[[1]])
  m.data$Cluster <- factor(m.data$Cluster, levels = 1:n)
  m.data$Transf <- factor(m.data$Transf)

  ## Plot one marker per page (all transformations)
  plots.l <- lapply(colnames(data_transf[[1]]), function(u) {
    lambda <- df$lambda[ df$marker == u ]
    cofactor <- df$cofactor[ df$marker == u ]
    maxc <- as.numeric(df$cluster[ df$marker == u ])
    m <- m.data[m.data$Marker == u & m.data$Cluster == maxc, ]
    p <- ggplot(m, aes(x = Value, color = Transf)) + 
      geom_density() + coord_cartesian(xlim = c(-5, 50)) +
      theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      ggtitle(paste0(u, " Cluster=", maxc, " Lambda=", lambda, " Cofactor=", cofactor))
  })
  
  plots <- marrangeGrob(plots.l, nrow=1, ncol=1)
  ggsave(file = filenames2[[j]], plots)
}


## Plot differences in F1 scores (single) relative to arcsinh(x/5)
## Violin plots + Points

plots.l <- list()
for (j in 1:3){
  diffs <- lapply(output[[j]],function(u){
    diff <- u$F1 - output[[j]][[2]]$F1
  })
  names(diffs) <- c("Raw", "Asinh5", "Power", "AsinhCof", "BoxCox", 
                    "PowerSc", "AsinhCofSc", "BoxCoxSc")
  labels <- data[[j]][, "label"]  # true labels
  
  # m.data <- melt(diffs[-2])
  m.data <- melt(diffs[3:8])
  colnames(m.data) <- c("Cluster", "Run", "F1_Deviation", "Transf")
  
  m.data <- cbind(m.data, "Scale" = 1)
  m.data$Scale[grep("Sc", m.data$Transf)] <- "Scaled"
  m.data$Scale[m.data$Scale == 1] <- "NotScaled"
  m.data$Transf <- sub("Sc", "", m.data$Transf)
    
  m.data$Cluster <- factor(m.data$Cluster)
  m.data$Transf <- factor(m.data$Transf)
  m.data$Scale <- factor(m.data$Scale)
  
  ## Ordering by cluster number (1,2,3...)
  # p <- ggplot(m.data, aes(x = Cluster, y = F1_Deviation, color = Scale, fill = Scale)) + 
  #   geom_violin() + geom_point(aes(color = Scale), 
  #                              position = position_dodge(width = 0.9),
  #                              size = 1, alpha = 0.6)  +
  #   scale_color_manual(values=c("tomato", "limegreen", "dodgerblue"),
  #                      labels=c("NotScaled", "Asinh5", "Scaled")) +
  #   facet_wrap(~ Transf, scales = "free") + coord_cartesian(ylim = c(-1, 1)) +
  #   geom_hline(aes(yintercept = 0, col = "red"), show.legend = FALSE) +
  #   theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  #   ggtitle(names[j]) + guides(fill=FALSE)
  
  # plots.l <- c(plots.l, list(p))
  
  
  ## To order by cluster size
  ord_clus <- order(table(labels)) # order by increasing cluster size
  m.data$Cluster <- factor(m.data$Cluster, levels(m.data$Cluster)[ord_clus])
  
  p2 <- ggplot(m.data, aes(x = Cluster, y = F1_Deviation, color = Scale, fill = Scale)) +
    geom_violin() + geom_point(aes(color = Scale), 
                               position = position_dodge(width = 0.9),
                               size = 1, alpha = 0.6)  +
    scale_color_manual(values=c("tomato", "limegreen", "dodgerblue"),
                       labels=c("NotScaled", "Asinh5", "Scaled")) +
    facet_wrap(~ Transf, scales = "free") + coord_cartesian(ylim = c(-1, 1)) +
    geom_hline(aes(yintercept = 0, col = "red"), show.legend = FALSE) +
    theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    guides(fill=FALSE) + 
    ggtitle(names[j])
    # ggtitle(paste0(names[j], " ordered by incr. Cluster Size"))
  
  plots.l <- c(plots.l, list(p2))
}

plots <- marrangeGrob(plots.l, nrow=1, ncol=1)
ggsave(file = "plots/clustering-performance-single.pdf", plots)


## Plot differences in F1 scores (averages) relative to arcsinh(x/5)
## Violin plots + Points

plots.l <- list()
for (i in 1:3){
  diffs <- sapply(output[[i]],function(u){
    diff <- apply(u$F1,1,mean) - apply(output[[i]][[2]]$F1,1,mean)
  })
  diffs <- diffs[, 2:8]
  colnames(diffs) <- c("Asinh5", "Power", "AsinhCof", "BoxCox", 
                       "PowerSc", "AsinhCofSc", "BoxCoxSc")
  labels <- data[[i]][, "label"]  # true labels
  row.names(diffs) <- as.character(1:nrow(diffs))
  
  m.diff <- melt(diffs[,-1])
  colnames(m.diff) <- c("Cluster", "Transf", "MeanF1_Deviation")
  m.diff$Cluster <- factor(m.diff$Cluster)
  m.diff$Transf <- factor(m.diff$Transf)
  
  ## Ordering by cluster number (1,2,3...)
  # p <- ggplot(m.diff, aes(x = Cluster, y = MeanF1_Deviation, group = Cluster)) + 
  #   geom_violin() + geom_point(aes(color = Transf), 
  #                              size = 2, alpha = 0.6)  +
  #   coord_cartesian(ylim = c(-1, 0.5)) +
  #   scale_color_manual(values=c("firebrick", "gold3", "forestgreen", "dodgerblue", 
  #                               "darkorchid", "hot pink", "black"),
  #                      labels=c("AsinhCof", "AsinhCofSc", "BoxCox", 
  #                               "BoxCoxSc", "Power", "PowerSc", "Asinh5")) +
  #   geom_hline(aes(yintercept = 0, col = "red"), show.legend = FALSE) +
  #   theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  #   ggtitle(names[i])
  
  # plots.l <- c(plots.l, list(p))
  
  ## To order by cluster size
  ord_clus <- order(table(labels)) # order by increasing cluster size
  m.diff$Cluster <- factor(m.diff$Cluster, levels(m.diff$Cluster)[ord_clus])
  
  p2 <- ggplot(m.diff, aes(x = Cluster, y = MeanF1_Deviation, group = Cluster)) + 
    geom_violin() + geom_point(aes(color = Transf), 
                               size = 2, alpha = 0.6)  +
    coord_cartesian(ylim = c(-1, 0.5)) +
    scale_color_manual(values=c("firebrick", "gold3", "forestgreen", "dodgerblue", 
                                "darkorchid", "hot pink", "black"),
                       labels=c("AsinhCof", "AsinhCofSc", "BoxCox", 
                                "BoxCoxSc", "Power", "PowerSc", "Asinh5")) +
    geom_hline(aes(yintercept = 0, col = "red"), show.legend = FALSE) +
    theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle(names[i])
    # ggtitle(paste0(names[i], " ordered by incr. Cluster Size"))
  
  plots.l <- c(plots.l, list(p2))
}

plots <- marrangeGrob(plots.l, nrow=1, ncol=1)
ggsave(file = "plots/clustering-performance-mean.pdf", plots)


# Heatmap of deviations from asinh5 reference

# plots.l <- list()

pdf(file = "plots/clustering-performance-heatmap.pdf")
for (i in 1:3){
  diffs <- sapply(output[[i]],function(u){
    diff <- apply(u$F1,1,mean) - apply(output[[i]][[2]]$F1,1,mean)
  })
  diffs <- diffs[,3:8]
  colnames(diffs) <- c("Power", "AsinhCof", "BoxCox", 
                       "PowerSc", "AsinhCofSc", "BoxCoxSc")
  rownames(diffs) <- 1:nrow(diffs)
  # m.diff <- melt(diffs)
  # colnames(m.diff) <- c("Cluster", "Transf", "MeanF1_Deviation")
  # m.diff$Cluster <- factor(m.diff$Cluster)
  
  ## Heatmap with ggplot
  # p <- ggplot(m.diff, aes(x = Transf, y = Cluster)) +
  #       geom_tile(aes(fill = MeanF1_Deviation)) + 
  #       scale_fill_gradient2(low = "dark blue", mid = "white", high = "red", 
  #                           limits = c(-1,1)) +
  #       geom_text(aes(label = round(MeanF1_Deviation, 2))) +
  #       theme(axis.text.x = element_text(angle = 90, size = 9)) +
  #       ggtitle(names[i])
  
  # plots.l <- c(plots.l, list(p))
  
  paletteLength <- 200
  myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)

  pheatmap(diffs, cluster_rows = TRUE, cluster_cols = TRUE, 
           color = myColor, breaks = seq(-1, 1, by = 0.01),
           display_numbers = TRUE, number_format = "%.3f",
           number_color = "black", fontsize = 8, fontsize_number = 6,
           main = names[i])
}
dev.off()

# plots <- marrangeGrob(plots.l, nrow=1, ncol=1)
# ggsave(file = "plots/clustering-performance-heatmap.pdf", plots)


## Plot of mean deviations from asinh5 reference

# plots.l <- list()

meandiffs <- matrix(nrow = 3, ncol = 6)
rownames(meandiffs) <- c("Levine32", "Levine13", "Samusik39")
colnames(meandiffs) <- c("Power", "AsinhCof", "BoxCox", 
                         "PowerSc", "AsinhCofSc", "BoxCoxSc")
for (i in 1:3){
  diffs <- sapply(output[[i]],function(u){
    diff <- apply(u$F1,1,mean) - apply(output[[i]][[2]]$F1,1,mean)
  })
  diffs <- diffs[,3:8]
  meandiffs[i,] <- apply(diffs, 2, mean)
}

m.data <- melt(meandiffs)
colnames(m.data) <- c("Dataset", "Transf", "Mean_Deviation")

p <- ggplot(m.data, aes(x = Transf, y = Mean_Deviation, col = Dataset)) + geom_point() +
        coord_cartesian(ylim = c(-0.4,0))

ggsave(file = "plots/clustering-performance-meandev-comparison.pdf", p)


#######################
# VLOG IMPLEMENTATION #
#######################

tmax <- function(a, b, xmax){
  tmax <- log( ( sqrt( xmax^2 + a*xmax + b^2 ) + xmax + a/2 ) / ( a/2 + b ) )
  return(tmax)
}

## VLog function
vlog <- function(x, a, b, zmax){
  xmax <- max(x)
  tmax <- tmax(a,b,xmax)
  z <- sign(x) * (zmax/tmax) * log( ( sqrt( x^2 + a*abs(x) + b^2 ) + abs(x) + a/2 ) / (a/2 + b) )
  return(z)
}

## Inverse VLog function
ivlog <- function(z, a, b, xmax){
  zmax <- max(z)
  tmax <- tmax(a,b,xmax)
  x <- sign(z) * a * sinh( (z/2)*(tmax/zmax) )^2 + b*sinh( z*(tmax/zmax) )
  return(x)
}

########################
# APPLY TRANSFORMATION #
########################

## Apply transformation with paramters a = 0, b = 10, zmax = asinhmax
t.data <- list()
for (j in 1:3){
  
  cof <- 5
  dataset <- data.s[[j]]          # raw dataset without labels
  ref <- asinh(dataset/cof)
  
  mdat <- apply(ref, 2, function(u){
    round(max(u),1)
  })
  
  t.dataset <- sapply(1:ncol(dataset),function(u){
    data <- dataset[,u]
    zmax <- mdat[u]
    vlog(x = data, a = 0, b = 10, zmax = zmax)
  })
  dimnames(t.dataset) <- dimnames(ref)
  
  t.data[[j]] <- list(asinh = ref, vlog = t.dataset)
}

## Plotting resulting distributions
## Facetting by cluster

filenames <- list(
  Levine32 = "plots/Levine32-vlog-marker.pdf", 
  Levine13 = "plots/Levine13-vlog-marker.pdf", 
  Samusik39 = "plots/Samusik39-vlog-marker.pdf"
)

plots.l <- list()
for (j in 1:3){
  data_transf <- t.data[[j]]      # get transformed dataset (zmax = custom)
  labels <- data[[j]][, "label"]  # get true gates for the given dataset
                                  # events in data_transf must be in same 
                                  # order as in labels
  
  ## For each transformed set, get all clusters
  cdata_transf <- lapply(data_transf, function(u){
    return_clusters(u, labels, orderit = FALSE) 
  })
  
  m.data <- melt(cdata_transf)
  colnames(m.data) <- c("Cell", "Marker", "Value", "Cluster", "Transf")
  n <- length(cdata_transf[[1]])
  m.data$Cluster <- factor(m.data$Cluster, levels = 1:n)
  m.data$Transf <- factor(m.data$Transf)
  
  ## Facetting by Cluster (one marker/page)
  plots.l <- lapply(colnames(data_transf[[1]]), function(u) {
    m <- m.data[m.data$Marker == u,]
    p <- ggplot(m, aes(x = Value, color = Transf, group = Transf)) + 
     geom_density() + facet_wrap(~ Cluster, scales = "free") +
     theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
     ggtitle(u)
  })
  
  plots <- marrangeGrob(plots.l, nrow=1, ncol=1)
  ggsave(file = filenames[[j]], plots)
}


## EVALUATE ON FLOWSOM CLUSTERING

# number of meta-clusters k
k <- 40

# number of repetitions of the analysis
nruns <- 10

output <- list()
for (j in 1:3){
  data_transf <- t.data[[j]]
  labels <- data[[j]][, "label"]  # true labels of dataset
  
  set.seed(123)
  res <- list()
  cat( system.time( res <- lapply(data_transf, function(x) { clus_and_eval(x, labels) }) ), '\n' )
  
  output[[j]] <- res
}

names(output) <- c("Levine32", "Levine13", "Samusik39")

save(output, file = paste0("vlog_z=c_",nruns,"run_",k,"k.RData"))

load("vlog_z=c_10run_40k.RData")


###################################
# COMPARE RESULTS ASINH5 VS. VLOG #
###################################

## Deviation values
names <- c("Levine32", "Levine13", "Samusik39")
plots.l <- list()
for (j in 1:3){

  diffs <- lapply(output[[j]],function(u){
    diff <- u$F1 - output[[j]][[1]]$F1
  })
  labels <- data[[j]][, "label"]  # true labels
  
  m.data <- melt(diffs[-1])
  colnames(m.data) <- c("Cluster", "Run", "F1", "Transf")
  m.data$Cluster <- factor(m.data$Cluster)
  m.data$Transf <- factor(m.data$Transf)
  

  ## to order by cluster size
  ord_clus <- order(table(labels)) # order by increasing cluster size
  m.data$Cluster <- factor(m.data$Cluster, levels(m.data$Cluster)[ord_clus])
  
  p <- ggplot(m.data, aes(x = Cluster, y = F1, group = Cluster)) + 
    geom_violin() + geom_point(size = 1.5, alpha = 0.6) +
    coord_cartesian(ylim = c(-1, 1)) +
    stat_summary(aes(group=Cluster), fun.y=mean, geom="point", 
                 size=1.5, color = "dodgerblue") + 
    geom_hline(aes(yintercept = 0, col = "red"), show.legend = FALSE) +
    theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle(names[j])
  
  plots.l <- c(plots.l, list(p))
}

plots <- marrangeGrob(plots.l, nrow=1, ncol=1)
ggsave(file = "plots/clustering-performance-vlog.pdf", plots)


## Single F1 scores
names <- c("Levine32", "Levine13", "Samusik39")
plots.l <- list()
for (j in 1:3){

  res <- lapply(output[[j]], function(u) u$F1)  # retrieve F1 scores
  labels <- data[[j]][, "label"]  # true labels
  
  m.data <- melt(res)
  colnames(m.data) <- c("Cluster", "Run", "F1", "Transf")
  m.data$Cluster <- factor(m.data$Cluster)
  m.data$Transf <- factor(m.data$Transf)

  ## to order by cluster size
  ord_clus <- order(table(labels)) # order by increasing cluster size
  m.data$Cluster <- factor(m.data$Cluster, levels(m.data$Cluster)[ord_clus])
  
  p <- ggplot(m.data, aes(x = Cluster, y = F1, fill = Transf, color = Transf)) + 
    geom_violin() + # geom_point(size = 1.5, alpha = 0.6) +
    coord_cartesian(ylim = c(0, 1)) +
    theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle(names[j])
  
  plots.l <- c(plots.l, list(p))
}

plots <- marrangeGrob(plots.l, nrow=1, ncol=1)
ggsave(file = "plots/clustering-performance-vlog-beta.pdf", plots)



## Heatmap
pdf(file = "plots/clustering-performance-vlog-heatmap.pdf")
mat <- matrix(data = NA, nrow = 24, ncol = 3)
for (j in 1:3){
  diffs <- sapply(output[[j]],function(u){
    diff <- apply(u$F1,1,mean) - apply(output[[j]][[1]]$F1,1,mean)
  })
  n <- length(diffs[,2])
  mat[1:n,j] <- diffs[,2]
}
colnames(mat) <- names(output)
rownames(mat) <- 1:24
paletteLength <- 200
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
  
pheatmap(mat, cluster_rows = FALSE, cluster_cols = FALSE, cellwidth = 50,
          color = myColor, breaks = seq(-1, 1, by = 0.01),
          display_numbers = TRUE, number_format = "%.3f",
          show_rownames = TRUE,
          number_color = "black", fontsize = 8, fontsize_number = 6)
dev.off()


## Plot differences in F1 scores (averages) relative to arcsinh(x/5)
## Mean points

m.data <- melt(mat, na.rm = TRUE)
colnames(m.data) <- c("Cluster", "Dataset", "Mean_Deviation")
m.data$Cluster <- factor(m.data$Cluster)
m.data$Dataset <- factor(m.data$Dataset)

labels <- data[[j]][, "label"]  # true labels
## to order by cluster size
ord_clus <- order(table(labels)) # order by increasing cluster size
m.data$Cluster <- factor(m.data$Cluster, levels(m.data$Cluster)[ord_clus])

p <- ggplot(m.data, aes(x = Cluster, y = Mean_Deviation)) + 
  geom_point(size = 1.5, alpha = 0.6) + 
  facet_grid(~Dataset, scales = "free") + coord_cartesian(ylim = c(-1, 1)) +
  geom_hline(aes(yintercept = 0, col = "red"), show.legend = FALSE) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(file = "plots/clustering-performance-vlog-mean.pdf", p)


## Average results across runs and clusters (unique F1 value per dataset)
mean.dat <- list()
for (j in 1:3){
  labels <- data[[j]][, "label"]  # true labels

  means <- lapply(output[[j]],function(u){
    t <- apply(u$F1, 1, mean) # mean across nruns
    m <- mean(t)   # normal mean across clusters 
    # wm <- weighted.mean(t, weights) # weighted mean across clusters 
  })
  
  mean.dat[[j]] <- means
}
names(mean.dat) <-c("Levine32", "Levine13", "Samusik39")


## Plot mean F1 results
# m.dat <- melt(mean.dat)
# colnames(m.dat) <- c("Mean_F1", "Transf", "Dataset")

# ggplot(m.dat, aes(x = Dataset, y = Mean_F1, col = Transf)) +
#   geom_point(size = 1.5, alpha = 0.6) +
#   coord_cartesian(ylim = c(0, 1)) +
#   ggtitle("Mean F1 score")


################################
# BIEXPONENTIAL IMPLEMENTATION #
################################

fcs_data <- lapply(files, function(u){
  read.FCS(u, transformation = FALSE, truncate_max_range = FALSE)
})
expr_data <- lapply(fcs_data, exprs)
labels_true <- lapply(expr_data, function(u){ u[, "label"] })

data.f <- list()
for (i in 1:3){
  data.f[[i]] <- expr_data[[i]][,marker_cols[[i]]]
}
names(data.f) <- c("Levine32", "Levine13", "Samusik39")
# data_transf <- lapply(data.f, function(u) asinh(u/5))

data_transf <- list()
biexp <- biexponentialTransform("biexp")
for (i in 1:3){
  ## ASINH5 
  trans.a <- asinh(data.f[[i]]/5)
  
  ## BIEXPONENTIAL
  # default parameters: selected to correspond to the hyperbolic sine
  #     (a = 0.5, b = 1, c = 0.5, d = 1, f = 0, w = 0, 
  #     tol = .Machine$double.eps^0.25, maxit = as.integer(5000))
  
  trans.b <- transform( fcs_data[[i]], transformList(colnames(fcs_data[[i]])[marker_cols[[i]]], biexp) )
  data_transf[[i]] <- c( list(asinh5 = trans.a), list( biexp_def = exprs(trans.b)[,marker_cols[[i]]] ) )
}
names(data_transf) <- c("Levine32", "Levine13", "Samusik39")


## Plotting resulting distributions
## Facetting by cluster

filenames <- list(
  Levine32 = "plots/Levine32-biexp-marker.pdf", 
  Levine13 = "plots/Levine13-biexp-marker.pdf", 
  Samusik39 = "plots/Samusik39-biexp-marker.pdf"
)

plots.l <- list()
for (j in 1:3){
  dataset <- data_transf[[j]]      # get transformed dataset (zmax = custom)
  labels <- data[[j]][, "label"]  # get true gates for the given dataset
  # events in data_transf must be in same 
  # order as in labels
  
  ## For each transformed set, get all clusters
  cdata_transf <- lapply(dataset, function(u){
    return_clusters(u, labels, orderit = FALSE) 
  })
  
  m.data <- melt(cdata_transf)
  colnames(m.data) <- c("Cell", "Marker", "Value", "Cluster", "Transf")
  n <- length(cdata_transf[[1]])
  m.data$Cluster <- factor(m.data$Cluster, levels = 1:n)
  m.data$Transf <- factor(m.data$Transf)
  
  ## Facetting by Cluster (one marker/page)
  plots.l <- lapply(colnames(dataset[[1]]), function(u) {
    m <- m.data[m.data$Marker == u,]
    p <- ggplot(m, aes(x = Value, color = Transf, group = Transf)) + 
      geom_density() + facet_wrap(~ Cluster, scales = "free") +
      theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      ggtitle(u)
  })
  
  plots <- marrangeGrob(plots.l, nrow=1, ncol=1)
  ggsave(file = filenames[[j]], plots)
}


## EVALUATE ON FLOWSOM CLUSTERING

# number of meta-clusters k
k <- 40

# number of repetitions of the analysis
nruns <- 10

output.b <- list()
for (j in 1:3){
  data <- data_transf[[j]]
  labels <- labels_true[[j]]  # true labels of dataset
  
  set.seed(123)
  res <- list()
  cat( system.time( res <- lapply(data, function(x) { clus_and_eval(x, labels) }) ), '\n' )
  
  output.b[[j]] <- res
}

names(output.b) <- c("Levine32", "Levine13", "Samusik39")

save(output.b, file = paste0("biexp",nruns,"run_",k,"k.RData"))
load("biexp10run_40k.RData")


##################################
# COMPARE RESULTS BIEXP VS. VLOG #
##################################

## Deviation values
names <- c("Levine32", "Levine13", "Samusik39")
plots.l <- list()
for (j in 1:3){
  
  diffs <- lapply(output.b[[j]],function(u){
    diff <- u$F1 - output.b[[j]][[1]]$F1
  })
  labels <- expr_data[[j]][, "label"]  # true labels
  
  m.data <- melt(diffs[[2]])
  colnames(m.data) <- c("Cluster", "Run", "F1_Deviation")
  m.data$Cluster <- factor(m.data$Cluster)

  ## to order by cluster size
  ord_clus <- order(table(labels)) # order by increasing cluster size
  m.data$Cluster <- factor(m.data$Cluster, levels(m.data$Cluster)[ord_clus])
  
  p <- ggplot(m.data, aes(x = Cluster, y = F1_Deviation, group = Cluster)) + 
    geom_violin() + geom_point(size = 1.5, alpha = 0.6) +
    coord_cartesian(ylim = c(-1, 1)) +
    stat_summary(aes(group=Cluster), fun.y=mean, geom="point", 
                 size=1.5, color = "dodgerblue") + 
    geom_hline(aes(yintercept = 0, col = "red"), show.legend = FALSE) +
    theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle(names[j])
  
  plots.l <- c(plots.l, list(p))
}

plots <- marrangeGrob(plots.l, nrow=1, ncol=1)
ggsave(file = "plots/clustering-performance-biexp.pdf", plots)

## Single F1 scores

plots.l <- list()
for (j in 1:3){
  
  res <- lapply(output.b[[j]], function(u) u$F1)  # retrieve F1 scores
  labels <- data[[j]][, "label"]  # true labels
  
  m.data <- melt(res)
  colnames(m.data) <- c("Cluster", "Run", "F1", "Transf")
  m.data$Cluster <- factor(m.data$Cluster)
  m.data$Transf <- factor(m.data$Transf)
  
  ## to order by cluster size
  ord_clus <- order(table(labels)) # order by increasing cluster size
  m.data$Cluster <- factor(m.data$Cluster, levels(m.data$Cluster)[ord_clus])
  
  p <- ggplot(m.data, aes(x = Cluster, y = F1, fill = Transf, color = Transf)) + 
    geom_violin() + # geom_point(size = 1.5, alpha = 0.6) +
    coord_cartesian(ylim = c(0, 1)) +
    theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle(names[j])
  
  plots.l <- c(plots.l, list(p))
}

plots <- marrangeGrob(plots.l, nrow=1, ncol=1)
ggsave(file = "plots/clustering-performance-biexp-beta.pdf", plots)


## Plot differences in F1 scores (averages) relative to arcsinh(x/5)
## Mean points

mat <- matrix(data = NA, nrow = 24, ncol = 3)
for (j in 1:3){
  diffs <- sapply(output.b[[j]],function(u){
    diff <- apply(u$F1,1,mean) - apply(output[[j]][[1]]$F1,1,mean)
  })
  n <- length(diffs[,2])
  mat[1:n,j] <- diffs[,2]
}
colnames(mat) <- names(output)

m.data <- melt(mat, na.rm = TRUE)
colnames(m.data) <- c("Cluster", "Dataset", "Mean_Deviation")
m.data$Cluster <- factor(m.data$Cluster)
m.data$Dataset <- factor(m.data$Dataset)

p <- ggplot(m.data, aes(x = Cluster, y = Mean_Deviation)) + 
  geom_point(size = 1.5, alpha = 0.6) + 
  facet_grid(~Dataset, scales = "free") + coord_cartesian(ylim = c(-1, 1)) +
  geom_hline(aes(yintercept = 0, col = "red"), show.legend = FALSE) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(file = "plots/clustering-performance-biexp-mean.pdf", p)


## Average results across runs and clusters (unique F1 value per dataset)

mean.dat <- list()
for (j in 1:3){
  labels <- data[[j]][, "label"]  # true labels
  
  means <- sapply(output.b[[j]],function(u){
    t <- apply(u$F1, 1, mean) # mean across nruns
    m <- mean(t)   # normal mean across clusters 
    # wm <- weighted.mean(t, weights) # weighted mean across clusters 
  })
  
  mean.dat[[j]] <- means
}
names(mean.dat) <-c("Levine32", "Levine13", "Samusik39")


## Plot mean F1 results
# m.dat <- melt(mean.dat)
# colnames(m.dat) <- c("Mean_F1", "Transf", "Dataset")

# ggplot(m.dat, aes(x = Dataset, y = Mean_F1, col = Transf)) +
#   geom_point(size = 1.5, alpha = 0.6) +
#   coord_cartesian(ylim = c(0, 1)) +
#   ggtitle("Mean F1 score")
