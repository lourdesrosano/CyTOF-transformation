# Function to run FlowSOM once with pre-transformed data

# inputs:
# data: data matrix, pre-transformed, containing protein marker columns only
# k: number of meta-clusters k, e.g. k <- 40

run_FlowSOM <- function(data, k) {
  # convert data to FlowFrame (required by FlowSOM)
  data <- flowCore::flowFrame(data)
  
  # run FlowSOM
  out <- FlowSOM::ReadInput(data, transform = FALSE, scale = FALSE)
  out <- FlowSOM::BuildSOM(out)
  # out <- FlowSOM::BuildMST(out)  ## not required here
  
  #FlowSOM::PlotStars(out)  ## visualization not required here
  
  # extract cluster labels prior to meta-clustering
  labels_pre <- out$map$mapping[, 1]
  
  # meta-clustering
  out <- ConsensusClusterPlus::ConsensusClusterPlus(t(out$map$codes), maxK = k, plot = 'pdf')
  out <- out[[k]]$consensusClass
  
  # extract cluster labels
  labels <- out[labels_pre]
  
  labels
}


