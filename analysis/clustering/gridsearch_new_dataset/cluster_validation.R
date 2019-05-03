########## Cluster validation #######


#### Libraries needed 



library('hdp')

#source('../../../src/tools.R')     # custom tools function
#source('../../../src/hdp_tools.R') # hdp related functions




#### Get data


# get cytogenetics data
dd_cyto <- read.table("../../../data/initial_dataset/dd_cyto_cut15.tsv", sep = '\t', stringsAsFactors = FALSE, header = TRUE)

# get mutation data
dd_mutation <- read.table("../../../data/initial_dataset/dd_mutation_hotspot_cut10.tsv", sep = '\t', stringsAsFactors = FALSE, header = TRUE)

# merge cytogenetics and mutation data
dd_all <- cbind(dd_mutation, dd_cyto)


# remove patients with no events
dd_all <- dd_all[rowSums(dd_all) !=0,]






#### Functions needed for down-sampling


# quality measure function
quality_measures2 <- function(hdp) {
  
  number_of_chains <- 3
  chain_list <- vector('list', number_of_chains)
  
  for (i in 1:number_of_chains) {
    seed <- i * 100
    
    activated_hdp <- dp_activate(hdp,
                                 dpindex = 1:numdp(hdp), # indices of the DP to activate: we choose all DP (1, ..., n + 1)
                                 initcc  = 16,           # number of starting clusters (every data item is randomly assigned to a cluster to start with)
                                 seed    = seed)
    
    hdp_output <- hdp_posterior(activated_hdp,
                                burnin = 5000, # number of burn-in iterations
                                n      = 400,      # number of posterior samples to collect
                                space  = 20,  # number of iterations between collected samples
                                # cpiter = 1,    # The number of iterations of concentration parameter sampling to perform after each iteration.
                                seed   = seed)
    
    chain_list[[i]] <- hdp_output
    
  }
  
  multi_output <- hdp_multi_chain(chain_list)
  multi_output <- hdp_extract_components(multi_output)
  
  
  dd_predicted <- data.frame(comp_dp_distn(multi_output)$mean[-1,])
  
  
  # change categories colnames
  colnames(dd_predicted) <- paste0('component_', 0:(ncol(dd_predicted)-1))
  components_colnames <- colnames(dd_predicted)
  
  
  # evaluate for each row the predicted component
  dd_predicted['predicted_component'] <- apply(dd_predicted, 1, function(x) {return(which.max(x)-1)})
  
  
  
  # evaluate for each row the maximum probability associated to the predicted component
  dd_predicted['max_proba'] <- apply(dd_predicted[,components_colnames], 1, function(x) {return(max(x))})
  
  
  # get second predicted component
  n_components <- ncol(dd_predicted) - 2 - 1 # without columns 'component_0', 'predicted_component' and 'max_proba'
  
  dd_predicted['second_max_proba'] <- apply(dd_predicted[,0:n_components + 1], 1, function(x) {return(sort(x,partial=n_components)[n_components])})
  
  dd_predicted['second_predicted_component'] <- apply(dd_predicted, 1, function(x) {return(which(x==x['second_max_proba'])[1]-1)})
  
  components_sd <- list()
  for( j in as.integer(levels(factor(dd_predicted$predicted_component)))){
    components_sd[[j+1]] <- sd(dd_predicted$max_proba[dd_predicted$predicted_component == j])
  }
  
  #components_sd <- apply(dd_predicted[, c('max_proba','predicted_component')], 2, function(x){ return(sd(x)) })
  
  # determine which patients are well assigned
  dd_predicted['assignment'] <- apply(dd_predicted[, c('max_proba', 'second_max_proba','second_predicted_component')], 1, function(x) { if(x['second_predicted_component']%in%dd_predicted$predicted_component)
    return( x['max_proba'] > x['second_max_proba']+components_sd[[x['second_predicted_component']+1]])
    else
      return(TRUE)})
  result <- list()
  result[[1]] <- c(n_components, nrow(dd_predicted[dd_predicted$predicted_component==0,]), nrow(dd_predicted[dd_predicted$assignment==TRUE,]))
  result[[2]] <- components_sd
  return(result)
}









# function for down-sampling
down_sampling2 <- function(dd_all, sample_sizes, n_experiment) {
  
  rlist <- list()
  result <- data.frame('sample_size' = integer(),
                       'n_components' = integer(),
                       'n_component_0' = integer(),
                       'assignment' = integer())
  
  compo_sd <- data.frame('sample_size'=integer(),
                         'component_0'=integer(),
                         'component_1'=integer(),
                         'component_2'=integer(),
                         'component_3'=integer(),
                         'component_4'=integer(),
                         'component_5'=integer(),
                         'component_6'=integer(),
                         'component_7'=integer(),
                         'component_8'=integer(),
                         'component_9'=integer(),
                         'component_10'=integer(),
                         'component_11'=integer(),
                         'component_12'=integer()
  )
  
  
  for (i in sample_sizes) {
    
    for (j in 1:n_experiment) {
      
      set.seed( i + j +5) 
      selected_rows <- sample(1:nrow(dd_all), size = i)
      dd_all_sampled <- dd_all[selected_rows,]
      
      # initialise hdp
      hdp <- initialise_hdp(dd_all_sampled)
      
      # get quality measures
      info <- quality_measures2(hdp)
      
      result[nrow(result)+1,] <- c(i, info[[1]])
      print(c(i,info[[1]]))
      
      temp <- c(i)
      for (k in 1:length(info[[2]])){
        temp<-c(temp,info[[2]][[k]])
      }
      temp<-c(temp,rep(-1,14-length(temp)))
      compo_sd[nrow(compo_sd)+1,]<-temp
    }
    
  }
  
  rlist[[1]] <- result
  rlist[[2]] <- compo_sd
  return(rlist)
}







#### Main



# down-sampling
result<-down_sampling2(dd_all,c(500,1000,2000,3000),5)

# write data.frame
write.table(result[[1]], 'cluster_validation_results_3.txt',row.names = FALSE, sep='\t',quote=F)
write.table(result[[2]], 'cluster_validation_results_3_sd.txt',row.names = FALSE, sep='\t',quote=F)
