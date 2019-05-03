
library('hdp')
library('clusterCrit')
library('grid')
library('gridExtra')
library('ggplot2')
library('ggrepel')

source('../../../src/tools.R')     # custom tools function
source('../../../src/hdp_tools.R') # hdp related functions
source('../../../src/hdp_tools_yanis.R')
theme_set(theme_minimal())

# set jupyer notebook parameters
options(repr.plot.res        = 100, # set a medium-definition resolution for the jupyter notebooks plots (DPI)
        repr.matrix.max.rows = 200, # set the maximum number of rows displayed
        repr.matrix.max.cols = 200) # set the maximum number of columns displayed

df_final <- read.table("../../../data/updated_dataset/final_clustering.csv",sep = ',' , header = T)
rownames(df_final)<- df_final$X
df_final <- df_final[,-1:-3]

###Binomial
num_cols = ncol(df_final)
bin <- function(x){
  set.seed(123)
  (rbinom(1, num_cols, mean(x))+1)/num_cols
}

###Normal

normal <- function(x){
  set.seed(123)
  abs(rnorm(1,mean(x),sd(x)))
}

###Poisson

poisson <- function(x){
  set.seed(123)
  (rpois(num_cols,1)+1)/num_cols
}

###Uniform equally over all columns

equally <- function(x){
  set.seed(123)
  1/num_cols
}

###Repet 1

repet <- function(x){
  set.seed(123)
  1
}

#sapply(df_gene,equally)
binomial <- unlist(sapply(df_final,bin))
gaussian <- unlist(sapply(df_final,normal))
pois <- as.numeric(unlist(sapply(df_final,poisson)))
unif <- unlist(sapply(df_final,equally))
repetition <- unlist(sapply(df_final,repet))

# define gridsearch function
test_assignment <- function(x, dd_predicted, components_sd){
  if(x['second_predicted_component']%in%dd_predicted$predicted_component) {
    return( x['max_proba'] > x['second_max_proba']+components_sd[[x['second_predicted_component']+1]])
  }
  else {return(TRUE)} 
}
hdp_gridsearch_yanis <- function( data, posterior_samples, initial_clusters, burnin, chains, space, base_dist, alphaa, alphab){  #+3 components
  
  result_table <- data.frame('chains' = integer(),
                             'inicc' = integer(),
                             'n' = integer(),
                             'burnin' = integer(),
                             'space' = integer(),
                             'alphaa' = integer(),
                             'alphab' = integer(),
                             'n_components' = integer(),
                             'component_0' = integer(),
                             'assignment' = numeric(),
                             'Davies_Bouldin' =numeric(),
                             'Silhouette' = numeric())
  
  
  for (a in chains) {
    for (b in initial_clusters) {
      for (c in posterior_samples) {
        for (d in burnin) {
          for (e in space) {
            for (f in base_dist){
              for (g in alphaa){
                for (h in alphab){
                  
                  cat(sprintf('### %d chains | %d initial clusters | %d posterior samples | %d burn-in | %d space ###\n', a, b, c, d, e))
                  flush.console()
                  
                  number_of_chains <- a
                  chain_list <- vector('list', number_of_chains)
                  hdp <- initialise_hdp_yanis(data=data,alphaa=g,alphab=h,hh=f)
                  for (i in 1:number_of_chains) {
                    seed <- i * 100
                    
                    activated_hdp <- dp_activate(hdp,
                                                 dpindex = 1:numdp(hdp), # indices of the DP to activate: we choose all DP (1, ..., n + 1)
                                                 initcc  = b,           # number of starting clusters (every data item is randomly assigned to a cluster to start with)
                                                 seed    = seed)
                    
                    hdp_output <- hdp_posterior(activated_hdp,
                                                burnin = d, # number of burn-in iterations
                                                n      = c,      # number of posterior samples to collect
                                                space  = e,  # number of iterations between collected samples
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
                  dd_predicted['predicted_component'] <- apply(dd_predicted, 1, function(x) { if (all(is.na(x)))
                    return(NaN)
                    else
                      return(which.max(x)-1)
                  })
                  
                  
                  
                  
                  # evaluate for each row the maximum probability associated to the predicted component
                  dd_predicted['max_proba'] <- apply(dd_predicted[,components_colnames], 1, function(x) { if (all(is.na(x)))
                    return(NaN)
                    else
                      return(max(x))})
                  
                  
                  # get second predicted component
                  n_components <- ncol(dd_predicted) - 2 - 1 # without columns 'component_0', 'predicted_component' and 'max_proba'
                  
                  dd_predicted['second_max_proba'] <- apply(dd_predicted[,0:n_components + 1], 1, function(x) { if (all(is.na(x)))
                    return(NaN)
                    else
                      return(sort(x,partial=n_components)[n_components])})
                  
                  dd_predicted['second_predicted_component'] <- apply(dd_predicted, 1, function(x) { if (all(is.na(x)))
                    return(NaN)
                    else
                      return(which(x==x['second_max_proba'])[1]-1)})
                  
                  
                  
                  # get standard deviations of components
                  components_sd <- list()
                  for( j in as.integer(levels(factor(dd_predicted$predicted_component)))){
                    if(! is.na(j)){
                      components_sd[[j+1]] <- sd(dd_predicted$max_proba[!is.na(dd_predicted$predicted_component) & dd_predicted$predicted_component == j])
                    }
                  }
                  
                  
                  # determine which patients are well assigned
                  dd_predicted['assignment'] <- apply(dd_predicted[, c('max_proba', 'second_max_proba','second_predicted_component')], 1, function(x) { if (all(is.na(x)))
                    return(NaN)
                    else
                      return(test_assignment(x, dd_predicted, components_sd))})
                  
                  # validation indexes
                  traj <- data.frame(data[rowSums(data) !=0,])
                  traj <- as.matrix(apply(traj, 2, function (x){as.numeric(x)}))
                  part <- as.integer(dd_predicted$predicted_component[!is.na(dd_predicted$predicted_component)])
                  intIdx <- intCriteria(traj, part, c("Davies_bouldin","Silhouette"))
                  
                  # result
                  result_table[nrow(result_table)+1,] <- c(a, b, c, d, e, g, h, numcomp(multi_output), nrow(dd_predicted[!is.na(dd_predicted$predicted_component) &dd_predicted$predicted_component==0,]), nrow(dd_predicted[!is.na(dd_predicted$assignment) & dd_predicted$assignment==TRUE,]),intIdx$davies_bouldin,intIdx$silhouette)
                  
                  cat('Done!\n')
                  flush.console()
                }
              }
            }
          }
        }
      }
    }
  }
  
  return(result_table)
  
}

# initialize the range of parameters used for the grid search
data <- df_final
posterior_samples <- 350
initial_clusters <- 13
burnin <- 5000
chains <- 3
space <- 20
base_dist <- list(binomial,gaussian)
alphaa <-2
alphab <-0.5
res <- hdp_gridsearch_yanis(data, posterior_samples, initial_clusters, burnin, chains, space,base_dist,alphaa,alphab)

write.table(res,'hdp_gridsearch.tsv',sep='\t',quote=F)