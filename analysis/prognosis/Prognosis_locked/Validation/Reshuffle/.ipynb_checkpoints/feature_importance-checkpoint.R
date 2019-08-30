# ---------------------------------- #
# FEATURE IMPORTANCE
# ---------------------------------- #

# This R script gathers the functions used to obtain the feature importance for any model and any algorithm


# ---------------------------------- #
# Algorithms and prediction
# ---------------------------------- #
                                   
# ---------------------------------- #
# Penalized Cox
# ---------------------------------- #
library(glmnet)

algo_Lasso <- function(designTrain, responseTrain, alpha=1, ninternalfolds=10) {
    cvfit = cv.glmnet(designTrain, responseTrain, family="cox", alpha=alpha, nfolds=ninternalfolds, grouped=TRUE)
    return(cvfit)
}

predictor_Lasso <- function(fit,designTest){
    risk.predict = predict(fit, newx=designTest, s="lambda.1se", type="response")
    risk.predict = as.vector(risk.predict[,1])
    return(risk.predict)
}


algo_Ridge <- function(designTrain, responseTrain, alpha=0, ninternalfolds=10) {
    cvfit = cv.glmnet(designTrain, responseTrain, family="cox", alpha=alpha, nfolds=ninternalfolds, grouped=TRUE)
    return(cvfit)
}

predictor_Ridge <- function(fit,designTest){
    risk.predict = predict(fit, newx=designTest, s="lambda.1se", type="response")
    risk.predict = as.vector(risk.predict[,1])
    return(risk.predict)
}


algo_Elastic_net <- function(designTrain, responseTrain, alpha=0.7, ninternalfolds=10) {
    cvfit = cv.glmnet(designTrain, responseTrain, family="cox", alpha=alpha, nfolds=ninternalfolds, grouped=TRUE)
    return(cvfit)
}

predictor_Elastic_net <- function(fit,designTest){
    risk.predict = predict(fit, newx=designTest, s="lambda.1se", type="response")
    risk.predict = as.vector(risk.predict[,1])
    return(risk.predict)
}


# ---------------------------------- #
# Random effects modelling
# ---------------------------------- #
library(CoxHD)

algo_RFX <- function(designTrain, responseTrain, max.iter = 50) {
    cvfit = CoxRFX(data.frame(designTrain), Surv(time=responseTrain[,1],event =responseTrain[,2]), max.iter =max.iter)
    cvfit$Z <- NULL
    return(cvfit)
}

predictor_RFX <- function(fit,designTest){
    risk.predict<-predict(fit,data.frame(designTest))
    
    return(risk.predict)
}
   

# ---------------------------------- #
# Complementary pairs stability selection
# ---------------------------------- #

algo_CPSS <- function(designTrain, responseTrain, bootstrap.samples = 50, nlambda = 100, mc.cores = 1) {
    cvfit = CoxCPSSInteractions(data.frame(designTrain), Surv(time=responseTrain[,1],event =responseTrain[,2]), bootstrap.samples = bootstrap.samples, nlambda = nlambda, mc.cores = mc.cores)
    return(cvfit)
}

predictor_CPSS <- function(fit,designTest){
     risk.predict = predict(fit, data.frame(designTest))
    
    return(risk.predict)
}
                                          
# ---------------------------------- #
# Random survival forests 
# ---------------------------------- #
library(randomForestSRC)
                                     
algo_RFS <- function(designTrain, responseTrain, ntree=1050, importance="none",nodesize=20) {
    cvfit = rfsrc(Surv(time, status) ~ ., data=data.frame(designTrain,responseTrain), ntree=ntree, importance=importance,nodesize=nodesize)
    return(cvfit)
}

predictor_RFS <- function(fit,designTest,importance="none"){
    risk.predict = predict(fit, data.frame(designTest), importance=importance)$predicted
    
    return(risk.predict)
}

# ---------------------------------- #
# Boosting 
# ---------------------------------- #

#### Warning: not working, error with newtime and new status for prediction
library(CoxBoost)

algo_BOOST <- function(designTrain, responseTrain) {
    cvfit<-CoxBoost(time=responseTrain[,1],
                  status=responseTrain[,2],
                  x=designTrain)
    return(cvfit)
}

predictor_BOOST <- function(fit,designTest){
     risk.predict<-predict(fit,designTest,newtime=responseTest[,1],newstatus=responseTest[,2],type='lp')
    
    return(risk.predict)
}
       

# ---------------------------------- #
# Cox
# ---------------------------------- #            
library('survival')                         
algo_Cox <- function(designTrain, responseTrain) {
    c <- coxph(Surv(time, status) ~ ., data=data.frame(designTrain,responseTrain))
    return(c)
}

predictor_Cox <- function(fit,designTest){
     risk.predict = predict(fit, data.frame(designTest))
    
    return(risk.predict)
}






# ---------------------------------- #
# Concordance index difference with test set
# ---------------------------------- #

get_CI_with_test <- function(designTrain, designTest, responseTrain, responseTest, features = colnames(designTrain),seed = 3456, npermutations, algorithm, predictor){
    # features: features to test for importance
    # algorithm: algorithm used to train, returns a fit
    # predictor: predicts the risk corresponding to the algorithm used for training
    # npermutations: number of permutations
    
    
    res <- data.frame('feature'=character(),
                      'ref_CI'=numeric(),
                      'permuted_CI'=numeric())
    
    
    # Train and compute ref CI on testing set
    fit <- algorithm(designTrain,responseTrain)
    risk_test_ref <- predictor(fit,designTest)
    ci_test_ref <- as.numeric(survConcordance(Surv(time,status) ~ risk_test_ref, as.data.frame(responseTest))$concordance) 
    
    # Permutations
    set.seed(seed)
    ntest <- nrow(responseTest)
    permutations <- list()
    for (j in seq(npermutations)) {
        permutations[[j]] <-sample(1:ntest)
    }
    
    # Compute CI for shuffled features
    for(feature in features){
        for(a in 1:length(permutations)){
            tmp_designTest <- designTest
            tmp_designTest[,feature]<- tmp_designTest[permutations[[a]],feature]
            risk_test_tmp <- predictor(fit,tmp_designTest)
            ci_test_tmp <- as.numeric(survConcordance(Surv(time,status) ~ risk_test_tmp, as.data.frame(responseTest))$concordance) 
            tmp <- data.frame('feature'=feature,
                              'ref_CI'=ci_test_ref,
                              'permuted_CI'=ci_test_tmp)
            res <- rbind(res,tmp)
        }
    }
    
    return(res)
}


                                         
# ---------------------------------- #
# Cross Validation
# ---------------------------------- #
library(parallel)
runCV_CI_with_test <- function(response, design, nfolds=5, nrepeats=10, seed=5396, mc.cores=1, features=colnames(design), npermutations=10, algorithm, predictor) {
    
    #  random number generator seed
    set.seed(seed)

    res <- data.frame('feature'=character(),
                      'ref_CI'=numeric(),
                      'permuted_CI'=numeric())
    
    # Make folds
    n = nrow(design)
    folds <- list()
    for (i in seq(nrepeats)) {
        folds <- c(folds, split(sample(seq(n)), rep(1:nfolds, length = n)))
    }
    nexp = length(folds) # the total number CV of experiments
    
    

    # Parallel CV
    print("start CV")
    rescv = mclapply(seq(nexp),
                   FUN=function(iexp) {
                       cat(".")
                       vTrain = design[-folds[[iexp]],,drop=F]
                       vTest = design[folds[[iexp]],,drop=F]
                       lTrain = response[-folds[[iexp]],]
                       lTest = response[folds[[iexp]],]
                       
                       tmp<- get_CI_with_test(vTrain, vTest, lTrain, lTest, features, seed=2*seed, npermutations, algorithm, predictor)
                       
                       return(tmp)
                   },
                   mc.cores=mc.cores
                   )

    for(i in 1:length(rescv)){
        res <- rbind(res,rescv[[i]])
    }
    
    return(res)

}

           




