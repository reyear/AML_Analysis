# ---------------------------------- #
# SURVIVAL PREDICTION COMPONENT 
# ---------------------------------- #

# ---------------------------------- #
# Design and Response matrix
# ---------------------------------- #
library(survival)
PrepareDesignMatrix <- function(data,
                                use.nuisance=TRUE,
                                response.name=c("os_diag_years","os_status"),
                                feature.category.name = c("ipssr","genenumber"),
                                list.features
                                ) {

    if (use.nuisance) {
        feature.category.name = unique(c("nuisance",feature.category.name))
    }
    myfeatures = as.vector(unlist(sapply(feature.category.name, function(x) list.features[x])))
    myformula = as.formula(paste("~ 0 +",paste(myfeatures,collapse="+")))

    design = model.matrix(myformula, data)
    rownames(design) = NULL

    response = as.matrix(data[,response.name])
    colnames(response) = c("time","status")

    return(list(design=design,response=response))

}

                                         
# ---------------------------------- #
# DEFINE FUNCTION FOR GROUPING
# ---------------------------------- #

get_category <- function(x,clin,wbc,mutations,geneint_double,geneint_triple,mutnumber,cytoaberrations,cytonumber){
    if (x %in% geneint_double | x %in% geneint_triple){return("GENE_INTERACTIONS")}
    if (x %in% c(mutations,mutnumber)){return("GENE")}
    if (x %in% c('AGE','SEXM')){return("DEMOGRAPHIC")}
    if (x %in% c(clin,wbc)){return('CLINICAL')}
    #if (x %in% c(mutnumber,cytonumber)){return("NUMBER")}
    if (x %in% c(cytoaberrations,cytonumber)){return('CYTO')}
    if (x %in% c('IPSSR_CALCULATEDINT','IPSSR_CALCULATEDGOOD','IPSSR_CALCULATEDPOOR','IPSSR_CALCULATEDVERY-GOOD','IPSSR_CALCULATEDVERY-POOR','CYTO_IPSSRINT','CYTO_IPSSRPOOR','CYTO_IPSSRVERY-POOR','CYTO_IPSSRGOOD','CYTO_IPSSVERY-GOOD')){return("SCORES")}
    if (x %in% c('componentComponent_5','componentComponent_0','componentComponent_4','componentComponent_NaN','componentComponent_1','componentComponent_2','componentComponent_3','componentComponent_9','componentComponent_7','componentComponent_6','componentComponent_8','componentComponent_10')){return("COMPONENT")}
    if (x %in% c('CENTER_MSK','CENTER_CCH','CENTER_CGM','CENTER_DUS','CENTER_DUTH','CENTER_FLO','CENTER_FUCE','CENTER_GESMD','CENTER_HIAE','CENTER_HMS','CENTER_ICO','CENTER_IHBT','CENTER_KI','CENTER_MUV','CENTER_PV','CENTER_REL','CENTER_RMCN','CENTER_ROM','CENTER_TUD','CENTER_UMG','CENTER_UOB','CENTER_UOXF','CENTER_VU')){return("CENTER")}
        
}                                         
                                         
# ---------------------------------- #
# Penalized Cox
# ---------------------------------- #
library(glmnet)
predictorLasso <- function(designTrain, designTest, responseTrain, alpha=1, ninternalfolds=10) {
    # alpha=1 --> l1 penalty
    # alpha=0 --> l2 penalty
    # alpha=1/2 --> elastic net

    # Train
    cvfit = cv.glmnet(designTrain, responseTrain, family="cox", alpha=alpha, nfolds=ninternalfolds, grouped=TRUE)
    # Predict
    risk.predict = predict(cvfit, newx=designTest, s="lambda.1se", type="response")
    risk.predict = as.vector(risk.predict[,1])

    return(risk.predict)
}

predictorRidge <- function(designTrain, designTest, responseTrain, alpha=0, ninternalfolds=10) {
    # alpha=1 --> l1 penalty
    # alpha=0 --> l2 penalty
    # alpha=1/2 --> elastic net

    # Train
    cvfit = cv.glmnet(designTrain, responseTrain, family="cox", alpha=alpha, nfolds=ninternalfolds, grouped=TRUE)
    # Predict
    risk.predict = predict(cvfit, newx=designTest, s="lambda.1se", type="response")
    risk.predict = as.vector(risk.predict[,1])

    return(risk.predict)
}    

predictorElastic_net <- function(designTrain, designTest, responseTrain, alpha=0.5, ninternalfolds=10) {
    # alpha=1 --> l1 penalty
    # alpha=0 --> l2 penalty
    # alpha=1/2 --> elastic net

    # Train
    cvfit = cv.glmnet(designTrain, responseTrain, family="cox", alpha=alpha, nfolds=ninternalfolds, grouped=TRUE)
    # Predict
    risk.predict = predict(cvfit, newx=designTest, s="lambda.1se", type="response")
    risk.predict = as.vector(risk.predict[,1])

    return(risk.predict)
}                                         
                                         
# ---------------------------------- #
# Random effects modelling
# ---------------------------------- #
library(CoxHD)
predictorRFX <- function(designTrain, designTest, responseTrain, max.iter = 50) {
    
    # Train
    cvfit = CoxRFX(data.frame(designTrain), Surv(time=responseTrain[,1],event =responseTrain[,2]) , max.iter =max.iter)
    cvfit$Z <- NULL
    # Predict
    risk.predict<-predict(cvfit,data.frame(designTest))
    
    return(risk.predict)
}
                                         
predictorRFX_groups <- function(designTrain, designTest, responseTrain, max.iter = 50,clin,wbc,mutations,geneint_double,geneint_triple,mutnumber,cytoaberrations,cytonumber) {
    
    # Remove one center 
    designTrain<-designTrain[,2:ncol(designTrain)]
    designTest<-designTest[,2:ncol(designTest)]
  
    # Groups
    groups<-sapply(colnames(designTrain),function(x){return(get_category(x,clin,wbc,mutations,geneint_double,geneint_triple,mutnumber,cytoaberrations,cytonumber))})
    groups<-factor(groups,levels=unique(groups))
    
    # Train
    cvfit = CoxRFX(data.frame(designTrain), Surv(time=responseTrain[,1],event =responseTrain[,2]) , max.iter =max.iter,groups=groups)
    cvfit$Z <- NULL
    # Predict
    risk.predict<-predict(cvfit,data.frame(designTest))
    
    return(risk.predict)
}                                         
   
# ---------------------------------- #
# Stepwise variable selection
# ---------------------------------- #            
library(survival)                                         
predictorAIC <- function(designTrain, designTest, responseTrain) {
    
    # Train
    c <- coxph(Surv(time, status) ~ ., data=data.frame(designTrain,responseTrain))
    scopeStep <- as.formula(paste("Surv(time,status) ~", paste(colnames(designTrain), collapse="+")))
    cvfit<-step(c, scope=scopeStep, k = 2, trace=0)
    # Predict
    risk.predict = predict(cvfit, data.frame(designTest))
    
    return(risk.predict)
}

# ---------------------------------- #
# Complementary pairs stability selection
# ---------------------------------- #
predictorCPSS <- function(designTrain, designTest, responseTrain, bootstrap.samples = 50, nlambda = 100, mc.cores = 1) {
    
    # Train
    cvfit = CoxCPSSInteractions(data.frame(designTrain), Surv(time=responseTrain[,1],event =responseTrain[,2]), bootstrap.samples = bootstrap.samples, nlambda = nlambda, mc.cores = mc.cores)
    # Predict
    risk.predict = predict(cvfit, data.frame(designTest))
    
    return(risk.predict)
}
                                          
# ---------------------------------- #
# Random survival forests 
# ---------------------------------- #
library(randomForestSRC)
predictorRF <- function(designTrain, designTest, responseTrain, ntree=500, importance="none") {
    
    # Train
    cvfit = rfsrc(Surv(time, status) ~ ., data=data.frame(designTrain,responseTrain), ntree=ntree, importance=importance)
    
    # Predict
    risk.predict = predict(cvfit, data.frame(designTest), importance=importance)$predicted
    
    return(risk.predict)
}                                     
                                         
predictorRFopti <- function(designTrain, designTest, responseTrain, ntree_range) {
  
  # Internal cross validation
  res <- c()
  for(ntree in ntree_range){
    res <- c(res,median(runCV(predictorRF, responseTrain, designTrain, nfolds=3, nrepeats=2, seed=1234, mc.cores=1, ntree=ntree)))
  }
  res <- data.frame(ntree=ntree_range,median=res)
  ntree_final <- res$ntree[which.max(res$median)]
  print(ntree_final)
  
  # Train
  cvfit = rfsrc(Surv(time, status) ~ ., data=data.frame(designTrain,responseTrain), ntree=ntree_final, importance='none')
  
  # Predict
  risk.predict = predict(cvfit, data.frame(designTest), importance='none')$predicted
  
  return(risk.predict)
}                                           

# ---------------------------------- #
# Boosting 
# ---------------------------------- #
library(CoxBoost)
predictorBOOST<-function(designTrain, designTest, responseTrain){
  cvfit<-CoxBoost(time=responseTrain[,1],
                  status=responseTrain[,2],
                  x=designTrain)
  
  risk.predict<-predict(cvfit,designTest,newtime=responseTest[,1],newstatus=responseTest[,2],type='lp')
  
  return(as.vector(risk.predict))
}

                                         
# ---------------------------------- #
# Cross Validation
# ---------------------------------- #
library(parallel)
runCV <- function(mypredictor, response, design, nfolds=5, nrepeats=10, seed=5396, mc.cores=1, ...) {
    # function that run "mypredictor" on a CV setting
    #
    # output a list of size the number of CV experiments (eg 50) (= nfolds x nrepeats)
    
    # "ref" contains the responses of the fold test set

    #  random number generator seed
    set.seed(seed)

    # Make folds
    n = nrow(design)
    folds <- list()
    for (i in seq(nrepeats)) {
        folds <- c(folds,split(sample(seq(n)), rep(1:nfolds, length = n)))
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
                       # Train and Predcit
                       predict.test = mypredictor(designTrain=vTrain, designTest=vTest, response=lTrain, ...)
                       # Evaluate CI on the test
                       ci.test = survConcordance(Surv(time,status) ~ predict.test, as.data.frame(lTest))
                       return(as.vector(ci.test$concordance))
                   },
                   mc.cores=mc.cores
                   )

    return(unlist(rescv))

}
