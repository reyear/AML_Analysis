library(glmnet)
library(doMC)
library(survival)
library(data.table)
library(mltools)
library(CoxBoost)
library(randomForestSRC)
library(CoxHD)

df_all <-read.table("df_prognosis.tsv",sep = '\t' , header = T) 

###    2094 rows 166 columns

df_all <-na.omit(df_all) # delete rows with na (161)
df_all <- df_all[df_all$os>0,] # delete when os is negative (2)

####

#Convert predicted_component to one hot encoder
df_all$new_eln<-factor(df_all$new_eln, levels = c("adverse","intermediate","favorable"), labels = 0:2, ordered = TRUE)  # convert categorical new_eln to numerical (0,1,2)
name <-rownames(df_all)
df_all$predicted_component <- as.factor(df_all$predicted_component)
df_final <- as.data.frame(one_hot(as.data.table(df_all),cols="predicted_component"))
rownames(df_final) <- name

### Predictors
predictorGLM <- function(designTrain, designTest, responseTrain, alpha, ninternalfolds=10) {
    # alpha=1 --> l1 penalty
    # alpha=0 --> l2 penalty
    # alpha=1/2 --> elastic net
    set.seed(1010)
    # Train
    cvfit = cv.glmnet(designTrain, responseTrain, family="cox", alpha=alpha, nfolds=ninternalfolds, grouped=TRUE)
    # Predict
    risk.predict = predict(cvfit, newx=designTest, s="lambda.1se", type="response")
    risk.predict = as.vector(risk.predict[,1])

    return(risk.predict)
}    

predictorBoost<-function(designTrain, designTest, responseTrain){
  set.seed(1010)
  cvfit<-CoxBoost(time=responseTrain[,1],
                  status=responseTrain[,2],
                  x=designTrain)
  
  risk.predict<-predict(cvfit,designTest,newtime=responseTest[,1],newstatus=responseTest[,2],type='lp')
  
  return(as.vector(risk.predict))
}
predictorRF <- function(designTrain, designTest, responseTrain, ntree=ntree, importance="none") {
    set.seed(1010)
    # Train
    cvfit = rfsrc(Surv(time, status) ~ ., data=data.frame(designTrain,responseTrain), ntree=ntree, importance=importance)
    
    # Predict
    risk.predict = predict(cvfit, data.frame(designTest), importance=importance)$predicted
    
    return(risk.predict)
} 
predictorAIC <- function(designTrain, designTest, responseTrain) {
    set.seed(1010)
    # Train
    c <- coxph(Surv(time, status) ~ ., data=data.frame(designTrain,responseTrain))
    scopeStep <- as.formula(paste("Surv(time,status) ~", paste(colnames(designTrain), collapse="+")))
    cvfit<-step(c, scope=scopeStep, k = 2, trace=0)
    # Predict
    risk.predict = predict(cvfit, data.frame(designTest))
    
    return(risk.predict)
}
predictorRFX <- function(designTrain, designTest, responseTrain, max.iter = 500) {
    set.seed(1010)
    # Train
    cvfit = CoxRFX(data.frame(designTrain), Surv(time=responseTrain[,1],event =responseTrain[,2]) , max.iter =max.iter)
    cvfit$Z <- NULL
    # Predict
    risk.predict<-predict(cvfit,data.frame(designTest))
    
    return(risk.predict)
}

### End Predictors

### Cross Validation

runCV <- function(mypredictor, response, design, nfolds=nfolds, nrepeats=nrepeats, seed=seed, mc.cores=mc.cores,alpha=alpha,use_alpha=FALSE,use_ntree=FALSE,ntree, ...) {
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
                       #predict.test = ifelse(use_alpha==TRUE,mypredictor(designTrain=vTrain, designTest=vTest, response=lTrain,alpha=alpha, ...),mypredictor(designTrain=vTrain, designTest=vTest, response=lTrain, ...))
                       if(use_alpha){
                           predict.test = mypredictor(designTrain=vTrain, designTest=vTest, response=lTrain,alpha=alpha, ...)
                       }else if(use_ntree) {
                           predict.test = mypredictor(designTrain=vTrain, designTest=vTest, response=lTrain,ntree=ntree, ...)
                       }else{
                           predict.test = mypredictor(designTrain=vTrain, designTest=vTest, response=lTrain, ...)
                       }
                       #predict.test = mypredictor(designTrain=vTrain, designTest=vTest, response=lTrain,alpha=alpha, ...)
                       # Evaluate CI on the test
                       ci.test = suppressWarnings(survConcordance(Surv(time,status) ~ predict.test, as.data.frame(lTest)))
                       return(as.vector(ci.test$concordance))
                   },
                   mc.cores=mc.cores
                   )

    return(unlist(rescv))

}

### End Cross Validation

x <- data.matrix(df_final[,1:84])
y <- data.matrix(df_final[,c("os","os_status")])
colnames(y) = c("time","status")

predictors <- c(rep(list(predictorGLM),10),rep(list(predictorRF),10),predictorBoost,predictorRFX)
str_predictors <-c(rep("CoxMod",10),rep("RFS",10),"CoxBoost","RFX")
res.CLIN_DEMO_CYTO <- c()
l_alpha <-seq(0.1,1,0.1)
l_ntree <- seq(100,1000,100)
i<-0
j<-0
k<-0
for(predictor in predictors){
    use_alpha<-ifelse(identical(predictorGLM,predictor),TRUE,FALSE)
    use_ntree<-ifelse(identical(predictorRF,predictor),TRUE,FALSE)
    i <- i+1
    j <-ifelse(use_alpha,j+1,j)
    k <-ifelse(use_ntree,k+1,k)
    alpha <- l_alpha[j]
    ntree <-l_ntree[k]
    tmp <- runCV(mypredictor=predictor,
          response=y, design=x,
          nfolds=5, nrepeats=10, seed=233,use_alpha=use_alpha,alpha=alpha,use_ntree=use_ntree,ntree=ntree, mc.cores=1)
    res.CLIN_DEMO_CYTO <- cbind(res.CLIN_DEMO_CYTO,tmp)
    colnames(res.CLIN_DEMO_CYTO) [i] <-paste(str_predictors[i],ifelse(use_alpha,alpha,
                                                                                   ifelse(use_ntree,ntree,"")),sep="_")
}

write.table(unlist(rescv),"test.tsv",quote=F,sep='\t')