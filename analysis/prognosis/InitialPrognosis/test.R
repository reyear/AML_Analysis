library(glmnet)
library(doMC)
library(survival)
library(data.table)
library(mltools)
library(CoxBoost)
library(randomForestSRC)

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

predictorRidge <- function(designTrain, designTest, responseTrain, alpha=0.9, ninternalfolds=10) {
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

predictorBOOST<-function(designTrain, designTest, responseTrain){
  cvfit<-CoxBoost(time=responseTrain[,1],
                  status=responseTrain[,2],
                  x=designTrain)
  
  risk.predict<-predict(cvfit,designTest,newtime=responseTest[,1],newstatus=responseTest[,2],type='lp')
  
  return(as.vector(risk.predict))
}
predictorRF <- function(designTrain, designTest, responseTrain, ntree=5, importance="none") {
    
    # Train
    cvfit = rfsrc(Surv(time, status) ~ ., data=data.frame(designTrain,responseTrain), ntree=ntree, importance=importance)
    
    # Predict
    risk.predict = predict(cvfit, data.frame(designTest), importance=importance)$predicted
    
    return(risk.predict)
} 
predictorAIC <- function(designTrain, designTest, responseTrain) {
    
    # Train
    c <- coxph(Surv(time, status) ~ ., data=data.frame(designTrain,responseTrain))
    scopeStep <- as.formula(paste("Surv(time,status) ~", paste(colnames(designTrain), collapse="+")))
    cvfit<-step(c, scope=scopeStep, k = 2, trace=0)
    # Predict
    risk.predict = predict(cvfit, data.frame(designTest))
    
    return(risk.predict)
}

set.seed(17)
x <- data.matrix(df_final[,1:177])
y <- data.matrix(df_final[,c("os","os_status")])
colnames(y) = c("time","status")
#y <- Surv(time = df_final$os, event = df_final$os_status)
nrepeats=1
nfolds=5 # to do 80% vs 20%
    # Make folds
n = nrow(x)
folds <- list()
# This splits the dataset into nfolds folds without repetition
for (i in seq(nrepeats)) {
    folds <- c(folds,split(sample(seq(n)), rep(1:nfolds, length = n)))
}

nexp=length(folds)
print("start CV")
rescv = mclapply(seq(nexp),
               FUN=function(iexp) {
                   cat(".")
                   vTrain = x[-folds[[iexp]],,drop=F]
                   vTest = x[folds[[iexp]],,drop=F]
                   lTrain = y[-folds[[iexp]],]
                   lTest = y[folds[[iexp]],]
                   # Train and Predcit
                   predict.test = predictorRidge(designTrain=vTrain, designTest=vTest, responseTrain=lTrain)
                   # Evaluate CI on the test
                   ci.test = survConcordance(Surv(time,status) ~ predict.test, as.data.frame(lTest))
                   print(as.vector(ci.test$concordance))
               },
               mc.cores=50
               )
write.table(unlist(rescv),"test.tsv",quote=F,sep='\t')