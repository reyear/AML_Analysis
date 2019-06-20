library('ggplot2')
library('reshape2')
library('ggpubr')
library(glmnet)
library(doMC)
library(survival)
library(data.table)
library(mltools)
library(CoxBoost)
library(randomForestSRC)
library(CoxHD)
library(Hmisc)
library(gridExtra)
library("survminer")
library(dplyr)
library(broom)
library(tidyr)
library(tidyverse)
source("../../../src/tools.R")
source('run_prognosis.R')

df_final <- read.table("prognosis_comp_final.tsv",sep='\t',header=T)
### Features that we can use
###-----------------------------------------------------------------------------
all_features <-c(1:180) #not used
clin_demo_comp <-c(155:180) #not used
clin_demo_cyto_gen_comp <- c(2:180) #not used
comp <- c(164:180) #not used
cyto_comp <-c(86:154,164:180) #not used
cyto_gen_comp <- c(2:154,164:180) #not used
eln_clin_demo_comp <- c(1,155:180) #not used
eln_cyto_comp <- c(1,86:154,164:180) #not used
eln_cyto_gen_comp <- c(1:154,164:180) #not used
eln_gen_comp <- c(1:85,164:180) #not used
gen_comp <- c(2:85,164:180) #not used
clin_comp <- c(155:161,164:180) #not used
clin_cyto_comp <- c(86:161,164:180) #not used
clin_gen_comp <- c(2:85,155:161,164:180) #not used
eln_clin_comp <- c(1,155:161,164:180) #not used
age <- c(163)
gen_age <- c(2:85,163)
eln_clin_gen <-  c(1:85,155:161)
eln_demo_gen <- c(1:85,162:163)
eln_clin_demo_cyto_gen <- c(1:163)
eln_clin_demo_cyto <- c(1,86:163)
eln_clin_demo_gen <- c(1:85,155:163)  ##START HERE
eln_clin_demo <- c(1,155:163)
eln_clin <- c(1,155:161)
eln_cyto_gen <- c(1:154)
clin_demo_cyto_gen <- c(2:163)
clin_demo_cyto <- c(86:163)
clin_demo_gen <- c(2:85,155:163)
clin_demo <- c(155:163)
cyto_gen <- c(2:154)
cyto <- c(86:154)
gen <- c(2:85)
clin_gen <- c(2:85,155:161)
clin_cyto <- c(86:161)
demo_gen <- c(2:85,162:163)
demo_cyto <- c(86:154,162:163)  

###Without age:
all_features_without_age <-c(1:162,164:180) #not used
clin_demo_comp_without_age <-c(155:162,164:180) #not used
clin_demo_cyto_gen_comp_without_age <- c(2:162,164:180) #not used
eln_clin_demo_comp_without_age <- c(1,155:162,164:180) #not used
eln_demo_gen_without_age <- c(1:85,162)
eln_clin_demo_cyto_gen_without_age <- c(1:162)
eln_clin_demo_cyto_without_age <- c(1,86:162)
eln_clin_demo_gen_without_age <- c(1:85,155:162)
eln_clin_demo_without_age <- c(1,155:162)
clin_demo_cyto_gen_without_age <- c(2:162)
clin_demo_cyto_without_age <- c(86:162)
clin_demo_gen_without_age <- c(2:85,155:162)
clin_demo_without_age <- c(155:162)
demo_gen_without_age <- c(2:85,162)
demo_cyto_without_age <- c(86:154,162) 




bootstrapping <- function(features=all_features,x,y,n_exp=100,alpha=0.7,mc.cores=50,model="glm"){
    set.seed(17)
    res_bootstrap <- data.frame('feature' = character(),
                      'coef' = numeric())
    design=x[,features]
    n = nrow(design)
    folds <- list()

    for (i in seq(n_exp)) {
        folds[[i]] <- sample(1:n, 0.8 * n, replace = TRUE)
    }

    nexp = length(folds)
    print("Start Bootstrapping")
    rescv = mclapply(seq(nexp),
                   FUN=function(iexp) {
                       set.seed(17)
                       cat(".")
                       x_sampling = design[folds[[iexp]],]
                       y_sampling = y[folds[[iexp]],]
                       if (model=="glm"){
                           cvfit <- cv.glmnet(x_sampling, y_sampling, family = 'cox', alpha=alpha, nfolds = 20, grouped = TRUE)
                           tmp <- as.data.frame(as.matrix(coef(cvfit, s = "lambda.min")))
                       } else if (model=="boost"){
                           cvfit<-CoxBoost(time=y_sampling[,1],status=y_sampling[,2],x=x_sampling)
                           tmp <- as.data.frame(as.matrix(coefficients(cvfit)))
                       } else if (model=="rfx"){
                           cvfit<-CoxRFX(data.frame(x_sampling),Surv(time=y_sampling[,1],event=y_sampling[,2]) , max.iter =50,tol=1e-3)
                           tmp <- as.data.frame(as.matrix(coef(cvfit)))
                       } else if (model=="rfs"){
                           cvfit <- rfsrc(Surv(time, status) ~ ., data=data.frame(x_sampling,y_sampling), ntree=1050, importance="TRUE",nodesize=20)
                           tmp <- as.data.frame(as.matrix(cvfit$importance))
                       }
                       colnames(tmp) <- 'coef'
                       tmp <- rownames_to_column(tmp, var = 'feature')


                   },
                   mc.cores=50
                   )

    for(i in 1:length(rescv)){
        res_bootstrap <- rbind(res_bootstrap,rescv[[i]])
    }


    res_bootstrap <- res_bootstrap[res_bootstrap$coef != 0,]
    return (res_bootstrap)
    }




x <- data.matrix(df_final)
y <- data.matrix(df_final[,c("os","os_status")])

colnames(y) = c("time","status")
response=y

prognosis_features<- list(all_features=all_features,clin_demo_comp=clin_demo_comp,clin_demo_cyto_gen_comp=clin_demo_cyto_gen_comp,comp=comp,cyto_comp=cyto_comp,cyto_gen_comp=cyto_gen_comp,eln_clin_demo_comp=eln_clin_demo_comp,
                         eln_cyto_comp=eln_cyto_comp,eln_cyto_gen_comp=eln_cyto_gen_comp,eln_gen_comp=eln_gen_comp,gen_comp=gen_comp,clin_comp=clin_comp,clin_cyto_comp=clin_cyto_comp,clin_gen_comp=clin_gen_comp,eln_clin_comp=eln_clin_comp,age=age)


predictors <- c(rep(list(predictorGLM),3),rep(list(predictorRF),1),predictorBoost,predictorRFX)
str_predictors <-c(rep("CoxGLM",3),rep("RFS",1),"CoxBoost","RFX")
l_alpha <-c(0,0.7,1)
l_ntree <- c(1050)
mc.cores <- 8
nodesize <- c(20)

