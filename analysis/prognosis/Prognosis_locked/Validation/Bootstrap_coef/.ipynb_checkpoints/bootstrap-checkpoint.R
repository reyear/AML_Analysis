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
source("../../../../../src/tools.R")
source('../../../tools_prognosis/run_prognosis.R')




### Features that we can use
###-----------------------------------------------------------------------------

df_final <- read.table("../full_data_validation.tsv")
eln <- c(113,114,115)
comp <- c(89:112)
age <- c(84)

all_gen <- c(1:57)
vect <- apply(X=df_final[,all_gen],2,FUN=function(x) 100*length(which(x==1))/dim(df_final)[1])
gen <- match(names(vect[vect>=2]),names(df_final))
gen_without <- setdiff(gen,grep("^NPM1$", colnames(df_final)))
gen_without <- setdiff(gen_without,grep("^CEBPA_bi$", colnames(df_final))) 
              
all_cyto <- c(58:80)
vect <- apply(X=df_final[,all_cyto],2,FUN=function(x) 100*length(which(x==1))/dim(df_final)[1])
cyto <- match(names(vect[vect>=2]),names(df_final))
cyto_without <- setdiff(cyto,grep("^inv_16$", colnames(df_final)))
cyto_without <- setdiff(cyto_without,grep("^t_8_21$", colnames(df_final)))      
cyto_without <- setdiff(cyto_without,grep("^t_v_11$", colnames(df_final)))  
              
clin <- c(85:87)
demo <- c(83,84)
demo_without_age <-c(83)

                          
                          
eln_comp <- c(eln,comp)          
eln_gen <- c(eln,gen)
eln_cyto <- c(eln,cyto_without)
eln_clin <- c(eln,clin)
eln_demo <- c(eln,demo)

# USEFUL FOR ELN COMPARISON
# with comp
eln_comp_gen <- c(eln_comp,gen_without)
eln_comp_cyto <- c(eln_comp,cyto_without)
eln_comp_clin <- c(eln_comp,clin)
eln_comp_demo <- c(eln_comp,demo)


eln_comp_gen_cyto <- c(eln_comp_gen,cyto_without)
eln_comp_gen_clin <- c(eln_comp_gen,clin)
eln_comp_gen_demo <- c(eln_comp_gen,demo)

eln_comp_cyto_clin <- c(eln_comp_cyto,clin)
eln_comp_cyto_demo <- c(eln_comp_cyto,demo)

eln_comp_clin_demo <- c(eln_comp_clin,demo)

eln_comp_gen_cyto_clin_demo <- c(eln_comp_gen_cyto,clin,demo)
eln_comp_gen_cyto_clin_demo_without_age <- c(eln_comp_gen_cyto,clin,demo_without_age)
              


# without comp


eln_gen_cyto <- c(eln_gen,cyto)
eln_clin_demo <- c(eln_clin,demo)
eln_clin <- c(eln,clin)
eln_demo <- c(eln,demo)

eln_gen_cyto_clin_demo <- c(eln_gen_cyto,clin,demo)

# USEFUL FOR COMP
comp_gen <- c(comp,gen_without)
comp_cyto <- c(comp,cyto_without)
comp_clin <- c(comp,clin)
comp_demo <- c(comp,demo)
comp_gen_cyto <- c(comp_gen,cyto_without)
comp_clin_demo <- c(comp_clin,demo)
comp_gen_cyto_clin_demo <- c(comp_gen_cyto,clin,demo)

#USEFUL FOR GEN
gen_cyto <- c(gen,cyto)
gen_clin <- c(gen,clin)
gen_demo <- c(gen,demo)
gen_clin_demo <- c(gen_clin,demo)
gen_cyto_clin_demo <- c(gen_cyto,clin,demo)

#USEFUL FOR CYTO 
cyto_clin <- c(cyto,clin)
cyto_demo <- c(cyto,demo)
cyto_clin_demo <- c(cyto_clin,demo)
cyto_gen_demo <- c(gen_cyto,demo)

clin_demo <-c(clin,demo)



y <- data.matrix(df_final[,c("OS","OS_Status")])

prognosis_features<- list(demo=demo,clin=clin,gen=gen,cyto=cyto,comp=comp,eln=eln,gen=gen,cyto=cyto,gen_cyto=gen_cyto,eln_gen_cyto=eln_gen_cyto,comp_gen_cyto=comp_gen_cyto,eln_comp=eln_comp,eln_clin_demo=eln_clin_demo,comp_clin_demo=comp_clin_demo,eln_comp_gen_cyto_clin_demo=eln_comp_gen_cyto_clin_demo)

              
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
y <- data.matrix(df_final[,c("OS","OS_Status")])

colnames(y) = c("time","status")
response=y
              

algos <-c("glm","rfs","boost","rfx")
alphas=c(0,0.7,1)
for (i in 1:length(prognosis_features)){
    for (algo in algos){
        if (algo=="glm"){
            for (alpha in alphas){
                print(alpha)
                print(algo)
                bootstrap <- bootstrapping(prognosis_features[[i]],x,y,100,alpha,30,algo)
                tmp_1 <- bootstrap  %>% group_by(feature) %>% summarise_all(sum)
                tmp_2 <- bootstrap  %>% group_by(feature) %>% count(feature)
                print(paste(paste(names(prognosis_features)[i],paste(algo,alpha,sep="_"),sep="_bootstrap_"),".tsv",sep=""))
                write.table(data.frame(merge(tmp_1,tmp_2,by='feature')),paste(paste(names(prognosis_features)[i],paste(algo,alpha,sep="_"),sep="_bootstrap_"),".tsv",sep=""),quote=F,sep='\t')

                if (alpha==0.7){
                    tmp_1_pos <- tmp_1[tmp_1$coef>0,]
                    tmp_1_neg <-  tmp_1[tmp_1$coef<0,]
                    features_reduced <- union(union(tmp_1_pos[tmp_1_pos$coef > quantile(tmp_1_pos$coef,0.90),]$feature,tmp_1_neg[tmp_1_neg$coef < quantile(tmp_1_neg$coef,0.15),]$feature),tmp_2[tmp_2$n > quantile(tmp_2$n,0.85),]$feature)
                    if (length(features_reduced)<2){features_reduced <- union(union(tmp_1_pos[tmp_1_pos$coef > quantile(tmp_1_pos$coef,0.90),]$feature,tmp_1_neg[tmp_1_neg$coef < quantile(tmp_1_neg$coef,0.15),]$feature)
                                                                              ,tmp_2[tmp_2$n > 0,]$feature)}
                    print(features_reduced)

                    predictors <- c(rep(list(predictorGLM),11),rep(list(predictorRF),1),predictorBoost,predictorRFX)
                    str_predictors <-c(rep("CoxGLM",11),"RFS","CoxBoost","RFX")
                    l_alpha <-seq(0,1,0.1)
                    l_ntree <- c(1050)
                    mc.cores <- 50
                    nodesize <- c(20)
                    print("DONE")
                    write.table(launch_prognosis(data.matrix(df_final[,features_reduced]),y=y,predictors=predictors,str_predictors=str_predictors,l_alpha=l_alpha,nrepeats=2,l_ntree=l_ntree,nodesize=nodesize,
                                mc.cores=mc.cores),paste(names(prognosis_features)[i],"_reduced.tsv",sep=""),quote=F,sep='\t')
                    print("DONE")
                }
            }
        } else {
                print(algo)
                if(algo=="rfs"){
                    bootstrap <- bootstrapping(prognosis_features[[i]],x,y,10,0.7,8,algo)
                }else {
                    bootstrap <- bootstrapping(prognosis_features[[i]],x,y,100,0.7,8,algo)
                    tmp_1 <- bootstrap  %>% group_by(feature) %>% summarise_all(sum)
                    tmp_2 <- bootstrap  %>% group_by(feature) %>% count(feature)
                    }

                write.table(data.frame(merge(tmp_1,tmp_2,by='feature')),paste(paste(names(prognosis_features)[i],algo,sep="_bootstrap_"),".tsv",sep=""),quote=F,sep='\t')
    print ('next')
        }
    }
}