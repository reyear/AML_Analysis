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
source("../../../../src/tools.R")
source('../../tools_prognosis/run_prognosis.R')

df_final <- read.table("../data_frame_final_prognosis.tsv")


### Features that we can use
###-----------------------------------------------------------------------------
eln <- c(1)
comp <- c(164:178)
age <- c(163)
all_gen <- c(2:85)
vect <- apply(X=df_final[,all_gen],2,FUN=function(x) 100*length(which(x==1))/dim(df_final)[1])
gen <- match(names(vect[vect>=2]),names(df_final))
cyto <- c(86:154)
clin <- c(155:161)
demo <- c(162:163)
demo_without_age <-c(162)
           
                          
                          
eln_comp <- c(eln,comp)
eln_age <- c(eln,age)
eln_gen <- c(eln,gen)
eln_cyto <- c(eln,cyto)
eln_clin <- c(eln,clin)
eln_demo <- c(eln,demo)
eln_demo_without_age <- c(eln,demo_without_age)

# USEFUL FOR ELN COMPARISON
# with comp
eln_comp_age <- c(eln_comp,age)
eln_comp_gen <- c(eln_comp,gen)
eln_comp_cyto <- c(eln_comp,cyto)
eln_comp_clin <- c(eln_comp,clin)
eln_comp_demo <- c(eln_comp,demo)
eln_comp_demo_without_age <- c(eln_comp,demo_without_age)

eln_comp_age_gen <- c(eln_comp_age,gen)
eln_comp_age_cyto <- c(eln_comp_age,cyto)
eln_comp_age_clin <- c(eln_comp_age,clin)

eln_comp_gen_cyto <- c(eln_comp_gen,cyto)
eln_comp_gen_clin <- c(eln_comp_gen,clin)
eln_comp_gen_demo <- c(eln_comp_gen,demo)
eln_comp_gen_demo_without_age <- c(eln_comp_gen,demo_without_age)

eln_comp_cyto_clin <- c(eln_comp_cyto,clin)
eln_comp_cyto_demo <- c(eln_comp_cyto,demo)
eln_comp_cyto_demo_without_age <- c(eln_comp_cyto,demo_without_age)

eln_comp_clin_demo <- c(eln_comp_clin,demo)
eln_comp_clin_demo_without_age <- c(eln_comp_clin,demo_without_age)

eln_comp_age_gen_cyto <- c(eln_comp_age_gen,cyto)
eln_comp_age_gen_clin <- c(eln_comp_age_gen,clin)
eln_comp_age_gen_demo <- c(eln_comp_age_gen,demo)
eln_comp_age_gen_demo_without_age <- c(eln_comp_age_gen,demo_without_age)

eln_comp_gen_cyto_clin_demo <- c(eln_comp_gen_cyto,clin,demo)

# without comp

eln_age_gen <- c(eln_age,gen)
eln_age_cyto <- c(eln_age,cyto)
eln_age_clin <- c(eln_age,clin)

eln_gen_cyto <- c(eln_gen,cyto)
eln_gen_clin <- c(eln_gen,clin)
eln_gen_demo <- c(eln_gen,demo)
eln_gen_demo_without_age <- c(eln_gen,demo_without_age)

eln_cyto_clin <- c(eln_cyto,clin)
eln_cyto_demo <- c(eln_cyto,demo)
eln_cyto_demo_without_age <- c(eln_cyto,demo_without_age)

eln_clin_demo <- c(eln_clin,demo)
eln_clin_demo_without_age <- c(eln_clin,demo_without_age)

eln_age_gen_cyto <- c(eln_age_gen,cyto)
eln_age_gen_clin <- c(eln_age_gen,clin)
eln_age_gen_demo <- c(eln_age_gen,demo)
eln_age_gen_demo_without_age <- c(eln_age_gen,demo_without_age)

eln_gen_cyto_clin_demo <- c(eln_gen_cyto,clin,demo)

# USEFUL FOR COMP
comp_age <- c(comp,age)
comp_gen <- c(comp,gen)
comp_cyto <- c(comp,cyto)
comp_clin <- c(comp,clin)
comp_demo <- c(comp,demo)
comp_demo_without_age <- c(comp,demo_without_age)
comp_gen_cyto <- c(comp_gen,cyto)
comp_clin_demo <- c(comp_clin,demo)
comp_gen_cyto_clin_demo <- c(comp_gen_cyto,clin,demo)

#USEFUL FOR GEN
gen_age <- c(gen,age)
gen_cyto <- c(gen,cyto)
gen_clin <- c(gen,clin)
gen_demo <- c(gen,demo)
gen_demo_without_age <- c(gen,demo_without_age)
gen_clin_demo <- c(gen_clin,demo)
gen_cyto_clin_demo <- c(gen_cyto,clin,demo)

#USEFUL FOR CYTO 
cyto_age <- c(cyto,age)
cyto_clin <- c(cyto,clin)
cyto_demo <- c(cyto,demo)
gen_demo_without_age <- c(gen,demo_without_age)
cyto_clin_demo <- c(cyto_clin,demo)
cyto_gen_demo <- c(gen_cyto,demo)


clin_age <-c(clin,age)
    
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
          
              
prognosis_features <- list(eln_age_gen=eln_age_gen,eln_age_cyto=eln_age_cyto,eln_age_clin=eln_age_clin,eln_gen_cyto=eln_gen_cyto,eln_gen_clin=eln_gen_clin,eln_gen_demo=eln_gen_demo,eln_gen_demo_without_age=eln_gen_demo_without_age,
                          eln_cyto_clin=eln_cyto_clin,eln_cyto_demo=eln_cyto_demo,eln_cyto_demo_without_age=eln_cyto_demo_without_age,eln_clin_demo=eln_clin_demo,eln_clin_demo_without_age=eln_clin_demo_without_age,
                          eln_age_gen_cyto=eln_age_gen_cyto,eln_age_gen_clin=eln_age_gen_clin,eln_age_gen_demo=eln_age_gen_demo,eln_age_gen_demo_without_age=eln_age_gen_demo_without_age,eln_gen_cyto_clin_demo=eln_gen_cyto_clin_demo)

algos <-c("glm","rfs","boost","rfx")
alphas=c(0,0.7,1)
for (i in 1:length(prognosis_features)){
    for (algo in algos){
        if (algo=="glm"){
            for (alpha in alphas){
                print(alpha)
                print(algo)
                bootstrap <- bootstrapping(prognosis_features[[i]],x,y,100,alpha,16,algo)
                tmp_1 <- bootstrap  %>% group_by(feature) %>% summarise_all(sum)
                tmp_2 <- bootstrap  %>% group_by(feature) %>% count(feature)
                print(paste(paste(names(prognosis_features)[i],paste(algo,alpha,sep="_"),sep="_bootstrap_"),".tsv",sep=""))
                write.table(data.frame(merge(tmp_1,tmp_2,by='feature')),paste(paste(names(prognosis_features)[i],paste(algo,alpha,sep="_"),sep="_bootstrap_"),".tsv",sep=""),quote=F,sep='\t')

                if (alpha==0.7){
                    tmp_1_pos <- tmp_1[tmp_1$coef>0,]
                    tmp_1_neg <-  tmp_1[tmp_1$coef<0,]
                    features_reduced <- union(union(tmp_1_pos[tmp_1_pos$coef > quantile(tmp_1_pos$coef,0.90),]$feature,tmp_1_neg[tmp_1_neg$coef < quantile(tmp_1_neg$coef,0.15),]$feature),tmp_2[tmp_2$n > quantile(tmp_2$n,0.85),]$feature)
                    if (length(features_reduced)<2){features_reduced <- union(union(tmp_1_pos[tmp_1_pos$coef > quantile(tmp_1_pos$coef,0.90),]$feature,tmp_1_neg[tmp_1_neg$coef < quantile(tmp_1_neg$coef,0.15),]$feature),tmp_2[tmp_2$n > 0,]$feature)}
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