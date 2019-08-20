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
df_final <- read.table("../../../../clustering/clustering_Final_1/df_final_full_component_ITD.tsv")
master <- read.table("../../../../../data/initial_dataset/Master_04_10_2019.tsv")
rownames(master) <- master$data_pd
df_final <- merge(df_final,master[,c("data_pd","intense")],by=0)
rownames(df_final) <- df_final$Row.names
df_final <- df_final[-1]
df_final$data_pd <- NULL

#Usual Features
eln <- c(2,3,4)
comp <- c(170:193)
age <- c(167)
intense <- c(196)

all_gen <- c(5:88)
vect <- apply(X=df_final[,all_gen],2,FUN=function(x) 100*length(which(x==1))/dim(df_final)[1])
gen <- match(names(vect[vect>=2]),names(df_final))
gen_without <- setdiff(gen,grep("^NPM1$", colnames(df_final)))
gen_without <- setdiff(gen_without,grep("^CEBPA_bi$", colnames(df_final))) 

all_cyto <- c(89:158)
vect <- apply(X=df_final[,all_cyto],2,FUN=function(x) 100*length(which(x==1))/dim(df_final)[1])
cyto <- match(names(vect[vect>=2]),names(df_final))
cyto_without <- setdiff(cyto,grep("^inv_16$", colnames(df_final)))
cyto_without <- setdiff(cyto_without,grep("^t_8_21$", colnames(df_final)))      
cyto_without <- setdiff(cyto_without,grep("^t_v_11$", colnames(df_final))) 

clin <- c(159:165)
demo <- c(166:167)
demo_without_age <-c(166)

name_genes <- colnames(df_final[,gen])
name_cyto <- colnames(df_final[,cyto])
name_comp <- colnames(df_final[,comp])
name_eln <- colnames(df_final[,eln])
for (col in c(name_eln,name_genes,name_cyto,name_comp)){
df_final[,paste(col,"intense",sep="_")] <- df_final[,col]*df_final$intense
}    
#Features with intensification

eln_intense <- c(197:199)
comp_intense <- c(255:278) 
gen_intense <- c(200:234)    

gen_without_intense <- setdiff(gen_intense,grep("^NPM1_intense$", colnames(df_final)))
gen_without_intense <- setdiff(gen_without_intense,grep("^CEBPA_bi_intense$", colnames(df_final)))
cyto_intense <- c(235:254)
cyto_without_intense <- setdiff(cyto_intense,grep("^inv_16_intense$", colnames(df_final)))
cyto_without_intense <- setdiff(cyto_without_intense,grep("^t_8_21_intense$", colnames(df_final)))      
cyto_without_intense <- setdiff(cyto_without_intense,grep("^t_v_11_intense$", colnames(df_final))) 



### Models to try
comp_comp_intense <- c (comp, comp_intense)
gen_gen_intense <- c (gen, gen_intense)
cyto_cyto_intense <- c (cyto, cyto_intense)
gen_cyto_gen_intense_cyto_intense <- c(gen,cyto,gen_intense,cyto_intense)
comp_gen_cyto_comp_intense_gen_intense_cyto_intense <- c(comp,gen_without,cyto_without,comp_intense,gen_without_intense,cyto_without_intense)
gen_cyto_clin_demo_gen_intense_cyto_intense <- c(gen_cyto_gen_intense_cyto_intense,clin,demo)
comp_clin_demo_comp_intense <- c(comp_comp_intense,clin,demo)
comp_gen_cyto_clin_demo_comp_intense_gen_intense_cyto_intense <- c(comp_gen_cyto_comp_intense_gen_intense_cyto_intense,clin,demo)
eln_eln_intense <- c(eln,eln_intense)
eln_comp_eln_intense_comp_intense <- c(eln_eln_intense,comp_comp_intense)
eln_gen_cyto_eln_intense_gen_intense_cyto_intense <- c(eln_eln_intense,gen_cyto_gen_intense_cyto_intense)
eln_comp_gen_cyto_eln_intense_comp_intense_gen_intense_cyto_intense <- c(eln_eln_intense,comp_gen_cyto_comp_intense_gen_intense_cyto_intense)
eln_comp_gen_cyto_clin_demo_eln_intense_comp_intense_gen_intense_cyto_intense <- c(eln_comp_gen_cyto_eln_intense_comp_intense_gen_intense_cyto_intense,clin,demo)
comp_intense_gen_intense_cyto_intense <- c(comp_intense,gen_without_intense,cyto_without_intense)
eln_intense_comp_intense_gen_intense_cyto_intense_clin_demo <- c(eln_intense,comp_intense_gen_intense_cyto_intense,clin,demo)
comp_intense_gen_intense <- c(comp_intense,gen_without_intense)
comp_intense_cyto_intense <- c(comp_intense,cyto_without_intense)
eln_intense_gen_intense <- c(eln_intense,gen_intense)
eln_intense_cyto_intense <- c(eln_intense,cyto_intense)
eln_intense_gen_intense_cyto_intense <- c(eln_intense,gen_intense,cyto_intense)

prognosis_features<- list(comp_comp_intense=c(comp_comp_intense,intense))

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
              

algos <-c("glm","rfs","boost","rfx")
alphas=c(0,0.7,1)
for (i in 1:length(prognosis_features)){
    for (algo in algos){
        if (algo=="glm"){
            for (alpha in alphas){
                print(alpha)
                print(algo)
                bootstrap <- bootstrapping(prognosis_features[[i]],x,y,10,alpha,30,algo)
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
                    l_ntree <- c(10)
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

                write.table(data.frame(merge(tmp_1,tmp_2,by='feature')),"test.tsv"),quote=F,sep='\t')
    print ('next')
        }
    }
}