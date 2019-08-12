library(glmnet)
library(doMC)
library(survival)
library(data.table)
library(mltools)
library(CoxBoost)
library(randomForestSRC)
library(CoxHD)
source('../../../tools_prognosis/run_prognosis.R')
source("feature_importance.R")




### Features that we can use
###-----------------------------------------------------------------------------
df_final <- read.table("../../../../clustering/clustering_Final_1/df_final_full_component_ITD.tsv")
master <- read.table("../../../../../data/initial_dataset/Master_04_10_2019.csv",sep= ",",header=T)
rownames(master) <- master$data_pd
df_final <- merge(df_final,master[,c("data_pd","intense")],by=0)
rownames(df_final) <- df_final$Row.names
df_final <- df_final[-1]
df_final$data_pd <- NULL

#Usual Features
eln <- c(2,3,4)
comp <- c(170:193)
age <- c(167)

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
eln_intense_comp_intense_gen_intense_cyto_intense_clin_demo <- c(eln_intense,comp_intense_gen_intense_cyto_intense,clin,demo )
comp_intense_gen_intense <- c(comp_intense,gen_without_intense)
comp_intense_cyto_intense <- c(comp_intense,cyto_without_intense)
eln_intense_gen_intense <- c(eln_intense,gen_intense)
eln_intense_cyto_intense <- c(eln_intense,cyto_intense)
eln_intense_gen_intense_cyto_intense <- c(eln_intense,gen_intense,cyto_intense)
              
y <- data.matrix(df_final[,c("os","os_status")])

prognosis_features<- list(comp_comp_intense=comp_comp_intense,gen_gen_intense=gen_gen_intense,cyto_cyto_intense=cyto_cyto_intense,gen_cyto_gen_intense_cyto_intense=gen_cyto_gen_intense_cyto_intense,
                         comp_gen_cyto_comp_intense_gen_intense_cyto_intense=comp_gen_cyto_comp_intense_gen_intense_cyto_intense,gen_cyto_clin_demo_gen_intense_cyto_intense=gen_cyto_clin_demo_gen_intense_cyto_intense,
                         comp_clin_demo_comp_intense=comp_clin_demo_comp_intense,comp_gen_cyto_clin_demo_comp_intense_gen_intense_cyto_intense=comp_gen_cyto_clin_demo_comp_intense_gen_intense_cyto_intense,
                         eln_eln_intense=eln_eln_intense,eln_comp_eln_intense_comp_intense=eln_comp_eln_intense_comp_intense,eln_gen_cyto_eln_intense_gen_intense_cyto_intense=eln_gen_cyto_eln_intense_gen_intense_cyto_intense,
                         eln_comp_gen_cyto_eln_intense_comp_intense_gen_intense_cyto_intense=eln_comp_gen_cyto_eln_intense_comp_intense_gen_intense_cyto_intense,
                         eln_comp_gen_cyto_clin_demo_eln_intense_comp_intense_gen_intense_cyto_intense=eln_comp_gen_cyto_clin_demo_eln_intense_comp_intense_gen_intense_cyto_intense,
                         comp_intense_gen_intense_cyto_intense=comp_intense_gen_intense_cyto_intense,eln_intense_comp_intense_gen_intense_cyto_intense_clin_demo=eln_intense_comp_intense_gen_intense_cyto_intense_clin_demo,
                         comp_intense_gen_intense=comp_intense_gen_intense,comp_intense_cyto_intense=comp_intense_cyto_intense,eln_intense_gen_intense=eln_intense_gen_intense,eln_intense_cyto_intense=eln_intense_cyto_intense,
                         eln_intense_gen_intense_cyto_intense=eln_intense_gen_intense_cyto_intense)

##---------------------------------------------------------------------------------PREPARING MODELS and ALGOS
                         

nrepeats=5
seed=1234
mc.cores=30
npermutations=4
nfolds=5

algorithms<-c(algo_Lasso, algo_Ridge, algo_Elastic_net,  algo_RFX, algo_RFS, algo_Cox)
predictors<-c(predictor_Lasso, predictor_Ridge, predictor_Elastic_net,  predictor_RFX, predictor_RFS,  predictor_Cox)
algo_names<-c('Lasso','Ridge','Elastic_net','RFX','RFS','Cox')

response <- data.matrix(df_final[,c("os","os_status")])
colnames(response) <- c("time","status")


##---------------------------------------------------------------------------------PREPARING MODELS and ALGOS

for (j in 1:length(prognosis_features)){
    print(names(prognosis_features[j]))
    res_data <- data.frame('feature'=character(),'ref_CI'=numeric(),'permuted_CI'=numeric(),'algo'=character(),'model'=character())
    for(i in 1:length(algorithms)){
        design <- data.matrix(data.frame(df_final[,prognosis_features[[j]]]))      
        tmp <- runCV_CI_with_test(response=response, design=design,
              nfolds=nfolds, nrepeats=nrepeats, seed=seed, mc.cores=mc.cores, features=colnames(design), npermutations=npermutations, 
                                  algorithm=algorithms[i][[1]], predictor=predictors[i][[1]])
        tmp$algo<-algo_names[i]
        tmp$model <- names(prognosis_features[j])
        res_data <- rbind(res_data,tmp)
    }
    write.table(res_data,paste(names(prognosis_features)[j],".tsv",sep="_reshuffle_importance"),quote=F,sep='\t')
}