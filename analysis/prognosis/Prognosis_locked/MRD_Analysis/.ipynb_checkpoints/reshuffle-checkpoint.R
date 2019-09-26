library(glmnet)
library(doMC)
library(survival)
library(data.table)
library(mltools)
library(CoxBoost)
library(randomForestSRC)
library(CoxHD)
source('../../tools_prognosis/run_prognosis.R')
source("feature_importance.R")



### Features that we can use
###-----------------------------------------------------------------------------
tmp <-read.table('../handovercompiled.Yanis.080919.csv',sep=",",header=T)
rownames(tmp) <- tmp$data_pd
cols_to_keep <- colnames(tmp)
master <- read.table('../../../../data/initial_dataset/Master_04_10_2019.csv',sep=",",header=T)
rownames(master) <- master$data_pd
df <- read.table('../../../clustering/clustering_Final_1/df_final_full_component_ITD.tsv')
df_merge <- merge(df,master[,cols_to_keep],by=0)
rownames(df_merge) <- df_merge$Row.names
df_merge <- df_merge[-1]
df_merge$AMLID <- as.character(df_merge$AMLID)
df_merge$MRD1 <- as.character(df_merge$MRD1)
df_merge <- df_merge[grep("17-", df_merge$AMLID),]   ### keep only AML 17
df_merge$CR_MRD_neg <- ifelse(df_merge$MRD1=="CR,MRD-",1,0)
df_merge$CR_MRD_pos <- ifelse(df_merge$MRD1=="CR,MRD+",1,0)
df_merge$all_others <- ifelse(df_merge$MRD1!="CR,MRD+" & df_merge$MRD1!="CR,MRD-",1,0)
df_merge <- df_merge[!is.na(df_merge$OS_CR),]
df_merge <- df_merge[!is.na(df_merge$MRD1),]
df_merge <- df_merge[df_merge$os >0 & df_merge$OS_CR >0 & df_merge$RFSyears>0,]


eln <- c(2,3,4)
comp <- c(170:193)
age <- c(167)

all_gen <- c(5:88)
vect <- apply(X=df_merge[,all_gen],2,FUN=function(x) 100*length(which(x==1))/dim(df_merge)[1])
gen <- match(names(vect[vect>=2]),names(df_merge))
gen_without <- setdiff(gen,grep("^NPM1$", colnames(df_merge)))
gen_without <- setdiff(gen_without,grep("^CEBPA_bi$", colnames(df_merge))) 
              
all_cyto <- c(89:158)
vect <- apply(X=df_merge[,all_cyto],2,FUN=function(x) 100*length(which(x==1))/dim(df_merge)[1])
cyto <- match(names(vect[vect>=2]),names(df_merge))
cyto_without <- setdiff(cyto,grep("^inv_16$", colnames(df_merge)))
cyto_without <- setdiff(cyto_without,grep("^t_8_21$", colnames(df_merge)))      
cyto_without <- setdiff(cyto_without,grep("^t_v_11$", colnames(df_merge)))        
              
clin <- c(159:165)
demo <- c(166:167)
demo_without_age <-c(166)
mrd <- c(234,235,236)
           
                          
                          
eln <- eln                         
eln_mrd <- c(eln,mrd)
              
comp <- comp               
comp_mrd <- c(comp,mrd)
              
gen <- gen              
gen_mrd <- c(gen,mrd)
              
cyto <- cyto              
cyto_mrd <- c(cyto,mrd)

clin <- clin
clin_mrd <- c(clin,mrd)

demo <- demo
demo_mrd <- c(demo,mrd)
              
clin_demo <- c(clin,demo)
clin_demo_mrd <- c(clin,demo,mrd)
              
gen_cyto <- c(gen,cyto)
gen_cyto_mrd <- c(gen,cyto,mrd)
              
gen_cyto_clin_demo <- c(gen,cyto,clin,demo)              
gen_cyto_clin_demo_mrd <- c(gen,cyto,clin,demo,mrd)
              
comp_clin_demo <- c(comp,clin,demo)              
comp_clin_demo_mrd <- c(comp,clin,demo,mrd)

eln_clin_demo <- c(eln,clin,demo)
eln_clin_demo_mrd <- c(eln,clin,demo,mrd)




prognosis_features<- list(eln_mrd=eln_mrd,comp_mrd=comp_mrd,gen_mrd=gen_mrd,cyto_mrd=cyto_mrd,clin_mrd=clin_mrd,demo_mrd=demo_mrd,clin_demo_mrd=clin_demo_mrd,gen_cyto_mrd=gen_cyto_mrd,gen_cyto_clin_demo_mrd=gen_cyto_clin_demo_mrd,comp_clin_demo_mrd=comp_clin_demo_mrd,eln_clin_demo_mrd=eln_clin_demo_mrd,
                         eln=eln,comp=comp,gen=gen,cyto=cyto,clin=clin,demo=demo,clin_demo=clin_demo,gen_cyto=gen_cyto,gen_cyto_clin_demo=gen_cyto_clin_demo,comp_clin_demo=comp_clin_demo,eln_clin_demo=eln_clin_demo)
                          
              
##---------------------------------------------------------------------------------PREPARING MODELS and ALGOS
                         

nrepeats=5
seed=1234
mc.cores=30
npermutations=4
nfolds=5

algorithms<-c(algo_Lasso, algo_Ridge, algo_Elastic_net,  algo_RFX, algo_RFS, algo_Cox)
predictors<-c(predictor_Lasso, predictor_Ridge, predictor_Elastic_net,  predictor_RFX, predictor_RFS,  predictor_Cox)
algo_names<-c('Lasso','Ridge','Elastic_net','RFX','RFS','Cox')

response_os <- data.matrix(df_merge[,c("os","os_status")])
response_RFS <- data.matrix(df_merge[,c("RFSyears","RFSStatus")])
response_OS_CR <- data.matrix(df_merge[,c("OS_CR","DEADStatus")])
              
colnames(response_os) <- c("time","status")
colnames(response_RFS) <- c("time","status")
colnames(response_OS_CR) <- c("time","status")
              
              
##---------------------------------------------------------------------------------PREPARING MODELS and ALGOS

for (j in 1:length(prognosis_features)){
    print(names(prognosis_features[j]))
    for (k in 1:3){
        res_data <- data.frame('feature'=character(),'ref_CI'=numeric(),'permuted_CI'=numeric(),'algo'=character(),'model'=character())
        for(i in 1:length(algorithms)){
            design <- data.matrix(data.frame(df_merge[,prognosis_features[[j]]]))  
            if(k==1){
                response <- response_os
            }else if(k==2){
                response <- response_RFS
            }else{
                response <- response_OS_CR
            }
            tmp <- runCV_CI_with_test(response=response, design=design,
                  nfolds=nfolds, nrepeats=nrepeats, seed=seed, mc.cores=mc.cores, features=colnames(design), npermutations=npermutations, 
                                      algorithm=algorithms[i][[1]], predictor=predictors[i][[1]])
            tmp$algo<-algo_names[i]
            tmp$model <- names(prognosis_features[j])
            res_data <- rbind(res_data,tmp)
        }
        write.table(res_data,paste(names(prognosis_features)[j],ifelse(k==1,"_os.tsv",ifelse(k==2,"_RFS.tsv","_OS_CR.tsv")),sep="_reshuffle_importance"),quote=F,sep='\t')
    }
}
              
