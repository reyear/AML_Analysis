library(glmnet)
library(doMC)
library(survival)
library(data.table)
library(mltools)
library(CoxBoost)
library(randomForestSRC)
library(CoxHD)
source('../InitialPrognosis/run_prognosis.R')
source("feature_importance.R")

##---------------------------------------------------------------------------------MODELS TO TRY

all_features <-c(1:180) 
clin_demo_comp <-c(155:180) 
clin_demo_cyto_gen_comp <- c(2:180) 
comp <- c(164:180) 
cyto_comp <-c(86:154,164:180) 
cyto_gen_comp <- c(2:154,164:180) 
eln_clin_demo_comp <- c(1,155:180) 
eln_cyto_comp <- c(1,86:154,164:180) 
eln_cyto_gen_comp <- c(1:154,164:180) 
eln_gen_comp <- c(1:85,164:180) 
gen_comp <- c(2:85,164:180) 
clin_comp <- c(155:161,164:180) 
clin_cyto_comp <- c(86:161,164:180) 
clin_gen_comp <- c(2:85,155:161,164:180) 
eln_clin_comp <- c(1,155:161,164:180) 

#Without age
all_features_without_age <-c(1:162,164:180) 
clin_demo_comp_without_age <-c(155:162,164:180) 
clin_demo_cyto_gen_comp_without_age <- c(2:162,164:180) 
eln_clin_demo_comp_without_age <- c(1,155:162,164:180) 




###With age:
eln_clin_gen <-  c(1:85,155:161)
eln_demo_gen <- c(1:85,162:163)
eln_clin_demo_cyto_gen <- c(1:163)
eln_clin_demo_cyto <- c(1,86:163)

eln_clin_demo_gen <- c(1:85,155:163)  
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
gen_age <- c(2:85,163)

eln_comp <-c(1,164:180)
eln_age <- c(1,163)
eln_gen <- c(1:85)
eln_cyto <- c(1,86:154)

##---------------------------------------------------------------------------------PREPARING MODELS and ALGOS
df_final <- read.table("../InitialPrognosis/prognosis_comp_final.tsv",sep='\t',header=T)

prognosis_features<-list(eln_clin_demo_cyto_gen_without_age=eln_clin_demo_cyto_gen_without_age,eln_clin_demo_cyto_without_age=eln_clin_demo_cyto_without_age,
eln_clin_demo_gen_without_age=eln_clin_demo_gen_without_age,eln_clin_demo_without_age=eln_clin_demo_without_age,
clin_demo_cyto_gen_without_age=clin_demo_cyto_gen_without_age,clin_demo_cyto_without_age=clin_demo_cyto_without_age,
clin_demo_gen_without_age=clin_demo_gen_without_age,clin_demo_without_age=clin_demo_without_age,demo_gen_without_age=demo_gen_without_age,
demo_cyto_without_age=demo_cyto_without_age,gen_age=gen_age,eln_comp=eln_comp,eln_age=eln_age,eln_gen=eln_gen,eln_cyto=eln_cyto)

                         
 

### PARAMETERS OF ANALYSIS:
nrepeats=5
seed=1234
mc.cores=30
npermutations=4
nfolds=5

algorithms<-c(algo_Lasso, algo_Ridge, algo_Elastic_net,  algo_RFX, algo_RFS, algo_Cox)
predictors<-c(predictor_Lasso, predictor_Ridge, predictor_Elastic_net,  predictor_RFX, predictor_RFS, predictor_BOOST, predictor_Cox)
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