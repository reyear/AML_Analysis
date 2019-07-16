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



df_final <- read.table("../data_frame_final_prognosis.tsv")
##---------------------------------------------------------------------------------MODELS TO TRY


eln <- c(1)
comp <- c(164:178)
age <- c(163)
all_gen <- c(2:85)
vect <- apply(X=df_final[,all_gen],2,FUN=function(x) 100*length(which(x==1))/dim(df_final)[1])
gen <- match(names(vect[vect>=2]),names(df_final))
all_cyto <- c(86:154)
vect <- apply(X=df_final[,all_cyto],2,FUN=function(x) 100*length(which(x==1))/dim(df_final)[1])
cyto <- match(names(vect[vect>=2]),names(df_final))
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

eln_comp_gen_cyto_clin_demo <- c(eln_comp_gen_cyto,clin,demo)
eln_comp_gen_cyto_clin_demo_without_age <- c(eln_comp_gen_cyto,clin,demo_without_age)
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
##---------------------------------------------------------------------------------PREPARING MODELS and ALGOS


prognosis_features<-list(eln_comp_gen_cyto=eln_comp_gen_cyto)


                         
 

### PARAMETERS OF ANALYSIS:
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