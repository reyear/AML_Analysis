library(glmnet)
library(doMC)
library(survival)
library(data.table)
library(mltools)
library(CoxBoost)
library(randomForestSRC)
library(CoxHD)
source('run_prognosis.R')
# RFS=105=15*7 for str but only right rep(list(predictorRF),length(l_ntree)) because there is an inside for loop

df_final <- read.table("prognosis_comp_final.tsv",sep='\t',header=T)
df_final$eln_favorable <- ifelse(df_final$eln_2017_ratio==3,1,0)
df_final$eln_intermediate <- ifelse(df_final$eln_2017_ratio==2,1,0)
df_final$eln_adverse <- ifelse(df_final$eln_2017_ratio==1,1,0)
df_final$NC_chr_splicing <- df_final$final_component_NC1 + df_final$final_component_NC2 + df_final$final_component_NC3 +df_final$final_component_NC6

eln <- c(184,185,186)
eln_comp <- c(164:180,184,185,186)
eln_age <- c(184,185,186,163)
eln_gen <- c(2:85,184,185,186)
eln_cyto <- c(184,185,186,86:154)
eln_comp_merged <- c(164:173,177,179,180,184,185,186,187)


prognosis_features<- list(eln=eln,eln_comp=eln_comp,eln_age=eln_age,eln_gen=eln_gen,eln_cyto=eln_cyto,eln_comp_merged=eln_comp_merged)

predictors <- c(rep(list(predictorGLM),11),rep(list(predictorRF),7),predictorBoost,predictorRFX)
str_predictors <-c(rep("CoxGLM",11),rep("RFS",28),"CoxBoost","RFX")
l_alpha <-seq(0,1,0.1)
l_ntree <- seq(300,1200,150)
mc.cores <- 50
nodesize <- c(5,10,20,30)
for (i in 1:length(prognosis_features)){
    print("DONE")
    x <- data.matrix(data.frame(df_final[,prognosis_features[[i]]]))
    write.table(launch_prognosis(x=x,y=y,predictors=predictors,str_predictors=str_predictors,l_alpha=l_alpha,nrepeats=5,
                l_ntree=l_ntree,mc.cores=mc.cores,nodesize=nodesize),paste(names(prognosis_features)[i],".tsv",sep=""),quote=F,sep='\t')
    print("DONE")
    }