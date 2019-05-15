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

df_final <- read.table("df_prognosis_features_ready.tsv",sep = '\t' , header = T)

### Different combinations to try ###
all_features <- c(1:177)
eln_clin<-c(1,169:175)
eln_clin_demo<-c(1,169:177)
eln_clin_demo_cyto <-c(1,85:153,169:177)
eln_clin_demo_gen <-c(1:84,169:177)
eln_clin_demo_cyto_gen <-c(1:153,169:177)
eln_clin_demo_comp <-c(1,154:177)
eln_cyto_gen<-c(1:153)
eln_cyto_gen_comp <-c(1:168)
eln_cyto_comp <-c(1,85:168)
eln_gen_comp <-c(1:84,154:168)

clin_demo<-c(169:177)
clin_demo_cyto <-c(85:153,169:177)
clin_demo_gen <-c(2:84,169:177)
clin_demo_cyto_gen <-c(2:153,169:177)
clin_demo_comp <-c(154:177)
cyto_gen<-c(2:153)
cyto_gen_comp <-c(2:168)
cyto_comp <-c(85:168)
gen_comp <-c(2:84,154:168)
clin_demo_cyto_gen_comp<-c(2:177)
gen<-c(2:84)
cyto<-c(85:153)
comp<-c(154:168)

#    eln_clin_demo_cyto_gen=eln_clin_demo_cyto_gen,eln_clin_demo_comp=eln_clin_demo_comp,
##    eln_cyto_gen=eln_cyto_gen,eln_cyto_gen_comp=eln_cyto_gen_comp,eln_cyto_comp=eln_cyto_comp,
#    eln_gen_comp=eln_gen_comp,clin_demo=clin_demo,clin_demo_cyto=clin_demo_cyto,clin_demo_gen=clin_demo_gen,
#    clin_demo_cyto_gen=clin_demo_cyto_gen,clin_demo_comp=clin_demo_comp,cyto_gen=cyto_gen,cyto_gen_comp=cyto_gen_comp,
#    cyto_comp=cyto_comp,gen_comp=gen_comp,clin_demo_cyto_gen_comp=clin_demo_cyto_gen_comp,gen=gen,cyto=cyto,comp=comp
prognosis_features<- list(all_features=all_features,eln_clin=eln_clin,eln_clin_demo=eln_clin_demo,
    eln_clin_demo_cyto=eln_clin_demo_cyto,eln_clin_demo_gen=eln_clin_demo_gen)
###--------------------------------------------------
y <- data.matrix(df_final[,c("os","os_status")])

predictors <- c(rep(list(predictorGLM),20),rep(list(predictorRF),15),predictorBoost,predictorRFX)
str_predictors <-c(rep("CoxGLM",20),rep("RFS",105),"CoxBoost","RFX")
l_alpha <-seq(0.05,1,0.05)
l_ntree <- seq(100,1500,100)
mc.cores <- 8
nodesize <- c(5,10,15,20,30,40,50)
for (i in 1:length(prognosis_features)){
    x <- data.matrix(df_final[,prognosis_features[[i]]])
    write.table(launch_prognosis(x=x,y=y,predictors=predictors,str_predictors=str_predictors,l_alpha=l_alpha,
                l_ntree=l_ntree,mc.cores=mc.cores,nodesize=nodesize),paste(names(prognosis_features)[i],".tsv",sep=""),quote=F,sep='\t')
    print("DONE")
    }
