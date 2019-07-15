library(glmnet)
library(doMC)
library(survival)
library(data.table)
library(mltools)
library(CoxBoost)
library(randomForestSRC)
library(CoxHD)
source('../../tools_prognosis/run_prognosis.R')

df_final <- read.table("../data_frame_final_prognosis.tsv")


df_final$eln_adverse <- ifelse(df_final$eln_2017_ratio==1,1,0)
df_final$eln_intermediate <- ifelse(df_final$eln_2017_ratio==2,1,0)
df_final$eln_favorable <- ifelse(df_final$eln_2017_ratio==3,1,0)

eln <- c(181:183)
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


y <- data.matrix(df_final[,c("os","os_status")])

prognosis_features<- list(eln_age=eln_age,eln_clin=eln_clin,eln_demo=eln_demo,eln_demo_without_age=eln_demo_without_age,eln_comp_age=eln_comp_age,eln_comp_clin=eln_comp_clin,eln_comp_demo=eln_comp_demo,eln_comp_demo_without_age=eln_comp_demo_without_age,
     eln_comp_age_gen=eln_comp_age_gen,eln_comp_age_cyto=eln_comp_age_cyto,eln_comp_age_clin=eln_comp_age_clin,eln_comp_gen_clin=eln_comp_gen_clin,eln_comp_gen_demo=eln_comp_gen_demo,
     eln_comp_gen_demo_without_age=eln_comp_gen_demo_without_age,eln_comp_cyto_clin=eln_comp_cyto_clin,eln_comp_cyto_demo=eln_comp_cyto_demo,eln_comp_cyto_demo_without_age=eln_comp_cyto_demo_without_age,
     eln_comp_clin_demo=eln_comp_clin_demo,eln_comp_clin_demo_without_age=eln_comp_clin_demo_without_age,eln_comp_cyto_clin=eln_comp_cyto_clin,eln_comp_cyto_demo=eln_comp_cyto_demo,eln_comp_cyto_demo_without_age=eln_comp_cyto_demo_without_age,
     eln_comp_clin_demo=eln_comp_clin_demo,eln_comp_clin_demo_without_age=eln_comp_clin_demo_without_age,eln_comp_age_gen_cyto=eln_comp_age_gen_cyto,eln_comp_age_gen_clin=eln_comp_age_gen_clin,eln_comp_age_gen_demo=eln_comp_age_gen_demo,
     eln_comp_gen_cyto_clin_demo=eln_comp_gen_cyto_clin_demo,eln_comp_gen_cyto_clin_demo_without_age=eln_comp_gen_cyto_clin_demo_without_age,eln_age_gen=eln_age_gen,eln_age_cyto=eln_age_cyto,eln_age_clin=eln_age_clin,
     eln_gen_clin=eln_gen_clin,eln_gen_demo=eln_gen_demo,eln_gen_demo_without_age=eln_gen_demo_without_age,eln_cyto_clin=eln_cyto_clin,eln_cyto_demo=eln_cyto_demo,eln_cyto_demo_without_age=eln_cyto_demo_without_age,eln_clin_demo=eln_clin_demo,
     eln_clin_demo_without_age=eln_clin_demo_without_age,eln_age_gen_cyto=eln_age_gen_cyto,eln_age_gen_clin=eln_age_gen_clin,eln_age_gen_demo=eln_age_gen_demo,eln_gen_cyto_clin_demo=eln_gen_cyto_clin_demo,                                           eln=eln,eln_comp=eln_comp,eln_gen=eln_gen,eln_cyto=eln_cyto,eln_gen_cyto=eln_gen_cyto,eln_comp_gen=eln_comp_gen,eln_comp_cyto=eln_comp_cyto,eln_comp_gen_cyto=eln_comp_gen_cyto,
                          eln_comp_gen_cyto_clin_demo=eln_comp_gen_cyto_clin_demo,
                        eln_comp_gen_cyto_clin_demo_without_age=eln_comp_gen_cyto_clin_demo_without_age)

predictors <- c(rep(list(predictorGLM),6),rep(list(predictorRF),1),predictorBoost,predictorRFX)
str_predictors <-c(rep("CoxGLM",6),rep("RFS",1),"CoxBoost","RFX")
l_alpha <-seq(0,1,0.2)
l_ntree <- c(1050)
mc.cores <- 30
nodesize <- c(20)
for (i in 1:length(prognosis_features)){
    print("DONE")
    x <- data.matrix(df_final[,prognosis_features[[i]]])
    write.table(launch_prognosis(x=x,y=y,predictors=predictors,str_predictors=str_predictors,l_alpha=l_alpha,nrepeats=5,
                l_ntree=l_ntree,mc.cores=mc.cores,nodesize=nodesize),paste(names(prognosis_features)[i],".tsv",sep=""),quote=F,sep='\t')
    print("DONE")
    }
