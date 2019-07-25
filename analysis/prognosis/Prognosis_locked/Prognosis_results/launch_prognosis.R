library(glmnet)
library(doMC)
library(survival)
library(data.table)
library(mltools)
library(CoxBoost)
library(randomForestSRC)
library(CoxHD)
source('../../tools_prognosis/run_prognosis.R')

df_final <- read.table("../prognosis_dataframe.tsv")


eln <- c(2,3,4)
comp <- c(184,198:213)
comp_overlap <- c(168:197)
age <- c(167)

all_gen <- c(5:88)
vect <- apply(X=df_final[,all_gen],2,FUN=function(x) 100*length(which(x==1))/dim(df_final)[1])
gen <- match(names(vect[vect>=2]),names(df_final))
              
all_cyto <- c(89:158)
vect <- apply(X=df_final[,all_cyto],2,FUN=function(x) 100*length(which(x==1))/dim(df_final)[1])
cyto <- match(names(vect[vect>=2]),names(df_final))
              
clin <- c(159:165)
demo <- c(166:167)
demo_without_age <-c(166)
           
                          
                          
eln_comp <- c(eln,comp)          
eln_gen <- c(eln,gen)
eln_cyto <- c(eln,cyto)
eln_clin <- c(eln,clin)
eln_demo <- c(eln,demo)

# USEFUL FOR ELN COMPARISON
# with comp
eln_comp_gen <- c(eln_comp,gen)
eln_comp_cyto <- c(eln_comp,cyto)
eln_comp_clin <- c(eln_comp,clin)
eln_comp_demo <- c(eln_comp,demo)


eln_comp_gen_cyto <- c(eln_comp_gen,cyto)
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
comp_gen <- c(comp,gen)
comp_cyto <- c(comp,cyto)
comp_clin <- c(comp,clin)
comp_demo <- c(comp,demo)
comp_gen_cyto <- c(comp_gen,cyto)
comp_clin_demo <- c(comp_clin,demo)
comp_clin <- c(comp,clin)
comp_demo <- c(comp,demo)
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



y <- data.matrix(df_final[,c("os","os_status")])

prognosis_features<- list(eln=eln,comp=comp,gen=gen,cyto=cyto,clin=clin,demo=demo,eln_gen=eln_gen,eln_cyto=eln_cyto,eln_clin=eln_clin,
                          eln_demo=eln_demo,eln_comp_cyto=eln_comp_cyto,eln_comp_clin=eln_comp_clin,eln_comp_demo=eln_comp_demo,
                          eln_comp_gen_cyto=eln_comp_gen_cyto,eln_comp_gen_clin=eln_comp_gen_clin,eln_comp_gen_demo=eln_comp_gen_demo,
                          eln_comp_cyto_clin=eln_comp_cyto_clin,eln_comp_cyto_demo=eln_comp_cyto_demo,eln_comp_clin_demo=eln_comp_clin_demo,
                          eln_comp_gen_cyto_clin_demo=eln_comp_gen_cyto_clin_demo,eln_comp_gen_cyto_clin_demo_without_age=eln_comp_gen_cyto_clin_demo_without_age,
                          eln_gen_cyto=eln_gen_cyto,eln_clin_demo=eln_clin_demo,eln_clin=eln_clin,eln_demo=eln_demo,eln_gen_cyto_clin_demo=eln_gen_cyto_clin_demo,
                          comp_gen=comp_gen,comp_cyto=comp_cyto,comp_clin=comp_clin,comp_demo=comp_demo,comp_gen_cyto=comp_gen_cyto,comp_clin_demo=comp_clin_demo,
                          comp_clin=comp_clin,comp_demo=comp_demo,comp_gen_cyto=comp_gen_cyto,comp_clin_demo=comp_clin_demo,comp_clin=comp_clin,comp_demo=comp_demo,
                          comp_gen_cyto_clin_demo=comp_gen_cyto_clin_demo,gen_cyto=gen_cyto,gen_clin=gen_clin,gen_demo=gen_demo,gen_clin_demo=gen_clin_demo,
                          gen_cyto_clin_demo=gen_cyto_clin_demo,cyto_clin=cyto_clin,cyto_demo=cyto_demo,cyto_clin_demo=cyto_clin_demo,cyto_gen_demo=cyto_gen_demo)
                          
                          

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
