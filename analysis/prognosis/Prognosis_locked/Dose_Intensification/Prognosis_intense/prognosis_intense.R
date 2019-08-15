library(glmnet)
library(doMC)
library(survival)
library(data.table)
library(mltools)
library(CoxBoost)
library(randomForestSRC)
library(CoxHD)
source('../../../tools_prognosis/run_prognosis.R')


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

              
all_cyto <- c(89:158)
vect <- apply(X=df_final[,all_cyto],2,FUN=function(x) 100*length(which(x==1))/dim(df_final)[1])
cyto <- match(names(vect[vect>=2]),names(df_final))

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

eln <- c(197:199)
comp <- c(255:278) 
gen <- c(200:234)                 
cyto <- c(235:254)
              
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
prognosis_features<- list(eln_gen=c(eln_gen,intense), eln_cyto=c(eln_cyto,intense),eln_clin=c(eln_clin,intense),eln_demo=c(eln_demo,intense),eln_gen_cyto=c(eln_gen_cyto,intense),eln_clin_demo=c(eln_clin_demo,intense),
                          eln_clin=c(eln_clin,intense),eln_demo=c(eln_demo,intense),eln_gen_cyto_clin_demo=c(eln_gen_cyto_clin_demo,intense),
                          gen=c(gen,intense) ,cyto =c(cyto,intense),gen_cyto=c(gen_cyto,intense),gen_clin=c(gen_clin,intense),gen_demo=c(gen_demo,intense),gen_clin_demo=c(gen_clin_demo,intense),gen_cyto_clin_demo=c(gen_cyto_clin_demo,intense),
                          cyto_clin=c(cyto_clin,intense),cyto_demo=c(cyto_demo,intense),cyto_clin_demo=c(cyto_clin_demo,intense),cyto_gen_demo=c(cyto_gen_demo,intense),
                          eln_comp_gen=c(eln_comp_gen,intense),eln_comp=c(eln_comp,intense),comp=c(comp,intense),eln_comp_cyto=c(eln_comp_cyto,intense),eln_comp_clin=c(eln_comp_clin,intense),eln_comp_demo=c(eln_comp_demo,intense),
                          eln_comp_gen_cyto=c(eln_comp_gen_cyto,intense),eln_comp_gen_clin=c(eln_comp_gen_clin,intense),eln_comp_gen_demo=c(eln_comp_gen_demo,intense),
                          eln_comp_cyto_clin=c(eln_comp_cyto_clin,intense),eln_comp_cyto_demo=c(eln_comp_cyto_demo,intense),eln_comp_clin_demo=c(eln_comp_clin_demo,intense),
                          eln_comp_gen_cyto_clin_demo=c(eln_comp_gen_cyto_clin_demo,intense),eln_comp_gen_cyto_clin_demo_without_age=c(eln_comp_gen_cyto_clin_demo_without_age,intense),
                          comp_gen=c(comp_gen,intense),comp_cyto=c(comp_cyto,intense),comp_clin=c(comp_clin,intense),comp_demo=c(comp_demo,intense),comp_gen_cyto=c(comp_gen_cyto,intense),comp_clin_demo=c(comp_clin_demo,intense),
                          comp_clin=c(comp_clin,intense),comp_demo=c(comp_demo,intense),comp_gen_cyto=c(comp_gen_cyto,intense),comp_clin_demo=c(comp_clin_demo,intense),comp_clin=c(comp_clin,intense),comp_demo=c(comp_demo,intense),
                          comp_gen_cyto_clin_demo=c(comp_gen_cyto_clin_demo,intense))
              
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
