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

correlates <- c("ahd","perf_status","bm_blasts","secondary","wbc","hb","plt","gender","age","os","os_status")

df_all_components <- read.table("../../../data/updated_dataset/refined_components_updated.tsv",sep = '\t' , header = T)
df_initial <- read.table("../../../data/initial_dataset/Master_04_10_2019.csv",sep = ',' , header = T)
df_initial <- read.table("../../../data/initial_dataset/Master_04_10_2019.csv",sep = ',' , header = T)
rownames(df_initial) <- df_initial$data_pd
df_initial <- df_initial[,-1:-4]
df_eln <- read.table("../../../data/updated_dataset/eln_ratio.tsv",sep = '\t' , header = T)
rownames(df_eln) <- df_eln$X
df_eln <- df_eln[-1]
df <- merge(df_eln,df_all_components,by=0)
rownames(df) <- df$Row.names
df <- df[-1]
df <- merge(df,df_initial[,correlates],by=0)
rownames(df) <- df$Row.names
df <- df[-1]
df <- df[,c(1:151,157,156,166:176)]
df <- na.omit(df)
df <- df[df$os>0,]

### To try when component availables: ###

all_features <-c(1:100) #not used
clin_demo_comp <-c(1:100) #not used
clin_demo_cyto_gen_comp <- c(1:100) #not used
comp <- c(1:100) #not used
cyto_comp <-c(1:100) #not used
cyto_gen_comp <- c(1:100) #not used
eln_clin_demo_comp <- c(1:100) #not used
eln_cyto_comp <- c(1:100) #not used
eln_cyto_gen_comp <- c(1:100) #not used
eln_gen_comp <- c(1:100) #not used
gen_comp <- c(1:100) #not used
clin_comp <- c(1:100) #not used
clin_cyto_comp <- c(1:100) #not used
clin_gen_comp <- c(1:100) #not used
eln_clin_comp <- c(1:100) #not used
### To try now: ###

###With age:
eln_clin_gen <-  c(1:84,154:160)
eln_demo_gen <- c(1:84,161:162)
eln_clin_demo_cyto_gen <- c(1:162)
eln_clin_demo_cyto <- c(1,85:162)
eln_clin_demo_gen <- c(1:84,154:162)
eln_clin_demo <- c(1,154:162)
eln_clin <- c(1,154:160)
eln_cyto_gen <- c(1:153)
clin_demo_cyto_gen <- c(2:162)
clin_demo_cyto <- c(85:162)
clin_demo_gen <- c(2:84,154:162)
clin_demo <- c(154:162)
cyto_gen <- c(2:153)
cyto <- c(85:153)
gen <- c(2:84)
clin_gen <- c(2:84,154:160)
clin_cyto <- c(85:160)
demo_gen <- c(2:84,161:162)
demo_cyto <- c(85,153,161:162)

###Without age:
eln_demo_gen_without_age <- c(1:84,161)
eln_clin_demo_cyto_gen_without_age <- c(1:161)
eln_clin_demo_cyto_without_age <- c(1,85:161)
eln_clin_demo_gen_without_age <- c(1:84,154:161)
eln_clin_demo_without_age <- c(1,154:161)
clin_demo_cyto_gen_without_age <- c(2:161)
clin_demo_cyto_without_age <- c(85:161)
clin_demo_gen_without_age <- c(2:84,154:161)
clin_demo_without_age <- c(154:161)
demo_gen_without_age <- c(2:84,161)
demo_cyto_without_age <- c(85,153,161)


y <- data.matrix(df_final[,c("os","os_status")])

prognosis_features<- list(eln_clin_gen=eln_clin_gen,eln_demo_gen=eln_demo_gen,eln_clin_demo_cyto_gen=eln_clin_demo_cyto_gen,eln_clin_demo_cyto=eln_clin_demo_cyto,
                          eln_clin_demo_gen=eln_clin_demo_gen,eln_clin_demo=eln_clin_demo,eln_clin=eln_clin,eln_cyto_gen=eln_cyto_gen,clin_demo_cyto_gen=clin_demo_cyto_gen,
                          clin_demo_cyto=clin_demo_cyto,clin_demo_gen=clin_demo_gen,clin_demo=clin_demo,cyto_gen=cyto_gen,cyto=cyto,gen=gen,clin_gen=clin_gen,clin_cyto=clin_cyto,
                          demo_gen=demo_gen,demo_cyto=demo_cyto,eln_demo_gen_without_age=eln_demo_gen_without_age,eln_clin_demo_cyto_gen_without_age=eln_clin_demo_cyto_gen_without_age,
                          eln_clin_demo_cyto_without_age=eln_clin_demo_cyto_without_age,eln_clin_demo_gen_without_age=eln_clin_demo_gen_without_age,
                          eln_clin_demo_without_age=eln_clin_demo_without_age,clin_demo_cyto_gen_without_age=clin_demo_cyto_gen_without_age,clin_demo_cyto_without_age=clin_demo_cyto_without_age,
                          clin_demo_gen_without_age=clin_demo_gen_without_age,clin_demo_without_age=clin_demo_without_age,demo_gen_without_age=demo_gen_without_age,demo_cyto_without_age=demo_cyto_without_age)

predictors <- c(rep(list(predictorGLM),11),rep(list(predictorRF),7),predictorBoost,predictorRFX)
str_predictors <-c(rep("CoxGLM",11),rep("RFS",28),"CoxBoost","RFX")
l_alpha <-seq(0,1,0.1)
l_ntree <- seq(300,1200,150)
mc.cores <- 8
nodesize <- c(5,10,20,30)
for (i in 1:length(prognosis_features)){
    print("DONE")
    x <- data.matrix(df_final[,prognosis_features[[i]]])
    write.table(launch_prognosis(x=x,y=y,predictors=predictors,str_predictors=str_predictors,l_alpha=l_alpha,nrepeats=2,
                l_ntree=l_ntree,mc.cores=mc.cores,nodesize=nodesize),paste(names(prognosis_features)[i],".tsv",sep=""),quote=F,sep='\t')
    print("DONE")
    }
