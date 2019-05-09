# GET READY
source("../../../src/merge_df.R")
source("./my_survival.R")
library('reshape2')
library('glmnet')
library('tibble')
library('grid')
library('gridExtra')

# -------------- DATA -------------- #
# CLINICAL DF
ddc = read.table("../../../data/updated_dataset/df_clinical_cyto_selected_process_ippsr.tsv",stringsAsFactors=F,header=T,sep="\t")
# MUTATIONS MAF
ddmaf = read.table("../../../data/updated_dataset/df_maf_driver_selected.tsv",stringsAsFactors=F,header=T,sep="\t")
# CYTO DF
ddcyto = read.table("../../../data/updated_dataset/df_cyto_binary_impute.tsv",stringsAsFactors=F,header=T,sep="\t")
# COMPONENTS DF
ddcomponents = read.table("../../../data/updated_dataset/df_components.tsv",stringsAsFactors=F,header=T,sep="\t")
# GENE INTERACTION DF
ddint_double = read.table("../../../data/updated_dataset/datagene_interactions_duplicate_cut1_v2.tsv",stringsAsFactors=F,header=T,sep="\t")
ddint_triple = read.table("../../../data/updated_dataset/datagene_interactions_triplicate_cut2_v2.tsv",stringsAsFactors=F,header=T,sep="\t")
# ---------------------------------- #
if(identical(ddc$LEUKID,ddcyto$LEUKID)) print("all good")
# CREATE MUTATIONS DF
res = merge_clinical_mutation(dd_clinical=ddc, dd_maf=ddmaf, binary=FALSE, col_field="EVENT")
#dd = res$ddmerge # merge
#ddmut = res$ddmut # binary mutation matrix
#if(identical(ddc$LEUKID,rownames(ddmut))) print("all good ^ 2")
# ---------------------------------- #
# ADD METRIC
ddc$AGE = round(ddc$AGE_AT_SAMPLE_TIME)
# -> number of cyto aberrations
ddcyto.reduce = ddcyto[,-which(colnames(ddcyto)%in%c("LEUKID","del4q24"))]
ddc$num_cyto_aberrations = apply(ddcyto.reduce,1,sum)
# -> correct for r_1_7
j17 = which(ddcyto[,"r_1_7"]==1 & ddcyto[,"plus1q"]==1 & ddcyto[,"del7q"]==1)
ddc$num_cyto_aberrations[j17] = ddc$num_cyto_aberrations[j17] - 1
# -> correct for iso17q
jiso17 = which(ddcyto[,"iso17q"]==1 & ddcyto[,"del17p"]==1 & ddcyto[,"plus17q"]==1)
ddc$num_cyto_aberrations[jiso17] = ddc$num_cyto_aberrations[jiso17] - 1
# -> total number of mutations + del4q24
ddc$num_mutations = apply(res$ddmut,1,sum) + ddcyto[,"del4q24"]
# ---------------------------------- #
# GET READY DATAFRAMES FOR PREDICTIONS
# ddc is clinical df
# ddmutation
ddmut5 = read.table("../../../data/updated_dataset/df_mut_hotspot_cut5.tsv", stringsAsFactors=F, header=T, sep="\t")
# cyto with features NK and CK
ddcyto5 = read.table("../../../data/updated_dataset/df_cyto_binary_impute_cut5.tsv", stringsAsFactors=F, header=T, sep="\t")
id <-ddcyto5[,1]
rownames(ddcyto5)<-id
ddcyto5<-ddcyto5[ddc$LEUKID,]
ddcyto5 = ddcyto5[,-1]
ddcyto5$NK = 0
ddcyto5$NK[ddc$num_cyto_aberrations==0] = 1
ddcyto5$CK = 0
ddcyto5$CK[ddc$num_cyto_aberrations>=3] = 1
# ---------------------------------- #
# WITH TP53 ALLELIC STATUS
ddtp53 = read.table("../../../data/updated_dataset/LEUKID_P53_status.tsv",sep="\t",header=T,stringsAsFactors=F)
rownames(ddtp53) = ddtp53$LEUKID ; ddtp53 = ddtp53[,-1]
ddtp53 = ddtp53[ddc$LEUKID,]
ddmut5[,"TP53"] = NULL
rownames(ddmut5)<-id
ddmut5<-ddmut5[ddc$LEUKID,]
ddmut5 = cbind(ddmut5,ddtp53)
ddmut5[,"U2AF1"]=ddmut5[,"U2AF1_157"]+ddmut5[,"U2AF1_34"]
ddmut5[,"U2AF1_157"]=NULL
ddmut5[,"U2AF1_34"]=NULL

rec<-colSums(ddcyto5)
rec<-rec[1:(length(rec)-2)]
rec<-sort(rec,decreasing = T)
rec<-rec[1:50]
rec<-names(rec)
# ---------------------------------- #

# ---------------------------------- #
# DEFINE FEATURE GROUPS
# ---------------------------------- #
nuisance = c("CENTER")
ipssr = c("AGE","HB","PLT","ANC","BM_BLAST","CYTO_IPSSR","IPSSR_CALCULATED")
clinicalother = c("SEX")
mutnumber = c("num_mutations")
cytonumber = c("num_cyto_aberrations")
mutations = colnames(ddmut5)
cytoaberrations = colnames(ddcyto5)
age=c("AGE")
ipssr.score=c("IPSSR_CALCULATED")
ipssr.wt.score = c("AGE","HB","PLT","ANC","BM_BLAST","CYTO_IPSSR")
wbc=c("WBC")
geneint_double=colnames(ddint_double)
geneint_triple=colnames(ddint_triple)
clin=c("HB","PLT","ANC","BM_BLAST")
component=c("component")
topcyto=c(rec,"CK","NK")
blast=c("BM_BLAST")

list.features = list(nuisance=nuisance,
                  ipssr=ipssr,
                  clinicalother=clinicalother,
                  mutnumber=mutnumber,
                  cytonumber=cytonumber,
                  mutations=mutations,
                  cytoaberrations=cytoaberrations,
                  age=age,
                  ipssr.score=ipssr.score,
                  ipssr.wt.score=ipssr.wt.score,
                  wbc=wbc,
                  geneint_double=geneint_double,
                  geneint_triple=geneint_triple,
                  clin=clin,
                  component=component,
                  topcyto=topcyto,
                  blast=blast)

# ---------------------------------- #
# COMPLETE CASE FOR OS
# ---------------------------------- #
# COMPLETE CASE WITH OS RESPONSE
i.complete = which(ddc$ipssr_complete_case=="complete")
response.name = c("os_diag_years","os_status")
tmp = ddc[i.complete,response.name]
i.survival = i.complete[-which(tmp[,1]<=0 | is.na(tmp[,1]) | is.na(tmp[,2]))]
ddcomponents = ddcomponents[as.character(ddc$LEUKID),]
ddgo = cbind(ddc,ddmut5,ddcyto5)
ddgo$component=ddcomponents
ddgo = ddgo[i.survival,]
ddgo = cbind(ddgo,ddint_double,ddint_triple)
ddgo = ddgo[,c(unlist(list.features),response.name)] # FINAL DF FOR OS PREDICTION STUDY
ddgo[is.na(ddgo[, "WBC"]), "WBC"] <- mean(ddgo[, "WBC"], na.rm = TRUE) # MEAN IMPUTATION FOR WBC MISSING
ddgo$CENTER <- paste0('_', ddgo$CENTER)




# ---------------------------------- #
# ANALYSIS
# ---------------------------------- #
nfolds=5
nrepeats=50
seed=1346
mc.cores=20

# ---------------- #
# MODEL COMPARISON
# ---------------- #
predictors <- c(predictorRidge,predictorElastic_net,predictorLasso,predictorCPSS,predictorRFX,predictorRF,predictorBOOST)


# 1) CLIN DEMO CYTO
prep.CLIN_DEMO_CYTO = PrepareDesignMatrix(data=ddgo,
                    use.nuisance=TRUE,
                    response.name=c("os_diag_years","os_status"),
                    feature.category.name = c("clin","wbc","age","clinicalother","cytoaberrations"),
                    list.features=list.features)
res.CLIN_DEMO_CYTO <- c()
for(predictor in predictors[6:7]){
    tmp <- runCV(mypredictor=predictor,
          response=prep.CLIN_DEMO_CYTO$response, design=prep.CLIN_DEMO_CYTO$design,
          nfolds=nfolds, nrepeats=nrepeats, seed=seed, mc.cores=mc.cores)
    res.CLIN_DEMO_CYTO <- cbind(res.CLIN_DEMO_CYTO,tmp)
}
colnames(res.CLIN_DEMO_CYTO) <- c('Ridge','Elastic_net','Lasso','CPSS','RFX','RFS','Boost')[6:7]
write.table(res.CLIN_DEMO_CYTO,"predictive_algo_CLIN_DEMO_CYTO_v2.tsv",quote=F,sep='\t')

# 2) CLIN DEMO CYTO GENE
prep.CLIN_DEMO_CYTO_GENE = PrepareDesignMatrix(data=ddgo,
                    use.nuisance=TRUE,
                    response.name=c("os_diag_years","os_status"),
                    feature.category.name = c("clin","wbc","age","clinicalother","cytoaberrations","mutations"),
                    list.features=list.features)
res.CLIN_DEMO_CYTO_GENE <- c()
for(predictor in predictors[6:7]){
    tmp <- runCV(mypredictor=predictor,
          response=prep.CLIN_DEMO_CYTO_GENE$response, design=prep.CLIN_DEMO_CYTO_GENE$design,
          nfolds=nfolds, nrepeats=nrepeats, seed=seed, mc.cores=mc.cores)
    res.CLIN_DEMO_CYTO_GENE <- cbind(res.CLIN_DEMO_CYTO_GENE,tmp)
}
colnames(res.CLIN_DEMO_CYTO_GENE) <- c('Ridge','Elastic_net','Lasso','CPSS','RFX','RFS','Boost')[6:7]
write.table(res.CLIN_DEMO_CYTO_GENE,"predictive_algo_CLIN_DEMO_CYTO_GENE_v2.tsv",quote=F,sep='\t')

# 3) CLIN DEMO CYTO GENE COMP
prep.CLIN_DEMO_CYTO_GENE_COMP = PrepareDesignMatrix(data=ddgo,
                    use.nuisance=TRUE,
                    response.name=c("os_diag_years","os_status"),
                    feature.category.name = c("clin","wbc","age","clinicalother","cytoaberrations","mutations","component"),
                    list.features=list.features)
res.CLIN_DEMO_CYTO_GENE_COMP <- c()
for(predictor in predictors){
    tmp <- runCV(mypredictor=predictor,
          response=prep.CLIN_DEMO_CYTO_GENE_COMP$response, design=prep.CLIN_DEMO_CYTO_GENE_COMP$design,
          nfolds=nfolds, nrepeats=nrepeats, seed=seed, mc.cores=mc.cores)
    res.CLIN_DEMO_CYTO_GENE_COMP <- cbind(res.CLIN_DEMO_CYTO_GENE_COMP,tmp)
}
colnames(res.CLIN_DEMO_CYTO_GENE_COMP) <- c('Ridge','Elastic_net','Lasso','CPSS','RFX','RFS','Boost')
write.table(res.CLIN_DEMO_CYTO_GENE_COMP,"predictive_algo_CLIN_DEMO_CYTO_GENE_COMP_v2.tsv",quote=F,sep='\t')

# 4) CLIN DEMO CYTO GENE NMUT
prep.CLIN_DEMO_CYTO_GENE_NMUT = PrepareDesignMatrix(data=ddgo,
                    use.nuisance=TRUE,
                    response.name=c("os_diag_years","os_status"),
                    feature.category.name = c("clin","wbc","age","clinicalother","cytoaberrations","mutations","mutnumber"),
                    list.features=list.features)
res.CLIN_DEMO_CYTO_GENE_NMUT <- c()
for(predictor in predictors[6:7]){
    tmp <- runCV(mypredictor=predictor, 
          response=prep.CLIN_DEMO_CYTO_GENE_NMUT$response, design=prep.CLIN_DEMO_CYTO_GENE_NMUT$design,
          nfolds=nfolds, nrepeats=nrepeats, seed=seed, mc.cores=mc.cores)
    res.CLIN_DEMO_CYTO_GENE_NMUT <- cbind(res.CLIN_DEMO_CYTO_GENE_NMUT,tmp)
}
colnames(res.CLIN_DEMO_CYTO_GENE_NMUT) <- c('Ridge','Elastic_net','Lasso','CPSS','RFX','RFS','Boost')[6:7]
write.table(res.CLIN_DEMO_CYTO_GENE_NMUT,"predictive_algo_CLIN_DEMO_CYTO_GENE_NMUT_v2.tsv",quote=F,sep='\t')

# 5) CLIN DEMO CYTO GENE NUMT COMP
prep.CLIN_DEMO_CYTO_GENE_NMUT_COMP = PrepareDesignMatrix(data=ddgo,
                    use.nuisance=TRUE,
                    response.name=c("os_diag_years","os_status"),
                    feature.category.name = c("clin","wbc","age","clinicalother","cytoaberrations","mutations","component"),
                    list.features=list.features)
res.CLIN_DEMO_CYTO_GENE_NMUT_COMP <- c()
for(predictor in predictors){
    tmp <- runCV(mypredictor=predictor,
          response=prep.CLIN_DEMO_CYTO_GENE_NMUT_COMP$response, design=prep.CLIN_DEMO_CYTO_GENE_NMUT_COMP$design,
          nfolds=nfolds, nrepeats=nrepeats, seed=seed, mc.cores=mc.cores)
    res.CLIN_DEMO_CYTO_GENE_NMUT_COMP <- cbind(res.CLIN_DEMO_CYTO_GENE_NMUT_COMP,tmp)
}
colnames(res.CLIN_DEMO_CYTO_GENE_NMUT_COMP) <- c('Ridge','Elastic_net','Lasso','CPSS','RFX','RFS','Boost')
write.table(res.CLIN_DEMO_CYTO_GENE_NMUT_COMP,"predictive_algo_CLIN_DEMO_CYTO_GENE_NMUT_COMP_v2.tsv",quote=F,sep='\t')


# 6) CLIN DEMO CYTO GENE GINT
prep.CLIN_DEMO_CYTO_GENE_GINT = PrepareDesignMatrix(data=ddgo,
                    use.nuisance=TRUE,
                    response.name=c("os_diag_years","os_status"),
                    feature.category.name = c("clin","wbc","age","clinicalother","cytoaberrations","mutations","geneint_double","geneint_triple"),
                    list.features=list.features)
res.CLIN_DEMO_CYTO_GENE_GINT <- c()
for(predictor in predictors[6:7]){
    tmp <- runCV(mypredictor=predictor, 
          response=prep.CLIN_DEMO_CYTO_GENE_GINT$response, design=prep.CLIN_DEMO_CYTO_GENE_GINT$design,
          nfolds=nfolds, nrepeats=nrepeats, seed=seed, mc.cores=mc.cores)
    res.CLIN_DEMO_CYTO_GENE_GINT <- cbind(res.CLIN_DEMO_CYTO_GENE_GINT,tmp)
}
colnames(res.CLIN_DEMO_CYTO_GENE_GINT) <- c('Ridge','Elastic_net','Lasso','CPSS','RFX','RFS','Boost')[6:7]
write.table(res.CLIN_DEMO_CYTO_GENE_GINT,"predictive_algo_CLIN_DEMO_CYTO_GENE_GINT_v2.tsv",quote=F,sep='\t')

# 7) CLIN DEMO CYTO GENE GINT NMUT
prep.CLIN_DEMO_CYTO_GENE_GINT_NMUT = PrepareDesignMatrix(data=ddgo,
                    use.nuisance=TRUE,
                    response.name=c("os_diag_years","os_status"),
                    feature.category.name = c("clin","wbc","age","clinicalother","cytoaberrations","mutations","geneint_double","geneint_triple","mutnumber"),
                    list.features=list.features)
res.CLIN_DEMO_CYTO_GENE_GINT_NMUT <- c()
for(predictor in predictors[6:7]){
    tmp <- runCV(mypredictor=predictor, 
          response=prep.CLIN_DEMO_CYTO_GENE_GINT_NMUT$response, design=prep.CLIN_DEMO_CYTO_GENE_GINT_NMUT$design,
          nfolds=nfolds, nrepeats=nrepeats, seed=seed, mc.cores=mc.cores)
    res.CLIN_DEMO_CYTO_GENE_GINT_NMUT <- cbind(res.CLIN_DEMO_CYTO_GENE_GINT_NMUT,tmp)
}
colnames(res.CLIN_DEMO_CYTO_GENE_GINT_NMUT) <- c('Ridge','Elastic_net','Lasso','CPSS','RFX','RFS','Boost')[6:7]
write.table(res.CLIN_DEMO_CYTO_GENE_GINT_NMUT,"predictive_algo_CLIN_DEMO_CYTO_GENE_GINT_NMUT_v2.tsv",quote=F,sep='\t')


# 8) CLIN DEMO TOPCYTO GENE NMUT 
prep.CLIN_DEMO_TOPCYTO_GENE_NMUT = PrepareDesignMatrix(data=ddgo,
                    use.nuisance=TRUE,
                    response.name=c("os_diag_years","os_status"),
                    feature.category.name = c("clin","wbc","age","clinicalother","topcyto","mutations","mutnumber"),
                    list.features=list.features)
res.CLIN_DEMO_TOPCYTO_GENE_NMUT <- c()
for(predictor in predictors){
    tmp <- runCV(mypredictor=predictor, 
          response=prep.CLIN_DEMO_TOPCYTO_GENE_NMUT$response, design=prep.CLIN_DEMO_TOPCYTO_GENE_NMUT$design,
          nfolds=nfolds, nrepeats=nrepeats, seed=seed, mc.cores=mc.cores)
    res.CLIN_DEMO_TOPCYTO_GENE_NMUT <- cbind(res.CLIN_DEMO_TOPCYTO_GENE_NMUT,tmp)
}
colnames(res.CLIN_DEMO_TOPCYTO_GENE_NMUT) <- c('Ridge','Elastic_net','Lasso','CPSS','RFX','RFS','Boost')
write.table(res.CLIN_DEMO_TOPCYTO_GENE_NMUT,"predictive_algo_CLIN_DEMO_TOPCYTO_GENE_NMUT_v2.tsv",quote=F,sep='\t')










