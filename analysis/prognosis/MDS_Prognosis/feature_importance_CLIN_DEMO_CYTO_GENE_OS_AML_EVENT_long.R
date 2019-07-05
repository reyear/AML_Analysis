# GET READY
source("../../../../src/merge_df.R")
source("../../../../src/my_survival.R")
source("../../../../src/feature_importance.R")
library('reshape2')
library('tibble')
library('grid')
library('gridExtra')

# -------------- DATA -------------- #
# CLINICAL DF
ddc = read.table("../../../../data/updated_dataset/df_clinical_cyto_selected_process_ippsr.tsv",stringsAsFactors=F,header=T,sep="\t")
# MUTATIONS MAF
ddmaf = read.table("../../../../data/updated_dataset/df_maf_driver_selected.tsv",stringsAsFactors=F,header=T,sep="\t")
# CYTO DF
ddcyto = read.table("../../../../data/updated_dataset/df_cyto_binary_impute.tsv",stringsAsFactors=F,header=T,sep="\t")
# COMPONENTS DF
#ddcomponents = read.table("../../data/updated_dataset/df_components.tsv",stringsAsFactors=F,header=T,sep="\t")
# GENE INTERACTION DF
#ddint_double = read.table("../../data/updated_dataset/datagene_interactions_duplicate_cut1_v2.tsv",stringsAsFactors=F,header=T,sep="\t")
#ddint_triple = read.table("../../data/updated_dataset/datagene_interactions_triplicate_cut2_v2.tsv",stringsAsFactors=F,header=T,sep="\t")
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
## event free survival
ddc$time<-with(ddc,ifelse(aml_status==0,os_diag_years,aml_diag_years))
ddc$status<-with(ddc,ifelse(aml_status==0,os_status,1))
# ---------------------------------- #
# GET READY DATAFRAMES FOR PREDICTIONS
# ddc is clinical df
# ddmutation
ddmut5 = read.table("../../../../data/updated_dataset/df_mut_hotspot_cut5.tsv", stringsAsFactors=F, header=T, sep="\t")
# cyto with features NK and CK
ddcyto5 = read.table("../../../../data/updated_dataset/df_cyto_binary_impute_cut5.tsv", stringsAsFactors=F, header=T, sep="\t")
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
ddtp53 = read.table("../../../../data/updated_dataset/LEUKID_P53_status.tsv",sep="\t",header=T,stringsAsFactors=F)
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
rec<-rec[1:20]
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
#geneint_double=colnames(ddint_double)
#geneint_triple=colnames(ddint_triple)
clin=c("HB","PLT","ANC","BM_BLAST")
#component=c("component")
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
                  #geneint_double=geneint_double,
                  #geneint_triple=geneint_triple,
                  clin=clin,
                  #component=component,
                  topcyto=topcyto,
                  blast=blast)

# ---------------------------------- #
# COMPLETE CASE FOR OS
# ---------------------------------- #
# COMPLETE CASE WITH OS AND AML RESPONSE 
i.complete = which(ddc$ipssr_complete_case=="complete")
response.name1 = c("os_diag_years","os_status")
tmp = ddc[i.complete,response.name1]
i.survival = i.complete[-which(tmp[,1]<=0 | is.na(tmp[,1]) | is.na(tmp[,2]))]
#ddcomponents = ddcomponents[as.character(ddc$LEUKID),]
ddgo = cbind(ddc,ddmut5,ddcyto5)
#ddgo$component=ddcomponents
ddgo = ddgo[i.survival,]
#ddgo = cbind(ddgo,ddint_double,ddint_triple)


response.name2 = c("aml_diag_years","aml_status")
i.not.aml = which(ddgo$WHO_2016_SIMPLIFY_2!="AML")
tmp2 = ddgo[i.not.aml,response.name2]
i.event.free = i.not.aml[-which(tmp2[,1]<=0 | is.na(tmp2[,1]) | is.na(tmp2[,2]))]

ddgo = ddgo[i.event.free,]

response.name3 = c("time","status")

ddgo = ddgo[,c(unlist(list.features),response.name1,response.name2,response.name3)] 

ddgo[is.na(ddgo[, "WBC"]), "WBC"] <- mean(ddgo[, "WBC"], na.rm = TRUE) # MEAN IMPUTATION FOR WBC MISSING
ddgo$CENTER <- paste0('_', ddgo$CENTER)









# ---------------------------------- #
# ANALYSIS
# ---------------------------------- #

nrepeats=5
seed=1234
mc.cores=30
npermutations=20
nfolds=5


# ---------------------------------- #
# MAIN
# ---------------------------------- #

algorithms<-c(algo_Lasso, algo_Ridge, algo_Elastic_net, algo_CPSS, algo_RFX, algo_RFS, algo_BOOST, algo_Cox)
predictors<-c(predictor_Lasso, predictor_Ridge, predictor_Elastic_net, predictor_CPSS, predictor_RFX, predictor_RFS, predictor_BOOST, predictor_Cox)
algo_names<-c('Lasso','Ridge','Elastic_net','CPSS','RFX','RFS','Boost','Cox')


# 1) CLIN DEMO CYTO GENE FOR OS
prep.CLIN_DEMO_CYTO_GENE = PrepareDesignMatrix(data=ddgo,
                    use.nuisance=TRUE,
                    response.name=c("os_diag_years","os_status"),
                    feature.category.name = c("clin","wbc","age","clinicalother","mutations","cytoaberrations"),
                    list.features=list.features)

res.CLIN_DEMO_CYTO_GENE <- data.frame('feature'=character(),
                                 'ref_CI'=numeric(),
                                 'permuted_CI'=numeric(),
                                 'algo'=character())

for(i in c(1,5,6)){
    tmp <- runCV_CI_with_test(response=prep.CLIN_DEMO_CYTO_GENE$response, design=prep.CLIN_DEMO_CYTO_GENE$design,
          nfolds=nfolds, nrepeats=nrepeats, seed=seed, mc.cores=mc.cores, features=colnames(prep.CLIN_DEMO_CYTO_GENE$design), npermutations=npermutations, 
                              algorithm=algorithms[i][[1]], predictor=predictors[i][[1]])
    tmp$algo<-algo_names[i]
    res.CLIN_DEMO_CYTO_GENE <- rbind(res.CLIN_DEMO_CYTO_GENE,tmp)
}

write.table(res.CLIN_DEMO_CYTO_GENE,"../../../../data/updated_dataset/feature_importance/feature_importance_CLIN_DEMO_CYTO_GENE_OS_long.tsv",quote=F,sep='\t')



# 2) CLIN DEMO CYTO GENE NMUT FOR AML TRANSFORMATION
prep.CLIN_DEMO_CYTO_GENE = PrepareDesignMatrix(data=ddgo,
                    use.nuisance=TRUE,
                    response.name=c("aml_diag_years","aml_status"),
                    feature.category.name = c("clin","wbc","age","clinicalother","mutations","cytoaberrations"),
                    list.features=list.features)

res.CLIN_DEMO_CYTO_GENE <- data.frame('feature'=character(),
                                 'ref_CI'=numeric(),
                                 'permuted_CI'=numeric(),
                                 'algo'=character())

for(i in c(1,5,6)){
    tmp <- runCV_CI_with_test(response=prep.CLIN_DEMO_CYTO_GENE$response, design=prep.CLIN_DEMO_CYTO_GENE$design,
          nfolds=nfolds, nrepeats=nrepeats, seed=seed, mc.cores=mc.cores, features=colnames(prep.CLIN_DEMO_CYTO_GENE$design), npermutations=npermutations, 
                              algorithm=algorithms[i][[1]], predictor=predictors[i][[1]])
    tmp$algo<-algo_names[i]
    res.CLIN_DEMO_CYTO_GENE <- rbind(res.CLIN_DEMO_CYTO_GENE,tmp)
}

write.table(res.CLIN_DEMO_CYTO_GENE,"../../../../data/updated_dataset/feature_importance/feature_importance_CLIN_DEMO_CYTO_GENE_AML_long.tsv",quote=F,sep='\t')





# 3) CLIN DEMO CYTO GENE FOR EVENT FREE
prep.CLIN_DEMO_CYTO_GENE = PrepareDesignMatrix(data=ddgo,
                    use.nuisance=TRUE,
                    response.name=c("time","status"),
                    feature.category.name = c("clin","wbc","age","clinicalother","mutations","cytoaberrations"),
                    list.features=list.features)

res.CLIN_DEMO_CYTO_GENE <- data.frame('feature'=character(),
                                 'ref_CI'=numeric(),
                                 'permuted_CI'=numeric(),
                                 'algo'=character())

for(i in c(1,5,6)){
    tmp <- runCV_CI_with_test(response=prep.CLIN_DEMO_CYTO_GENE$response, design=prep.CLIN_DEMO_CYTO_GENE$design,
          nfolds=nfolds, nrepeats=nrepeats, seed=seed, mc.cores=mc.cores, features=colnames(prep.CLIN_DEMO_CYTO_GENE$design), npermutations=npermutations, 
                              algorithm=algorithms[i][[1]], predictor=predictors[i][[1]])
    tmp$algo<-algo_names[i]
    res.CLIN_DEMO_CYTO_GENE <- rbind(res.CLIN_DEMO_CYTO_GENE,tmp)
}

write.table(res.CLIN_DEMO_CYTO_GENE,"../../../../data/updated_dataset/feature_importance/feature_importance_CLIN_DEMO_CYTO_GENE_EVENT_FREE_long.tsv",quote=F,sep='\t')


