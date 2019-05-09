# GET READY
source("../../src/merge_df.R")
source("./scripts/my_survival.R")
library('reshape2')
library('glmnet')
library('tibble')
library('grid')
library('gridExtra')
library('crrp')

# -------------- DATA -------------- #
# CLINICAL DF
ddc = read.table("../../data/updated_dataset/df_clinical_cyto_selected_process_ippsr.tsv",stringsAsFactors=F,header=T,sep="\t")
# MUTATIONS MAF
ddmaf = read.table("../../data/updated_dataset/df_maf_driver_selected.tsv",stringsAsFactors=F,header=T,sep="\t")
# CYTO DF
ddcyto = read.table("../../data/updated_dataset/df_cyto_binary_impute.tsv",stringsAsFactors=F,header=T,sep="\t")
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
## competing risks
ddc$time<-with(ddc,ifelse(aml_status==0,os_diag_years,aml_diag_years))
ddc$fstatus<-with(ddc,ifelse(aml_status==0,2*os_status,1))
# ---------------------------------- #
# GET READY DATAFRAMES FOR PREDICTIONS
# ddc is clinical df
# ddmutation
ddmut5 = read.table("../../data/updated_dataset/df_mut_hotspot_cut5.tsv", stringsAsFactors=F, header=T, sep="\t")
# cyto with features NK and CK
ddcyto5 = read.table("../../data/updated_dataset/df_cyto_binary_impute_cut5.tsv", stringsAsFactors=F, header=T, sep="\t")
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
ddtp53 = read.table("../../data/updated_dataset/LEUKID_P53_status.tsv",sep="\t",header=T,stringsAsFactors=F)
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
rec<-rec[1:10]
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
i.competing.risks = i.not.aml[-which(tmp2[,1]<=0 | is.na(tmp2[,1]) | is.na(tmp2[,2]))]

ddgo = ddgo[i.competing.risks,]

response.name3 = c("time","fstatus")

ddgo = ddgo[,c(unlist(list.features),response.name1,response.name2,response.name3)] 

ddgo[is.na(ddgo[, "WBC"]), "WBC"] <- mean(ddgo[, "WBC"], na.rm = TRUE) # MEAN IMPUTATION FOR WBC MISSING
ddgo$CENTER <- paste0('_', ddgo$CENTER)










# ---------------------------------- #
# COMPETING RISKS FUNCTION
# ---------------------------------- #

crrp2 <- function(time, fstatus, X, failcode=1, cencode=0,
                 penalty=c("MCP", "SCAD", "LASSO"), gamma=switch(penalty, SCAD=3.7, 2.7), 
                 alpha=1, lambda.min=.001, nlambda=50, lambda, eps=.001, 
                 max.iter=1000, penalty.factor=rep(1, ncol(X)), weighted=FALSE){
  
  dpenalty <- function(beta, gamma=switch(penalty, SCAD=3.7, 2.7), lambda, penalty=c("LASSO", "SCAD", "MCP"),
                       penalty.factor=rep(1, length(beta))){
    penalty <- match.arg(penalty)
    m <- penalty.factor
    ans <- rep(NA, length(beta))
    if (penalty=="LASSO"){
      ans[abs(beta) > 0] = (abs(beta[abs(beta) > 0])*m[abs(beta)>0])^{-1}
      ans[beta==0]=0
      w <- diag(as.vector(ans))
      ww <- w*lambda
      return(ww)
    }
    else if (penalty=="SCAD"){
      ans[abs(beta) > 0]=(lambda*(abs(beta[abs(beta) > 0])<=lambda)+
                            pmax(gamma*lambda-abs(beta[abs(beta) > 0]), 0)/(gamma-1)*(abs(beta[abs(beta) > 0])>lambda))/abs(beta[abs(beta) > 0])
      ans[beta==0]=0
      return(diag(as.vector(ans)))
    }
    else if (penalty=="MCP"){
      ans[abs(beta)>0]=1*(abs(beta[abs(beta) > 0])<=gamma*lambda)*(lambda-abs(beta[abs(beta) > 0])/gamma)/abs(beta[abs(beta) > 0])
      ans[beta==0]=0
      return(diag(as.vector(ans)))
      
    }
  }
  
  ginv = function(X, tol = sqrt(.Machine$double.eps)){
    s = svd(X)
    nz = s$d > tol * s$d[1]
    if(any(nz)) s$v[,nz] %*% (t(s$u[,nz])/s$d[nz])
    else X*0
  }
  
  ##sort time
  n <- length(time)
  p <- ncol(X)
  d <- data.frame(time=time, fstatus=fstatus)
  if (!missing(X)) d$X <- as.matrix(X)
  d <- d[order(d$time),]
  ftime <- d$time
  cenind <- ifelse(d$fstatus==cencode,1,0)
  fstatus <- ifelse(d$fstatus==failcode, 1,2*(1-cenind))
  X <- d$X
  u <- do.call('survfit',list(formula=Surv(ftime,cenind)~1,data=data.frame(ftime,cenind)))  
  #### uuu is weight function
  u <- approx(c(0,u$time,max(u$time)*(1+10*.Machine$double.eps)),c(1,u$surv,0),
              xout=ftime*(1-100*.Machine$double.eps),method='constant',f=0,rule=2)
  uuu <- u$y
  #end of data preparation##################################################
  #Output time, fstatus, uuu, X
  
  ## Error checking
  if (gamma <= 1 & penalty=="MCP") stop("gamma must be greater than 1 for the MC penalty")
  if (gamma <= 2 & penalty=="SCAD") stop("gamma must be greater than 2 for the SCAD penalty")
  if (nlambda < 2) stop("nlambda must be at least 2")
  if (alpha <= 0) stop("alpha must be greater than 0; choose a small positive number instead")
  if (length(penalty.factor)!=ncol(X)) stop("penalty.factor does not match up with X")
  
  ## Set up XX, yy, lambda
  std <- .Call("standardize", X, PACKAGE="crrp")
  XX <- std[[1]]
  center <- std[[2]]
  scale <- std[[3]]
  
  if (weighted) penalty.factor <- penalty.factor*scale
  
  nz <- which(scale > 1e-6)
  if (length(nz) != ncol(XX)) XX <- XX[ ,nz, drop=FALSE]
  penalty.factor <- penalty.factor[nz]
  
  # initial value
  #set value of lambda
  if (missing(lambda)) {
    eta0 <- rep(0,n)  
    sw <- .C("scorehessian", as.double(ftime),as.integer(fstatus), as.double(XX),
             as.integer(p), as.integer(n),  as.double(uuu), as.double(eta0), 
             double(n), double(n), double(1), PACKAGE="crrp")
    score0 <- sw[[8]]
    w0 <- sw[[9]]
    r0 <- ifelse(w0==0, 0, score0/w0)
    z <- eta0+r0
    l.max <- max(t(w0*z)%*%XX)/n
    l.min <- lambda.min
    lambda=exp(seq(log(l.max),log(l.min*l.max),len=nlambda)) 
    user.lambda <- FALSE
  } else {
    nlambda <- length(lambda)
    user.lambda <- TRUE
  }
  
  ## Fit
  
  res <- .Call("cdfit_psh", XX, as.numeric(ftime), as.integer(fstatus), uuu, penalty, lambda, 
               eps, as.integer(max.iter), gamma, penalty.factor, alpha, PACKAGE="crrp")
  b <- matrix(res[[1]], p, (nlambda))
  dev <- res[[2]]
  iter <- res[[3]]  
  score <- matrix(res[[5]], n, nlambda)
  hessian <- matrix(res[[6]], n, nlambda)  
  ## Unstandardize
  beta <- b/scale  
  ## Names
  if (is.null(colnames(X))) varnames <- paste("V",1:ncol(X),sep="")
  else varnames <- colnames(X)
  colnames(beta) <-  round(lambda,digits=4)
  rownames(beta) <- varnames
  
  #calculate GCV, BIC, and standard error
  sd <- matrix(0, p, nlambda)
  bic <- gcv <- ps <- rep(0, nlambda)
  for ( l in 1:nlambda){
    #zwz <- t(XX)%*%diag(hessian[,l])%*%XX
    #pp <- dpenalty(beta=b[,l ],lambda=lambda[l], penalty=penalty, penalty.factor=penalty.factor)
    #sd[,l ] <- sqrt(diag(ginv(zwz+n*pp)%*%zwz%*%ginv(zwz+n*pp)))/scale
    
    ll <-dev[l]/-2
    bic[l] <- -2*ll+sum(1-(b[,l ]==0))*log(n)
    
    #ps[l] <- sum(diag(XX%*%ginv(zwz+n*pp)%*%t(XX)%*%diag(hessian[,l])))-sum(b[,l ] == 0)
    #gcv[l] = -ll/(n*(1-ps[l]/n)^2)
    sd[l]<-NA
    gcv[l]<-NA
  }
  
  ## Output
  val <- structure(list(beta = beta,
                        iter = iter,
                        lambda = lambda,
                        penalty = penalty,
                        gamma = gamma,
                        alpha = alpha,
                        loglik = ll,
                        GCV = gcv,
                        BIC = bic,
                        SE=sd,
                        call=sys.call()),
                   class = "crrp")
  
  val
  
}





# ---------------------------------- #
# BOOTSTRAPPING FUNCTIONS
# ---------------------------------- #

bootstrapping_OS_AML <- function(model, ddgo, n_experiment = 100, prediction = 'os', mc.cores = 1, seed =1234){
    set.seed(seed)
    
    res_bootstrap <- data.frame('feature' = character(),
                  'coef' = numeric())
    
    formula <- paste(models[[model]], collapse = '+')
    formula <- paste0('~ 0 +', formula)
    
    x <- model.matrix(as.formula(formula), ddgo[, models[[model]]])
    if(prediction=='os'){
        y <- ddgo[, c('os_diag_years', 'os_status')]
    } else if (prediction=='aml'){
        y <- ddgo[, c("aml_diag_years","aml_status")]
    }
    
    colnames(y) <- c('time', 'status')
    y <- as.matrix(y)
    
    
    n = nrow(ddgo)
    folds <- list()
    for (i in seq(n_experiment)) {
        folds[[i]] <- sample(1:n, 0.8 * n, replace = TRUE)
    }
    nexp = length(folds)
    
    ###
    print("Start bootstrap")
    rescv = mclapply(seq(nexp),
                   FUN=function(iexp) {
                       cat(".")
                       x_sampling = x[folds[[iexp]],]
                       y_sampling = y[folds[[iexp]],]
                       
                       cvfit <- cv.glmnet(x_sampling, y_sampling, family = 'cox', alpha=1, nfolds = 20, grouped = TRUE)

                       tmp <- as.data.frame(as.matrix(coef(cvfit, s = "lambda.1se")))
                       colnames(tmp) <- 'coef'
                       tmp <- rownames_to_column(tmp, var = 'feature')
                       
                   },
                   mc.cores=mc.cores
                   )

    for(i in 1:length(rescv)){
        res_bootstrap <- rbind(res_bootstrap,rescv[[i]])
    }
    
    ###
    res_bootstrap <- res_bootstrap[res_bootstrap$coef != 0,]
    return(res_bootstrap)

}











# bootstrapping_CR <- function(model, ddgo, n_experiment = 100, mc.cores = 1, seed =1234){
#     set.seed(seed)
#     
#     res_bootstrap <- data.frame('feature' = character(),
#                   'coef' = numeric())
#     
#     formula <- paste(models[[model]], collapse = '+')
#     formula <- paste0('~ 0 +', formula)
#     
#     x <- model.matrix(as.formula(formula), ddgo[, models[[model]]])
#     y <- ddgo[, c("time","fstatus")]
#     
#     colnames(y) <- c('time', 'status')
#     y <- as.matrix(y)
#     
#     # 1 run on full model to determine best lambda using BIC
#     
#     fit=crrp2(y[,1],
#               y[,2],
#               x,
#               penalty = 'LASSO',
#               nlambda=50,
#               max.iter=1000)
#     lambda = fit$lambda[which.min(fit$BIC)]
#     
#     
#     n = nrow(ddgo)
#     folds <- list()
#     for (i in seq(n_experiment)) {
#         folds[[i]] <- sample(1:n, 0.8 * n, replace = TRUE)
#     }
#     nexp = length(folds)
#     
#     ###
#     print("Start bootstrap")
#     rescv = mclapply(seq(nexp),
#                    FUN=function(iexp) {
#                        cat(".")
#                        x_sampling = x[folds[[iexp]],]
#                        y_sampling = y[folds[[iexp]],]
#                        
#                        cvfit <- crrp2(y_sampling[,1],
#                                       y_sampling[,2],
#                                       x_sampling,
#                                       penalty = 'LASSO',
#                                       nlambda=20,
#                                       lambda=c(lambda,0),   #need at least 2 lambda
#                                       max.iter=100)
#                        
#                        tmp <- as.data.frame(cvfit$beta[,1])
#                        colnames(tmp) <- 'coef'
#                        tmp <- rownames_to_column(tmp, var = 'feature')
#                        
#                    },
#                    mc.cores=mc.cores
#                    )
# 
#     for(i in 1:length(rescv)){
#         res_bootstrap <- rbind(res_bootstrap,rescv[[i]])
#     }
#     
#     ###
#     res_bootstrap <- res_bootstrap[res_bootstrap$coef != 0,]
#     return(res_bootstrap)
# 
# }

bootstrapping_CR <- function(model, ddgo, n_experiment = 100, mc.cores = 1, seed =1234){
  set.seed(seed)
  
  res_bootstrap <- data.frame('feature' = character(),
                              'coef' = numeric())
  
  formula <- paste(models[[model]], collapse = '+')
  formula <- paste0('~ 0 +', formula)
  
  x <- model.matrix(as.formula(formula), ddgo[, models[[model]]])
  y <- ddgo[, c("time","fstatus")]
  
  colnames(y) <- c('time', 'status')
  y <- as.matrix(y)
  
  
  n = nrow(ddgo)
  folds <- list()
  for (i in seq(n_experiment)) {
    folds[[i]] <- sample(1:n, 0.8 * n, replace = TRUE)
  }
  nexp = length(folds)
  
  ###
  print("Start bootstrap")
  rescv = mclapply(seq(nexp),
                   FUN=function(iexp) {
                     cat(".")
                     x_sampling = x[folds[[iexp]],]
                     y_sampling = y[folds[[iexp]],]
                     
                     cvfit <- crrp2(y_sampling[,1],
                                    y_sampling[,2],
                                    x_sampling,
                                    penalty = 'LASSO',
                                    nlambda=10,
                                    max.iter=100)
                     
                     
                     tmp <- as.data.frame(cvfit$beta[,which.min(cvfit$BIC)])
                     colnames(tmp) <- 'coef'
                     tmp <- rownames_to_column(tmp, var = 'feature')
                     
                   },
                   mc.cores=mc.cores
  )
  
  for(i in 1:length(rescv)){
    res_bootstrap <- rbind(res_bootstrap,rescv[[i]])
  }
  
  ###
  res_bootstrap <- res_bootstrap[res_bootstrap$coef != 0,]
  return(res_bootstrap)
  
}


# ---------------------------------- #
# ANALYSIS
# ---------------------------------- #
n_experiment=100
seed=1234
mc.cores=20



# ---------------------------------- #
# MODELS
# ---------------------------------- #

models <- list()



models[['CLIN_DEMO_CYTO']] <- c(nuisance,clinicalother,wbc,age,clin,cytoaberrations)
models[['CLIN_DEMO_CYTO_GENE']] <- c(nuisance,clinicalother,wbc,age,clin,mutations,cytoaberrations)





# ---------------------------------- #
# MAIN
# ---------------------------------- #




# 1) CLIN_DEMO_CYTO
CLIN_DEMO_CYTO_OS <- bootstrapping_OS_AML('CLIN_DEMO_CYTO',
                                          ddgo, 
                                          n_experiment = n_experiment, 
                                          prediction = 'os', 
                                          mc.cores = mc.cores, 
                                          seed =seed)
CLIN_DEMO_CYTO_OS$pred <- 'OS'

CLIN_DEMO_CYTO_AML <- bootstrapping_OS_AML('CLIN_DEMO_CYTO',
                                          ddgo, 
                                          n_experiment = n_experiment, 
                                          prediction = 'aml', 
                                          mc.cores = mc.cores, 
                                          seed =seed)
CLIN_DEMO_CYTO_AML$pred <- 'AML'

CLIN_DEMO_CYTO_CR <- bootstrapping_CR('CLIN_DEMO_CYTO',
                                          ddgo, 
                                          n_experiment = n_experiment,
                                          mc.cores = mc.cores, 
                                          seed =seed)
CLIN_DEMO_CYTO_CR$pred <- 'Competing_risks'

CLIN_DEMO_CYTO<-rbind(CLIN_DEMO_CYTO_OS,CLIN_DEMO_CYTO_AML,CLIN_DEMO_CYTO_CR)
write.table(CLIN_DEMO_CYTO,"CLIN_DEMO_CYTO_bootstrapping_OS_AML_CR.tsv",quote=F,sep='\t')



# 2) CLIN_DEMO_CYTO_GENE
CLIN_DEMO_CYTO_GENE_OS <- bootstrapping_OS_AML('CLIN_DEMO_CYTO_GENE',
                                          ddgo, 
                                          n_experiment = n_experiment, 
                                          prediction = 'os', 
                                          mc.cores = mc.cores, 
                                          seed =seed)
CLIN_DEMO_CYTO_GENE_OS$pred <- 'OS'

CLIN_DEMO_CYTO_GENE_AML <- bootstrapping_OS_AML('CLIN_DEMO_CYTO_GENE',
                                          ddgo, 
                                          n_experiment = n_experiment, 
                                          prediction = 'aml', 
                                          mc.cores = mc.cores, 
                                          seed =seed)
CLIN_DEMO_CYTO_GENE_AML$pred <- 'AML'

CLIN_DEMO_CYTO_GENE_CR <- bootstrapping_CR('CLIN_DEMO_CYTO_GENE',
                                          ddgo, 
                                          n_experiment = n_experiment,
                                          mc.cores = mc.cores, 
                                          seed =seed)
CLIN_DEMO_CYTO_GENE_CR$pred <- 'Competing_risks'

CLIN_DEMO_CYTO_GENE<-rbind(CLIN_DEMO_CYTO_GENE_OS,CLIN_DEMO_CYTO_GENE_AML,CLIN_DEMO_CYTO_GENE_CR)
write.table(CLIN_DEMO_CYTO_GENE,"CLIN_DEMO_CYTO_GENE_bootstrapping_OS_AML_CR.tsv",quote=F,sep='\t')
