library(reshape2)
library(DescTools)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(survival)
library(survminer)
library(stats)

# --------------------- #
# DF
# --------------------- #
extract_df_hotspot <- function(gene="SF3B1",
                               type="aa",aa.range1=c(700,700),aa.range2=c(600,666),
                               dd_maf, dd_clinical) {
  
  if (type=="aa") {
    ia1 = which(dd_maf$GENE==gene & dd_maf$aa_number>=aa.range1[1] & dd_maf$aa_number<=aa.range1[2])
    name1 = paste(gene,paste(unique(aa.range1),collapse="_"),sep="_")
    dd_maf$GENE[ia1] = name1
    if (identical(aa.range2,"other")) {
      itmp = which(dd_maf$GENE==gene)
      ia2 = itmp[!itmp %in% ia1]
      name2 = paste(gene,"other",sep="_")
      dd_maf$GENE[ia2] = name2
    } else {
      ia2 = which(dd_maf$GENE==gene & dd_maf$aa_number>=aa.range2[1] & dd_maf$aa_number<=aa.range2[2])
      name2 = paste(gene,paste(unique(aa.range2),collapse="_"),sep="_")
      dd_maf$GENE[ia2] = name2
    }
  }
  if (type=="effect") {
    itrunc = which(dd_maf$GENE==gene & dd_maf$IS_TRUNCATED=="truncated")
    inottrunc = which(dd_maf$GENE==gene & dd_maf$IS_TRUNCATED=="not-truncated")
    name1 = paste(gene,"truncated",sep="_")
    name2 = paste(gene,"missense",sep="_")
    dd_maf$GENE[itrunc] = name1       
    dd_maf$GENE[inottrunc] = name2
  }
  ddh = merge_clinical_mutation(dd_clinical=dd_clinical, dd_maf=dd_maf)$ddmut
  return(ddh[,c(name1,name2)])
}

# --------------------- #
# BARPLOT AND VOLCANO
# --------------------- #
barplot_volcano_hotspot <- function(dd_hotspot,ddalt,gene="SF3B1",loose.cutoff=10,qval_cutoff=0.05,...) {
  
  ddtmp = ddalt
  i = apply(ddtmp, 2, function(x) sum(((x+ddalt[,gene])==2))>=loose.cutoff)
  ddtmp = ddtmp[,i]
  ddtmp[,gene] = NA
  ddtmp[dd_hotspot[,1]==1,gene] = names(dd_hotspot)[1]
  ddtmp[dd_hotspot[,2]==1,gene] = names(dd_hotspot)[2]
  mygtmp = colnames(ddtmp)[!colnames(ddtmp)%in%gene]
  
  ddres = fisher_gene_group(mydata=ddtmp, mygenes=mygtmp, group=gene, mylevels=NULL, method.adjust="BH")
  
  ddres_print = ddres[ddres$qval<=qval_cutoff,c(1,2,3,4,5,7,10)]
  rownames(ddres_print) = NULL
  
  fpb = fisher_plot_bar(ddres=ddres, qval_cutoff=qval_cutoff, cnum=loose.cutoff,...)$plot
  
  fpv = fisher_plot_volcano(ddres,
                            baselog=10,
                            label_qval_cutoff=qval_cutoff
  )$volcano
  
  if (nrow(ddres_print)==0) { ddres_print = ddres[sort(ddres$qval,index.r=T)$ix[1:3],c(1,2,3,4,5,7,10)] }
  
  gplots = arrangeGrob(tableGrob(ddres_print),arrangeGrob(fpb,fpv,ncol=2),heights=c(1,2))
  
  return(list(ddres=ddres,ddresprint=ddres_print,barplot=fpb,volcano=fpv,plots=gplots))
}

# --------------------- #
# ODD RATIO PANEL PLOT
# --------------------- #
comutation_hotspot <- function(dd_hotspot, dd_top, OR.range=c(-3,3), method="BH") {
  myhotspots = colnames(dd_hotspot)
  nhotspots = as.integer(apply(dd_hotspot,2,sum))
  topgenes = colnames(dd_top)
  labels = paste(myhotspots,"\n",paste0("(N=",nhotspots,")"))
  num = 2
  # odds per 2 hotspots with all topgenes
  # test of homogeneity of odds
  resh = lapply(topgenes, function(gg) {
    if (!grepl(gg,myhotspots[1])) {
      # tables
      tt_1 = table(dd_hotspot[,myhotspots[1]],dd_top[,gg])
      tt_2 = table(dd_hotspot[,myhotspots[2]],dd_top[,gg])
      # fisher test
      f1 = fisher.test(tt_1)
      f2 = fisher.test(tt_2)
      vec.odds = c(f1$estimate, f2$estimate)
      vec.pval = c(f1$p.val, f2$p.val)
      # non-significant odds
      vec.odds.simp = vec.odds
      vec.odds.simp[vec.pval>0.05] = 1
      vec.odds.simp[vec.odds.simp<10^OR.range[1]] = 10^OR.range[1]
      vec.odds.simp[vec.odds.simp>10^OR.range[2]] = 10^OR.range[2]
      # homogeneity test
      aa = array(c(tt_1,tt_2),dim=c(2,2,2))
      bdt = BreslowDayTest(aa)
      bdt.pval = bdt$p.val
      # result
      dr = data.frame(hotspots=myhotspots,labels=labels,gene=c(gg,gg),odds=vec.odds,odds.simp=vec.odds.simp,f.pval=vec.pval, bdt.pval=c(bdt.pval,bdt.pval))
      return(dr)
    } else {
      return(data.frame(hotspots=myhotspots,labels=labels,gene=c(gg,gg),odds=c(NA,NA),odds.simp=c(NA,NA),f.pval=c(NA,NA), bdt.pval=c(NA,NA)))
    } } )
  # dataframe of results
  ddres = do.call("rbind",resh)
  ddres$bdt.qval = rep(p.adjust(ddres$bdt.pval[seq(1,nrow(ddres),num)],method=method),each=num)
  # for plot
  label.OR = as.character(10^(OR.range[1]:OR.range[2]))
  #cc = cut(log10(ddres$odds.simp), breaks=c(OR.range[1]:OR.range[2]),labels=FALSE,include.lowest=TRUE)
  # cut with a specific interval around odds of 1 for low coloring
  cc = cut(log10(ddres$odds.simp), breaks=c(OR.range[1]:0-.Machine$double.eps,0:OR.range[2]),labels=FALSE,include.lowest=TRUE)
  ddres$cut = cc
  ddres$myx = rep(1:(nrow(ddres)/num),each=num)
  # start plot
  gg = ggplot(ddres) + geom_tile(aes(x=gene,y=labels,fill=cut), colour="grey20") + theme_minimal() +
    scale_fill_gradientn(colours=brewer.pal(length(label.OR),"RdYlBu"),
                         na.value = "grey80",
                         name="Odds ratio",
                         labels=label.OR,
                         breaks=1:length(label.OR),
                         limits=c(1,length(label.OR))) + noxtitle + noytitle + bold15 + angle45
  # add 0.01 significance
  i001 = which(ddres$bdt.qval<0.01)
  if (length(i001)>0) {
    dd001 = ddres[i001[seq(1,length(i001),num)],,drop=F]
    gg = gg + geom_point(data=dd001,aes(x=gene,y=1.5),shape=8,size=5,colour="black") + geom_rect(data=dd001,aes(xmin=myx-0.5,xmax=myx+0.5,ymin=0.5,ymax=2.5),fill=NA,colour="black",size=0.8)
  }
  # add 0.05 significance
  i005 = which(ddres$bdt.qval<0.05 & ddres$bdt.qval>0.01)
  if (length(i005)>0) {
    dd005 = ddres[i005[seq(1,length(i005),num)],,drop=F]
    gg = gg + geom_point(data=dd005,aes(x=gene,y=1.5),shape="+",size=5,colour="black") + geom_rect(data=dd005,aes(xmin=myx-0.5,xmax=myx+0.5,ymin=0.5,ymax=2.5),fill=NA,colour="black",size=0.5)
  }
  ggnude = gg + noxlabel + noleg
  return(list(ddres=ddres, gg=gg, ggnude=ggnude))
}

# --------------------- #
# ALL COMUT PLOTS
# --------------------- #
barplot_volcano_oddr_hotspot <- function(dd_hotspot,ddalt,ddtop,gene="SF3B1",losse.cutoff=20,qval_cutoff=0.05,...) {
  
  res1 = barplot_volcano_hotspot(dd_hotspot,ddalt=ddalt,gene=gene,losse.cutoff=losse.cutoff,qval_cutoff=qval_cutoff,...)
  
  p2 = comutation_hotspot(dd_hotspot=dd_hotspot, dd_top=ddtop, OR.range=c(-3,3), method="BH")$gg
  
  plots = arrangeGrob(tableGrob(res1$ddresprint),
                      arrangeGrob(res1$barplot,res1$volcano,ncol=2),
                      p2,
                      heights=c(1,2,1))
  
  return(plots)
}

# --------------------- #
# CLINICAL CORRELATES
# --------------------- #
clinical_hotspot <- function(dd_hotspot, dd_mutation, dd_clinical, mycol=c("#bdbdbd","#af8dc3","#7fbf7b"), myclinical=c("AGE","HB","PLT","WBC","ANC","MONOCYTES","PB_BLAST","BM_BLAST","RINGED_SIDEROBLASTS"),myg=NULL) {
  
  # FEATURES
  if (is.null(myg)) {
    myg = strsplit(colnames(dd_hotspot)[1],split="_")[[1]][1]
  }
  iwt = which(dd_mutation[,myg]==0)
  samples.wt = rownames(dd_mutation)[iwt]
  ih1 = which(dd_hotspot[,1]==1)
  samples.h1 = rownames(dd_hotspot)[ih1]
  ih2 = which(dd_hotspot[,2]==1)
  samples.h2 = rownames(dd_hotspot)[ih2]
  
  dd_clinical$cat = NA
  dd_clinical$cat[dd_clinical$LEUKID%in%samples.wt] = "WT"
  dd_clinical$cat[dd_clinical$LEUKID%in%samples.h1] = colnames(dd_hotspot)[1]
  dd_clinical$cat[dd_clinical$LEUKID%in%samples.h2] = colnames(dd_hotspot)[2]
  dd_clinical$AGE = dd_clinical$AGE_AT_SAMPLE_TIME
  dd_clinical = dd_clinical[!is.na(dd_clinical$cat),]
  
  ddm = melt(dd_clinical[,c(myclinical,"LEUKID","cat")], id.vars=c("LEUKID","cat"))
  ddm$cat = factor(ddm$cat, levels=c("WT",colnames(dd_hotspot)[1],colnames(dd_hotspot)[2]))
  
  my_comparisons = list(levels(ddm$cat)[2:3])
  
  gg1 = ggplot(ddm) + geom_boxplot(aes(x=cat,y=value,fill=cat),notch=T) + facet_wrap(~variable,scales="free") + theme1 + scale_y_sqrt() + topleg + scale_fill_manual(values=mycol) + nolegtitle
  gg2 = ggplot(ddm,aes(x=cat,y=value,fill=cat)) + geom_boxplot() +  stat_compare_means(comparisons=list(c(2, 3))) + facet_wrap(~variable,scales="free") + theme1 + topleg + scale_fill_manual(values=mycol) + nolegtitle + scale_y_sqrt()
  gg3 = ggplot(ddm) + geom_violin(aes(x=cat,y=value,fill=cat)) + facet_wrap(~variable,scales="free") + theme1 + scale_y_sqrt() + topleg + scale_fill_manual(values=mycol) + nolegtitle
  
  return(list(ddmelt=ddm, gg1=gg1, gg2=gg2, gg3=gg3))
  
}

# --------------------- #
# OUTCOME
# --------------------- #
outcome_hotspot <- function(dd_hotspot, dd_mutation, dd_clinical, mycol=c("#bdbdbd","#af8dc3","#7fbf7b"),conf.int=FALSE,myg=NULL) {
  
  hhc = dd_clinical
  
  # Ready Df
  if (is.null(myg)) {
    myg = strsplit(colnames(dd_hotspot)[1],split="_")[[1]][1]
  }
  iwt = which(dd_mutation[,myg]==0)
  samples.wt = rownames(dd_mutation)[iwt]
  ih1 = which(dd_hotspot[,1]==1)
  samples.h1 = rownames(dd_hotspot)[ih1]
  ih2 = which(dd_hotspot[,2]==1)
  samples.h2 = rownames(dd_hotspot)[ih2]
  hhc$cat = NA
  hhc$cat[hhc$LEUKID%in%samples.wt] = "WT"
  hhc$cat[hhc$LEUKID%in%samples.h1] = colnames(dd_hotspot)[1]
  hhc$cat[hhc$LEUKID%in%samples.h2] = colnames(dd_hotspot)[2]
  hhc = hhc[!is.na(hhc$cat),]
  hhc$cat = factor(hhc$cat, levels=c("WT",colnames(dd_hotspot)[1],colnames(dd_hotspot)[2])) 
  
  # Overall Survival
  ff = as.formula(Surv(os_diag_years,os_status)~cat)
  kmfit = survfit(Surv(os_diag_years,os_status)~cat,data=hhc)
  myleg = paste0(gsub("cat=","",names(kmfit$strata))," (N=",kmfit$strata,")")
  ggs = ggsurvplot(kmfit, data = hhc, legend.title="",legend.labs=myleg,conf.int=conf.int)$plot + scale_color_manual(values=mycol) + xlab("Years")
  #print(ggs)
  value_1_name = colnames(dd_hotspot)[1]
  value_2_name = colnames(dd_hotspot)[2]
  pt = pairwise_survdiff(formula=ff, data=hhc[hhc[,"cat"] %in% c(value_1_name, value_2_name),], p.adjust.method = "BH")
  ggOS = ggs + annotate("text",label=paste0("p-value ",value_1_name," vs ",value_2_name," =  ",signif(pt$p.value[1,1],3)),x=20,y=1)
  # AML incidence
  ff = as.formula(Surv(aml_diag_years,aml_status)~cat)
  kmfit = survfit(Surv(aml_diag_years,aml_status)~cat,data=hhc)
  myleg = paste0(gsub("cat=","",names(kmfit$strata))," (N=",kmfit$strata,")")
  ggs = ggsurvplot(kmfit, data = hhc, fun="event", legend.title="",legend.labs=myleg,conf.int=conf.int)$plot + scale_color_manual(values=mycol) + xlab("Years")
  pt = pairwise_survdiff(formula=ff, data=hhc[hhc[,"cat"] %in% c(value_1_name, value_2_name),], p.adjust.method = "BH")
  ggAML = ggs + annotate("text",label=paste0("p-value ",value_1_name," vs ",value_2_name," = ",signif(pt$p.value[1,1],3)),x=20,y=1)
  
  gg = arrangeGrob(ggOS,ggAML,ncol=2)
  
  return(list(kmfit=kmfit,gg=gg,data=hhc))
  
  
}