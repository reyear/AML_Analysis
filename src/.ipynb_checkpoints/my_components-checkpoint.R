ComponentEventFreq <- function(data, mycomponent="comp_10",molecular.features,mycol="orange", mycut=1, gopairs=T) {
  
  # ----- #
  data$category = "other"
  data$category[data$predicted_component==mycomponent] = mycomponent
  # ----- #
  # Barplot of events for the component
  icomp = which(data$category==mycomponent)
  tmp = data[icomp,molecular.features]
  tmp = tmp[,apply(tmp,2,sum)>0]
  tt = rev(sort(apply(tmp,2,sum)))
  # df of count and within class freq
  dtmp = data.frame(event=names(tt), count=as.vector(tt))
  dtmp$freq = 100 * dtmp$count / nrow(tmp)
  dtmp$event = factor(dtmp$event, levels=dtmp$event)
  
  gbar = ggplot(dtmp[dtmp$count>mycut,], aes(x=event,y=count)) + geom_bar(fill=mycol,alpha=.8,stat="identity") + theme1 + angle45 + noxtitle + ylab(paste("count in",mycomponent)) + scale_y_continuous(sec.axis = sec_axis(~100*./nrow(tmp), name = "frequency"))
  # return dtmp and barplot
  # ----- #
  
  # ----- #
  # Pairwise Frequency
  # ----- #
  gfreq = NA
  if (gopairs) {
    ttmp = tmp[,as.vector(dtmp$event[dtmp$count>mycut])]
    pairs = sapply(1:ncol(ttmp), function(i) colMeans(ttmp * ttmp[,i], na.rm=TRUE))
    colnames(pairs) = rownames(pairs)
    pairs = 100*pairs
    pairs = as.data.frame(pairs)
    aa = pairs
    aa$id = rownames(pairs)
    meltpairs = melt(aa,id="id")
    colnames(meltpairs) = c("x","y","frequency")
    meltpairs$x = factor(meltpairs$x, levels=aa$id)
    meltpairs$y = factor(meltpairs$y, levels=aa$id)
    gfreq = ggplot(meltpairs, aes(x=x,y=y,fill=frequency)) + geom_tile() + scale_fill_continuous(low="#deebf7", high="#e31a1c", na.value="#cccccc") + angle45 + noxtitle + noytitle
  }
  
  return(list(dcount=dtmp,gbar=gbar,gfreq=gfreq))
  
}

OverallComponentLandscape <- function (data, components, gene.features, cyto.features, maxcutfreq=3, col.component=c(brewer.pal(n=12, "Paired")[-11],"orchid3","darkcyan")) {
  
  cna.features = cyto.features[!grepl("r_",cyto.features)]
  rearr.features = cyto.features[grepl("r_",cyto.features)]
  all.features = c(gene.features,cna.features,rearr.features)
  
  thecomponents = components[!grepl("NaN",components)]
  lltmp = lapply(thecomponents, function(mycomponent) {
    htmp = ComponentEventFreq(data=data,mycomponent=mycomponent,molecular.features=all.features,gopairs=F)$dcount
    htmp$event = as.vector(htmp$event)
    htmp$component = mycomponent
    return(htmp)
  })
  datatmp = do.call('rbind',lltmp)
  datatmp$type = "Genes"
  datatmp$type[datatmp$event%in%cna.features] = "CNAs"
  datatmp$type[datatmp$event%in%rearr.features] = "Fusions"
  
  keep.event = unique(datatmp$event)[sapply(unique(datatmp$event), function(x) max(datatmp[datatmp$event==x,"freq"])>=maxcutfreq)]
  
  datatmp = datatmp[datatmp$event %in% keep.event,]
  datatmp$event = factor(datatmp$event, levels=c(all.features))
  datatmp$type = factor(datatmp$type, levels=c("Genes","CNAs","Fusions"))
  datatmp$component = factor(datatmp$component, levels=thecomponents)
  tt = table(data$predicted_component)[thecomponents]
  datatmp$label = paste(datatmp$component,"\n","N=",tt[datatmp$component])
  datatmp$label = factor(datatmp$label, levels=unique(datatmp$label))
  
  #ggplot(datatmp, aes(x=event,y=freq,fill=component)) + geom_bar(alpha=0.8,stat="identity") + facet_grid(component ~ type,scales="free_x", space="free_x") + theme_light() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + noxtitle + ylab("frequency within components") + theme(strip.text.y = element_text(angle=0,face="bold",size=10), strip.text.x = element_text(face="bold",size=10)) + scale_fill_manual(values=c(brewer.pal(n=length(thecomponents), "Set3"))) + noleg
  
  gg = ggplot(datatmp, aes(x=event,y=freq,fill=component)) + geom_bar(stat="identity") + facet_grid(label ~ type,scales="free_x", space="free_x") + theme_light() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + noxtitle + ylab("frequency within components") + theme(strip.text.y = element_text(angle=0,face="bold",size=10), strip.text.x = element_text(face="bold",size=10)) + scale_fill_manual(values=col.component) + noleg
  
  return(gg)
  
}

ComponentContinuousClinical <- function(data, mycomponent="comp_10", mycol="orange",myclinical=c("AGE","HB","PLT","WBC","ANC","MONOCYTES","PB_BLAST","BM_BLAST","RINGED_SIDEROBLASTS")) {
  
  data$category = "other"
  data$category[data$predicted_component==mycomponent] = mycomponent
  data$AGE = data$AGE_AT_SAMPLE_TIME
  hh = melt(data[,c(myclinical,"LEUKID","category")], id.vars=c("LEUKID","category"))
  hh$category = factor(hh$category, levels=c(mycomponent,"other"))
  hh$label = paste0(hh$category," (N=",table(data$category)[as.vector(hh$category)],")")
  # ggplot(hh,aes(x=category,y=value,fill=label)) + geom_boxplot(alpha=0.8) + facet_wrap(~variable,scales="free") + theme1 + topleg + scale_fill_manual(values=c(mycol,"grey")) + nolegtitle + noxtitle + stat_compare_means() + scale_y_sqrt()
  gg = ggplot(hh,aes(x=category,y=value,fill=label)) + geom_violin(alpha=0.8) + geom_boxplot(width=0.1) + facet_wrap(~variable,scales="free") + theme1 + topleg + scale_fill_manual(values=c(mycol,"grey")) + nolegtitle + noxtitle + stat_compare_means() + scale_y_sqrt()
  
  return(gg)
}

ComponentDiscreteClinical <- function(data, mycomponent="comp_10",mydiscrete="SEX",col.discrete=NULL) {
  
  data$category = "other"
  data$category[data$predicted_component==mycomponent] = mycomponent
  
  data$label = paste0(data$category,"\n","(N=",table(data$category)[as.vector(data$category)],")")
  data$value = data[,mydiscrete]
  
  gg = ggplot(data) + geom_bar(aes(x=label,fill=value),position="fill",alpha=.8) + theme1 + topleg + noxtitle + ylab("proportion") + nolegtitle
  
  if (!is.null(col.discrete)) { gg = gg + scale_fill_manual(values=col.discrete,na.value="grey") }
  
  return(gg)
}

ComponentWHO <- function(data, mycomponent="comp_10",col.who=c("#e34a33","#33a02c","#b2df8a","#fa9fb5","#b2abd2","#6a3d9a","orchid4","#1f78b4","#a6cee3","#a8ddb5")) {
  
  data$category = "other"
  data$category[data$predicted_component==mycomponent] = mycomponent
  
  data$WHO_2016_SIMPLIFY_2[which(data$WHO_2016_SIMPLIFY_2=="other")]=NA
  data$mywho = data$WHO_2016_SIMPLIFY_2
  data$mywho = factor(data$mywho, levels=c("MDS-del5q","MDS-SLD/MLD","MDS-RS-SLD/MLD","MDS-EB1/2","CMML","AML","aCML","MDS-U","MDS/MPN-RS-T","MDS/MPN-U"))
  
  data$label = paste0(data$category,"\n","(N=",table(data$category)[as.vector(data$category)],")")
  
  gg = ggplot(data) + geom_bar(aes(x=label,fill=mywho),position="fill",alpha=.8) + theme1 + noxtitle + ylab("proportion") + scale_fill_manual(values=col.who,na.value="grey") + nolegtitle
  
  return(gg)
  
}


ComponentOutcome <- function(data, mycomponent="comp_10", mycol="orange") {
  
  data$category = "other"
  data$category[data$predicted_component==mycomponent] = mycomponent
  # Overall Survival
  ff = as.formula(paste("Surv(os_diag_years,os_status)~category"))
  kmfit = survfit(ff,data=data)
  kmfit$call$formula <- ff
  myleg = paste0(gsub("category=","",names(kmfit$strata))," (N=",kmfit$strata,")")
  ggs = ggsurvplot(kmfit, data = data, legend.title="",legend.labs=myleg,conf.int=F)$plot + scale_color_manual(values=c(mycol,"grey")) + xlab("Years")
  pt = pairwise_survdiff(formula=ff, data=data, p.adjust.method = "BH")
  ggOS = ggs + annotate("text",label=paste0("p-value ",signif(pt$p.value[1,1],3)),x=15,y=1)
  # AML incidence
  ff = as.formula(paste("Surv(aml_diag_years,aml_status)~category"))
  kmfit = survfit(ff,data=data)
  kmfit$call$formula <- ff
  myleg = paste0(gsub("category=","",names(kmfit$strata))," (N=",kmfit$strata,")")
  ggs = ggsurvplot(kmfit, data = data, legend.title="",legend.labs=myleg,conf.int=F,fun="event")$plot + scale_color_manual(values=c(mycol,"grey")) + xlab("Years")
  pt = pairwise_survdiff(formula=ff, data=data, p.adjust.method = "BH")
  ggAML = ggs + annotate("text",label=paste0("p-value ",signif(pt$p.value[1,1],3)),x=15,y=1)
  
  return(list(ggOS=ggOS,ggAML=ggAML))
  
}

ComponentOutcomeStratified <- function(data,mycomponent="comp_10",col.who=c("#e34a33","#33a02c","#b2df8a","#fa9fb5","#b2abd2","#6a3d9a","orchid4","#1f78b4","#a6cee3","#a8ddb5"),col.ipssr=c("#2ca25f","#99d8c9","#bcbddc","#fdbb84","#e34a33")) {
  
  data$category = "other"
  data$category[data$predicted_component==mycomponent] = mycomponent
  data$WHO_2016_SIMPLIFY_2[which(data$WHO_2016_SIMPLIFY_2=="other")]=NA
  data$mywho = data$WHO_2016_SIMPLIFY_2
  levels.who = c("MDS-del5q","MDS-SLD/MLD","MDS-RS-SLD/MLD","MDS-EB1/2","CMML","AML","aCML","MDS-U","MDS/MPN-RS-T","MDS/MPN-U")
  data$mywho = factor(data$mywho, levels=levels.who)
  
  kmfit = survfit(Surv(os_diag_years,os_status)~mywho,data=data[data$category==mycomponent,])
  myleg = paste0(gsub("mywho=","",names(kmfit$strata))," (N=",kmfit$strata,")")
  mycolwho = col.who[levels.who%in%gsub("mywho=","",names(kmfit$strata))]
  ggwho = ggsurvplot(kmfit, data = data, legend.title="",legend.labs=myleg,conf.int=F)$plot + scale_color_manual(values=mycolwho) + xlab("Years") + guides(colour=guide_legend(ncol=2))
  
  levels.ipssr = c("VERY-GOOD","GOOD","INT","POOR","VERY-POOR")
  data$IPSSR_CALCULATED = factor(data$IPSSR_CALCULATED, levels=levels.ipssr)
  kmfit = survfit(Surv(os_diag_years,os_status)~IPSSR_CALCULATED,data=data[data$category==mycomponent,])
  myleg = paste0(gsub("IPSSR_CALCULATED=","",names(kmfit$strata))," (N=",kmfit$strata,")")
  mycolipssr = col.ipssr[levels.ipssr%in%gsub("IPSSR_CALCULATED=","",names(kmfit$strata))]
  ggipssr = ggsurvplot(kmfit, data = data, legend.title="",legend.labs=myleg,conf.int=F)$plot + scale_color_manual(values=mycolipssr) + xlab("Years") + guides(colour=guide_legend(ncol=2))
  
  return(list(ggwho=ggwho,ggipssr=ggipssr))
  
}

ComponentOutcomeAMLStratified <- function(data,mycomponent="comp_10",col.who=c("#e34a33","#33a02c","#b2df8a","#fa9fb5","#b2abd2","#6a3d9a","orchid4","#1f78b4","#a6cee3","#a8ddb5"),col.ipssr=c("#2ca25f","#99d8c9","#bcbddc","#fdbb84","#e34a33")) {
  
  data$category = "other"
  data$category[data$predicted_component==mycomponent] = mycomponent
  data$WHO_2016_SIMPLIFY_2[which(data$WHO_2016_SIMPLIFY_2=="other")]=NA
  data$mywho = data$WHO_2016_SIMPLIFY_2
  levels.who = c("MDS-del5q","MDS-SLD/MLD","MDS-RS-SLD/MLD","MDS-EB1/2","CMML","AML","aCML","MDS-U","MDS/MPN-RS-T","MDS/MPN-U")
  data$mywho = factor(data$mywho, levels=levels.who)
  
  kmfit = survfit(Surv(os_diag_years,os_status)~mywho,data=data[data$category==mycomponent,])
  myleg = paste0(gsub("mywho=","",names(kmfit$strata))," (N=",kmfit$strata,")")
  mycolwho = col.who[levels.who%in%gsub("mywho=","",names(kmfit$strata))]
  ggwho = ggsurvplot(kmfit, data = data, legend.title="",legend.labs=myleg,conf.int=F,fun="event")$plot + scale_color_manual(values=mycolwho) + xlab("Years") + guides(colour=guide_legend(ncol=2))
  
  levels.ipssr = c("VERY-GOOD","GOOD","INT","POOR","VERY-POOR")
  data$IPSSR_CALCULATED = factor(data$IPSSR_CALCULATED, levels=levels.ipssr)
  kmfit = survfit(Surv(os_diag_years,os_status)~IPSSR_CALCULATED,data=data[data$category==mycomponent,])
  myleg = paste0(gsub("IPSSR_CALCULATED=","",names(kmfit$strata))," (N=",kmfit$strata,")")
  mycolipssr = col.ipssr[levels.ipssr%in%gsub("IPSSR_CALCULATED=","",names(kmfit$strata))]
  ggipssr = ggsurvplot(kmfit, data = data, legend.title="",legend.labs=myleg,conf.int=F,fun="event")$plot + scale_color_manual(values=mycolipssr) + xlab("Years") + guides(colour=guide_legend(ncol=2))
  
  return(list(ggwho=ggwho,ggipssr=ggipssr))
  
}


MiscComponentOutcome <- function(data, mycomponent="comp_10", mycol="orange", covariate="ETNK1", addcol = c("brown","coral","darkcyan","aquamarine3")) {
  
  data$category = "other"
  data$category[data$predicted_component==mycomponent] = mycomponent
  # Overall Survival
  ff = as.formula(paste("Surv(os_diag_years,os_status)~category +",covariate))
  kmfit = survfit(ff,data=data)
  kmfit$call$formula <- ff
  myleg = paste0(gsub("category=","",names(kmfit$strata))," (N=",kmfit$strata,")")
  ggOS = ggsurvplot(kmfit, data = data, legend.title="",legend.labs=myleg,conf.int=F)$plot + scale_color_manual(values=c(addcol)) + xlab("Years") + guides(colour=guide_legend(ncol=2))
  #pt = pairwise_survdiff(formula=ff, data=data, p.adjust.method = "BH")
  #ggOS = ggs + annotate("text",label=paste0("p-value ",signif(pt$p.value[1,1],3)),x=15,y=1)
  # AML incidence
  ff = as.formula(paste("Surv(aml_diag_years,aml_status)~category +",covariate))
  kmfit = survfit(ff,data=data)
  kmfit$call$formula <- ff
  myleg = paste0(gsub("category=","",names(kmfit$strata))," (N=",kmfit$strata,")")
  ggAML = ggsurvplot(kmfit, data = data, legend.title="",legend.labs=myleg,conf.int=F,fun="event")$plot + scale_color_manual(values=addcol) + xlab("Years") + guides(colour=guide_legend(ncol=2))
  #pt = pairwise_survdiff(formula=ff, data=data, p.adjust.method = "BH")
  #ggAML = ggs + annotate("text",label=paste0("p-value ",signif(pt$p.value[1,1],3)),x=15,y=1)
  
  return(list(ggOS=ggOS,ggAML=ggAML))
  
}

MiscLenComponent6Outcome <- function(data, mycomponent="comp_6", mycol="orange") {
  
  data$category = "other"
  data$category[data$predicted_component==mycomponent] = mycomponent
  # Overall Survival
  ff = as.formula(Surv(os_diag_years,os_status)~category + CSNK1A1 + lenalidomid)
  kmfit = survfit(Surv(os_diag_years,os_status)~category + CSNK1A1 + lenalidomid,data=data)
  myleg = paste0(gsub("category=","",names(kmfit$strata))," (N=",kmfit$strata,")")
  ggOS = ggsurvplot(kmfit, data = data, legend.title="",legend.labs=myleg,conf.int=F)$plot + scale_color_manual(values=c(mycol,"orchid4","grey","pink","steelblue","green","brown","cyan")) + xlab("Years") + guides(colour=guide_legend(ncol=2))
  #pt = pairwise_survdiff(formula=ff, data=data, p.adjust.method = "BH")
  #ggOS = ggs + annotate("text",label=paste0("p-value ",signif(pt$p.value[1,1],3)),x=15,y=1)
  # AML incidence
  ff = as.formula(Surv(aml_diag_years,aml_status)~category + CSNK1A1 + lenalidomid)
  kmfit = survfit(Surv(aml_diag_years,aml_status)~category + CSNK1A1 + lenalidomid,data=data)
  myleg = paste0(gsub("category=","",names(kmfit$strata))," (N=",kmfit$strata,")")
  ggAML = ggsurvplot(kmfit, data = data, legend.title="",legend.labs=myleg,conf.int=F,fun="event")$plot + scale_color_manual(values=c(mycol,"orchid4","grey","pink","steelblue","green","brown","cyan")) + xlab("Years") + guides(colour=guide_legend(ncol=2))
  #pt = pairwise_survdiff(formula=ff, data=data, p.adjust.method = "BH")
  #ggAML = ggs + annotate("text",label=paste0("p-value ",signif(pt$p.value[1,1],3)),x=15,y=1)
  
  return(list(ggOS=ggOS,ggAML=ggAML))
  
}

