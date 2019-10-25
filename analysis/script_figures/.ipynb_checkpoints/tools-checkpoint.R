####################################### Libraries #######################################
library(doMC)
library(dplyr)
library(ggplot2)
####################################### End Libraries #######################################





####################################### Prepare the colors #######################################
df_final <- read.table("/Users/taziy/AML_analysis/analysis/clustering/clustering_Final_1/df_final_full_component.tsv")
df_w_correlates <- read.table("/Users/taziy/AML_analysis/analysis/clustering/clustering_Final_1/df_final_full_component.tsv")
df_final_itd <- read.table("/Users/taziy/AML_analysis/analysis/clustering/clustering_Final_1/df_final_full_component_ITD.tsv")
cols_component <- colnames(df_final %>% dplyr:: select(starts_with("full_component_")))
tmp <- NULL
for (co in cols_component){
    tmp1 <- df_w_correlates[df_w_correlates[,co]==1,]
    tmp1$comparison <- str_remove(co,"full_component_")
    tmp <- rbind(tmp,tmp1)
    }
tmp$fill_pal <- factor(tmp$comparison)
color_values <- c("grey45", "#e79f00", "#009E73","#0072B2",  "#CC79A7","#9ad0f3", "#D55E00", "lightgoldenrod","lightskyblue","#F0E442",
                  "firebrick3","#000000","#870C14","#a6bddb","mistyrose4","#fdbb84","gray","deeppink","darkblue","darkred","darkgreen",
                  "purple","forestgreen")
names(color_values) <- levels(factor(tmp$fill_pal))
color_values["overlap"] <- "purple"

eln <- colnames(df_final_itd[,c(2,3,4)])
comp <- colnames(df_final_itd[,c(170:193)])
age <- c(167)

gen<- colnames(df_final_itd[,c(5:88)])

cyto <- colnames(df_final_itd[,c(89:158)])       
              
clin <- colnames(df_final_itd[,c(159:165)])
demo <- colnames(df_final_itd[,c(166:167)])
mrd <- c("CR_MRD_neg","CR_MRD_pos","all_others")

colors_analysis <- c(eln="#EE9937",comp="#59A08B",gen="#BFBFBF",cyto="#2b8cbe",clin="#870C14",demo="#a6bddb",age="#a6bddb",gen_cyto="pink",eln_gen_cyto="#fdbb84",comp_gen_cyto="lightgoldenrod",all="lightsalmon",mrd="pink")

val <- c("#EE9937","#5C5C5C","#870C14","#BFBFBF","#59A08B","#2b8cbe","#a6bddb","#fdbb84","#e79f00","#000000","darkseagreen","lightskyblue","#0072B2","pink","blue","green")

####################################### End Prepare the colors #######################################





####################################### Survival Plots #######################################

plot_surv_curves <-function(fit,submain="",vals=val,legend="top",risk.tab=F,y="",linetype=1,size=4,pval=TRUE,pval.coord=c(0,0.05),font.legend=28,xlab="Time (years)",legend.title="",...){
    
    # Create the vector to use for leg.labs
    vec <- NULL
    for (i in 1:length(fit$strata)){
    vec <- cbind(vec,paste(str_remove(names(fit$strata[i]),"comparison="),fit$n[i],sep=" , N="))
    }
    
    ggsurvplot(fit,  pval = pval,main = "Survival curve",risk.table=risk.tab,submain = submain,palette=vals,legend=legend,pval.size=8,pval.coord=pval.coord,risk.table.fontsize=10,xlab=xlab,legend.labs=vec,ylab=y,legend.title=legend.title,...,
               linetype=linetype,size=size,
                  ggtheme = theme_survminer(
                 font.main = c(30, "plain", "black"),
                 font.submain = c(30, "plain", "black"),
                 font.legend=font.legend,
                 font.caption = c(30, "plain", "black"),
                 font.x = c(24, "plain", "black"),
                 font.y = c(24, "plain", "black"),
                 font.tickslab = c(26, "plain", "black")))

}

####################################### End Survival Plots #######################################



####################################### Overall Analysis #######################################

heatmaps <- function(df_final, type=c("not_ordered","ordered"),cols_to_plot=c(genes,cytos)){
    
    cols_to_keep <- c(cols_to_plot,"final_component_numeric")
    splicing <-  c("ZRSR2", "U2AF1_p.S34","U2AF1_p.Q157", "SRSF2", "SF3B1", "SF1", "NF1", "CUX1")
    chromatin <- c("ASXL1", "STAG2", "BCOR", "MLL", "EZH2", "PHF6")
    transcriptor <- c("RUNX1","SETBP1") 
    factors <- c(splicing,chromatin,transcriptor)
    features_cols <- colnames(df_final[,cols_to_plot])
    tmp_all <-df_final[,cols_to_keep]
    for (ord in type){
    #     png(paste(ord,"_heatmap.png"),width=3000,height=3500,res=200)
        if (ord=="ordered"){
            tmp_all <- df_final[,cols_to_keep][order(df_final[,cols_to_keep]$final_component_numeric),]
        }

        transp_df_all <-(as.data.frame(t(tmp_all)))
        transp_df_all$features <- rownames(transp_df_all)
        transp_df_all <- transp_df_all[! row.names(transp_df_all) == "final_component_numeric",]
        col_order <- c("t_15_17","inv_16","t_8_21","t_v_11","t_9_11","t_6_9","inv_3","add_8","add_13","add_21","add_22","add_11","TP53","complex","NPM1","CEBPA_bi","CEBPA_mono","DNMT3A","IDH1","IDH2_p.R140","IDH2_p.R172","WT1","ITD","TET2",factors,
                           "FLT3_TKD","FLT3_other",
                           "del_5","del_7","del_17","del_9","del_13","del_20","del_18","del_16","del_12","del_3","minusy")
        col_order <- c(col_order,setdiff(features_cols,col_order))
        col_order <- rev(col_order)
        transp_df_all <- transp_df_all[col_order,]
        melt.data_all<-melt(transp_df_all,id.vars ="features", variable_name=colnames(transp_df_all))
        melt.data_all$features <- factor(melt.data_all$features,levels=col_order)  # rev because otherwise the first are at the bottom


        pal <- ifelse(is.element(melt.data_all$features,genes),"#2b8cbe","#BFBFBF") 
    #     options(repr.plot.res = 200)
        col_sep <- tabulate(df_final[,cols_to_keep]$final_component_numeric)   ## to separate by component
        if(ord=="ordered"){
            cat("Order of Components :","\n")
            cat(paste(unique(df_final[order(df_final[,cols_to_keep]$final_component_numeric),]$final_component)," ,"))
        }
            # plot the heatmap
        p <- ggplot(melt.data_all, aes(variable,features)) +  
        geom_raster(aes(fill = factor(value)), show.legend = FALSE) +
        scale_fill_manual(values = c("0" = "grey45", "1" = "#e79f00")) +
        theme(axis.text.y = element_text(colour = pal,size=10,face="bold"))+
        labs(x = 'patients')+ylab("")
        if(ord=="ordered"){
            p <- p+geom_vline(xintercept = cumsum(col_sep)[c(1,2,3,4,6,7,8,9,10)] + 0.5, col = "#009E73", linetype = 1, size = 1) +
            geom_vline(xintercept = cumsum(col_sep)[c(5,11,12,13,14)] + 0.5, col = "#CC79A7", linetype = 2, size = 1) +
            geom_vline(xintercept = cumsum(col_sep)[c(15,16)] + 0.5, col = "white", linetype = 3, size = 1) 
        } else {
            p <- p+theme(axis.text.y=element_blank())
        }
        
        return(p) 
        #     dev.off()
        }
    
} 

unique_genotypes <- function(df_final,features,top=20){   
    
    
    ### features: Features to select from which features we want the genotype (example : genes and cytos)
    
    df_final <- df_final[,features]
    df_genotypes <- cbind.data.frame(df_final,overall = apply(df_final, 1, function(x)paste(colnames(df_final)[x[1:length(x)] == 1], collapse = ", ")))
                                               
    levels(df_genotypes$overall)[levels(df_genotypes$overall) == ""] <- "no events"
                                               
# png("unique_genotypes.png",width=4500,height=2500,res=200)
                                               
    df_genotypes <- data.frame(table(df_genotypes$overall))    
    df_genotypes <-df_genotypes[order(df_genotypes$Freq,decreasing=T),]                                               
    ggplot(df_genotypes[1:top,],aes(x=reorder(Var1,Freq),y=Freq))+geom_bar(stat="identity")+coord_flip()+
    theme(plot.title = element_text(hjust = 0.5,size=25),axis.text=element_text(size=22),axis.title=element_text(size=34,face="bold"),
    legend.position="none")+xlab("Top 20 Unique genotypes") + ylab("Number of Patients")
}

                                                              
overall_strat_frequency <- function(df_w_correlates,top=70){
    tmp <- df_w_correlates
    df_to_plot <- tmp[,c(genes,cytos,"eln_2017")] %>% group_by(eln_2017) %>% summarise_all(funs(sum))
    df_to_plot <- df_to_plot[,names(sort(colSums(df_to_plot), decreasing = TRUE))]
    df_to_plot <- cbind(eln = c("adverse","intermediate","favorable"),df_to_plot[,1:top])
    df_to_plot [1,2:ncol(df_to_plot)] <- 100 * df_to_plot [1,2:ncol(df_to_plot)] /nrow(tmp)
    df_to_plot [2,2:ncol(df_to_plot)] <- 100 * df_to_plot [2,2:ncol(df_to_plot)] / nrow(tmp)
    df_to_plot [3,2:ncol(df_to_plot)] <- 100 * df_to_plot [3,2:ncol(df_to_plot)] / nrow(tmp)
    gg <- melt(df_to_plot,id="eln",value.name="Frequency", variable.name="Variable")

    ggplot(gg, aes(x = Variable, y = Frequency,fill=eln)) +
        geom_bar(stat='identity')+scale_fill_manual(values=c("#EE9937","#2b8cbe","#59A08B"))+theme(axis.title.x=element_text(hjust = 0.5,size=30),axis.title.y=element_text(hjust = 0.5,size=30),
                                                                                                   axis.text.x=element_text(angle=90, hjust=1),plot.title = element_text(hjust = 0.5,size=55),
                                                                                                   axis.text=element_text(size=30,face="bold"),axis.title=element_text(size=14,face="bold"),
                                                                                                  legend.text=element_text(size=25),legend.title=element_text(size=45))+
        labs(fill="ELN 2017")
    }
####################################### End Heatmaps #######################################




####################################### Component Plots #######################################
                                                              
comp_features_frequency <- function(df_final="",genes="",cytos="",comp="ss",cols_to_keep=c(gen,cyto,"inv_3","t_6_9","t_9_11","t_15_17"),facet_type=T,colors=""){
#     png("comp_features_frequency.png",width=4500,height=5500))    


    set_notebook_plot_size(60,10*length(comp))

    columns <- c(cols_to_keep,comp)
    tmp_cat <- as.data.frame(df_final[,columns])
    categories_repartition1 <- data.frame(category = colnames(tmp_cat))

    for (i in comp)
        categories_repartition1[sprintf('%s', i)] <- apply(categories_repartition1, 1, function(s) (sum(tmp_cat[tmp_cat[,i] == 1, s['category']])/nrow(tmp_cat[tmp_cat[,i] ==1,])))
    categories_repartition1 <- categories_repartition1[!grepl(c("full|overlap"), categories_repartition1$category),]

    test <- t(categories_repartition1)
    colnames(test) = test[1, ] 
    # test = test[-1, ] 
    names <- rownames(test)
    rownames(test) <- NULL
    data <- cbind(names,test)
    
    data <- data.frame(data)
    data <- data[-1,]
    data$names <- str_replace(data$names,"full_component_","")
    data <- data[!(data$names=="no_events"),]
    melt_data <- melt(data,id="names",name="Frequency", variable.name="Variable")
    melt_data$value <- as.numeric(melt_data$value)
    melt_data$fact <- factor(melt_data$names, levels=c(data$names))
    melt_data$separation <- ifelse(melt_data$Variable %in% genes,"genes",
                                    ifelse(melt_data$Variable %in% c("t_15_17","t_v_11","t_9_11",'t_8_21','t_6_9','t_8_21',"others_transloc","inv_3","inv_16"),"fusions","CNAs"))
    melt_data$separation <- factor(melt_data$separation, levels=c("genes","CNAs","fusions"))

    p <- ggplot(melt_data, aes(x=Variable, y=value, fill=fact))+
      geom_bar(stat="identity")
    if(facet_type==T){
        p <- p + facet_grid(fact~separation,scales = "free_x", space = "free_x")
    }else{
        p <- p + facet_grid(fact~.)
    }                                                           
    
    p+theme(strip.text = element_text(face="bold", size=40,lineheight=5.0),
    strip.background = element_rect( colour="black",size=1),legend.position="none")+ylab("Frequency")+xlab("")+theme(plot.title = element_text(hjust = 0.5,size=120,face="bold"),axis.text.x = element_text(angle = 90, hjust = 1, size=ifelse(length(comp)>1,80,30),face="bold"),axis.text.y = element_text(size=50,face="bold"))+
         theme(axis.title.y = element_text(size = 40, angle = 90, vjust = 0.25,face="bold"))+
            scale_fill_manual(values=colors)+theme(legend.position = "none")+ylab("Frequency")+xlab("")

# dev.off()
}

comp_repartition <- function(df_final,cols_component,cols_order){

    df_tmp <- setNames(data.frame(matrix(ncol = length(cols_order), nrow = length(cols_order))), cols_order)
    rownames(df_tmp) <- gsub("full_component_","",cols_order)
    for (col in cols_component){
        for (col_bis in cols_component){
            df_tmp[gsub("full_component_","",col),gsub("full_component_","",col_bis)] <- ifelse(is.null(dim(df_final[df_final[,col]==1 & df_final[,col_bis]==1,])[1]),0,dim(df_final[df_final[,col]==1 & df_final[,col_bis]==1,])[1])
            }}
    
#     get_lower_tri<-function(cormat){
#     cormat[upper.tri(cormat)] <- NA
#     return(cormat)
#     }
    
    # Get upper triangle of the correlation matrix
    get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
    }
    upper_tri <- get_upper_tri(as.matrix(df_tmp))
    melted_cormat <- melt(upper_tri, na.rm = TRUE)
    melted_cormat <- cbind(melted_cormat, value_bis=0)
    for (col in unique(melted_cormat$Var1)){
        melted_cormat[melted_cormat$Var1==col,]$value_bis <- 100*melted_cormat[melted_cormat$Var1==col,]$value / melted_cormat[melted_cormat$Var1==col & melted_cormat$Var2==col,]$value
        }
  
    p <- ggplot(melted_cormat, aes(Var2, Var1, fill = value_bis))+
        geom_tile(color = "white")+
        scale_fill_gradient2(low = "#0072B2", high = "#009E73", mid = "#BFBFBF",space = "Lab", name="Proportion of overlapping patients (%)") +
        coord_fixed()+
        xlab("")+
        ylab("")+
        theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 22, hjust = 1,face="bold"),
        axis.text.y = element_text( vjust = 1, size = 22, hjust = 1,face="bold"), panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        legend.justification = c(2, 0),
        legend.direction = "horizontal",
        legend.text=element_text(size=15),
        legend.title=element_text(size=15),
        plot.title=element_text(size=25,face="bold",hjust=0.5))+guides(fill = guide_colorbar(barwidth = 15, barheight = 4,title.position = "top", title.hjust = 0.5),
        legend.title=element_text(size=22))+
        geom_text(aes(Var2, Var1, label = value), color = "black", size =10,face="bold") +
        ggtitle("Component Repartition")
              
    p
}

comp_frequency <- function(df_final,comp){

    df_tmp <- df_final[,comp]
    r <- data.frame(freq = 100*colSums(df_tmp)/nrow(df_final))
    r$names <- rownames(r)
    r$names <- str_replace(r$names,"full_component_","")
    p <- ggplot(r,aes(reorder(names,freq),freq))+
            geom_bar(stat="identity")+
            coord_flip()+
            theme(plot.title = element_text(hjust = 0.5,size=25),axis.text=element_text(size=22),axis.title=element_text(size=34,face="bold"),
            legend.position="none")+xlab("Components") + ylab("Frequency over the cohort (%)")
    
    p
}


comp_frequency_by_age <- function(df_final,comp){
    
    df_tmp <- df_final[,c(comp,"age")]

    df_tmp_young <- df_tmp[df_tmp$age<60,]
    young <- data.frame(freq = 100*colSums(df_tmp_young[,comp])/nrow(df_final),age="<60")
    young$names <- rownames(young)

    df_tmp_old <- df_tmp[df_tmp$age>=60,]
    old <- data.frame(freq = 100*colSums(df_tmp_old[,comp])/nrow(df_final),age=">=60")
    old$names <- rownames(old)

    r <- rbind(young,old)
    r
    r$names <- str_replace(r$names,"full_component_","")
    p <- ggplot(r,aes(reorder(names,freq),freq,fill=age))+
            geom_bar(stat="identity",aes(fill=age))+
            coord_flip()+
            xlab("Components") + ylab("Frequency over the cohort (%)")+scale_fill_manual(values=c("#e79f00", "#009E73"))+
    scale_y_continuous(breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1)))))+
    labs(fill='Age stratification') +ggtitle("% Repartition of components representation in dataset across age groups")+
    theme(axis.text.x=element_text(size=22,face="bold"),axis.title.x = element_text(size=25),axis.text.y=element_text(size=22,face="bold"),axis.title.y = element_text(size=25),legend.text = element_text(size=25),legend.title = element_text(size=30),
          plot.title = element_text(size=25,face="bold",hjust=0.5))+
                       ylab("Frequency (%)")+xlab("Components")
                       
    p
}
                       
                       
comp_stacked_frequency <- function(df_final,comp,colors=color_values){     ### color_values : dictionary with name of comp and color values
                       
    df_tmp <- df_final[,comp]
    r <- data.frame(freq = colSums(df_tmp))
    r$freq <- 100*r$freq/sum(r$freq)
    r$names <- rownames(r)
    r$names <- str_replace(r$names,"full_component_","")
    r$test <- "try"
    p <- ggplot(r,aes(test,freq,fill=reorder(names,freq),label=names))+
            geom_bar(stat="identity")+#coord_flip()+
    theme(plot.title = element_text(hjust = 0.5,size=25),axis.text.x=element_blank(),axis.text=element_text(size=22),axis.title=element_text(size=34,face="bold"),
    legend.position="none")+xlab("Components") + ylab("Frequency over the cohort (%)")+scale_fill_manual(values=colors)+
            geom_text(size = 5, position = position_stack(vjust = 0.5))

    p
    
}
                       

comp_clinical_correlates <- function(df_w_correlates,cols_component,binary_correlates,continuous_correlates,colors=color_values){

set_notebook_plot_size(30,5)

v=c("0" = "grey45", "1" = "#e79f00", "2" = "#009E73", "3" ="#0072B2", "4"="#CC79A7")
    
    #### Binary Correlates
tmp <- NULL
for (comp in cols_component){
    tmp1 <- df_w_correlates[df_w_correlates[,comp]==1,]
    tmp1$comparison <- str_remove(comp,"full_component_")
    tmp <- rbind(tmp,tmp1)
    }
tmp$fill_pal <- factor(tmp$comparison)
    
for (col in binary_correlates){
#     png(paste(col,"_clinical.png",sep=""),width=4500,height=1500,res=200)
    p <- ggplot(tmp, aes(reorder(factor(comparison),-tmp[,col]), fill = factor(tmp[,col]))) + 
    geom_bar(position = "fill") +
    scale_y_continuous(labels = scales::percent)+
    xlab("Full Components")+scale_fill_manual(values = v)+theme(plot.title = element_text(hjust = 0.5,size=25),axis.text=element_text(size=10,face="bold"),axis.title=element_text(size=14,face="bold"))+
    labs(fill=col)
    plot(p)
#     dev.off()
}
    
    #### Continuous Correlates

for (col in continuous_correlates){
#     png(paste(col,"_clinical.png",sep=""),width=4500,height=2500,res=200)
    p <- ggplot(tmp,aes(x=reorder(comparison,-tmp[,col],FUN=median),y=tmp[,col],fill=fill_pal)) +
    geom_boxplot() + 
    theme1+
    ylab(col)+xlab("Full Components")+theme(plot.title = element_text(hjust = 0.5,size=25),axis.text=element_text(size=10,face="bold"),axis.title=element_text(size=14,face="bold"))+
    scale_y_continuous(limits = quantile(tmp[,col], c(0.01, 0.99),na.rm = T))+
    scale_fill_manual(values=colors)+theme(legend.position = "none")

    
    suppressWarnings(plot(p))
    }
}
                       
                       
gene_gene_interactions <- function(genomicData){
    genomicData <- df_final[,c(genes,cytos)]
    genomicData <- (sapply(unique(sub("(_ITD)|(_TKD)|(_other)|(_mono)|(_bi)|(_p172)|(_p140)","",colnames(genomicData) )), function(x) rowSums(genomicData[grep(paste(x,"($|_.+)",sep=""), colnames(genomicData))])) > 0)+0
    genomicData <- genomicData[,colSums(genomicData, na.rm=TRUE)>=8]
    genomicGroups <- factor(grepl("^[a-z]", colnames(genomicData)) + grepl("t_*[0-9M]", colnames(genomicData) ), labels=c("Genetics","CNA","Fusions"))
    dim(genomicData)

    set_notebook_plot_size(20,20)
#     png("gene_gene_interaction.png",width=3000,height=3000,res=250)
    logPInt <- sapply(1:ncol(genomicData), function(i) sapply(1:ncol(genomicData), function(j) {
        f<- try(fisher.test(genomicData[,i], genomicData[,j]), silent=TRUE)
        if(class(f)=="try-error"){
            0
        } else{
            ifelse(f$estimate>1, -log10(f$p.val),log10(f$p.val))}}
                                                             ))

    odds <- sapply(1:ncol(genomicData), function(i) sapply(1:ncol(genomicData), function(j) {f<- try(fisher.test(table(genomicData[,i], genomicData[,j])), silent=TRUE); if(class(f)=="try-error") f=NA else f$estimate} ))
    pairs <- sapply(1:ncol(genomicData), function(i) colMeans(genomicData * genomicData[,i], na.rm=TRUE))
    diag(logPInt) <- 0
    diag(odds) <- 1
    colnames(odds) <- rownames(odds) <- colnames(logPInt) <- rownames(logPInt) <- colnames(genomicData) 
    odds[odds<1e-3] = 1e-4
    odds[odds>1e3] = 1e4
    odds[10^-abs(logPInt) > 0.05] = 1
    logOdds=log10(odds)
    diag(logPInt) <- NA
    par(bty="n", mgp = c(2,.5,0), mar=c(4,4,4,4)+.1, las=2, tcl=-.33) 
    ix = TRUE#colnames(interactions) %in% colnames(all_genotypes)
    d <- dist(t(genomicData[,ix]) + 10*as.numeric(genomicGroups))
    h = hclust(d, method="average")
    o = order(genomicGroups,-colSums(genomicData, na.rm=TRUE))#order(cmdscale(d, k=1))#h$order #c(h$order,(length(h$o rder) +1):ncol(interactions))
    M <- matrix( NA, ncol=ncol(odds), nrow=nrow(odds))
    M[lower.tri(M)] <- cut(logOdds[o,o][lower.tri(M)], breaks = c(-4:0-.Machine$double.eps,0:4), include.lowest=TRUE) 
    M[upper.tri(M, diag=TRUE)] <- as.numeric(cut(pairs[o,o][upper.tri(M, diag=TRUE)]*nrow(genomicData), breaks=c(-1,0 ,5,10,20,50,100,200,600))) + 9
    image(x=1:ncol(logPInt), y=1:nrow(logPInt), M, col=c(brewer.pal(9,"BrBG"), c("white",brewer.pal(7,"Reds"))), breaks=0:max(M,na.rm=TRUE), xaxt="n", yaxt="n", xlab="",ylab="", xlim=c(0, ncol(logPInt)+3), ylim=c(0, ncol(logPInt)+3))
    l <- colnames(logPInt)[o]
    mtext(side=1, at=1:ncol(logPInt), l, cex=.66, font=ifelse(grepl("^[A-Z]",l),3,1))
    mtext(side=2, at=1:ncol(logPInt), l, cex=.66, font=ifelse(grepl("^[A-Z]",l),3,1))
    abline(h=0:ncol(logPInt)+.5, col="white", lwd=.5)
    abline(v=0:ncol(logPInt)+.5, col="white", lwd=.5)
    P <- 10^-abs(logPInt[o,o])
    P[upper.tri(P)] <- NA
    w = arrayInd(which(p.adjust(P, method="BH") < .1), rep(nrow(logPInt),2))
    points(w, pch=".", col="black")
    w = arrayInd(which(p.adjust(P) < .05), rep(nrow(logPInt),2))
    points(w, pch="*", col="black")
    image(y = 1:9 +18, x=rep(ncol(logPInt),2)+c(2,3), z=matrix(c(1:8), nrow=1), col=c("white",brewer.pal(7,"Reds")),
    add=TRUE)
    axis(side = 4, at = seq(1,7) + 19, cex.axis=.66, tcl=-.15, label=c(1,5,10,20,50,100,200), las=1, lwd=.5)
    mtext(side=4, at=28, "Frequency", las=2, line=-1,cex=.66)
    image(y = 1:8 +5, x=rep(ncol(logPInt),2)+c(2,3), z=matrix(c(1:8), nrow=1), col=brewer.pal(8,"BrBG"), add=TRUE)
    axis(side = 4, at = seq(1,7) + 5.5, cex.axis=.66, tcl=-.15, label=10^seq(-3,3), las=1, lwd=.5)
    mtext(side=4, at=14, "Odds ratio", las=2, line=-1,cex=.66)
    mtext(side=4, at=4, "Significance", las=2, line=-1,cex=.66)
    points(x=rep(ncol(logPInt),2)+2.5, y=1:2, pch=c("*","."))
    image(x=rep(ncol(logPInt),2)+c(2,3), y=(2:3) +0.5, z=matrix(1), col=brewer.pal(3,"BrBG"), add=TRUE)
    mtext(side=4, at=3:1, c("P > 0.05", "FDR < 0.1", "FWER < 0.05"), cex=.66, line=0.2)
#     dev.off()
}
                    
                    
num_events_comp <- function(df_final,cols_component){
    tmp <- NULL
    for (co in cols_component){
        tmp1 <- df_w_correlates[df_w_correlates[,co]==1,]
        tmp1$comparison <- str_replace(co,"full_component_","")
        tmp <- rbind(tmp,tmp1)
        }
    tmp$genes <- rowSums(tmp[,c(genes)])
    tmp$cytos <- rowSums(tmp[,c(cyto)])
    tmp$total <- rowSums(tmp[,c(genes,cyto)])
    p <- ggplot(tmp, aes(x=reorder(comparison,total), y=total,fill=comparison)) + 
      geom_boxplot()+coord_flip()+scale_fill_manual(values=color_values)+xlab("Component")+ylab("Number of Total alterations")+
    theme(axis.title = element_text(hjust = 0.5,size=25,face="bold"),axis.text.x = element_text( hjust = 1, size=20,face="bold"),axis.text.y = element_text(size=20,face="bold"))+guides(fill=F)+ylim(c(0,15))

    q <- ggplot(tmp, aes(x=reorder(comparison,total), y=genes,fill=comparison)) + 
      geom_boxplot()+coord_flip()+scale_fill_manual(values=color_values)+xlab("")+ylab("Number of Gene alterations")+
    theme(axis.title = element_text(hjust = 0.5,size=25,face="bold"),axis.text.y=element_blank(),axis.text.x = element_text( hjust = 1, size=20,face="bold"))+guides(fill=F)+ylim(c(0,15))

    r <- ggplot(tmp, aes(x=reorder(comparison,total), y=cytos,fill=comparison)) + 
      geom_boxplot()+coord_flip()+scale_fill_manual(values=color_values)+xlab("")+ylab("Number of Cyto alterations")+
    theme(axis.title = element_text(hjust = 0.5,size=25,face="bold"),axis.text.y=element_blank(),axis.text.x = element_text( hjust = 1, size=20,face="bold"))+guides(fill=F)+ylim(c(0,15))
    
    grid.arrange(p,q,r,ncol=3)
}
                    
####################################### End Component Plots #######################################               
                    
                    
                    
                    
####################################### Regression Tools #######################################
         
multivariate_regression <- function(df_final,target,features=c(gen,cyto),fam="gaussian",num_iterations=100,lambda=NULL,top_n=10,size_title=25,size_axis=12){                  
    registerDoMC(cores=50)


    data <- df_final
    df_multi <- NULL
    i <- 1
    for (i in c(1:num_iterations)){
        res1 <- cv.glmnet(data.matrix(data[,features]), data[,target], family=fam,alpha=1,nfolds=10,parallel =T , lambda=lambda)
        df_multi <- cbind(df_multi,as.matrix(coef(res1,s="lambda.min")))
        i <- i+1
        }
    l <- data.frame(coef = rowSums(df_multi)/num_iterations)
    l$names <- rownames(l)
    l <- l[-1,]
    l$abs_coef <- abs(l$coef)
    l$Model <- ifelse(l$names %in% gen,"gen",
                  ifelse(l$names %in% cyto,"cyto",
                        ifelse(l$names %in% clin, "clin",
                              ifelse(l$names %in% demo, "demo",
                                    ifelse(l$names %in% eln, "eln",
                                           ifelse(l$names %in% mrd,"mrd","comp"))))))
    l <- l %>% top_n(top_n,abs_coef)
    t <- ggplot(l[l$coef!=0,],aes(x=reorder(names,abs_coef),y=coef,fill=Model))+geom_bar(stat="identity")+coord_flip()+
    theme(plot.title = element_text(hjust = 0.5,size=25),axis.text=element_text(size=size_axis),axis.title=element_text(size=size_title,face="bold"),
         legend.position="none")+
    scale_fill_manual(values=colors_analysis)+ylab(paste(paste("Lasso Regression Coefficients for ",target,sep="")," (averaged over 100 repetitions)",sep="")) + xlab("Top 10 Features")
    #     png(paste(target,"_lasso.png",sep=""),width=4500,height=2500,res=200)
    return(t)
    #     dev.off()
}
                    
                    
univariate_volcano <- function(df_final,target,features=c(gen,cyto),type="continuous",quantile=c(0,1),p_value_threshold=0.05,size_title=25,legend.position="top"){
    
    df <- data.frame(beta = double(),pvalue = double(),Frequency = double()) 
    
    #     png(paste(target,"_univariate.png",sep=""),width=4500,height=2500,res=300)
    data <- df_final
    for (col in features){
        if (type=="continuous"){
            fit <- lm(as.formula(paste(paste(target," ~ ",sep= ""),col,sep="")), data=data)
            df[col,1:3] <- c(coef(fit)[[col]],glance(fit)[["p.value"]],100*sum(data[,col])/dim(data)[1])
        }else{
            fit <- logistf(formula = as.formula(paste(paste(target," ~ ",sep= ""),col,sep="")), data = data)
            df[col,1:3] <- c(fit$coef[[2]],fit$prob[[2]],100*sum(data[,col])/dim(data)[1])
            }
        }
    

    df[,"adjusted_pvalue"] <- p.adjust(df$pvalue)
    df[,"-log10(adjusted_pvalue)"] <- -log(p.adjust(df$pvalue),10)
    df["names"] <- rownames(df)
    df$Model <- ifelse(df$names %in% gen,"gen",
                      ifelse(df$names %in% cyto,"cyto",
                            ifelse(df$names %in% clin, "clin",
                                  ifelse(df$names %in% demo, "demo",
                                        ifelse(df$names %in% eln, "eln",
                                               ifelse(df$names %in% mrd,"mrd","comp"))))))
    df$Frequency <- ifelse(df$names %in% clin, 5,
                                  ifelse(df$names %in% demo, 5, df$Frequency))

    df <- df[order(df$adjusted_pvalue),]
    s <- ggplot(df, aes(x=beta, y=-log10(adjusted_pvalue))) + #volcanoplot with log2Foldchange versus pvalue
            geom_point(aes(size=Frequency,col=Model)) + geom_text_repel(data=df[(df["pvalue"]<p_value_threshold) ,], aes(label=names,fontface=8))+scale_size_continuous(range = c(3,12)) + 
             scale_color_manual(values=colors_analysis,limits=names(colors_analysis))+  ## respect color in feature importance
            theme(plot.title = element_text(hjust = 0.5,size=25),axis.text=element_text(size=22),axis.title=element_text(size=size_title,face="bold"))  + theme(legend.position=legend.position) + scale_x_continuous(limits = quantile(df$beta, quantile,na.rm = T))
     #     png(paste(target,"_univariate.png",sep=""),width=4500,height=2500,res=300)
    return(s)
#     dev.off()
}
    
    
    
                    
                    
####################################### End Regression Tools #######################################
                    
                    
                    
####################################### Density Tools #######################################

shift_density <- function(tmp_final,cols_component){
    
    # png("tmp.png",width=4000,height=10500,res=200)
    df <- data.frame(mean_survival = double() )
    for (co in colnames(tmp_final[,cols_component])){
        df[co,] <- mean(tmp_final[tmp_final[,co]==1,]$os)
        df[co,"name"] <- co
    }
    components_ordered <- df[order(df$mean_survival,decreasing=T),]$name
    rows <- c('eln_2017_favorable','eln_2017_intermediate','eln_2017_adverse',components_ordered)
    set_notebook_plot_size(20,41.8)
    p <- list()
    # pdf(file="graphs/density_plots.pdf",width=15,height=3)
    for (co in c('eln_2017_favorable','eln_2017_intermediate','eln_2017_adverse',cols_component)){
        if (is.element(co,c("full_component_t_8_21","full_component_t_6_9","full_component_NPM1","full_component_CEBPA_bi","full_component_DNMT3A_IDH1_2","full_component_WT1","full_component_not_assigned"))){    
            tmp <- tmp_final[tmp_final[,co]==1 &  tmp_final$ITD==0,]
            fit <- kphaz.fit(tmp$os,tmp$os_status,q=1,method="nelson")
            df <- fit%>%
            as.data.frame()
            df[,"category"] <- co
            tmp <- tmp_final[tmp_final[,co]==1 &  tmp_final$ITD==1,]
            fit <- kphaz.fit(tmp$os,tmp$os_status,q=1,method="nelson")
            df1 <- fit%>%
            as.data.frame()
            df1[,"category"] <- paste(co,"ITD",sep="_")        
            df <- rbind(df,df1)

         }else{
            tmp <- tmp_final[tmp_final[,co]==1 ,]
            fit <- kphaz.fit(tmp$os,tmp$os_status,q=1,method="nelson")
            df <- fit%>%
            as.data.frame()
            df[,"category"] <- co
            }   
        p[[co]] <- ggplot(df, aes(x=log(haz,10),color=category,fill=category)) + 
      geom_density(aes(y=..count..),alpha=1)+guides(fill=FALSE)+guides(fill=F,color=F)+ggtitle(co)+xlim(c(-2,1)) +theme(plot.title = element_text(hjust = 0.5,size=45),axis.text=element_text(size=18),
            axis.title=element_text(size=22,face="bold"))
        if(co=="eln_2017_favorable"){
            p[[co]]<- p[[co]]+ scale_fill_manual(values=c("#2b8cbe"))
        }
        if(co=="eln_2017_intermediate"){
            p[[co]]<- p[[co]]+ scale_fill_manual(values=c("#59A08B"))
        }
        if(co=="eln_2017_adverse"){
            p[[co]]<- p[[co]]+ scale_fill_manual(values=c("#EE9937"))
        }
        if(is.element(co,components_ordered)){
            p[[co]] <- p[[co]]+ scale_fill_manual(values=c("#870C14","#a6bddb","#a6bddb","pink","#fdbb84"))

            }

        }
    do.call(grid.arrange,c(p,nrow=length(rows)))
    # dev.off()
}
                    
####################################### End Density Tools #######################################