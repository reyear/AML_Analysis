#-------------#
# utils
#-------------#
# pvalue to code
pval_to_signif_code <- function(p_vals) {
  return(ifelse(p_vals < 0.001, "***",
                ifelse(p_vals < 0.01, "**",
                       ifelse(p_vals < 0.05, "*", ""))))
}

#-------------#
# univariate
#-------------#

# univariate
# perform fisher association test
# for frequency of mutated genes across of group with 2 levels
fisher_gene_group <- function(mydata, mygenes, group="SEX", mylevels=NULL, method.adjust="BH") {
  # INPUT:
  # - mydata: NxP matrix with patients in row and genes as columns as well as the group
  # - gene mutations should be binary here
  
  if (is.null(mylevels)) {
    mylevels = sort(unique(mydata[,group]))
  }
  # one test per gene
  lres = lapply(mygenes, function(gg) {
    tt = table(mydata[,gg],mydata[,group])
    freq = 100* tt[2,] / apply(tt,2,sum)
    resf = fisher.test(tt)
    mydatares = data.frame(gene=gg, count=sum(tt[2,]), frequency=100*sum(tt[2,])/nrow(mydata),
                           t(freq), 
                           pval=resf$p.value, OR=resf$estimate, 
                           confLow=resf$conf.int[1], confUp=resf$conf.int[2])
    rownames(mydatares) = NULL
    return(mydatares) })
  ddres = do.call("rbind",lres)
  
  # multiple testing
  ddres$qval = p.adjust(ddres$pval,method=method.adjust)
  ddres$code = pval_to_signif_code(ddres$qval)
  ddres$label = paste(ddres$gene,ddres$code)
  
  return(ddres)
  
}

# univariate
# frequency plot
# with fisher significance
# inverted barplot
library(reshape2)
fisher_plot_bar <- function(ddres, qval_cutoff=0.01, cnum=NULL,
                            col_for_groups = c("#f1a340","#998ec3"),
                            col_significant_axis = "#a50f15",
                            col_regular_axis = "#636363"
) {
  # INPUT:
  # - ddres: result of fisher_gene_group()
  
  mylevels = colnames(ddres)[4:5]
  ddres$gene= as.vector(ddres$gene)
  # cutoff
  if (is.null(cnum)) {
    cnum = min(ddres$count[which(ddres$qval<qval_cutoff)])
  }
  ddres = ddres[sort(ddres$count,decreasing=T,index.r=T)$ix,]
  ddres = ddres[ddres$count>=cnum,]
  
  # melt
  mres = melt(ddres, id.vars=colnames(ddres)[!colnames(ddres)%in%mylevels])
  colnames(mres)[ncol(mres)] = "group_frequency"
  colnames(mres)[ncol(mres)-1] = "group"
  
  
  # inverted barplot
  mres$group_frequency[mres$group==mylevels[1]] = -mres$group_frequency[mres$group==mylevels[1]]
  ii = sort(ddres[,mylevels[1]],index.r=T)$ix
  mres$gene = factor(mres$gene, levels=as.vector(ddres$gene[ii]))
  mres$label = factor(mres$label, levels=as.vector(ddres$label[ii]))
  
  forcol = rep(col_regular_axis,length(levels(mres$gene)))
  mygf = as.vector(ddres$gene[ddres$code!=""])
  forcol[levels(mres$gene)%in%mygf] = col_significant_axis
  
  mymax = round(max(abs(mres$group_frequency)))
  mylim = round(mymax/10) * 10
  
  gg = ggplot(mres) + geom_bar(aes(x=label,y=group_frequency,fill=group),position="stack",stat="identity") + coord_flip() + theme_classic() + bold15 + topleg + nolegtitle + theme(axis.text.y = element_text(color=forcol)) + theme(axis.text.y = element_text(hjust=0)) + noytitle + scale_fill_manual(values=col_for_groups) + theme(axis.text.x=element_text(color=col_regular_axis)) + scale_y_continuous(breaks=seq(-mylim,mylim,10),limits=c(-(mymax+1),mymax+1)) + ylab("frequency within group")
  
  return(list(melted=mres, plot=gg))
  
}

# univariate
# frequency plot
# with fisher significance
# volcano plot
fisher_plot_volcano <- function(ddres,
                                baselog=10,
                                label_qval_cutoff = 0.001
) {
  
  # INPUT:
  # - ddres: result of fisher_gene_group()
  
  # minus log q value
  ddres$logq = -log(ddres$qval,base=baselog)
  
  # for text annotation
  ddres$toannotate = ""
  ddres$toannotate[ddres$qval<label_qval_cutoff] = as.vector(ddres$gene[ddres$qval<label_qval_cutoff])
  
  gg = ggplot(ddres) + geom_point(aes(x=OR,y=-log(ddres$qval,base=baselog),size=frequency)) + theme_classic() + bold20 + xlab("Odds Ratio") + geom_vline(xintercept=1,linetype="dotted",size=0.5,color="#636363") + geom_hline(yintercept=-log(label_qval_cutoff,base=baselog),linetype="dotted",size=0.5,color="#636363") + scale_x_continuous(trans=paste0("log",baselog)) + topleg + geom_text(aes(x=OR,y=-log(ddres$qval,base=baselog),label=toannotate),hjust=0, nudge_x=0.03)
  
  if (baselog==10) {
    gg = gg + ylab(expression(paste(-log[10], ,Qvalue)))
  }
  if (baselog==2) {
    gg = gg + ylab(expression(paste(-log[2], ,Qvalue)))
  }
  
  return(list(volcano=gg))
  
}

#-------------#
# KM curve
#-------------#
print_and_flush <- function(string) {
  # Print a string in stdout and flush the result (useful to print in Jupyter Notebook in for loops)
  # → Arguments
  #     - string: string to print
  
  cat(string)
  flush.console()
}

print_p_value_from_survival_curve <- function(data, column_name, value_1_name, value_2_name) {
  # work in progress...
  stmp <- survdiff(as.formula(sprintf('Surv(os_diag_years, os_status) ~ %s', column_name)), data = data[data[,column_name] %in% c(value_1_name, value_2_name),])
  pval <- 1 - pchisq(stmp$chisq, length(stmp$n) - 1)
  
  print_and_flush(sprintf('%-21s and %21s: p = %.3f\n', value_1_name, value_2_name, pval))
}

print_p_value_from_boxplot <- function(data, column_name, value_1_name, value_2_name, clinical_name) {
  # work in progress...
  
  vec1 = data[ which(data[,column_name] %in% value_1_name), clinical_name]
  vec2 = data[ which(data[,column_name] %in% value_2_name), clinical_name]
  
  pval <- t.test(vec1,vec2)$p.val 
  
  print_and_flush(sprintf('%-21s and %21s: p = %.3f\n', value_1_name, value_2_name, pval))
}