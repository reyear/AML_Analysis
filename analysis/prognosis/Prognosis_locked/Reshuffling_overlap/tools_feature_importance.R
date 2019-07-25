# ---------------------------------- #
# tools_feature_importance.R
# ---------------------------------- #


# This file gathers some useful R tool functions for studying feature importance



library('ggplot2')
library('reshape2')
library('grid')
library('gridExtra')
library('ggthemes')
library('stringr')



# get summary of ratio for one feature and for one couple (algo, endpoint)
get_summary<-function(data,feature,algo,endpoint){
    # data: dataframe of ratios of CI for different features, endpoints and algorithms
    # feature: character of the feature of interest
    # algo: character, name of the algorithm of interest
    # endpoint: character, name of the endpoint of interest
    # returns summary table of ratio of that feature for that algorithm and that endpoint
    
    tmp <- data[data$feature==feature & data$algo==algo & data$endpoint==endpoint,]
    return(summary(tmp$ratio))
}





# get max median ratio for all features and all couples (algorithm,endpoint) in a given dataset
get_max_median<-function(data){
    # data: dataframe of ratios of CI for different features, endpoints and algorithms
    # returns dataframe of the max median ratio of all endpoints and algorithms for each feature, used to classify features
    
    # get all couples (algorithm,endpoint)
    couples <- expand.grid(unique(data$algo),unique(data$endpoint))

    feature_median_max <- data.frame(feature=unique(data$feature))

    feature_median_max$max_median <- sapply(feature_median_max$feature, function(x){
        med<-apply(couples,1,function(y){return(get_summary(data,x,y[1],y[2])[3])})
        return(max(med))})
    
    return(feature_median_max)
    
}






# get limits of axis for all features and all couples in a given dataset for a particular category
get_limits <- function(data) {
    # data: dataframe of ratios of CI for different features, endpoints and algorithms
    # returns vector (min,max) to get limits of the plots, min corresponding to the smallest 1st quantile for all features from one category and for all couples (endpoint,algorithm)
    
    couples <- expand.grid(unique(data$algo),unique(data$endpoint))

    feature_range <- data.frame(feature=unique(data$feature))

    feature_range$q1 <- sapply(feature_range$feature, function(x){
        q1<-apply(couples,1,function(y){return(get_summary(data,x,y[1],y[2])[2])})
        return(min(q1))})
    
    feature_range$q3 <- sapply(feature_range$feature, function(x){
        q3<-apply(couples,1,function(y){return(get_summary(data,x,y[1],y[2])[5])})
        return(max(q3))})
    
    limits<-c(min(feature_range$q1),max(feature_range$q3))
    
    return(limits)
}







# plots features grouped by categories for all endpoints and all algorithms
caterogies_boxplot<-function(data,center=TRUE){
    # data: dataframe of ratios of CI for different features, endpoints and algorithms
    # center: boolean, whether to display the center features
    
    
    categories <- unique(data$category)
    
    if (!center){
        categories<-categories[which(categories!='CENTER')]
    }
    
    
    glist<- gList()
    for(cat in categories){
        tmp<- data[data$category==cat,]
        
        # order features by importance
        feature_median_max<-get_max_median(tmp)
        feature_median_max<-feature_median_max[order(feature_median_max$max_median),]
        
        tmp$feature <- factor(tmp$feature, levels=feature_median_max$feature)
        
        # get limits for the graph, only 2nd and 3rd quantiles are shown for a better visualization
        ylimits<-get_limits(tmp)
        
        g <- ggplot(tmp, aes(feature, ratio,fill=endpoint)) + 
                    scale_fill_manual(values = c('#66c2a5','#fc8d62','#8da0cb','#e78ac3')) +
                    facet_wrap(~algo,ncol = length(unique(tmp$algo))) +
                    coord_flip(ylim=ylimits) + # only use this function to set the limits of the plot, otherwise boxplots will be recalculated according to the data between the set limits
                    geom_hline(yintercept = 1,linetype = 3)+
                    geom_boxplot(outlier.shape = NA,coef=0)+
                    theme(legend.position = 'top',plot.title = element_text(size=20, face="bold"))+
                    ggtitle(cat)
        
        glist[[cat]]<-g
    }
    
    # custom sizing of plots for each category
    size<-data.frame(category=categories)
    size$n_features<-sapply(size$category, function(x){return(length(unique(data$feature[data$category==x])))})
    
    return(grid.arrange(grobs = glist, ncol=1,heights=sqrt(size$n_features)))
}








# get p-value of comparison of permuted CI distribution and not permuted CI distribution for one feature and for one couple (algo, endpoint)
get_pval <- function(data,feature,algo,endpoint){
    # data: dataframe of ratios of CI for different features, endpoints and algorithms
    # feature: character of the feature of interest
    # algo: character, name of the algorithm of interest
    # endpoint: character, name of the endpoint of interest
    # returns p-value for that feature for that algorithm and that endpoint
    tmp <- data[data$feature==feature & data$algo==algo & data$endpoint==endpoint,]
    
    return(wilcox.test(tmp$ref_CI,tmp$permuted_CI))
    
}









# get p-value ratio for all features and all couples (algorithm,endpoint) in a given dataset
get_table_pval<-function(data){
    # data: dataframe of ratios of CI for different features, endpoints and algorithms
    # returns dataframe of the p-value for each endpoint and algorithm for each feature
    
    features_pval <- c()

    
    for( endpoint in unique(data$endpoint)){
        tmp <- data.frame(feature=unique(data$feature))
        for (algo in unique(data$algo) ){
            tmp[,algo]<- sapply(tmp$feature, function(x){return(get_pval(data,as.character(x),algo,endpoint)$p.value)})
        }    
        tmp$endpoint<-endpoint
        features_pval<-rbind(features_pval,tmp)
    }
    
    
    return(features_pval)
    
}





# get standard deviation of ratio for one feature and for one couple (algo, endpoint)
get_sd<-function(data,feature,algo,endpoint){
    # data: dataframe of ratios of CI for different features, endpoints and algorithms
    # feature: character of the feature of interest
    # algo: character, name of the algorithm of interest
    # endpoint: character, name of the endpoint of interest
    # returns standard deviation of ratio of that feature for that algorithm and that endpoint
    tmp <- data[data$feature==feature & data$algo==algo & data$endpoint==endpoint,]
    return(sd(tmp$ratio))
}





# get max (median-sd) ratio for all features and all couples in a given dataset
get_max_median_sd<-function(data){
    # data: dataframe of ratios of CI for different features, endpoints and algorithms
    # returns dataframe of the max (median-sd) ratio of all endpoints and algorithms for each feature, used to classify features
    
    couples <- expand.grid(unique(data$algo),unique(data$endpoint))

    feature_median_sd_max <- data.frame(feature=unique(data$feature))

    feature_median_sd_max$max_median_sd <- sapply(feature_median_sd_max$feature, function(x){
        med<-apply(couples,1,function(y){return(get_summary(data,x,y[1],y[2])[3])})
        sd<-apply(couples,1,function(y){return(get_sd(data,x,y[1],y[2]))})
        return(max(med-sd))})
    
    return(feature_median_sd_max)
    
}




# get max of 1st quantile for all features and all couples in a given dataset
get_max_1stq<-function(data){
    # data: dataframe of ratios of CI for different features, endpoints and algorithms
    # returns dataframe of the max of 1st quantile of all endpoints and algorithms for each feature, used to classify features
 
    couples <- expand.grid(unique(data$algo),unique(data$endpoint))

    feature_max_1stq <- data.frame(feature=unique(data$feature))

    feature_max_1stq$max_1stq <- sapply(feature_max_1stq$feature, function(x){
        q1<-apply(couples,1,function(y){return(get_summary(data,x,y[1],y[2])[2])})
        return(max(q1))})
    
    return(feature_max_1stq)
    
}