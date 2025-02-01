#####################################################################################
# Required packages
library(ggrepel)
library(stringr)
library(ggplot2)
library(DescTools)
library(stringr)
library(dplyr)
library(gridExtra)
library(grid)
library(tidyr)

######################################################################################
#function: sampling random regions of size n
sample_loci <- function(regions, blacklisted=NA, size=NA, res=NA, n=NA, flag=FALSE, exp_name=NA){
    ori<-regions
    cat(paste('Sampling regions for background model from:',nrow(regions), '\n'))
   
    #Step 1: keep only regions with signal
    regions<-na.omit(regions)
    cat(paste('After removing NAs:', nrow(regions), '\n'))
   
    #Step 2: remove blacklisted regions accordingly
    if(flag){ #flag
        regions<-regions[which( !paste(regions$chr, regions$start)%in%paste(blacklisted$chr, blacklisted$start) ),]
        #regions<-anti_join(regions, blacklisted, by=c('chr','start')) #alt approach
        cat(paste('Removed blacklisted regions...\n', nrow(regions),'\n'))
        print(head(blacklisted))
        print(nrow(blacklisted))
        regions<-na.omit(regions[,c('chr','start')])
    }
    cat(paste('Removed NAs again after binding:', nrow(regions),'\n'))
    #sort dataframe
    regions<-regions[order(regions$chr,regions$start),]
   
    #Step 3: permute all possible non-empty contiguous regions of requested size
    a<-size/res #num of bins at resolution to get req size
    possible<-c()
    for (i in 1:(nrow(regions)-a)){
        if ((regions$chr[i+a-1]==regions$chr[i]) & (regions$start[i+a-1]==(regions$start[i]+size))){
            possible<-c(possible, paste(regions$chr[i], regions$start[i]))
            #print(paste0(i, '/', nrow(regions), length(possible)))
        }
    }
   
    #Step 4: sample with replacement and return chr and start coordinates from original regions df
    set.seed(123) #since sampling regions for background is random, set seed for reproducibility
    cat(paste('Possible regions for background creation, n=', length(which(paste(ori$chr,ori$start)%in%possible)),'\n'))
    coords<-with(ori, sample(x = which(paste(chr,start)%in%possible), size = n, replace = TRUE))
    cat('Created background.\n')

    print(length(coords))
    return(coords)
}

######################################################################################
#function: compute difference (different methods)
compute_diff <- function(x, method='single'){
    x$x<-(1:nrow(x))
    #rename columns per condition
    names(x)<-str_replace(names(x), pattern='c1.', replacement='a')
    names(x)<-str_replace(names(x), pattern='c2.', replacement='b')
    #print(names(x))
    #number of reps per condition
    ra<-max(na.omit(as.integer(str_extract(string = names(x), pattern = '^a([[:digit:]]+)$', group=TRUE))))
    rb<-max(na.omit(as.integer(str_extract(string = names(x), pattern = '^b([[:digit:]]+)$', group=TRUE))))
    #print(r)
    if(ra==1 | rb==1){
        method='single'
    }
    x<-na.omit(x)
    
    if(method=='single'){
        if(ra==1){a<-x[,grep(paste0('a1'), names(x))]}
        else{ a<-rowMeans(x[,grep(paste0('^a'), names(x))], na.rm=TRUE)}
        if(rb==1){ b<-x[,grep(paste0('b1'), names(x))] }
        else{ b<-rowMeans(x[,grep(paste0('^b'), names(x))], na.rm=TRUE) }
        return(AUC(x=x$x, y=a-b, absolutearea = FALSE, na.rm=TRUE))
    }
       
    cat('ERROR: unrecognized method name. Use: "single".\n') #built for future approach extensions, only relevant provided
}

######################################################################################
#function: build background
background <- function(regions, blacklisted=NA, size=NA, res=NA, n=NA, method=NA, flag=FALSE, exp_name=NA){
   
    #Step 1: simplify regions list
    #column signal will represent both conditions - used to eliminate NAs
    regions$signal<-rowMeans(regions[,grep(pattern = 'c[12]\\.[[:digit:]]', x = names(regions))],na.rm = FALSE)
   
    #Step 2: call sample_loci function to obtain coordinates of background
    bckgrnd<-sample_loci(regions=regions[,c('chr','start','signal')], blacklisted=blacklisted, size=size, res=res, n=n, flag=flag, exp_name=exp_name)
    
    #Step 3: iterate over background regions to obtain measurement of background given the chosen metric (method)
    a<-size/res
    conc<-c()
    for (i in bckgrnd){
        sub<-regions[i:(i+a-1),]
        diff_tmp<-compute_diff(x=sub, method=method)
       
        conc<-c(conc, diff_tmp)
    }
   if(length(which(is.na(conc))>0)){cat('Warning: NAs in background\n')}
    
    #Step 4: return background
    return(conc)
}

######################################################################################
# function: calculate empirical p value
empirical_pvalue<- function(query, bckgrnd_model, exp_name, type_test){
    #type_test can be 'chg', 'del', 'adv'
    # query + background = universe
    bckgrnd_model<-c(query, bckgrnd_model)
    vals<-c()
    #different questions can be asked
    if (type_test=='del'){
        for (i in 1:length(query)){
            vals[i]<-length(which(bckgrnd_model<=query[i]))/length(bckgrnd_model)
        }
        return(vals)
    } else if (type_test=='adv'){
        for (i in 1:length(query)){
            vals[i]<-length(which(bckgrnd_model>=query[i]))/length(bckgrnd_model)
        }
        return(vals)
    } else if (type_test=='chg'){
        for (i in 1:length(query)){
            vals[i]<-length(which(abs(bckgrnd_model)>=abs(query[i])))/length(bckgrnd_model)
        }
        return(vals)
    }
    stop('Type of test (type_test) not recognized. Use: "chg", "del", or "adv"')
}

######################################################################################
# function: main
significance <- function(query='DPPA', regions, blacklisted=NA, exp_name='run1', type_test='del',
                         size=2000000, res=50000, n=10000, exclude_query=TRUE, method='area'){
    # Step 1: query input
    # fill in query for DPPA coordinates
    if(is.data.frame(query)){
        cat('Received query dataframe.\n')
        head(query)
    }
    else if(query=='DPPA'){
        query<-data.frame(chr='chr16', start=47500000, end=49500000)
    }
    else if(query=='DPPA_right'){
        query<-data.frame(chr='chr16', start=48550000, end=49500000)
        size<-query$end[1]-query$start[1]
    }
    else if(query=='DPPA_left'){
        query<-data.frame(chr='chr16', start=47500000, end=48550000)
        size<-query$end[1]-query$start[1]
    }
    # read df from user-provided queries in bed file
    else {
        cat(paste('Reading queries from file:', query,'\n'))
        query<-read.table(query, sep='\t', header=FALSE, comment.char ='#')[,1:3]
        names(query)<-c('chr','start','end')
        cat(paste('Read bed file with', nrow(query), 'queries.'))
    }
    cat(paste('Using background size=', size))
   
    # Step 2: process regions (.bg) data
    regions<-regions[order(regions$chr,regions$start),]
    cat('Succesfully ordered df.\n')
    print(dim(regions))
   
    # Step 3: Handle blacklisted regions
    black_flag<-FALSE
    if( !is.na(blacklisted) ){
        black_flag<-TRUE
        cat('Warning: Assuming that blacklisted regions share coordinates with resolution of regions.bg\t')
        blacklisted<-read.table(blacklisted, header = FALSE, sep='\t', comment.char = '#')[,1:3]
        names(blacklisted)<-c('chr','start','end')
        print(head(blacklisted))
        if(exclude_query){
            blacklisted<-rbind(blacklisted, query[,c('chr','start','end')])
        }
    } else if(exclude_query) {
        blacklisted<-query
    } 
       
    #Step 4: build background and obtain dist
    bckgrnd_model<-background(regions=regions, blacklisted=blacklisted, size=size, res=res, n=n, method=method, flag=black_flag, exp_name=exp_name)
    #report background model stats
    cat(paste0('Background model. Mean=',mean(bckgrnd_model),'. sd=', sd(bckgrnd_model)))
   
    #Step 5: calc measure for each query
    for (i in 1:nrow(query)){
        sub<-query[i,]
        sub<-regions[which(regions$chr==sub$chr & regions$start>=sub$start & regions$end <=sub$end),]
        #print(sub)
        query$diff[i]<-compute_diff(x=sub, method=method)
    }
   
    #Step 6: calc empirical pval
    query$name<-exp_name
    query$pval_empirical<- empirical_pvalue(query=query$diff, bckgrnd_model=bckgrnd_model, exp_name=exp_name, type_test=type_test)
    #adjust p value (only when mult. queries)
    if(nrow(query)>1){ #if multiple queries are tested against the same background need to correct pvals
    query$padj_empirical<-p.adjust(query$pval_empirical, method="BH")
    }
    #return results
    return(query)
}


##################################PLOTS

#for plotting, receives dataframe, plots everything depending on received columns,
#column name will be used as condition, colorG indicates color for each condition
plot_experiment2 <- function(data, experiment_name, colorG, Fdppa=TRUE, pos_L='bottom', order=FALSE, levels=NA, subT='') {
    #use only DPPA domain
    if(Fdppa){
        data<-subset(data, chr=="chr16" &start>=46.5e6&start<=50.5e6)
    }

    #plotting aesthetics - global
    XLAB <- "Chr16 (Mb)"
    XLIM <- c(46.5, 50.5)
    YLIM <- c(-2.2, 2.2)
    YLAB <- " Log2 E/L RT"
    XTICKS_breaks <- c(4.7e+07, 4.8e+07 , 4.9e+07 , 5.0e+07)
    XTICKS_labs <- c("47","48", "49", "50")
    
    data <- data %>%
    pivot_longer(cols = -c(chr, start, end), 
               names_to = "exp", 
               values_to = "rt")
    print(head(data))
    
    if(order){
        data$exp<-factor(data$exp, levels = levels)
    }
    
    p <- ggplot(data, aes(x = start, y=rt, color=exp)) +   geom_hline(yintercept=0, linetype=2, color='lightgrey') +
    geom_line(linewidth=0.5) +
        labs( x = XLAB, y = YLAB, color='') + ggtitle(label = experiment_name, subtitle=subT ) +
        theme(axis.line = element_line(linewidth =0.75, color='black'), 
              plot.title = element_text(hjust = 0.5, size=8.5),
              plot.subtitle = element_text(hjust=0.5, size=6, color='blue'),
              panel.background = element_blank(),
             legend.position = pos_L) +
        scale_x_continuous(breaks=XTICKS_breaks, labels=XTICKS_labs) +
        ylim(YLIM) + 
        scale_color_manual(values = colorG)
    
    return(p)
}
