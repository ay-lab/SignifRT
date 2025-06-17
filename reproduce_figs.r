library(optparse)
library(stringr)

# Define the command-line arguments
option_list <- list(
  make_option(c("-i", "--inputfile"), type="character", default=NULL,
              help="Path to the input file", metavar="character"),
  make_option(c("-t", "--inputtable"), type="character", default=NULL,
              help="Path to the significance table", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=NULL,
              help="Path to the output folder and output file prefix", metavar="character"),
  make_option(c("-s", "--sign"), type="character", default="./Significance.R",
              help="Path to Significance.R", metavar="character")
)

# Parse the command-line arguments
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

dppa<- opt$inputfile
sig_table<-opt$inputtable
out<-opt$out
functions_file<-opt$sign

cat(paste('input:', dppa,'\t\t'))
cat(paste('table:', sig_table,'\t\t'))
cat(paste('out:', out,'\t\t'))

source(functions_file)

############copied from plotting from table nb
results3<-read.table(sig_table, sep='\t', header=TRUE)
head(results3)

#read and process data from master table
dppa <- read.table(dppa, header=TRUE)
colnames(dppa)[1] <- c("chr")
print(dim(dppa))

# these are all the deletions we have:
views<-names(dppa)[-which(names(dppa) %in% c("start","chr",'end'))]
length(views)
#unique(str_extract(pattern = '^([^_]+)_clone.+', string = views,  group = TRUE))
count<-as.data.frame(table(str_extract(pattern = '^([^_]+)_.+', string = views,  group = TRUE)))
sum(count$Freq)

#get prefixes for all names - skip coord columns from df
prefixes <- unique(sapply(strsplit(names(dppa)[4:ncol(dppa)], "_"), `[`, 1))

#dels
dels<-c('X100kdel','X145kdel','X245kdel','X435kdel','X680kdel','A.B.C', #old, from Jiao's paper
        'A.B','A.B.COSN123','A.BOSN12.BTSS.COSN123','A.BOSN12.COSN123','A.BOSN1.COSN123','A.BTSS.C',
        'A.C','A','A.COSN123','A.COSN12','AOSN.B.COSN123','AOSN.BOSN12.COSN123',
        'AOSN.BOSN1.COSN123','AOSN.COSN123','AOSN.COSN12','ATSS.BTSS.CTSS','B.C','B',
        'BOSN12.COSN123','BOSN1.COSN123','BTSS','BTSS.CTSS','C','CL20k','COSN123',
    'COSN12','COSN1','CR11.8k','CR13.4k','CR14.6k','CTSS','Z.AOSN.COSN123','Z.B.C', 'ATSS'
)

##################################### extra columns and column renaming for internal iterations
#1. one column with the average of all controls (WT) - per genome
dppa$All_avg_genome1_CON<-rowMeans(dppa[,grep('.+genome1_CON', names(dppa))], na.rm = TRUE)
dppa$All_avg_genome2_CON<-rowMeans(dppa[,grep('.+genome2_CON', names(dppa))], na.rm = TRUE)
cols_CON<-grep('All_avg_genome*_CON',names(dppa))

#2. Rename ATSS columns from triple-TSS clones, genome2
find<-grep('ATSS.BTSS.CTSS_clone.+_genome2_ATSS',names(dppa))
names(dppa)[find]<-str_replace(string = str_replace(string = names(dppa)[find], pattern = '(ATSS.BTSS.CTSS)', replacement = 'ATSS'),
            pattern= '(genome2_ATSS)', replacement = 'genome2_DEL')
cat("\nAdjusted ATSS alts.\n")
cat(names(dppa)[find])

#finally, one BTSS.CTSS clone genome 2 is also BTSS
find<-grep('BTSS.CTSS_clone.+_genome2_BTSS',names(dppa))
names(dppa)[find]<-str_replace(string = str_replace(string = names(dppa)[find], pattern = '(BTSS.CTSS)', replacement = 'BTSS'),
           pattern= '(genome2_BTSS)', replacement = 'genome2_DEL')
cat("\nAdjusted BTSS alts. p2\n")
cat(names(dppa)[find])

num=0.05
results3$pass<-ifelse(results3$pval_empirical<=num,TRUE,FALSE)
results3$name<-factor(results3$name, levels=rev(sort(results3$name)))

###########figure 1 specific stats panel
plotf='fig1_area_between_the_curves.pdf'
cairo_pdf(paste0(out,plotf), height =3 , width =3.75 )

chosen<-c('CL20k', 'CR14.6k', 'CR13.4k','CR11.8k', 'CTSS', 'COSN1', 'COSN12', 'COSN123','C', 'A.B.COSN123','A.B.C')
tmp<-results3[results3$name%in%chosen,]
tmp$name<-factor(tmp$name, levels=chosen)
fig1_signS<-ggplot(tmp, aes(y=name, x=abs(diff), color=pass, size=-log10(pval_empirical))) + 
geom_point() + theme_classic(base_size = 10) +
ggtitle(label = paste('')) +
labs(size='-log10(pval)', color=paste0('p-value<',num), x='|\u0394 area|', y='Deletion') +
scale_color_manual(values = c('maroon','navyblue'))
fig1_signS
dev.off()

#supplementary fig S1: stats panel
plotf='figS1_area_between_the_curves.pdf'
cairo_pdf(paste0(out,plotf), height =3 , width =5 )

chosen<-c('A', 'B', 'C', 'A.B', 'B.C','A.C', 'A.B.C')
tmp<-results3[results3$name%in%chosen,]
tmp$name<-factor(tmp$name, levels=chosen)
tmp$pass_label<-ifelse(tmp$pass, paste("p-value ≤",num),'not significant')
figS1_signS<-ggplot(tmp, aes(y=name, x=abs(diff), color=pass_label, size=-log10(pval_empirical))) + 
geom_point() + theme_classic(base_size = 10) +
ggtitle(label = paste('')) +
labs(size='-log10(pval)', color='', x='|\u0394 area|', y='Deletion') +
scale_color_manual(values = c('maroon','navyblue'))
figS1_signS
dev.off()

###########figure 2 specific stats panel
plotf='fig2_area_between_the_curves.pdf'
cairo_pdf(paste0(out,plotf), height =3 , width =5 )

chosen<-c('A','Z.AOSN.COSN123','AOSN.COSN123','COSN123','A.COSN123','B.C', 'Z.B.C','A.C')
tmp<-results3[results3$name%in%chosen,]
tmp$name<-factor(tmp$name, levels=chosen)
fig1_signS<-ggplot(tmp, aes(y=name, x=abs(diff), color=pass, size=-log10(pval_empirical))) + 
geom_point() + theme_classic(base_size = 10) +
ggtitle(label = paste('')) +
labs(size='-log10(pval)', color=paste0('p-value<',num), x='|\u0394 area|', y='Deletion') +
scale_color_manual(values = c('maroon','navyblue'))
fig1_signS
dev.off()

###########figure 3 specific stats panel
plotf='fig3_area_between_the_curves.pdf'
cairo_pdf(paste0(out,plotf), height =3 , width =5)

chosen<-c('AOSN.BOSN1.COSN123',
          'BOSN1.COSN123', 'COSN123', 'A.BOSN1.COSN123', 'A.COSN123', 
          'BOSN12.COSN123','AOSN.BOSN12.COSN123','B.C','A.BOSN12.COSN123', 'A.B.C')

tmp<-results3[results3$name%in%chosen,]
tmp$name<-factor(tmp$name, levels=chosen)
fig1_signS<-ggplot(tmp, aes(y=name, x=abs(diff), color=pass, size=as.character(round(-log10(pval_empirical), 0) ))) + 
geom_point() + theme_classic(base_size = 10) +
ggtitle(label = paste('')) +
labs(size='-log10(pval)', color=paste0('p-value<',num), x='|\u0394 area|', y='Deletion') +
scale_color_manual(values = c('FALSE'='maroon','TRUE'='navyblue')) +
scale_size_manual(values= c('1'=1,'2'=2,'3'=3,'4'=4) )
fig1_signS
dev.off()

###########figure 4 specific stats panel
plotf='fig4_area_between_the_curves.pdf'
cairo_pdf(paste0(out,plotf), height =2.5 , width =5 )

chosen<-c('ATSS', 'BTSS', 'ATSS.BTSS.CTSS', 'BTSS.CTSS', 'B.C', 
'A.BOSN12.BTSS.COSN123', 'A.BOSN12.COSN123','A.BTSS.C', 'A.B.C')

tmp<-results3[results3$name%in%chosen,]
tmp$name<-factor(tmp$name, levels=chosen)
fig1_signS<-ggplot(tmp, aes(y=name, x=abs(diff), color=pass, size=-log10(pval_empirical))) + 
geom_point() + theme_classic(base_size = 10) +
ggtitle(label = paste('')) +
labs(size='-log10(pval)', color=paste0('p-value<',num), x='|\u0394 area|', y='Deletion') +
scale_color_manual(values = c('maroon','navyblue'))
fig1_signS
dev.off()


###########
#2. Supplementary Fig 2 - all deletions RT panels with replicates¶
###########

#iterate over all target deletions, plot avg for each
plots<-list()

#iterate over prefixes of deletions
for (p in 1:length(dels)) {
#for (p in c(21)) { #run just one to test plotting
    prefix<-dels[[p]]
    
    ##########find columns with prefix match
    cols <- grep(paste0("^", prefix, "_"), names(dppa))
    #check
    if(length(cols) == 0) {
        cat(paste('Warning: no columns match for:', prefix, '\n'))
        next 
    } 
        
    ###########Find and classify columns
    #find prefix matches to ID genome
    del_cols<-grep(paste0("^", prefix, "_.+DEL"), names(dppa))
    #ID genome
    gen=Mode(rep(str_extract(pattern = 'genome(.)_', string = names(dppa)[del_cols], group=TRUE),2))[[1]]
    cat(paste('Using genome:', gen, '\n'))
    
    #find gen deletion columns (redundant but good to have)
    #and build names as per condition{genome}.{replicate} 
    #del_cols<-grep(paste0("^", prefix, "_clone.+genome",gen,"_DEL"), names(dppa))
    del_cols<-grep(paste0("^", prefix, "_clone.+_DEL"), names(dppa))
    del_col_names<-paste0("c1.",1:length(del_cols))

    #find control or obtain 'master'
    con_cols<-grep(paste0("^", prefix, "_clone.+_CON"), names(dppa))
    control_color<-'black'
    #use master control for corresponding genome if no control found
    if(length(con_cols)<1){
        #use 'master' control
        cat('Warning: No internal control found, defaulting to master\n')
        gen2<-ifelse(gen=='1','2','1')
        con_cols<-grep(paste0('All_avg_genome',gen2,'_CON'),names(dppa))
        control_color<-'grey35'

    }
    con_col_names<-paste0("c2.",1:length(con_cols))
    
    
    #fetch coordinates, dels and controls
    regs<- dppa[,c(1:3,del_cols,con_cols)]
    cat("Done selecting columns. using:\n")
    print(names(regs))
    #name adjustments
    names(regs)<-c('chr','start','end', del_col_names, con_col_names)
    cat("Adjusted column names:\n")
    print(names(regs))
    
    #adjust some names of all deletions when displayed at top of plot
    deprecated_names<-c('X100kdel','X145kdel','X245kdel','X435kdel','X680kdel')
    if(prefix%in%deprecated_names){
      deletion_display_name<-paste0('DEL',str_extract(prefix, pattern='^X(.{3}k)del', group=TRUE))
    }
    else{deletion_display_name<-paste0('\u0394',prefix)}
    
    ############# plot
    colorG<-rep(c('red',control_color), c(length(del_cols),length(con_cols)))
    names(colorG)<-names(regs)[-c(1:3)]
    p<- plot_experiment2(data=regs, experiment_name=deletion_display_name, colorG=colorG, 
                         subT= paste0('|\u0394 area|=', round(abs(results3$diff[results3$name==prefix]), digits = 2),
                                      ', pval=', round(results3$pval_empirical[results3$name==prefix],digits = 5)),
                         pos_L='none', order=TRUE, levels=rev(names(colorG)))
    plots[[prefix]] <- p 
 }

###########save
cairo_pdf(paste0(out,'FigS2.pdf'), height =16 , width =12.5 ) # RULE: width=2.5*n, height=2*m

subset<-c('A.B.C', 'A.B','A.B.COSN123','A.BOSN12.BTSS.COSN123','A.BOSN12.COSN123','A.BOSN1.COSN123','A.BTSS.C',
        'A.C','A','A.COSN123','A.COSN12','AOSN.B.COSN123','AOSN.BOSN12.COSN123',
        'AOSN.BOSN1.COSN123','AOSN.COSN123','AOSN.COSN12','ATSS.BTSS.CTSS','B.C','B',
        'BOSN12.COSN123','BOSN1.COSN123','BTSS','BTSS.CTSS','C','CL20k','COSN123',#'COSN1.2.8k',
    'COSN12','COSN1','CR11.8k','CR13.4k','CR14.6k','CTSS','Z.AOSN.COSN123','Z.B.C', 'ATSS',
    'X100kdel','X145kdel','X245kdel','X435kdel','X680kdel')

#tmp<-plots[sort(names(plots))]
tmp<-plots[sort(subset)]
tmp2<-lapply(seq_along(tmp), function(i) add_roman_annotation(tmp[[i]], i))

layout <- rbind(c(1:5),
                c(6:10),
                c(11:15),
                c(16:20),
                c(21:25),
                c(26:30),
                c(31:35),
                c(36:40))

grid.arrange(grobs =tmp2 , layout_matrix = layout)

dev.off()


########################## panel B - with Jiao's
cairo_pdf(paste0(out,'FigS2B.pdf'), height =2 , width =12.5 ) # RULE: width=2.5*n, height=2*m

subset<- c('X100kdel','X145kdel','X245kdel','X435kdel','X680kdel')

#tmp<-plots[sort(names(plots))]
tmp<-plots[sort(subset)]
tmp2<-lapply(seq_along(tmp), function(i) add_roman_annotation(tmp[[i]], i))

layout <- rbind(c(1:5))

grid.arrange(grobs =tmp2 , layout_matrix = layout)

dev.off()

#############
#3. Figure-specific RepliSeq panels: avg of replicates
#############
# generate figures with replicates averaged into a single line

#iterate over all target deletions, find all replicates, average replicates and plot
plots2<-list()

#iterate over prefixes of deletions
for (p in 1:length(dels)) {
    prefix<-dels[[p]]
    print(prefix)
        
    ###########Find and classify columns
    #find prefix matches to ID genome
    del_cols<-grep(paste0("^", prefix, "_.+DEL"), names(dppa))
    #ID genome
    gen=Mode(rep(str_extract(pattern = 'genome(.)_', string = names(dppa)[del_cols], group=TRUE),2))[[1]]
    
    #find gen deletion columns (redundant but good to have, higher specificity)
    #and build names as per condition{genome}.{replicate} 
    del_cols<-grep(paste0("^", prefix, "_clone.+genome",gen,"_DEL"), names(dppa))

    
    #find control or obtain 'master'
    con_cols<-grep(paste0("^", prefix, "_clone.+_CON"), names(dppa))
    control_color<-'black'
    #use master control for corresponding genome if no control found
    if(length(con_cols)<1){
        #use 'master' control
        cat('Warning: No internal control found, defaulting to master\n')
        gen2<-ifelse(gen=='1','2','1')
        control_color<-'grey35'
        con_cols<-grep(paste0('All_avg_genome',gen2,'_CON'),names(dppa))
    }
    
    #fetch coordinates
    regs<- dppa[,c(1:3)]
    if(length(del_cols)<2){regs$del<-dppa[,del_cols]}
    else {regs$del<-rowMeans(dppa[,del_cols], na.rm=TRUE)}
    if(length(con_cols)<2){regs$con<-dppa[,con_cols]}
    else {regs$con<-rowMeans(dppa[,con_cols], na.rm=TRUE)}
    
    print(con_cols)
    print(del_cols)
    
    #plot
    colorG<-c('del'='red','con'=control_color)
    p<- plot_experiment2(data=regs, experiment_name=paste0('\u0394',prefix), colorG=colorG, 
                         subT= paste0('|\u0394 area|=', round(abs(results3$diff[results3$name==prefix]), digits = 2),
                                      ', pval=', round(results3$pval_empirical[results3$name==prefix],digits = 5)),
                         pos_L='none')
    plots2[[prefix]] <- p + theme(legend.position=c(0.9,0.8),
                                  legend.key.size = unit(0.3, "cm"),
                                  legend.text = element_text(size = 6),
                                  legend.title = element_blank(),
                                  legend.background = element_blank()) +
    labs(color='Allele')

 }

###########save
cairo_pdf(paste0(out,'AllRepliSeq-Avg.pdf'), height =16 , width =10 )

layout <- rbind(c(1:5),
                c(6:10),
                c(11:15),
                c(16:20),
                c(21:25),
                c(26:30),
                c(31:35),
                c(36:40))

grid.arrange(grobs = plots2[sort(names(plots2))], layout_matrix = layout)

dev.off()

####bonus
#NPC and ESC
prefix='ESC_NPC'
regs<-dppa[,c(1:3)]
regs$ESC<-rowMeans(dppa[,grep('^F121.9.4DN_ESC_rep._genome1_quantile_scaled_loess', names(dppa))], na.rm=TRUE)
regs$NPC<-rowMeans(dppa[,grep('^F121.9.4DN_NPC_rep._genome1_quantile_scaled_loess', names(dppa))], na.rm=TRUE)

#############
#Supplementary Fig 1
#############

#plot, add to plot list
colorG<-c('NPC'='maroon','ESC'='black')
p<- plot_experiment2(data=regs, experiment_name='ESC → NPC', colorG=colorG, 
                         subT= '', pos_L='none')
plots2[[prefix]] <- p + theme(legend.position=c(0.9,0.8),
                                  legend.key.size = unit(0.3, "cm"),
                                  legend.text = element_text(size = 6),
                                  legend.title = element_blank(),
                                  legend.background = element_blank())
                    labs(color='Allele')

#supplementary fig 1: A.B.C combos
cairo_pdf(paste0(out,'_Sfig1_abc.pdf'), height =4 , width =10 )

chosen<-c('ESC_NPC','A', 'B', 'C', 'A.B', 'B.C','A.C', 'A.B.C')
tmp<-plots2[chosen]
tmp2<-lapply(seq_along(tmp), function(i) add_roman_annotation(tmp[[i]], i))
layout <- rbind(c(1:4),
                c(5:8))
             
grid.arrange(grobs = tmp2, layout_matrix = layout)
dev.off()

##########
# Fig 1
##########

#figure 1 independent panels
cairo_pdf(paste0(out,'_fig1_indP.pdf'), height =6 , width =7.5 )

chosen<-c('CL20k', 'CR14.6k', 'CTSS', 'COSN1', 'CR11.8k','COSN12', 'COSN123', 'CR13.4k', 'A.B.COSN123')
tmp<-plots2[chosen]
tmp2<-lapply(seq_along(tmp), function(i) add_roman_annotation(tmp[[i]], i))
layout <- rbind(c(1:3),
                c(4:6),
               c(7:9))
grid.arrange(grobs = tmp2, layout_matrix = layout)
dev.off()

############
##figure 1 comparison panel: 1. 
exp='A.B.COSN123'
colsDel <- grep(paste0('^',exp,'_clone.+DEL'), names(dppa))
colsCon <- grep(paste0('^A.B.C_clone.+CON'), names(dppa))
colsC3<-grep(paste0('^A.B.C_clone.+DEL'), names(dppa)) 
colsC2<-grep(paste0('^A.B_clone.+DEL'), names(dppa)) 
colsC<-grep(paste0('^C_clone.+DEL'), names(dppa))

sub<-dppa[,c(1:3)]

sub$control<-rowMeans(dppa[,colsCon],na.rm = TRUE)
sub$`ΔA.B.COSN123`<-rowMeans(dppa[,colsDel],na.rm = TRUE)

#add Jiao's A.B.C deletion
sub[,6]<-rowMeans(dppa[,colsC3],na.rm = TRUE)
names(sub)[6]<-'\u0394A.B.C'

#add A.B
sub[,7]<-rowMeans(dppa[,colsC2],na.rm = TRUE)
names(sub)[7]<-'\u0394A.B'

#add C
sub[,8]<-rowMeans(dppa[,colsC],na.rm = TRUE)
names(sub)[8]<-'\u0394C'

colorG<-c('control'='black',
          '\u0394A.B.COSN123'='#F33B7F',
          '\u0394A.B.C'='#0004B7',
          '\u0394A.B'='#1A85FF',
          '\u0394C'='grey')

colorG<-c('control'='black', 
          '\u0394A.B.COSN123'='red',
          '\u0394A.B.C'='#59a89c',
          '\u0394A.B'='#f0c571',
          '\u0394C'='darkblue')

fig1_abc<-plot_experiment2(data = sub, experiment_name = '\u0394A.B.COSN123 \u2245 \u0394A.B.C', pos_L = 'bottom',
                            colorG = colorG)

#fig1_abc
#figure 1 comparison panel: 2
exp='COSN123'
colsDel <- grep(paste0('^',exp,'_clone.+DEL'), names(dppa))
colsCon <- grep(paste0('^C_clone.+CON'), names(dppa))
cols1<-grep(paste0('^COSN1_clone.+DEL'), names(dppa)) 
cols12<-grep(paste0('^COSN12_clone.+DEL'), names(dppa)) 
colsC<-grep(paste0('^C_clone.+DEL'), names(dppa))

sub<-dppa[,c(1:3)]

#add COSN123 and control
sub$control<-rowMeans(dppa[,colsCon],na.rm = TRUE)
sub$`ΔCOSN123`<-rowMeans(dppa[,colsDel],na.rm = TRUE)

#add COSN1 deletion
sub[,6]<-rowMeans(dppa[,cols1],na.rm = TRUE)
names(sub)[6]<-'\u0394COSN1'

#add COSN12
sub[,7]<-rowMeans(dppa[,cols12],na.rm = TRUE)
names(sub)[7]<-'\u0394COSN12'

#add C
sub[,8]<-rowMeans(dppa[,colsC],na.rm = TRUE)
names(sub)[8]<-'\u0394C'

colorG<-c('control'='black', 
          '\u0394COSN123'='red',
          '\u0394COSN12'='#59a89c',
          '\u0394COSN1'='#f0c571',
          '\u0394C'='darkblue')

fig1_osn<-plot_experiment2(data = sub, experiment_name = '\u0394COSN123 \u2245 \u0394C', pos_L = 'bottom',
                            colorG = colorG)

#fig1_osn
############

#merge comparison panels
cairo_pdf(paste0(out,'_fig1_compP.pdf'), height =8 , width =5 )

grid.arrange(grobs = list(fig1_abc, fig1_osn), nrow=2)
dev.off()

################
# figure 2 independent panels
################
cairo_pdf(paste0(out,'_Oldfig2_indP.pdf'), height =4 , width =12.5 )

chosen<-c('AOSN.COSN123', 'A.COSN123', 'Z.AOSN.COSN123', 'Z.B.C', 'BOSN1.COSN123', 'A.BOSN1.COSN123', 
          'BOSN12.COSN123', 'A.BOSN12.COSN123', 'AOSN.BOSN1.COSN123', 'AOSN.BOSN12.COSN123')

tmp<-plots2[chosen]
tmp2<-lapply(seq_along(tmp), function(i) add_roman_annotation(tmp[[i]], i))

layout <- rbind(c(1:5),
                c(6:10))

grid.arrange(grobs = tmp2, layout_matrix = layout)
dev.off()
################
# Fig2 independent panels
################

cairo_pdf(paste0(out,'_fig2_indP.pdf'), height =4 , width =5 )

chosen<-c('AOSN.COSN123', 'A.COSN123', 'Z.AOSN.COSN123', 'Z.B.C')

tmp<-plots2[chosen]
tmp2<-lapply(seq_along(tmp), function(i) add_roman_annotation(tmp[[i]], i))

layout <- rbind(c(1,2),
                c(3,4))

grid.arrange(grobs = tmp2, layout_matrix = layout)
dev.off()

################
# Fig3 independent panels
################

cairo_pdf(paste0(out,'_fig3_indP.pdf'), height =6 , width =5 )

chosen<-c('BOSN1.COSN123', 'A.BOSN1.COSN123', 
          'BOSN12.COSN123', 'A.BOSN12.COSN123', 
          'AOSN.BOSN1.COSN123', 'AOSN.BOSN12.COSN123')

tmp<-plots2[chosen]
tmp2<-lapply(seq_along(tmp), function(i) add_roman_annotation(tmp[[i]], i))

layout <- rbind(c(1,2),
                c(3,4),
                c(5,6))

grid.arrange(grobs = tmp2, layout_matrix = layout)
dev.off()

#############
# figure 4
#############

# figure 4 independent panels
cairo_pdf(paste0(out,'_fig4_indP.pdf'), height =6 , width =5 )
chosen<-c('BTSS', 'BTSS.CTSS', 'ATSS', 'ATSS.BTSS.CTSS', 'A.BTSS.C', 'A.BOSN12.BTSS.COSN123')
tmp<-plots2[chosen]
tmp2<-lapply(seq_along(tmp), function(i) add_roman_annotation(tmp[[i]], i))

layout <- rbind(c(1,2),
                c(3,4),
                c(5,6))
grid.arrange(grobs = tmp2, layout_matrix = layout)
dev.off()

##############
# figure 2-3 all putative comparison panels
#############
#Plot 1: ∆AC, ∆COSN123, ∆ACOSN123, ∆AOSNCOSN123
#figure 2 comparison panel: 1. 
exp=c('A.C','COSN123','A.COSN123', 'AOSN.COSN123')
panel_name='Fig2.p1'

sub<-dppa[,c(1:3)]

#add control from first
colsCon <- grep(paste0('^',exp[1],'_clone.+CON'), names(dppa))
sub$control<-rowMeans(dppa[,colsCon],na.rm = TRUE)

#add avg del allele from all others
for(prefix in exp){
  cols<-grep(paste0('^',prefix,'_clone.+DEL'), names(dppa))
  sub[,ncol(sub)+1]<-rowMeans(dppa[,cols],na.rm = TRUE)
}
print(dim(sub))
#add names
Fnames<-c('control',paste0('\u0394',exp))
names(sub) <- c('chr','start','end',Fnames)

colorG<-c('black','red','#59a89c','#f0c571','darkblue')
names(colorG)<-c(Fnames)

fig2p1<-plot_experiment2(data = sub, experiment_name = panel_name, pos_L = 'bottom',
                         colorG = colorG)

#Plot 2: ΔAOSN.COSN123, ΔB.C, ΔZ.AOSN.COSN123, ΔZ.B.C
#figure 2 comparison panel: 2. 
exp=c('AOSN.COSN123', 'B.C', 'Z.AOSN.COSN123', 'Z.B.C')
panel_name='Fig2.p2'

sub<-dppa[,c(1:3)]

#add control from first
colsCon <- grep(paste0('^',exp[1],'_clone.+CON'), names(dppa))
sub$control<-rowMeans(dppa[,colsCon],na.rm = TRUE)

#add avg del allele from all others
for(prefix in exp){
  cols<-grep(paste0('^',prefix,'_clone.+DEL'), names(dppa))
  sub[,ncol(sub)+1]<-rowMeans(dppa[,cols],na.rm = TRUE)
}
print(dim(sub))
#add names
Fnames<-c('control',paste0('\u0394',exp))
names(sub) <- c('chr','start','end',Fnames)

colorG<-c('black','red','#59a89c','#f0c571','darkblue')
names(colorG)<-c(Fnames)

fig2p2<-plot_experiment2(data = sub, experiment_name = panel_name, pos_L = 'bottom',
                         colorG = colorG)

#Plot 3: ΔBOSN12.COSN123, ΔBOSN1.COSN123, ΔCOSN123, ΔB.C
#Plot 3
#figure 2 comparison panel: 3. 
exp=c('B.C','BOSN1.COSN123','BOSN12.COSN123', 'COSN123')
panel_name='Fig2.p3'

sub<-dppa[,c(1:3)]

#add control from first
colsCon <- grep(paste0('^',exp[1],'_clone.+CON'), names(dppa))
sub$control<-rowMeans(dppa[,colsCon],na.rm = TRUE)

#add avg del allele from all others
for(prefix in exp){
  cols<-grep(paste0('^',prefix,'_clone.+DEL'), names(dppa))
  if(length(cols)>1){sub[,ncol(sub)+1]<-rowMeans(dppa[,cols],na.rm = TRUE)}
  else{sub[,ncol(sub)+1]<-dppa[,cols]}
}
print(dim(sub))
#add names
Fnames<-c('control',paste0('\u0394',exp))
names(sub) <- c('chr','start','end',Fnames)

colorG<-c('black','red','#59a89c','#f0c571','darkblue')
names(colorG)<-c(Fnames)

fig2p3<-plot_experiment2(data = sub, experiment_name = panel_name, pos_L = 'bottom',
                         colorG = colorG)

#Plot 4: ΔA.BOSN12.COSN123, ΔA.BOSN1.COSN123, ΔA.COSN123, ΔA.B.COSN123
#figure 2 comparison panel: 4. 
exp=c('A.B.COSN123', 'A.COSN123', 'A.BOSN1.COSN123', 'A.BOSN12.COSN123')
panel_name='Fig2.p4'

sub<-dppa[,c(1:3)]

#add control from first
colsCon <- grep(paste0('^',exp[1],'_clone.+CON'), names(dppa))
sub$control<-rowMeans(dppa[,colsCon],na.rm = TRUE)

#add avg del allele from all others
for(prefix in exp){
  cols<-grep(paste0('^',prefix,'_clone.+DEL'), names(dppa))
  if(length(cols)>1){sub[,ncol(sub)+1]<-rowMeans(dppa[,cols],na.rm = TRUE)}
  else{sub[,ncol(sub)+1]<-dppa[,cols]}
}
print(dim(sub))
#add names
Fnames<-c('control',paste0('\u0394',exp))
names(sub) <- c('chr','start','end',Fnames)

colorG<-c('black','red','#59a89c','#f0c571','darkblue')
names(colorG)<-c(Fnames)

fig2p4<-plot_experiment2(data = sub, experiment_name = panel_name, pos_L = 'bottom',
                         colorG = colorG)

#Plot 5: ΔA.BOSN12.COSN123, ΔA.C, ΔB.C, ΔA.C, ΔA.B.C
#figure 2 comparison panel: 5. 
exp=c('A.B.C','A.C','A.B','B.C','A.BOSN12.COSN123')
panel_name='Fig2.p5'

sub<-dppa[,c(1:3)]

#add control from first
colsCon <- grep(paste0('^',exp[1],'_clone.+CON'), names(dppa))
sub$control<-rowMeans(dppa[,colsCon],na.rm = TRUE)

#add avg del allele from all others
for(prefix in exp){
  cols<-grep(paste0('^',prefix,'_clone.+DEL'), names(dppa))
  if(length(cols)>1){sub[,ncol(sub)+1]<-rowMeans(dppa[,cols],na.rm = TRUE)}
  else{sub[,ncol(sub)+1]<-dppa[,cols]}
}
print(dim(sub))
#add names
Fnames<-c('control',paste0('\u0394',exp))
names(sub) <- c('chr','start','end',Fnames)

colorG<-c('black','red','#59a89c','#f0c571','darkblue','lightblue')
names(colorG)<-Fnames

fig2p5<-plot_experiment2(data = sub, experiment_name = panel_name, pos_L = 'bottom',
                         colorG = colorG)


####merge all fig2 comparison panels
cairo_pdf(paste0(out,'_fig2_compP.pdf'), height =12.6 , width =12 )

grid.arrange(grobs = list(fig2p1,fig2p2,fig2p3,fig2p4,fig2p5), ncol=2)
dev.off()

##########################


#############
# figure S3 - all stats comp
#############
num=0.05

subset<-c('A.B.C', 'A.B','A.B.COSN123','A.BOSN12.BTSS.COSN123','A.BOSN12.COSN123','A.BOSN1.COSN123','A.BTSS.C',
        'A.C','A','A.COSN123','A.COSN12','AOSN.B.COSN123','AOSN.BOSN12.COSN123',
        'AOSN.BOSN1.COSN123','AOSN.COSN123','AOSN.COSN12','ATSS.BTSS.CTSS','B.C','B',
        'BOSN12.COSN123','BOSN1.COSN123','BTSS','BTSS.CTSS','C','CL20k','COSN123',
    'COSN12','COSN1','CR11.8k','CR13.4k','CR14.6k','CTSS','Z.AOSN.COSN123','Z.B.C', 'ATSS')

results3$pass<-ifelse(results3$pval_empirical<=0.1,'<=0.1','NS')
results3$pass[which(results3$pval_empirical<=0.05)]<-'<=0.05'
results3$pass[which(results3$pval_empirical<=0.01)]<-'<=0.01'

results3$name<-factor(results3$name, levels=rev(sort(results3$name)))

plotf=paste0(out, '_FigS3_summaryStats.pdf')
cairo_pdf(plotf, height=8 , width=8 )

ggplot(results3[which(results3$name%in%subset),], aes(y=name, x=abs(diff), color=pass, size=-log10(pval_empirical))) + 
geom_point() + theme_classic(base_size = 15) +
ggtitle(label = paste('')) +
labs(size='-log10(pval)', 
     #color=paste0('p-value<=',num), 
     color='p-value',
     x='| \u0394 area |', y='Deletion') +
scale_color_manual(values = c('maroon','navyblue', 'gold','darkgreen'))

dev.off()


######################
#function: create comparison plot using deltas 

compPanel <- function(dppa, colorG){
  deletion_names<-names(colorG)
  sub<-dppa[,c(1:3)]

  #find avg del and control alleles and calc delta
  for(prefix in deletion_names){
    #del cols
    cols<-grep(paste0('^',prefix,'_clone.+DEL'), names(dppa))
    if(length(cols)>1){sub[,ncol(sub)+1]<-rowMeans(dppa[,cols],na.rm = TRUE)}
    else{sub[,ncol(sub)+1]<-dppa[,cols]}

    #control cols
    cols<-grep(paste0('^',prefix,'_clone.+CON'), names(dppa))
    if(length(cols)>1){sub[,ncol(sub)]<- sub[,ncol(sub)] - rowMeans(dppa[,cols],na.rm = TRUE)}
    else{sub[,ncol(sub)]<- sub[,ncol(sub)] - dppa[,cols]}

  }
  print(dim(sub))
  head(sub)
  
  #add names
  Fnames<-paste0('\u0394',deletion_names)
  names(sub) <- c('chr','start','end',Fnames)
  names(colorG)<-Fnames

  p<-plot_experiment2(data = sub, experiment_name = '', pos_L = 'bottom',
                          colorG = colorG) + ylab('\u0394 RT (DEL - CON)')
  return(p)
}

#####################

#comparison panels using function

#fig1-1
plotf=paste0(out, '_fig1CompP1.pdf')
cairo_pdf(plotf, height=3, width=5)
colorG<-c('A.B.COSN123'='red',
          'A.B.C'='#59a89c',
          'A.B'='#f0c571',
          'C'='darkblue')
compPanel(dppa=dppa, colorG=colorG)
dev.off()

#fig1-2
plotf=paste0(out, '_fig1CompP2.pdf')
cairo_pdf(plotf, height=3, width=5)
colorG<-c('COSN123'='red',
          'COSN12'='#59a89c',
          'COSN1'='#f0c571',
          'C'='darkblue')
compPanel(dppa=dppa, colorG=colorG)
dev.off()


#fig2-1
plotf=paste0(out, '_fig2CompP1.pdf')
cairo_pdf(plotf, height=3, width=5)
colorG<-c('red','#59a89c','#f0c571','darkblue')
names(colorG)<-c('A.C','COSN123','A.COSN123', 'AOSN.COSN123')
compPanel(dppa=dppa, colorG=colorG)
dev.off()

#fig2-2
plotf=paste0(out, '_fig2CompP2.pdf')
cairo_pdf(plotf, height=3, width=5)
colorG<-c('maroon','#59a89c','#f0c571','darkblue')
names(colorG)<-c('AOSN.COSN123', 'B.C', 'Z.AOSN.COSN123', 'Z.B.C')
compPanel(dppa=dppa, colorG=colorG)
dev.off()

#fig3
plotf=paste0(out, '_fig3CompP.pdf')
cairo_pdf(plotf, height=2, width=3)
colorG<-c('maroon','#59a89c','#f0c571','darkblue','lightblue')
names(colorG)<-c('A.B.C','A.C','A.B','B.C','A.BOSN12.COSN123')
compPanel(dppa=dppa, colorG=colorG) + theme(legend.position='none') #run once without theme adj to get legend
dev.off()

#fig4 putative
plotf=paste0(out, '_fig4CompP.pdf')
cairo_pdf(plotf, height=3, width=5)
colorG<-c('maroon','#f0c571','darkblue')
names(colorG)<-c('A.B.C','A.BOSN12.COSN123','A.BOSN12.BTSS.COSN123')
compPanel(dppa=dppa, colorG=colorG) #+ theme(legend.position='none') #run once without theme adj to get legend
dev.off()


########################

results3$pass<-ifelse(results3$pval_empirical<=0.05,'<=0.05','NS')
#results3$pass[which(results3$pval_empirical<=0.05)]<-'<=0.05'
#results3$pass[which(results3$pval_empirical<=0.01)]<-'<=0.01'

results3$name<-factor(results3$name, levels=rev(sort(results3$name)))

plotf=paste0(out, '_FigS3_summaryStats.pdf')
cairo_pdf(plotf, height=8 , width=8 )

ggplot(results3, aes(y=name, x=abs(diff), color=pass, size=-log10(pval_empirical))) + 
geom_point() + theme_classic(base_size = 15) +
ggtitle(label = paste('')) +
labs(size='-log10(pval)', 
     #color=paste0('p-value<=',num), 
     color='p-value',
     x='| \u0394 area |', y='Deletion') +
     scale_color_manual(values = c('<=0.05'='navyblue','NS'='maroon')) 
#scale_color_manual(values = c('maroon','navyblue', 'gold','darkgreen')) #LH 01.17.25 changed for consistency

dev.off()