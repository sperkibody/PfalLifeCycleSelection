library(ggplot2)
library(tidyverse)
library(dplyr)
###############MAIN DIVERSITY STATISTICS##########################
setwd('~')
path="" # insert path to data directory

breadth_labels <- read.csv(paste(path,'breadth_final.csv',sep="")) %>% dplyr::rename(c("GENE"="gene"))
props=read.csv(paste(path,'props_adj_breadth.txt',sep = ""), sep='\t')
props$X <- gsub("location=", "", props$COORD)
props$X <- sapply(props$X, function(x) substr(x, 0, nchar(x)-3))
props <- props %>% arrange(GENE, desc(TOTAL_CODING_LENGTH))
ags <- read.csv(paste(path,'helb_seropos.txt',sep=""))[[1]]
#props %>% filter(value <= 0.7 | ag) %>% nrow()
# retain only the longest transcript 
props <- props[!duplicated(dplyr::select(props, GENE)),]
# exclude any gene outside of the core region or with all variants filtered at variant quality level 
all_core <- read.csv(paste(path,'all.bed',sep=""),
                     sep="\t",header=FALSE) 
min_expressed_assay <- read.table(paste(path,'min_expression_2.5.txt',sep=""),sep='\n')
sum(!props$GENE %in% all_core$V6)
sum(!props$GENE %in% min_expressed_assay$V1)
qc_variants_filtered <- c("PF3D7_0115000","PF3D7_0201900") # variants not found in any pop VCFs after variant QC filtering step 
props <- props[props$GENE %in% all_core$V6,]
props <- props[props$GENE %in% min_expressed_assay$V1,]
props <- props[!props$GENE %in% qc_variants_filtered,]
nrow(props)
###############PASS RATE CHECKS##########################
# first: full expressed core genome, including antigens 
hist(props$value)
histprops <- ggplot(props) + geom_histogram(aes(x=value)) + geom_vline(xintercept=0.7) + 
  ggtitle('Per-gene SNP call pass rate in Pf7 callset') + 
  xlab('Proportion of passing SNP calls') + ylab('Number of genes')

# pass rate vs coding length
propscorr <- ggplot(props, mapping=aes(x=log10(TOTAL_CODING_LENGTH), y=value)) + 
  geom_point() + geom_smooth(method = "lm") +
  ggtitle('Per-gene pass rate vs. coding length') + 
  xlab(expression(log10("Total Coding Length"))) + ylab('Proportion of passing SNP calls')
cor.test(props$value, props$TOTAL_CODING_LENGTH, method="spearman")
cor.test(filter(props, value>0.7)$value, filter(props, value>0.7)$TOTAL_CODING_LENGTH)

props$ag <- props$GENE %in% ags
stage_labels <- read.csv(paste(path,'m3_combined.csv',sep="")) %>% dplyr::rename(c("GENE"="gene"))
propls <- left_join(props, stage_labels, by = 'GENE') %>% left_join(breadth_labels, by='GENE') %>% filter(value>0.7, !ag)
# check how many genes are excluded by pass rate 
left_join(props, stage_labels, by = 'GENE') %>% left_join(breadth_labels, by='GENE') %>% filter(!ag) %>% 
  filter((breadth==1 | is.na(breadth)), !is.na(stage)) %>% mutate(pass = value > 0.7) %>%
  group_by(pass) %>% summarize(n())

stage_labels50 <- read.csv(paste(path,'stage_specific.csv',sep="")) %>% dplyr::rename(c("GENE"="gene"))
propls50 <- left_join(props, stage_labels50, by = 'GENE') %>% left_join(breadth_labels, by='GENE') %>% filter(value>0.7, !ag)
# gene set sizes before augmentation with additional sporozoite set 
propls50 %>% filter(!ag) %>% 
  filter((breadth==1 | is.na(breadth)), !is.na(stage)) %>% mutate(pass = value > 0.7) %>%
  group_by(stage) %>% summarize(n())
spz_full <- read.delim(file = paste(path,"sporozoite_full.txt",sep=""), 
                       header = FALSE)$V1
propls50$stage[propls50$GENE %in% spz_full] <- 'sporozoite'
propls$stage <- factor(propls$stage, levels=c('sporozoite', 'ring', 'trophozoite', 'schizont', 
                                              'gametocyte', 'ookinete'))
propls50$stage <- factor(propls50$stage, levels=c('sporozoite', 'ring', 'trophozoite', 'schizont', 
                                                  'gametocyte', 'ookinete'))
propls50 %>% group_by(stage) %>% summarize(n())


passrate_by_stagem3 <- ggplot(data = filter(propls, (breadth==1 | is.na(breadth)), !is.na(stage))) + geom_point(aes(y=value, x = stage, fill=stage), position = position_jitterdodge(), alpha=0.5) +  geom_boxplot(alpha=0.5, outlier.color = NA, aes(y=value, x = stage, fill=stage))  + xlab("") + 
  ylab("Proportion of passing SNP calls") + theme_bw() + 
  theme(axis.text = element_text(size=8),   
        axis.title.x=element_text(size=16), axis.title.y=element_text(size=12), 
        axis.text.x = element_text(size=12, angle = 45, hjust=1),
        plot.title = element_text(size=30), panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
        legend.text=element_text(size=15), legend.title=element_text(size=15), legend.position = "none") + 
  scale_fill_viridis(discrete=TRUE, labels=c('sporozoite', 'ring', 'trophozoite', 'schizont', 
                                             'gametocyte', 'ookinete'), name='Life stage') + 
  scale_color_viridis(discrete=TRUE, labels=c('sporozoite', 'ring', 'trophozoite', 'schizont', 
                                              'gametocyte', 'ookinete'), name='Life stage')

passrate_by_stage50 <- ggplot(data = filter(propls50, (breadth==1 | is.na(breadth)), !is.na(stage))) + geom_point(aes(y=value, x = stage, fill=stage), position = position_jitterdodge(), alpha=0.5) +  geom_boxplot(alpha=0.5, outlier.color = NA, aes(y=value, x = stage, fill=stage))  + xlab("") + 
  ylab("Proportion of passing SNP calls") + theme_bw() + 
  theme(axis.text = element_text(size=8),   
        axis.title.x=element_text(size=16), axis.title.y=element_text(size=12), 
        axis.text.x = element_text(size=12, angle = 45, hjust=1),
        plot.title = element_text(size=30), panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
        legend.text=element_text(size=15), legend.title=element_text(size=15), legend.position = "none") + 
  scale_fill_viridis(discrete=TRUE, labels=c('sporozoite', 'ring', 'trophozoite', 'schizont', 
                                             'gametocyte', 'ookinete'), name='Life stage') + 
  scale_color_viridis(discrete=TRUE, labels=c('sporozoite', 'ring', 'trophozoite', 'schizont', 
                                              'gametocyte', 'ookinete'), name='Life stage')

###############FIGURE S5##########################
w = 8
h = 8
rat = h/w
scaled = w/5.5
ggarrange(histprops, propscorr, passrate_by_stagem3, passrate_by_stage50, 
          font.label = list(size=15), 
          labels = c("A", "B", "C", "D"), ncol = 2, nrow = 2)
ggsave(filename=paste(path,"figS5.png",sep=""), dpi=300,
       width = 5.5, height=5.5*rat, scale = scaled, units="in")

cor.test(props$TOTAL_CODING_LENGTH, log10(props$value), method="kendall")
kruskal.test(propls$value, propls$stage,p.adjust.method = "bonf")
kruskal.test(propls50$value, propls50$stage,p.adjust.method = "bonf")


###############FIGURE S13--FFD SITES BY STAGE##########################
ggplot(data = filter(propls50, (breadth==1 | is.na(breadth)), !is.na(stage))) + 
  geom_point(aes(y=PROP_FFD, x = stage, fill=stage), position = position_jitterdodge(), alpha=0.5) +  
  geom_boxplot(alpha=0.5, outlier.color = NA, aes(y=PROP_FFD, x = stage, fill=stage))  + xlab("") + 
  ylab("FFD sites/total sites") + theme_bw() + 
  theme(axis.text = element_text(size=8),   
        axis.title.x=element_text(size=16), axis.title.y=element_text(size=12), 
        axis.text.x = element_text(size=12, angle = 45, hjust=1),
        plot.title = element_text(size=30), panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
        legend.text=element_text(size=15), legend.title=element_text(size=15), legend.position = "none") + 
  scale_fill_viridis(discrete=TRUE, labels=c('sporozoite', 'ring', 'trophozoite', 'schizont', 
                                             'gametocyte', 'ookinete'), name='Life stage') + 
  scale_color_viridis(discrete=TRUE, labels=c('sporozoite', 'ring', 'trophozoite', 'schizont', 
                                              'gametocyte', 'ookinete'), name='Life stage')


w = 8
h = 8
rat = h/w
scaled = w/5.5
ggsave(filename=paste(path,"figS13.png",sep=""), dpi=300,
       width = 5.5, height=5.5*rat, scale = scaled, units="in")

pmat <- as.matrix(summary(glht(model, linfct = mcp(stage = "Tukey")), test = adjusted("none"))$test$pvalues)
propsls_filtered <- filter(propls, (breadth==1 | is.na(breadth)), !is.na(stage))
kruskal.test(propsls_filtered$PROP_FFD, propsls_filtered$stage,p.adjust.method = "bonf")

model <- lm(PROP_FFD~stage+log10(TOTAL_CODING_LENGTH), data=propsls_filtered)
anova(model)
as.matrix(summary(glht(model, linfct = mcp(stage = "Tukey")), 
                          test = adjusted("BH"))$test$pvalues)

