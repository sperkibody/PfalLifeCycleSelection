library(ggplot2)
library(tidyverse)
library(dplyr)
library(ppcor)
library(PResiduals)
library(multcomp)

###############################################################################
# DN/DS BY STAGE 
################################################################################
setwd('~')
path = ''
gs="m3_combined"
fn_gs = paste(path,gs,'.csv',sep='')
fn_out = paste(path,'lsfigs/fig2.png',sep='')
props=read.csv(paste(path,'props_adj_breadth.txt',sep=''), sep='\t')
props$X <- gsub("location=", "", props$COORD)
props$X <- sapply(props$X, function(x) substr(x, 0, nchar(x)-3))
# breadth labels from combined approach: tebben 50% + DE marker ID 
breadth_labels <- read.csv(paste(path,'breadth_final.csv',sep='')) %>% dplyr::rename(c("GENE"="gene"))
ags <- read.table(paste(path,'helb_seropos.txt',sep=''))[[1]]
stage_labels <- read.csv(fn_gs) %>% dplyr::rename(c("GENE"="gene"))
# get stage-specific dN/dS data
get_dnds_stage_df <- function(species="preich"){
  fn = paste(path,'breadth_divergence_by_gene_',species,".csv",sep="")
  dnds_dat <- read.csv(fn) %>% dplyr::select(-stage)
  # calculate dN and dS from site counts 
  dnds_dat <- dnds_dat %>% mutate(dN1 = ns/ns_sites, dS2 = s/s_sites, dn_ds = dN1/dS2)
  dnds_dat$TRANS <- dnds_dat$Pfal_gene
  dnds_dat <- left_join(dnds_dat, props, by='TRANS') 
  dnds_dat <- left_join(dnds_dat, breadth_labels, by='GENE') 
  dnds_dat$ag <- dnds_dat$GENE %in% ags 
  dnds_dat <- dnds_dat %>% filter(ag==FALSE)
  dnds_stage <- left_join(dnds_dat, stage_labels, by='GENE') %>% filter(!is.na(stage), (breadth < 2 | is.na(breadth)))
  dnds_stage$stage <- factor(dnds_stage$stage, levels=c('sporozoite', 'ring', 'trophozoite', 'schizont', 
                                                                            'gametocyte', 'ookinete'))
  return(dnds_stage)
}

# P. reich results 
dnds_stage <- get_dnds_stage_df()
# get small values to transform data above 0 for log-transformation
dn_shift <- round(min(filter(dnds_stage,dN1>0)$dN1),3)
ds_shift <- round(min(filter(dnds_stage,dS2>0)$dS2),2)
dnds_shift <- round(min(filter(dnds_stage,dn_ds>0)$dn_ds),2)
### isolate effects of dN and dS rates 
model <- lm(log10(dN1 + dn_shift)~stage + log10(TOTAL_CODING_LENGTH), data = dnds_stage)
print(anova(model))
pmat <- as.matrix(summary(glht(model, linfct = mcp(stage = "Tukey")), test = adjusted("BH"))$test$pvalues)
ps <- pmat 
colnames(ps) <- "p.value"
comp_names <- names(summary(glht(model, linfct = mcp(stage = "Tukey")), test = adjusted("BH"))$test$pvalues)
group1 <- unname(sapply(comp_names, function(x) strsplit(x," - ")[[1]][1]))
group2 <- unname(sapply(comp_names, function(x) strsplit(x," - ")[[1]][2]))
a <- data.frame(group1, group2, p.value=ps) %>% filter(p.value < 0.05)
a$p.signif <- sapply(a$p.value, function(x) ifelse(x < 0.05, ifelse(x < 0.01, ifelse(x < 0.001, '***', '**'), '*'), ' '))
a$group2 <- factor(a$group2, levels=c('sporozoite', 'ring', 'trophozoite', 'schizont', 
                                      'gametocyte', 'ookinete'))
#a[5,] <- a[5,c(2,1,3,4)]
#a <- arrange(a,by_group=group2)
a$y.position <- seq(0.06, 0.08, length.out=nrow(a))
pdn <- ggplot(dnds_stage, mapping = aes(y=dN1, x = stage)) + 
  geom_jitter(alpha=0.3) + 
  geom_boxplot(alpha=0.7, outlier.color = NA, mapping = aes(y=dN1, x = stage, fill=stage)) + 
  xlab('Life Stage') + 
  ylab(expression(dN)) + theme_bw() + 
  theme(axis.text = element_text(size=12), axis.text.x = element_text(size=10),  
        axis.title.x=element_text(size=16), axis.title.y=element_text(size=30), 
        plot.title = element_text(size=30), panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
        legend.position = "none") + 
  scale_fill_viridis(discrete=TRUE, labels=c('sporozoite', 'ring', 'trophozoite', 'schizont', 
                                             'gametocyte', 'ookinete'), name='Life stage') + 
  scale_color_viridis(discrete=TRUE, labels=c('sporozoite', 'ring', 'trophozoite', 'schizont', 
                                              'gametocyte', 'ookinete'), name='Life stage') 

pdn <- pdn + stat_pvalue_manual(data.frame(a), label = "p.signif", size=5, 
                                bracket.size = 0.5, bracket.nudge.y = 0, 
                                tip.length = 0.01)

model <- lm(log10(dS2 + ds_shift)~stage + log10(TOTAL_CODING_LENGTH), data = dnds_stage)
print(anova(model))
pmat <- as.matrix(summary(glht(model, linfct = mcp(stage = "Tukey")), test = adjusted("BH"))$test$pvalues)
sum(pmat[,1]<0.05)
# no significant differences in dS 
pds <- ggplot(dnds_stage, mapping = aes(y=dS2, x = stage, fill=stage)) + 
  geom_jitter(alpha=0.3) + 
  geom_boxplot(alpha=0.7, outlier.color = NA) + 
  xlab('Life Stage') + 
  ylab(expression(dS)) + theme_bw() + 
  theme(axis.text = element_text(size=12), axis.text.x = element_text(size=10),  
        axis.title.x=element_text(size=16), axis.title.y=element_text(size=30), 
        plot.title = element_text(size=30), panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
        legend.position = "none") + 
  scale_fill_viridis(discrete=TRUE, labels=c('sporozoite', 'ring', 'trophozoite', 'schizont', 
                                             'gametocyte', 'ookinete'), name='Life stage') + 
  scale_color_viridis(discrete=TRUE, labels=c('sporozoite', 'ring', 'trophozoite', 'schizont', 
                                              'gametocyte', 'ookinete'), name='Life stage') 


model <- lm(log10(dn_ds + dnds_shift)~stage + log10(TOTAL_CODING_LENGTH), data = filter(dnds_stage, !is.infinite(dnds_stage$dn_ds)))
print(anova(model))
pmat <- as.matrix(summary(glht(model, linfct = mcp(stage = "Tukey")), test = adjusted("BH"))$test$pvalues)
ps <- pmat 
colnames(ps) <- "p.value"
comp_names <- names(summary(glht(model, linfct = mcp(stage = "Tukey")), test = adjusted("BH"))$test$pvalues)
group1 <- unname(sapply(comp_names, function(x) strsplit(x," - ")[[1]][1]))
group2 <- unname(sapply(comp_names, function(x) strsplit(x," - ")[[1]][2]))
a <- data.frame(group1, group2, p.value=ps) %>% filter(p.value < 0.05)
a$p.signif <- sapply(a$p.value, function(x) ifelse(x < 0.05, ifelse(x < 0.01, ifelse(x < 0.001, '***', '**'), '*'), ' '))
a$group2 <- factor(a$group2, levels=c('sporozoite', 'ring', 'trophozoite', 'schizont', 
                                      'gametocyte', 'ookinete'))
a$y.position <- seq(0.82, 1.2, length.out=nrow(a))
dnds_stage$stage <- factor(dnds_stage$stage, levels=c('sporozoite', 'ring', 'trophozoite', 'schizont', 
                                                                  'gametocyte', 'ookinete'))
pstage <- ggplot(filter(dnds_stage, !is.infinite(dnds_stage$dn_ds)), 
                 mapping = aes(y=dn_ds, x = stage)) + 
  geom_jitter(alpha=0.3) + 
  geom_boxplot(mapping = aes(y=dn_ds, x = stage, fill=stage), 
               alpha=0.7, outlier.color = NA) + 
  xlab('Life Stage') + 
  ylab(expression(dN/dS)) + theme_bw() + 
  theme(axis.text = element_text(size=12), axis.text.x = element_text(size=10),  
        axis.title.x=element_text(size=16), axis.title.y=element_text(size=30), 
        plot.title = element_text(size=30), panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
        legend.position = "none") + 
  scale_fill_viridis(discrete=TRUE, labels=c('sporozoite', 'ring', 'trophozoite', 'schizont', 
                                             'gametocyte', 'ookinete'), name='Life stage') + 
  scale_color_viridis(discrete=TRUE, labels=c('sporozoite', 'ring', 'trophozoite', 'schizont', 
                                              'gametocyte', 'ookinete'), name='Life stage') + coord_cartesian(ylim = c(0,1.25)) #+ ggtitle(title)
pstage <- pstage + stat_pvalue_manual(data.frame(a), label = "p.signif", size=5, 
                                      bracket.size = 0.5, bracket.nudge.y = 0, 
                                      tip.length = 0.01)
pstage

# examine missing and Inf value
dnds_stage[dnds_stage$dn_ds > 1.25,]

# examine medians 
m1 <- summary(dnds_stage$dn_ds[dnds_stage$stage=="sporozoite"])[3]
m2 <- summary(dnds_stage$dn_ds[dnds_stage$stage!="sporozoite"])[3]
m1
m2
m1/m2
###############FIGURE 2##########################
w = 8
h = 10
rat = h/w
scaled = w/5.5
ggarrange(pdn, pstage, pds, labels = c("A", "B", "C"), 
          font.label = list(size = 20),
          ncol = 1, nrow = 3)
ggsave(filename=fn_out, dpi=300, width = 5.5, height=5.5*rat, scale = scaled, 
       units="in")

