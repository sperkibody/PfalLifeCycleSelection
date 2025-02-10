library(ggplot2)
library(tidyverse)
library(dplyr)
library(ppcor)
library(PResiduals)

###############################################################################
# BREADTH - DN/DS PREICH
################################################################################
setwd('~')
path = '' # insert path to data directory with diversity output 
props=read.csv(paste(path,'props_adj_breadth.txt',sep=""), sep='\t')
props$X <- gsub("location=", "", props$COORD)
props$X <- sapply(props$X, function(x) substr(x, 0, nchar(x)-3))
# breadth labels from combined approach: 50% + DE marker ID 
breadth_labels <- read.csv(paste(path,'breadth_final.csv',sep="")) %>% dplyr::rename(c("GENE"="gene"))
ags <- read.table(paste(path,'helb_seropos.txt',sep=""))[[1]]

# calculate dN/dS
breadth_partial <- function(species="preich"){
  fn = paste(path,'breadth_divergence_by_gene_',species,".csv",sep="")
  dnds_dat <- read.csv(fn) %>% dplyr::select(-stage)
  # calculate dN and dS from site counts 
  dnds_dat <- dnds_dat %>% mutate(dN1 = ns/ns_sites, dS2 = s/s_sites, dn_ds = dN1/dS2)
  dnds_dat$TRANS <- dnds_dat$Pfal_gene
  dnds_dat <- left_join(dnds_dat, props, by='TRANS') 
  dnds_dat <- left_join(dnds_dat, breadth_labels, by='GENE') 
  dnds_dat <- dnds_dat[!is.na(dnds_dat$breadth),]
  dnds_dat$ag <- dnds_dat$GENE %in% ags 
  dnds_dat <- dnds_dat %>% filter(ag==FALSE)
  mpartial <- partial_Spearman(breadth|log10(dn_ds + 0.001)~log10(TOTAL_CODING_LENGTH), data=dnds_dat)
  print(sum(is.infinite(dnds_dat$dn_ds)))
  mpartial
  return(mpartial)
}

# calculate dN/dS
breadth_plot <- function(species="preich"){
  fn = paste(path,'breadth_divergence_by_gene_',species,".csv",sep="")
  dnds_dat <- read.csv(fn) %>% dplyr::select(-stage)
  # calculate dN and dS from site counts 
  dnds_dat <- dnds_dat %>% mutate(dN1 = ns/ns_sites, dS2 = s/s_sites, dn_ds = dN1/dS2)
  dnds_dat$TRANS <- dnds_dat$Pfal_gene
  dnds_dat <- left_join(dnds_dat, props, by='TRANS') 
  dnds_dat <- left_join(dnds_dat, breadth_labels, by='GENE') 
  dnds_dat <- dnds_dat[!is.na(dnds_dat$breadth),]
  dnds_dat$ag <- dnds_dat$GENE %in% ags 
  dnds_dat <- dnds_dat %>% filter(ag==FALSE)
  print(sum(is.infinite(dnds_dat$dn_ds)))
  # plot 
  dnds_dat$breadth <- sapply(dnds_dat$breadth, function(x) ifelse(x>3, 4, x))
  cnames <- c("1", "2", "3", "4+")
  dnds_dat %>% filter(is.infinite(dn_ds)) %>% nrow()
  pdnds <- ggplot(data = filter(dnds_dat, !is.na(breadth), !is.infinite(dnds_dat$dn_ds)), aes(y=dn_ds, x = as.factor(breadth))) + geom_jitter(fill="grey",alpha=0.35,pch=21) + geom_boxplot(fill=NA, outlier.shape = NA, color='black') +  
    xlab('Number of Life Stage Associations') + 
    ylab(expression(dN/dS)) + theme_bw() + geom_smooth(method = 'lm') + ylim(0,1.25) + 
    theme(axis.text = element_text(size=12), axis.text.x = element_text(vjust = 0.5, hjust=1),   
          axis.title.x=element_text(size=16), axis.title.y=element_text(size=30), 
          plot.title = element_text(size=70), panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
          legend.text=element_text(size=15), legend.title=element_text(size=15), 
          strip.background = element_blank(), strip.text = element_text(size=25), legend.position="none") + 
    scale_x_discrete(labels=cnames) 
  print(pdnds)
  print('not shown: ')
  print(nrow(filter(dnds_dat, !is.na(breadth), is.finite(dn_ds), dn_ds > 1.25)))
  return(pdnds)
}

# P. reich results 
breadth_partial()
breadth_plot()
# P. praefal results 
breadth_partial("praefal")
breadth_plot("praefal")

###############################################################################
# READ IN DIV. DF
################################################################################
diversity_breadth <- read.csv(paste(path,'diversity_values.csv',sep=""))
diversity_breadth <- diversity_breadth[!is.na(diversity_breadth$breadth),]
diversity_breadth$pop <- as.factor(diversity_breadth$population)
diversity_breadth_overall <- diversity_breadth
###############################################################################
# BREADTH - pN/pS 
################################################################################
pops = c('cambodia', 'drc', 'tanzania', 'ghana')
diversity_breadth <- diversity_breadth %>% filter(!is.na(pnps))
diversity_breadth <- diversity_breadth_overall %>% filter(pop==pops[1])
mpartial1 <- partial_Spearman(breadth|log10(pnps+0.001)~log10(TOTAL_CODING_LENGTH), data=diversity_breadth)
diversity_breadth <- diversity_breadth_overall %>% filter(pop==pops[2])
mpartial2 <- partial_Spearman(breadth|log10(pnps+0.001)~log10(TOTAL_CODING_LENGTH), data=diversity_breadth)
diversity_breadth <- diversity_breadth_overall %>% filter(pop==pops[3])
mpartial3 <- partial_Spearman(breadth|log10(pnps+0.001)~log10(TOTAL_CODING_LENGTH), data=diversity_breadth)
diversity_breadth <- diversity_breadth_overall %>% filter(pop==pops[4])
mpartial4 <- partial_Spearman(breadth|log10(pnps + 0.001)~log10(TOTAL_CODING_LENGTH), data=diversity_breadth)

# print partial spearman results 
mpartial1
mpartial2
mpartial3
mpartial4

diversity_breadth_overall$breadth <- sapply(diversity_breadth_overall$breadth, function(x) ifelse(x>3, 4, x))
cnames <- c("1", "2", "3", "4+")
filter(diversity_breadth_overall, !is.na(breadth)) %>% filter(is.finite(pnps), pnps > 10) %>% nrow()
filter(diversity_breadth_overall, !is.na(breadth)) %>% filter(is.infinite(pnps)) %>% nrow()

pdiv <- ggplot(filter(diversity_breadth_overall, !is.na(breadth), !is.na(pop), is.finite(pnps)), aes(y=pnps, x = as.factor(breadth))) + geom_jitter(fill="grey",alpha=0.35,pch=21) + geom_boxplot(fill=NA, outlier.shape = NA, color='black') +  
  xlab('Number of Life Stage Associations') + facet_grid(~pop, 
                                                         labeller = labeller(pop = c("drc" = "DRC", "ghana" = "Ghana", "tanzania" = "Tanzania", "cambodia" = "Cambodia"))) +
  ylab(expression(pi[NS]/pi[S])) + theme_bw() + geom_smooth(method = 'lm') + coord_cartesian(ylim = c(0, 10)) + 
  theme(axis.text = element_text(size=12), axis.text.x = element_text(),   
        axis.title.x=element_text(size=16), axis.title.y=element_text(size=30), 
        plot.title = element_text(size=70), panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
        legend.text=element_text(size=15), legend.title=element_text(size=15), 
        strip.background = element_blank(), strip.text = element_text(size=25)) + scale_x_discrete(labels=cnames)


###############FIGURE S6##########################
w = 8
h = 10
rat = h/w
scaled = w/5.5

# add in pdnds w/ praefal as b 
ggarrange(breadth_plot(), breadth_plot("praefal"), pdiv, labels = c("A", "B", "C"), 
          font.label = list(size = 25), 
          ncol = 1, nrow = 3)
ggsave(filename=paste(path,"lsfigs/figS8.png",sep=""), 
       dpi=300, width = 5.5, height=5.5*rat, scale = scaled, units="in")

# save breadth above 3 genes as file to input into shinyGO for enrichment analysis 
breadth_above3 <- breadth_labels %>% filter(breadth > 3) %>% dplyr::select(GENE) 
write.table(breadth_above3, paste(path,"breadth_above3.txt",sep=""), 
            row.names = FALSE, sep = '\n', quote=FALSE, col.names = FALSE)
