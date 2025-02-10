library(ggplot2)
library(dplyr)
###############GENOME-WIDE CODING LENGTH CORRELATIONS for non-antigens##########################
###############DIVERSITY DATA##########################
path = '' # add path to data directory
fn = paste(path, 'diversity_values.csv', sep='')
diversity <- read.csv(fn)
# filter to just Ghana data 
country <- filter(diversity, pop=="ghana")
pi_ns_shift <- round(min(filter(country, pi_ns > 0)$pi_ns),5)
pi_s_shift <- round(min(filter(country, pi_s > 0)$pi_s), 5)
pnps_shift <- round(min(filter(country, pnps > 0)$pnps),3)
##### MAKE PLOTS 
pn_corr <- ggplot(country,mapping = 
                   aes(x=log10(TOTAL_CODING_LENGTH),y=log10(pi_ns+pi_ns_shift))) + 
  geom_point() + theme_bw() +  
  theme(axis.text = element_text(size=15),   
        axis.title.x=element_text(size=16), axis.title.y=element_text(size=20), 
        plot.title = element_text(size=30), panel.border = element_blank(), 
        axis.line = element_line(colour = "black")) + xlab(expression(log10("Total Coding Length"))) + 
  ylab(expression(log10(pi[NS] + 1E-5)))
ps_corr <- ggplot(country,mapping = 
                    aes(x=log10(TOTAL_CODING_LENGTH),y=log10(pi_s+pi_s_shift))) + 
  geom_point() + theme_bw() +  
  theme(axis.text = element_text(size=15),   
        axis.title.x=element_text(size=16), axis.title.y=element_text(size=20), 
        plot.title = element_text(size=30), panel.border = element_blank(), 
        axis.line = element_line(colour = "black")) + xlab(expression(log10("Total Coding Length"))) + 
  ylab(expression(log10(pi[S] + 1E-5)))
pnps_corr <- ggplot(country,mapping = 
                    aes(x=log10(TOTAL_CODING_LENGTH),y=log10(pnps+pnps_shift))) + 
  geom_point() + theme_bw() +
  theme(axis.text = element_text(size=15),   
        axis.title.x=element_text(size=16), axis.title.y=element_text(size=20), 
        plot.title = element_text(size=30), panel.border = element_blank(), 
        axis.line = element_line(colour = "black")) + xlab(expression(log10("Total Coding Length"))) + 
  ylab(expression(log10(pi[NS]/pi[S] + 1E-3)))
d_corr <- ggplot(country,mapping = 
                   aes(x=log10(TOTAL_CODING_LENGTH),y=tajimas_d)) + 
  geom_point() + theme_bw() + 
  theme(axis.text = element_text(size=15),   
        axis.title.x=element_text(size=16), axis.title.y=element_text(size=20), 
        plot.title = element_text(size=30), panel.border = element_blank(), 
        axis.line = element_line(colour = "black")) + xlab(expression(log10("Total Coding Length"))) + 
  ylab(expression(D))

###############FIGURE S4-supp corr 2##########################
w = 16
h=9
rat = h/w
scaled = w/5.5
ggarrange(pn_corr, ps_corr, pnps_corr, d_corr, 
          font.label = list(size = 25),
          labels = c("A", "B", "C", "D"), 
          ncol = 2, nrow = 2)
ggsave(filename=paste(path,"lsfigs/figS5.png",sep=''), 
       dpi=300, width = 5.5, height=5.5, scale = 2.5, units="in")

# test genome-wide correlations in Ghana data via Kendall's tau 
nrow(country)
cor.test(log10(as.numeric(country$TOTAL_CODING_LENGTH)), as.numeric(country$tajimas_d), 
         method = 'kendall')
cor.test(as.numeric(country$TOTAL_CODING_LENGTH), as.numeric(country$pnps), 
         method = 'kendall')
cor.test(as.numeric(country$TOTAL_CODING_LENGTH), as.numeric(country$pi_ns), 
         method = 'kendall')
cor.test(as.numeric(country$TOTAL_CODING_LENGTH), as.numeric(country$pi_s), 
         method = 'kendall')
cor.test(as.numeric(country$TOTAL_CODING_LENGTH), as.numeric(country$breadth), 
         method = 'kendall')
###############DIVERGENCE DATA##########################
setwd('~')
fn_props = paste(path, 'props_adj_breadth.txt', sep='')
props=read.csv(fn_props, sep='\t')
props$X <- gsub("location=", "", props$COORD)
props$X <- sapply(props$X, function(x) substr(x, 0, nchar(x)-3))
# breadth labels from combined approach: 50% + DE marker ID 
breadth_labels <- read.csv(paste(path,'de50_combined.csv',sep="")) %>% dplyr::rename(c("GENE"="value", "breadth"="variable"))
ags <- read.table(paste(path,'helb_seropos.txt',sep=""))[[1]]
# calculate dN/dS
fn = paste(path,'breadth_divergence_by_gene_preich.csv',sep="")
dnds_dat <- read.csv(fn) %>% dplyr::select(-stage)
# calculate dN and dS from site counts 
dnds_dat <- dnds_dat %>% mutate(dN1 = ns/ns_sites, dS2 = s/s_sites, dn_ds = dN1/dS2)
dnds_dat$TRANS <- dnds_dat$Pfal_gene
dnds_dat <- left_join(dnds_dat, props, by='TRANS') 
dnds_dat <- left_join(dnds_dat, breadth_labels, by='GENE') 
#dnds_dat <- dnds_dat[!is.na(dnds_dat$breadth),]
dnds_dat$ag <- dnds_dat$GENE %in% ags 
dnds_dat <- dnds_dat %>% filter(ag==FALSE)

min(filter(dnds_dat, dN1 > 0)$dN1)
dn_shift <- 1E-4
min(filter(dnds_dat, dS2 > 0)$dS2)
ds_shift <- 1E-3
min(filter(dnds_dat, dn_ds > 0)$dn_ds)
dnds_shift <- 1E-3
# genome-wide correlations 
dn_corr <- ggplot(dnds_dat,mapping = aes(x=log10(TOTAL_CODING_LENGTH),y=log10(dN1 + dn_shift))) + geom_point() + 
  theme_bw() + 
  theme(axis.text = element_text(size=15),   
        axis.title.x=element_text(size=16), axis.title.y=element_text(size=20), 
        plot.title = element_text(size=30), panel.border = element_blank(), 
        axis.line = element_line(colour = "black")) + xlab(expression(log10("Total Coding Length"))) + 
  ylab(expression(log10(dN + 1E-4)))
ds_corr <- ggplot(dnds_dat,mapping = aes(x=log10(TOTAL_CODING_LENGTH),y=log10(dS2 + ds_shift))) + geom_point() + 
  theme_bw() + 
  theme(axis.text = element_text(size=15),   
        axis.title.x=element_text(size=16), axis.title.y=element_text(size=20), 
        plot.title = element_text(size=30), panel.border = element_blank(), 
        axis.line = element_line(colour = "black")) + xlab(expression(log10("Total Coding Length"))) + 
  ylab(expression(log10(dS + 1E-3)))
dnds_corr <- ggplot(dnds_dat,mapping = aes(x=log10(TOTAL_CODING_LENGTH),y=log10(dn_ds + dnds_shift))) + geom_point() + 
  theme_bw() + 
  theme(axis.text = element_text(size=15),   
        axis.title.x=element_text(size=16), axis.title.y=element_text(size=20), 
        plot.title = element_text(size=30), panel.border = element_blank(), 
        axis.line = element_line(colour = "black")) + xlab(expression(log10("Total Coding Length"))) + 
  ylab(expression(log10(dN/dS + 1E-3)))
propns_corr <- ggplot(props,mapping = aes(x=log10(TOTAL_CODING_LENGTH),y=PROP_NS)) + geom_point() + 
  theme_bw() + 
  theme(axis.text = element_text(size=15),   
        axis.title.x=element_text(size=16), axis.title.y=element_text(size=20), 
        plot.title = element_text(size=30), panel.border = element_blank(), 
        axis.line = element_line(colour = "black")) + xlab(expression(log10("Total Coding Length"))) + 
  ylab(expression("NS sites"/"total sites"))
props_corr <- ggplot(props,mapping = aes(x=log10(TOTAL_CODING_LENGTH),y=PROP_SYN)) + geom_point() + 
  theme_bw() + 
  theme(axis.text = element_text(size=15),   
        axis.title.x=element_text(size=16), axis.title.y=element_text(size=20), 
        plot.title = element_text(size=30), panel.border = element_blank(), 
        axis.line = element_line(colour = "black")) + xlab(expression(log10("Total Coding Length"))) + 
  ylab(expression("S sites"/"total sites"))
propffd_corr <- ggplot(props,mapping = aes(x=log10(TOTAL_CODING_LENGTH),y=PROP_FFD)) + geom_point() + 
  theme_bw() + 
  theme(axis.text = element_text(size=15),   
        axis.title.x=element_text(size=16), axis.title.y=element_text(size=20), 
        plot.title = element_text(size=30), panel.border = element_blank(), 
        axis.line = element_line(colour = "black")) + xlab(expression(log10("Total Coding Length"))) + 
  ylab(expression("FFD sites"/"total sites"))

###############FIGURE S4-supp corr 1##########################
w = 16
h=9
rat = h/w
scaled = w/5.5
ggarrange(dn_corr, ds_corr, dnds_corr, propns_corr, props_corr, propffd_corr,
          font.label = list(size = 25),
          labels = c("A", "B", "C", "D", "E", "F"), ncol = 2, nrow = 3)
ggsave(filename=paste(path,"figS4.png",sep=''), 
       dpi=300, width = 5.5, height=8.25, scale = 2.5, units="in")

### examine correlations
cor.test(as.numeric(dnds_dat$TOTAL_CODING_LENGTH), as.numeric(dnds_dat$dn_ds), 
         method = 'kendall')
cor.test(as.numeric(dnds_dat$TOTAL_CODING_LENGTH), as.numeric(dnds_dat$dN1), 
         method = 'kendall')
cor.test(as.numeric(dnds_dat$TOTAL_CODING_LENGTH), as.numeric(dnds_dat$dS2), 
         method = 'kendall')
cor.test(props$PROP_NS,props$TOTAL_CODING_LENGTH,method="kendall")
cor.test(props$PROP_SYN,props$TOTAL_CODING_LENGTH,method="kendall")
cor.test(props$PROP_FFD,props$TOTAL_CODING_LENGTH,method="kendall")
nrow(dnds_dat)
