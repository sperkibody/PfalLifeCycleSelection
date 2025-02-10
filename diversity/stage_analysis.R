path = '' # add path to data directory 
fn = paste(path, 'diversity_values.csv', sep='')
diversity <- read.csv(fn)
gs = 'stage_specific'
gs="m3_combined"
fn_gs = paste(path,gs,'.csv',sep='')
fn_out = paste(path,'lsfigs/DivTable_',gs,'.png',sep='')

# add in stage labels 
# LS labels from M3Drop + DE marker ID approach. thus includes / overlaps w/ 50% except at some -- filter post hoc. 
if(gs=='stage_specific'){
  stage_labels <- read.csv(fn_gs) %>% dplyr::rename(c("GENE"="gene"))
  diversity <- full_join(diversity, stage_labels, by="GENE")
  diversity$stage <- factor(diversity$stage, levels=c('sporozoite', 'ring', 'trophozoite', 'schizont', 
                                                      'gametocyte', 'ookinete'))
  # add in sporozoite gold standard set 
  spz_full <- read.delim(file = paste(path,"sporozoite_full.txt",sep=''), 
                         header = FALSE)$V1
  diversity$stage[diversity$GENE %in% spz_full] <- as.factor('sporozoite')
  # set plot y-axis limits for data visualization
  pi_ns_lim <- 0.002
  pi_s_lim <- 0.002
  pnps_lim <- 10
  fn_out_plot = paste(path,'lsfigs/figS9.png',sep='')
} else{
  stage_labels <- read.csv(fn_gs) %>% dplyr::rename(c("GENE"="gene"))
  diversity <- full_join(diversity, stage_labels, by="GENE")
  diversity$stage <- factor(diversity$stage, levels=c('sporozoite', 'ring', 'trophozoite', 'schizont', 
                                                      'gametocyte', 'ookinete'))
  # set plot y-axis limits for data visualization
  pi_ns_lim <- 0.005
  pi_s_lim <- 0.002
  pnps_lim <- 20
  fn_out_plot = paste(path,'lsfigs/fig3.png',sep='')
}
diversity$pop <- factor(diversity$population, levels=c("ghana", "drc", "tanzania", "cambodia"))
pops=c("ghana", "drc", "tanzania", "cambodia")
diversity <- filter(diversity, !is.na(pop))
diversity %>% group_by(stage,population) %>% summarize(n()) %>% View()
###############META-ANALYSIS OF DIVERSITY STATISTICS##########################
# transform data as described 
metapval <- function(stat){
  for(n in 1:length(pops)){
    country <- diversity %>% filter(pop==pops[n], (breadth==1 | is.na(breadth)), !is.na(stage)) #%>% drop_na('pi_ns')
    if(stat=="pi_ns") {formula <- paste("log10(pi_ns+", min(country$pi_ns[country$pi_ns>0]), ") ~ stage + log10(TOTAL_CODING_LENGTH)")}
    if(stat=="pi_s") {formula <- paste("log10(pi_s+", min(country$pi_s[country$pi_s>0]), ") ~ stage + log10(TOTAL_CODING_LENGTH)")}
    if(stat=="pnps") {
      country <- country[is.finite(country$pnps) & !is.na(country$pnps),]
      formula <- paste("log10(pnps+", min(country$pnps[country$pnps>0]), ") ~ stage + log10(TOTAL_CODING_LENGTH)")
      print(formula)}
    else { formula <- paste(stat, "~ stage + log10(TOTAL_CODING_LENGTH)") }
    model <- lm(eval(formula), data = filter(country, !is.na(stage)))
    print(anova(model))
    pmat <- as.matrix(summary(glht(model, linfct = mcp(stage = "Tukey")), test = adjusted("none"))$test$pvalues)
    if(n==1) { ps <- pmat}
    else {ps <- cbind(ps,pmat)}
  }
  return(p.adjust(AWFisher_pvalue(ps)$pvalues, "BH"))
}

res_tajimas = metapval("tajimas_d")
res_pi_ns = metapval("pi_ns")
res_pi_s = metapval("pi_s")
res_pnps = metapval("pnps")
q_all <- unlist(c(res_pi_ns, res_pi_s, res_pnps, res_tajimas))
# example model to extract comparison names 
country <- filter(diversity,pop=="ghana")
pairwise.wilcox.test(country$TOTAL_CODING_LENGTH, country$stage,p.adjust.method = "bonf")
model <- lm(log(pi_ns + 1E-3) ~ stage + log(TOTAL_CODING_LENGTH), data = filter(country, (breadth==1 | is.na(breadth)), !is.na(stage)))
comp_names <- names(summary(glht(model, linfct = mcp(stage = "Tukey")), test = adjusted("none"))$test$pvalues)

# format output table
pdf <- data.frame(matrix(q_all, nrow=15, ncol=4))
colnames(pdf) <- c("pi_ns", "pi_s", "pnps", "D")
rownames(pdf) <- comp_names
setwd("~")
write.csv(pdf,paste(path,"stage_div_comparisons.csv", sep=""))

##########GENERATE PLOTS#####################
library(ggpubr)
cnames <- c('Ghana', 'DRC', 'Tanzania', 'Cambodia')
pi_ns_plot <- ggplot(data = filter(diversity, (breadth==1 | is.na(breadth)), !is.na(stage))) + 
  geom_point(aes(y=log10(pi_ns + 0.001), fill = stage, x=pop), 
             position = position_jitterdodge(dodge.width=0.9), alpha=0.3) +  
  geom_boxplot(alpha=0.7, width = 0.8, outlier.color = NA, 
               aes(y=log10(pi_ns + 0.001), fill = stage, x=pop), 
               position = position_dodge(width = 0.9))  + xlab("") + 
  ylab(expression(pi[NS])) + theme_bw() + 
  theme(axis.text = element_text(size=15), axis.text.x = element_text(size=20),   
        axis.title.x=element_text(size=30), axis.title.y = element_text(size=40), 
        plot.title = element_text(size=90), panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
        legend.text=element_text(size=20), legend.title=element_text(size=15), 
        strip.background = element_blank(), strip.text = element_text(size=30), legend.position = "none") + 
  scale_fill_viridis(discrete=TRUE, labels=c('sporozoite', 'ring', 'trophozoite', 'schizont', 
                                             'gametocyte', 'ookinete'), name='Life stage') + 
  scale_color_viridis(discrete=TRUE, labels=c('sporozoite', 'ring', 'trophozoite', 'schizont', 
                                              'gametocyte', 'ookinete'), name='Life stage') +
  scale_x_discrete(labels=cnames) #+ coord_cartesian(ylim = c(0, pi_ns_lim)) 

pnps_plot <- ggplot(data = filter(diversity, (breadth==1 | is.na(breadth)), !is.na(stage),
                                  is.finite(pnps))) + 
  geom_point(aes(y=log10(pnps + 0.05), fill = stage, x=pop), 
             position = position_jitterdodge(dodge.width=0.9), alpha=0.3) +  
  geom_boxplot(alpha=0.7, width = 0.8, outlier.color = NA, 
               aes(y=log10(pnps + 0.05), fill = stage, x=pop), 
               position = position_dodge(width = 0.9))  + xlab("") + 
  ylab(expression(pi[NS]/pi[S])) + theme_bw() + 
  theme(axis.text = element_text(size=15), axis.text.x = element_text(size=20),   
        axis.title.x=element_text(size=30), axis.title.y = element_text(size=40), 
        plot.title = element_text(size=90), panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
        legend.text=element_text(size=20), legend.title=element_text(size=15), 
        strip.background = element_blank(), strip.text = element_text(size=30), legend.position = "none") + 
  scale_fill_viridis(discrete=TRUE, labels=c('sporozoite', 'ring', 'trophozoite', 'schizont', 
                                             'gametocyte', 'ookinete'), name='Life stage') + 
  scale_color_viridis(discrete=TRUE, labels=c('sporozoite', 'ring', 'trophozoite', 'schizont', 
                                              'gametocyte', 'ookinete'), name='Life stage') +
  scale_x_discrete(labels=cnames) #+ coord_cartesian(ylim = c(0, pnps_lim)) 

d_plot <- ggplot(data = filter(diversity, (breadth==1 | is.na(breadth)), !is.na(stage))) + 
  geom_point(aes(y=tajimas_d, fill = stage, x=pop), 
             position = position_jitterdodge(dodge.width=0.9), alpha=0.3) +  
  geom_boxplot(alpha=0.7, width = 0.8, outlier.color = NA, 
               aes(y=tajimas_d, fill = stage, x=pop), 
               position = position_dodge(width = 0.9))  + xlab("") + 
  ylab(expression(D)) + theme_bw() + 
  theme(axis.text = element_text(size=15), axis.text.x = element_text(size=20),   
        axis.title.x=element_text(size=30), axis.title.y = element_text(size=40), 
        plot.title = element_text(size=90), panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
        legend.text=element_text(size=20), legend.title=element_text(size=15), 
        strip.background = element_blank(), strip.text = element_text(size=30), legend.position = "none") + 
  scale_fill_viridis(discrete=TRUE, labels=c('sporozoite', 'ring', 'trophozoite', 'schizont', 
                                             'gametocyte', 'ookinete'), name='Life stage') + 
  scale_color_viridis(discrete=TRUE, labels=c('sporozoite', 'ring', 'trophozoite', 'schizont', 
                                              'gametocyte', 'ookinete'), name='Life stage') +
  scale_x_discrete(labels=cnames) 

pi_s_plot <- ggplot(data = filter(diversity, (breadth==1 | is.na(breadth)), !is.na(stage))) + 
  geom_point(aes(y=log10(pi_s+0.001), fill = stage, x=pop), 
             position = position_jitterdodge(dodge.width=0.9), alpha=0.3) +  
  geom_boxplot(alpha=0.7, width = 0.8, outlier.color = NA, 
               aes(y=log10(pi_s+0.001), fill = stage, x=pop), 
               position = position_dodge(width = 0.9))  + xlab("") + 
  ylab(expression(pi[S])) + theme_bw() + 
  theme(axis.text = element_text(size=15), axis.text.x = element_text(size=20),   
        axis.title.x=element_text(size=30), axis.title.y = element_text(size=40), 
        plot.title = element_text(size=90), panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
        legend.text=element_text(size=20), legend.title=element_text(size=15), 
        strip.background = element_blank(), strip.text = element_text(size=30), legend.position = "none") + 
  scale_fill_viridis(discrete=TRUE, labels=c('sporozoite', 'ring', 'trophozoite', 'schizont', 
                                             'gametocyte', 'ookinete'), name='Life stage') + 
  scale_color_viridis(discrete=TRUE, labels=c('sporozoite', 'ring', 'trophozoite', 'schizont', 
                                              'gametocyte', 'ookinete'), name='Life stage') +
  scale_x_discrete(labels=cnames) #+ coord_cartesian(ylim = c(-3, pi_s_lim)) 

# count outliers not plotted for each stat 
filter(diversity, (breadth==1 | is.na(breadth)), !is.na(stage)) %>%
  filter(is.finite(pnps)) %>%
  filter(pnps > pnps_lim) %>% nrow()

filter(diversity, (breadth==1 | is.na(breadth)), !is.na(stage)) %>%
  filter(is.infinite(pnps)) %>% nrow()

filter(diversity, (breadth==1 | is.na(breadth)), !is.na(stage)) %>%
  filter(pi_ns > pi_ns_lim) %>% nrow()

filter(diversity, (breadth==1 | is.na(breadth)), !is.na(stage)) %>%
  filter(pi_s > pi_s_lim) %>% nrow()

###############FIGURE 3##########################
w = 12
h = 20
rat = h/w
scaled = w/5.5
scaled=2
ggarrange(pi_ns_plot, pnps_plot, d_plot, pi_s_plot, labels = c("A", "B", "C", "D"), 
          font.label = list(size = 25),
          ncol = 1, nrow = 4)
ggsave(filename=fn_out_plot, 
       dpi=300, width = 5.5, height=5.5*rat, scale = scaled, units="in")

