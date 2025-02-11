library(reshape2)
library(ggpubr)
library(viridis)
library(ggplot2)
library(dplyr)
#### READ IN DF AND ADD STAGE LABELS 
# add path to data directory 
path=""
fn = paste(path, 'diversity_values.csv', sep='')
diversity <- read.csv(fn)
gs="m3_combined"
fn_gs = paste(path,gs,'.csv',sep='')
fn_out = paste(path,'lsfigs/FstTable2.png',sep='')

# add in stage labels 
stage_labels <- read.csv(fn_gs) %>% dplyr::rename(c("GENE"="gene"))
diversity <- full_join(diversity, stage_labels, by="GENE")
diversity$stage <- factor(diversity$stage, levels=c('sporozoite', 'ring', 'trophozoite', 'schizont', 
                                                      'gametocyte', 'ookinete'))

fn_out_plot = paste(path,'lsfigs/fig4.png',sep='')

diversity$pop <- factor(diversity$population, levels=c("ghana", "drc", "tanzania", "cambodia"))
pops=c("ghana", "drc", "tanzania", "cambodia")

fst_all <- filter(diversity, (breadth==1 | is.na(breadth)), !is.na(stage)) %>% 
  dplyr::select(COORD,  contains("fst"), stage)
fst_all <- fst_all %>% melt(variable.name='comparison', value.name = 'Fst') 
fst_all$comparison <- recode_factor(fst_all$comparison, fst.drc.ghana='DRC/Ghana', 
                                    fst.drc.tanzania='DRC/Tanzania', fst.ghana.tanzania='Ghana/Tanzania', 
                                    fst.ghana.cambodia='Cambodia/Ghana',  fst.drc.cambodia='Cambodia/DRC', 
                                    fst.tanzania.cambodia='Cambodia/Tanzania')
fst_all <- left_join(fst_all, dplyr::select(diversity, c("COORD","TOTAL_CODING_LENGTH")))
fst_all <- fst_all[!duplicated(fst_all),]
fst_all[duplicated(fst_all),] %>% View()
# ensure negative values have been replaced with 0 (should have been already)
sum(fst_all$Fst<0,na.rm = TRUE)
fst_all$Fst[fst_all$Fst<0] <- 0
sum(fst_all$Fst<0,na.rm = TRUE)
fst_all_sub <- fst_all %>% filter(!grepl('Cambodia', comparison))
fst_all_cam <- fst_all %>% filter(grepl('Cambodia', comparison))
#### PLOTS MATCHING OTHER DIV PLOTS FOR FST 
fst_all %>% group_by(comparison, stage) %>% summarize(n()) %>% View()

fst1 <- ggplot(data = fst_all_sub, aes(y=Fst, x = comparison, fill=stage)) + 
  geom_point(position = position_jitterdodge(), alpha=0.3) +  
  geom_boxplot(alpha=0.7, outlier.color = NA) + xlab('') + 
  ylab(expression(F[ST])) + theme_bw() + 
  theme(axis.text = element_text(size=15), axis.text.x = element_text(size=17),   
        axis.title.x=element_text(size=30), axis.title.y=element_text(size=30), 
        plot.title = element_text(size=90), panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
        legend.text=element_text(size=20), legend.title=element_text(size=15), 
        strip.background = element_blank(), strip.text = element_text(size=30), legend.position = "none") + 
  scale_fill_viridis(discrete=TRUE, labels=c('sporozoite', 'ring', 'trophozoite', 'schizont', 
                                             'gametocyte', 'ookinete'), name='Life stage') + 
  scale_color_viridis(discrete=TRUE, labels=c('sporozoite', 'ring', 'trophozoite', 'schizont', 
                                              'gametocyte', 'ookinete'), name='Life stage') + coord_cartesian(ylim = c(0, 0.14))

fst2 <- ggplot(data = fst_all_cam, aes(y=Fst, x = comparison, fill=stage)) + 
  geom_point(position = position_jitterdodge(), alpha=0.3) +  
  geom_boxplot(alpha=0.7, outlier.color = NA) + xlab('') + 
  ylab(expression(F[ST])) + theme_bw() + 
  theme(axis.text = element_text(size=15), axis.text.x = element_text(size=17),   
        axis.title.x=element_text(size=30), axis.title.y=element_text(size=30), 
        plot.title = element_text(size=90), panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
        legend.text=element_text(size=20), legend.title=element_text(size=15), 
        strip.background = element_blank(), strip.text = element_text(size=30), legend.position = "none") + 
  scale_fill_viridis(discrete=TRUE, labels=c('sporozoite', 'ring', 'trophozoite', 'schizont', 
                                             'gametocyte', 'ookinete'), name='Life stage') + 
  scale_color_viridis(discrete=TRUE, labels=c('sporozoite', 'ring', 'trophozoite', 'schizont', 
                                              'gametocyte', 'ookinete'), name='Life stage') + 
  coord_cartesian(ylim = c(0, 0.87))

fst_all_sub %>% filter(Fst > 0.14) %>% nrow()
fst_all_cam %>% filter(Fst > 0.87) %>% nrow()
###############FIGURE 4##########################
w = 8
h = 10
rat = h/w
scaled = w/5.5
ggarrange(fst1, fst2, labels = c("A", "B"),  font.label = list(size = 25), ncol = 1, nrow = 2)
ggsave(filename=fn_out_plot, 
       dpi=300, width = 5.5, height=5.5*rat, scale = scaled, units="in")

### Fst statistical testing 
comparisons <- levels(fst_all$comparison)

metapvalfst <- function(){
  for(n in 1:length(comparisons)){
    country <- fst_all %>% filter(comparison==comparisons[n]) 
    formula <- "Fst ~ stage + log10(TOTAL_CODING_LENGTH)"
    model <- lm(eval(formula), data = filter(country, !is.na(stage)))
    print(anova(model))
    tukey_out <- summary(glht(model, linfct = mcp(stage = "Tukey")), test = adjusted("none"))$test$pvalues
    pmat <- as.matrix(tukey_out)
    comp_names <- names(tukey_out)
    if(n==1) { ps <- pmat}
    else {ps <- cbind(ps,pmat)}
  }
  return(list(p.adjust(AWFisher_pvalue(ps)$pvalues, "BH"), comp_names))
}

res_fst = metapvalfst()
comp_names = res_fst[[2]]
res_fst = res_fst[[1]]

### FORMAT INTO TABLE 
pdf <- data.frame(matrix(res_fst, nrow=15, ncol=1))
colnames(pdf) <- c("FST")
rownames(pdf) <- comp_names
fn_fst <- paste(path,"fst_results.csv",sep="")
write.csv(pdf,fn_fst)

