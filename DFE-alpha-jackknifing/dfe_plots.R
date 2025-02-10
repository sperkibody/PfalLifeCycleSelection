##### DFE-ALPHA RESULTS (JACKKNIFING)
setwd('~')
path = ""

read_dfe <- function(run_id,species="preich"){
  fn = paste(path, run_id, sep="")
  dfea <- read.csv(fn)
  # ensure including jackknife runs excluding genes with divergence or diversity data influencing dfe-alpha output
  fn = paste(path, 'diversity_values.csv', sep='')
  diversity <- read.csv(fn)
  fn = paste(path,'breadth_divergence_by_gene_',species,".csv",sep="")
  serre_dnds <- read.csv(fn) %>% dplyr::select(-stage)
  genes_found <- unique(c(diversity$GENE,serre_dnds$Gene.ID))
  # only keep runs excluding genes in the dfs
  dfea <- dfea[dfea$excluded_gene %in% genes_found,]
  # significance test: pairwise wilcoxon rank sum test
  dfea$stage <- factor(dfea$stage, levels=c('sporozoite', 'ring', 'trophozoite', 'schizont', 
                                            'gametocyte', 'ookinete'))
  return(dfea)
}

# print summary stats
dfe_stats <- function(dfea_ghana){
  print(kruskal.test(alpha~stage, data=dfea_ghana))
  print(pairwise.wilcox.test(dfea_ghana$alpha, dfea_ghana$stage,p.adjust.method = "bonf"))
  print("sporozoite medians (omega, Q, alpha): ")
  print(c(summary(filter(dfea_ghana,stage=='sporozoite')$omega)[3], 
          summary(filter(dfea_ghana,stage=='sporozoite')$Q)[3], 
          summary(filter(dfea_ghana,stage=='sporozoite')$alpha)[3]))
  print("non-sporozoite medians (omega, Q, alpha): ")
  print(c(summary(filter(dfea_ghana,stage!='sporozoite')$omega)[3], 
          summary(filter(dfea_ghana,stage!='sporozoite')$Q)[3], 
          summary(filter(dfea_ghana,stage!='sporozoite')$alpha)[3]))
  print("blood stages medians (omega, Q, alpha): ")
  mosq <- c("sporozoite", "ookinete","gametocyte")
  print(c(summary(filter(dfea_ghana,!(stage %in% mosq))$omega)[3], 
          summary(filter(dfea_ghana,!(stage %in% mosq))$Q)[3], 
          summary(filter(dfea_ghana,!(stage %in% mosq))$alpha)[3]))
  print("ookinete medians (omega, Q, alpha): ")
  print(c(summary(filter(dfea_ghana,stage=='ookinete')$omega)[3], 
          summary(filter(dfea_ghana,stage=='ookinete')$Q)[3], 
          summary(filter(dfea_ghana,stage=='ookinete')$alpha)[3]))
  print("Q gam / ook medians (Q):")
  print(c(summary(filter(dfea_ghana,stage!='gametocyte')$Q)[3], 
          summary(filter(dfea_ghana,stage!='ookinete')$Q)[3]))
  print("spz to blood stage Q ratio: ")
  print(c(summary(filter(dfea_ghana,stage=='sporozoite')$Q)[3]/summary(filter(dfea_ghana,!(stage %in% mosq))$Q)[3]))
  print("gam to blood stage Q ratio: ")
  print(c(summary(filter(dfea_ghana,stage=='gametocyte')$Q)[3]/summary(filter(dfea_ghana,!(stage %in% mosq))$Q)[3]))
  print("ook to blood stage Q ratio: ")
  print(c(summary(filter(dfea_ghana,stage=='ookinete')$Q)[3]/summary(filter(dfea_ghana,!(stage %in% mosq))$Q)[3]))
  print("trophozoite median alpha: ")
  print(summary(filter(dfea_ghana,stage=='trophozoite')$alpha[3]))
  }

# plotting 
dfe_plots <- function(pop, stat0, supfig=FALSE,shiftup=2, bars=TRUE,pop_title=TRUE) {
  # title and sign. bars 
  dfea_ghana <- filter(dfea, country==pop)
  xlab0 <- ""
  title0 <- ""
  if(stat0=="alpha"){expchar <- expression(alpha)}
  else if(stat0=="omega"){
    expchar <- expression(omega)
    xlab0 <- "Life stage"}
  else{
    expchar <- "Q"
    if(pop_title) {
      if(pop=='drc') {title0 <- toupper(pop)}
      else{title0 <- str_to_title(pop)}
    }
  }
  if(supfig){ if(pop!='drc'){expchar <- ""} }
  title <- expression(paste("Estimated ", expchar, " by life stage"))
  a <- pairwise.wilcox.test(dfea_ghana[,eval(stat0)], dfea_ghana$stage,p.adjust.method = "BH") %>% tidy() %>% rowwise() %>% filter(p.value < 0.05)
  shift <- range(dfea_ghana[,eval(stat0)])[2]-range(dfea_ghana[,eval(stat0)])[1]
  a$y.position <- seq(max(dfea_ghana[,eval(stat0)]) + shift*0.05, max(dfea_ghana[,eval(stat0)]) + shift*shiftup, length.out=nrow(a))
  a$p.signif <- sapply(a$p.value, function(x) ifelse(x < 0.05, ifelse(x < 0.01, ifelse(x < 0.001, '***', '**'), '*'), ' '))
  # plotting ftn 
  p <- ggplot(dfea_ghana) + geom_boxplot(aes(y=eval(parse(text=stat0)), x = stage, fill=stage), alpha=0.5, outlier.color = NA) + 
    geom_jitter(aes(y=eval(parse(text=stat0)), x = stage), position=position_jitter(0.2)) +  xlab(xlab0) +
    ylab(expchar) + theme_bw() + ggtitle(title0) + 
    theme(axis.text = element_text(size=15), axis.text.x = element_text(size=12, angle = 45, hjust=1),   
          axis.title.x=element_text(size=16), axis.title.y=element_text(size=40), 
          plot.title = element_text(size=20, hjust=0.5), panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
          legend.text=element_text(size=20), legend.title=element_text(size=15), 
          strip.background = element_blank(), strip.text = element_text(size=30), legend.position = "none") + 
    scale_fill_viridis(discrete=TRUE, labels=c('sporozoite', 'ring', 'trophozoite', 'schizont', 
                                               'gametocyte', 'ookinete'), name='Life stage') + 
    scale_color_viridis(discrete=TRUE, labels=c('sporozoite', 'ring', 'trophozoite', 'schizont', 
                                                'gametocyte', 'ookinete'), name='Life stage')  
  
  if(bars) { p <- p + stat_pvalue_manual(data.frame(a), label = "p.signif", size=5, bracket.size = 0.7) }
  return(p)
}

# fig 5
run_id_main = "dfe_files/DFE_results_m3_reich_syn_demoexp.csv" # main dataset
dfea <- read_dfe(run_id_main)
# plot main fig: Ghana data, synonymous sites
ggarrange(dfe_plots("ghana","Q",bars = FALSE,pop_title = FALSE), 
          dfe_plots("ghana","alpha",bars = FALSE,pop_title = FALSE), 
          dfe_plots("ghana","omega",bars = FALSE,pop_title = FALSE),
          labels = c("A", "B", "C"), font.label = list(size = 25), ncol=2,nrow=2)
ggsave(filename="fig5.png", dpi=300,
       width = 10, height=10)

dfe_stats(dfea)

# fig. S10
# plot FFD and praefal together 
run_id <- "dfe_files/DFE_results_m3_reich_ffd_demoexp.csv"
dfea <- read_dfe(run_id)
# significance test: pairwise wilcoxon rank sum test
q_ffd <- dfe_plots("ghana", "Q")
a_ffd <- dfe_plots("ghana", "alpha")
o_ffd <- dfe_plots("ghana", "omega")
dfe_stats(dfea)
run_id = "dfe_files/DFE_results_m3_praefal_syn_demoexp.csv" # praefal data 
dfea <- read_dfe(run_id,species="praefal")
q_pr <- dfe_plots("ghana", "Q")
a_pr <- dfe_plots("ghana", "alpha")
o_pr <- dfe_plots("ghana", "omega")
dfe_stats(dfea)
ggarrange(q_ffd, q_pr,
          a_ffd, a_pr,
          o_ffd, o_pr,
          labels = c("A", "B", 
                     "C", "D", 
                     "E", "F"), 
          ncol = 2, nrow = 3)
w = 8.5
h = 12
rat = h/w
scaled=2
ggsave(filename="lsfigs/figS10.png", 
       dpi=300, width = 5.5, height=5.5*rat, scale = scaled, units="in")

# fig S12
run_id = dfea <- "dfe_files/DFE_results_m3_final_reich.csv"
dfea <- read_dfe(run_id)
# plot main fig: Ghana data, synonymous sites
ggarrange(dfe_plots("ghana","Q",pop_title = FALSE), 
          dfe_plots("ghana","alpha",pop_title = FALSE), 
          dfe_plots("ghana","omega",pop_title = FALSE),
          labels = c("A", "B", "C"), font.label = list(size = 25), ncol=2,nrow=2)
ggsave(filename="lsfigs/figS12.png", dpi=300,
       width = 10, height=10)
dfe_stats(dfea)

# fig. S11
run_id = "dfe_files/DFE_results_m3_reich_syn_demoexp_suppCountries.csv" # alternate pop data 
dfea <- read_dfe(run_id)
w = 8.5
h = 12
rat = h/w
scaled=2
ggarrange(dfe_plots("drc", "Q"), dfe_plots("tanzania", "Q"), 
          dfe_plots("drc", "alpha"), dfe_plots("tanzania", "alpha"), 
          dfe_plots("drc", "omega"), dfe_plots("tanzania", "omega"), 
          labels = c("A", "B", "C", 
                     "D", "E", "F"), 
          ncol = 2, nrow = 3)

ggsave(filename="lsfigs/figS11.png", 
       dpi=300, width = 5.5, height=5.5*rat, scale = scaled, units="in")

dfe_stats(filter(dfea,country=="tanzania"))
dfe_stats(filter(dfea,country=="drc"))
