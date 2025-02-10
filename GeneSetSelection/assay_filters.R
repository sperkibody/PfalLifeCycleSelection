#***********************************************************************
#*File for filtering life stage specific gene sets, including: 
#*1) filtering for variable genes genes via M3 Drop (Feature selection): https://github.com/tallulandrews/M3Drop/blob/master/vignettes/M3Drop_Vignette.Rnw
#*2) defining stage-specific genes as described by Howick et al. and Tebben et al.--genes expressed in 50% of cells with a life stage label for only one life stage. 
#*3) filtering at assay level for genes with low level expression 
#*4) examining overlap between gene sets 
#*5) examing overlap between life stage expression patterns overall 
#***********************************************************************
# load libraries
suppressWarnings({
  suppressPackageStartupMessages({
    library(dplyr)
    library(Seurat)
    library(patchwork)
    #library(scAlign)
    library(gridExtra)
    library(ggplot2)
    library(M3Drop)
    library(stats)
    library(scran)
    library(scuttle)
    library(SingleCellExperiment)
    library(VennDiagram)
    library(reshape2)
  })
})

setwd('~')
path="" # insert path
ags <- read.csv(paste(path,'helb_seropos.txt',sep=""))[[1]]
# M3Drop feature selection 
fn_cmd = paste(path,'pf-ss2-set1/pf-ss2-set1-ss2-data.csv',sep="")
fn_expr = paste(path,'pf-ss2-set1/pf-ss2-set1-ss2-exp.csv',sep="")

cmd <- read.csv(fn_cmd)
labels <- cmd$CELL_ID
counts <- read.csv(fn_expr, row.names = 1)
# MAKE CUSTOM LABELS: LOWRES LABELS WITH GAMETOCYTE SPLIT INTO HIGH RES CATEGORIES
ls <- cmd$STAGE_LR
data.frame(ls) %>% group_by(ls) %>% summarize(n())
# M3Drop QC: remove cells with fewer than 2000 transcripts counted 
total_features <- colSums(counts >= 0)
counts <- counts[,total_features >= 2000]
labels <- labels[total_features >= 2000]

#***********************************************************************
##########*M3DROP FEATURE SELECTION STEP################################
#***********************************************************************
norm <- M3DropConvertData(counts, is.counts=TRUE)
# get M3Drop genes (feature selection step)
M3Drop_genes <- M3DropFeatureSelection(norm, mt_method="fdr", mt_threshold=0.01)
heat_out <- M3DropExpressionHeatmap(M3Drop_genes$Gene, norm, cell_labels = labels)
marker_genes <- M3DropGetMarkers(norm, ls)
marker_genes$fdr <- p.adjust(marker_genes$pval, method = 'BH')
marker_genes <- marker_genes %>% filter(fdr < 0.05)
marker_genes <- marker_genes[rownames(marker_genes) %in% M3Drop_genes$Gene,]
marker_genes %>% group_by(Group) %>% summarize(n = n(), p_max = max(pval))

# check counts with antigen removal
m3drop_nonag <- marker_genes[!rownames(marker_genes) %in% ags,] 
m3drop_nonag <- data.frame(gene=rownames(m3drop_nonag),stage=m3drop_nonag$Group)
m3drop_nonag %>% group_by(stage) %>% summarize(n = n())

fn_out_m3 = paste(path,"m3_combined.csv",sep="")
# write to output all m3 genes without breadth exclusions if desired 
# write.csv(m3drop_nonag, fn_out_m3, row.names = FALSE)

#***********************************************************************
##########*ASSAY CHARACTERISTICS#######################################
#***********************************************************************
# re-read in DF 
cmd <- read.csv(fn_cmd, row.names = 1)
cmd$cell_name <- rownames(cmd)
# filter as per seurat guidelines, removing cells with MIT and API 
rownames(cmd) <- cmd$cell_name
counts <- read.csv(fn_expr, row.names = 1)
pfss2 <- CreateSeuratObject(counts, project = "SeuratProject", assay = "RNA",  meta.data = cmd)
FetchData(object = pfss2, vars = c("STAGE_LR"))
pfss2[["percent.mt"]] <- PercentageFeatureSet(pfss2, pattern = "*-MIT*")
pfss2[["percent.api"]] <- PercentageFeatureSet(pfss2, pattern = "PF3D7-API")
# (min feature count tailored to this assay / parasite)
pfss2 <- subset(pfss2, subset = nFeature_RNA > 150 & nFeature_RNA < 2500 & percent.mt < 5 & percent.api < 5)
pfss2@meta.data %>% group_by(STAGE_LR) %>% summarize(n())
# graph features by stage
lr <- pfss2@meta.data$'STAGE_LR'
Idents(pfss2) <- factor(lr, levels=c('sporozoite', 'ring', 'trophozoite', 'schizont', 'gametocyte', 'ookinete'))
# visualize expression by lowres life stage annotation - save plots for fig. S1 (end of file)
ss2_counts <- VlnPlot(pfss2, features = 'nCount_RNA', log=TRUE, cols=scale_color_viridis(discrete=TRUE, labels=c('sporozoite', 'ring', 'trophozoite', 'schizont', 
                                                                                                                 'gametocyte', 'ookinete')) )
ss2_features <- VlnPlot(pfss2, log=TRUE, features = 'nFeature_RNA', cols=scale_color_viridis(discrete=TRUE, labels=c('sporozoite', 'ring', 'trophozoite', 'schizont', 
                                                                                                                     'gametocyte', 'ookinete')))
dev.off()

# save assay features to txt for diversity analysis 
all <- rownames(counts)
write.table(all, file = paste(path,"gene_sets_all/all.txt",sep=""), 
            quote = FALSE, row.names = FALSE, col.names = FALSE)

# extract counts 
sdf <- pfss2@assays[['RNA']]$counts
sdf <- data.frame(sdf)
# create matrix indicating binary detectable expression (0|1)
exprdf <- data.frame()
for(i in 1:nrow(sdf)){
  print(i)
  cells <- sapply(1:ncol(sdf), function(x) (sdf[i,x] > 0))
  exprdf <- rbind(exprdf, cells)
}
edf2 <- exprdf
colnames(edf2) <- unname(pfss2$STAGE_LR)
# check for expression in at least 2.5% of cells in one cell type 
expressed_all <- data.frame()
for(i in unique(colnames(edf2))) {
  sub <- edf2[,colnames(edf2)==i]
  min_cells <- round(0.025*ncol(sub))
  sums <- sapply(1:nrow(sub), function(x) sum(sub[x,])>min_cells)
  pass <- rownames(counts)[sums]
  if(nrow(expressed_all)==0){
    expressed_all <- data.frame(gene=pass,stage=rep(i, length(pass)))
  }
  else{
    expressed_all <- rbind(expressed_all, data.frame(gene=pass,stage=rep(i, length(pass))))
  } 
}
n_distinct(expressed_all$gene)
min_expression_2.5 <- unique(expressed_all$gene)
write.table(min_expression_2.5,paste(path,'min_expression_2.5.txt',sep=""),sep='\n',
            quote = FALSE, row.names = FALSE, 
            col.names = FALSE)

#***********************************************************************
##########*ALTERNATIVE GENE SETS/BREADTH################################
#***********************************************************************
#*repeat the above but with 50% threshold 
expressed_all <- data.frame()
for(i in unique(colnames(edf2))) {
  sub <- edf2[,colnames(edf2)==i]
  min_cells <- round(0.5*ncol(sub))
  sums <- sapply(1:nrow(sub), function(x) sum(sub[x,])>min_cells)
  pass <- rownames(counts)[sums]
  if(nrow(expressed_all)==0){
    expressed_all <- data.frame(gene=pass,stage=rep(i, length(pass)))
  }
  else{
    expressed_all <- rbind(expressed_all, data.frame(gene=pass,stage=rep(i, length(pass))))
  } 
}
length(min_expression_2.5)
length(min_expression_2.5) - n_distinct(expressed_all$gene)
tallies <- expressed_all %>% group_by(gene) %>% summarize(n=n()) %>% arrange(desc(n))
groupings_50_min <- expressed_all %>% group_by(gene) %>% summarize(cell_types = paste(sort(unique(stage)), collapse = ", "))
a <- tallies %>% group_by(n) %>% summarize(n = n())
sum(a$n)
ggplot(tallies, aes(x=n)) + geom_histogram()
sum(unique(props$GENE) %in% unique(expressed_all$gene))
sum(a$n)

# isolate stage-specific genes
ss <- tallies %>% filter(n==1) 
ss_labeled <- expressed_all[expressed_all$gene %in% ss$gene,]
head(ss_labeled)
ss_labeled %>% group_by(stage) %>% summarize(n())
ss_labeled[!ss_labeled$gene %in% ags,] %>% group_by(stage) %>% summarize(n())

library(data.table)
mdf <- data.frame()
# apply additional differential expression filter to these 
for(n in 1:5) {
  print(n)
  combos <- combn(unique(ls), n)
  print(paste(ncol(combos), 'unique combinations'))
  # update labels 
  for (i in 1:ncol(combos)){
    combo_string <- paste(sort(combos[,i]),  collapse = ", ")
    print(paste(i,'/', ncol(combos)))
    grouped_ls <- ls
    in_combo <- ls %in% combos[,i][[1]]
    # extract group as list 
    grouped_ls[in_combo] <- paste('g', n, i, sep = '_')
    # pairwise t-tests for DE 
    marker_genes <- M3DropGetMarkers(norm, grouped_ls)
    marker_genes$gene <- rownames(marker_genes)
    union_genes <- groupings_50_min$gene[groupings_50_min$cell_types==combo_string]
    marker_genes <- marker_genes[marker_genes$gene %in% union_genes,]
    if (!nrow(marker_genes)){next}
    marker_genes$breadth <- n
    # only keep results specific to tested group 
    marker_genes <- marker_genes %>% filter(Group==paste('g', n, i, sep = '_'))
    # concat uncorrected results to full df 
    if(!nrow(mdf)) { mdf <- marker_genes}
    else{mdf <- rbind(mdf, marker_genes)}
  }
}

# within final DF, adjust p-values and find minimum (or max?) breadth 
mdf$fdr <- p.adjust(mdf$pval, method = 'BH')
mdf <- mdf %>% filter(fdr < 0.01)
breadth_6 <- tallies[tallies$n==6,]$gene
breadth_6 <- genes[!(genes %in% mdf$gene) & (genes %in% filter(n_50_min, n==6)$gene)]
breadth_de <- mdf %>% group_by(gene) %>% summarize(breadth = max(breadth)) # %>% 
# filter(!gene %in% breadth_6)

# add in the final category 
breadth_6 <- data.frame(gene = breadth_6, breadth=rep(as.integer(6),length(breadth_6)))
breadth_de <- breadth_de %>% rbind(breadth_6) %>% filter(gene %in% min_expression_2.5)
breadth_de %>% group_by(breadth) %>% summarize(n())
# save breadth file without excluding antigens 
fn_out_breadth = paste(path,"breadth_final.csv",sep="",row.names=FALSE)
write.csv(breadth_de, fn_out_breadth)
4930 - nrow(breadth_de)

# also save this to txt files 
for(i in unique(breadth_de$breadth)){
  sub <- breadth_de[breadth_de$breadth==i,]
  fn_sub <- paste(path,'gene_sets_all/',i,'.txt',sep="")
  write.table(sub$gene, file = fn_sub, quote = FALSE, row.names = FALSE, 
              col.names = FALSE)
}
# check how many non-antigenic core genes got excluded from breadth estimation
props %>% mutate(ag = GENE %in% ags) %>% filter(!ag) %>% 
  left_join(dplyr::rename(breadth_de, c("GENE"="gene")),by="GENE") %>% group_by(breadth) %>% summarize(n=n()) 

# pull out stage labels from original 50% comparisons for genes with evidence 
# for differential expression in a single stage 
ss <- breadth_de %>% filter(breadth==1) 
ss_labeled <- expressed_all[expressed_all$gene %in% ss$gene,]
#ss_labeled <- ss_labeled[ss_labeled$gene %in% M3Drop_genes$Gene,]
ss_labeled %>% group_by(stage) %>% summarize(n())
# remove antigens 
ss_labeled <- ss_labeled %>% filter(!gene %in% ags)
ss_labeled %>% group_by(stage) %>% summarize(n())
# save stage-specific 50%/DE files 
fn_out_gs50 = paste(path,"stage_specific.csv",sep="")
# write to output all m3 genes without antigen/breadth exclusions
write.csv(ss_labeled, fn_out_gs50, row.names = FALSE)
genes_not_ss <- breadth_de %>% filter(breadth>1) %>% dplyr::select(gene)


#***********************************************************************
##################*SAVE MAIN GENE SETS AS TXT FILES#####################
#***********************************************************************
# remove m3 genes with breadth > 1 and check category sizes
m3drop_final <- m3drop_nonag %>% filter(!gene %in% genes_not_ss$gene)
m3drop_final %>% group_by(stage) %>% summarize(n = n())

# write to output all m3 genes without breadth exclusions if desired 
write.csv(m3drop_final, fn_out_m3, row.names = FALSE)

for(i in unique(m3drop_final$stage)){
  sub <- m3drop_final[m3drop_final$stage==i,]
  fn_sub <- paste(path,'gene_sets_m3_final/',i,'.txt',sep="")
  write.table(sub$gene, file = fn_sub, quote = FALSE, row.names = FALSE, 
              col.names = FALSE)
}

##############EXAMINE OVERLAP BETWEEN GENE SETS###################
stage_labels_m3 <- read.csv(fn_out_m3) %>% dplyr::rename(c("GENE"="gene"))
stage_labels_50 <- read.csv(fn_out_gs50) %>% dplyr::rename(c("GENE"="gene"))
# excluge breadth > 1 from M3 sets for visualization given its exclusion in analysis 
stage_labels_m3 <- stage_labels_m3 %>% filter(!GENE %in% genes_not_ss$gene)
stage_labels_m3 %>% group_by(stage) %>% summarize(n = n())

# overwrite gs50 with sporozoite gold standard set 
spz_full <- read.delim(file = paste(path,"sporozoite_full.txt",sep=""), 
                       header = FALSE)

# R venn diagram output: base formatting code from R graph gallery https://r-graph-gallery.com/14-venn-diagramm
myCol <- brewer.pal(2, "Pastel2")
overlap <- function(stagename){
  sub50 <- stage_labels_50[stage_labels_50$stage==stagename,]
  subm3 <- stage_labels_m3[stage_labels_m3$stage==stagename,]
  venn.diagram(x = list(sub50$GENE, subm3$GENE), category.names = c('', ''), 
               filename = paste(path,'gene_set_overlaps/', stagename, '.png', sep=""), 
               output = TRUE,
               # Output features
               imagetype="png" ,
               height = 480 , 
               width = 480 , 
               resolution = 300, 
               # Circles
               lwd = 2,
               lty = 'blank',
               fill = myCol[1:2],
               
               # Numbers
               cex = .6,
               fontface = "bold",
               fontfamily = "sans",
               
               # Set names
               cat.cex = 0.3,
               cat.fontface = "bold",
               cat.default.pos = "outer",
               cat.pos = c(-10, 15),
               cat.dist = c(0.055, 0.055),
               cat.fontfamily = "sans")
}

#overlap('sporozoite')
overlap('ookinete')
overlap('trophozoite')
overlap('gametocyte')
overlap('ring')
overlap('schizont')

# for sporozoite, also add in gold standard set for comparison 
# (since that's used later on in the analysis)
stagename='sporozoite'
sub50 <- stage_labels_50[stage_labels_50$stage==stagename,]
subm3 <- stage_labels_m3[stage_labels_m3$stage==stagename,]
spz_full_non_ags <- spz_full$V1[!spz_full$V1 %in% ags]
venn.diagram(x = list(sub50$GENE, subm3$GENE, spz_full_non_ags), category.names = c('', '', ''), 
             filename = paste(path,'gene_set_overlaps/sporozoite.png', sep=""), 
             output = TRUE,
             # Output features
             imagetype="png" ,
             height = 480 , 
             width = 480 , 
             resolution = 300, 
             # Circles
             lwd = 2,
             lty = 'blank',
             fill = myCol[1:3],
             
             # Numbers
             cex = .6,
             fontface = "bold",
             fontfamily = "sans",
             
             # Set names
             cat.cex = 0.3,
             cat.fontface = "bold",
             cat.default.pos = "outer",
             cat.pos = c(-10, 15, 2.5),
             cat.dist = c(0.055, 0.055, 0.005),
             cat.fontfamily = "sans")
#***********************************************************************
##########*MAKE PAIRWISE CORRELATION HEATMAP FOR STAGES#################
#***********************************************************************
# using all non-antigenic genes detected in at least one stage by this definition
cordat <- data.frame()
for(stage in unique(colnames(edf2))){
  print(stage)
  stagecol <- data.frame()
  for(gene in unique(expressed_all$gene)){
    detected <- gene %in% expressed_all$gene[expressed_all$stage==stage]
    newrow <- data.frame(detected)
    if(nrow(stagecol)==0){ stagecol <- newrow}
    else{stagecol <- rbind(stagecol,newrow) }
  }
  print(nrow(stagecol))
  if(ncol(cordat)==0) { cordat <- stagecol}
  else{ cordat <- cbind(cordat, stagecol) }
}
rownames(cordat) <- unique(expressed_all$gene)
colnames(cordat) <- unique(colnames(edf2))
cordat <- cordat[!rownames(cordat) %in% ags,] # remove antigens 
nrow(cordat)

library(Hmisc)

cor <- rcorr(as.matrix(cordat), type = "pearson")
cor_r <- cor$r
cor_p <- cor$P*15
signif_level <- ifelse(cor_p < 0.001, '***', 
                       ifelse(cor_p < 0.01, '**', 
                              ifelse(cor_p < 0.05, '*','')))

cor_full <- melt(cor_r)
signif_full <- melt(signif_level)
p_full <- melt(cor_p)
cor_full$signiflevel <- signif_full$value
cor_full$p <- p_full
order <- c("gametocyte","ookinete", "sporozoite","ring","trophozoite",
           "schizont")
# remove redundant comparisons 
cor_full <- cor_full %>% rowwise() %>% mutate(order= (which(order==Var1) > which(order==Var2))) %>% filter(order) %>% dplyr::select(-order)
cor_full$Var1 <- factor(cor_full$Var1, levels=order)
cor_full$Var2 <- factor(cor_full$Var2, levels=order)
ggplot(cor_full, mapping=aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() + theme_minimal() + 
  geom_text(aes(label = paste(round(value,3), signiflevel)), color="black", size=5) + 
  scale_fill_gradientn(colors=c("blue","white","red"),
                      limit=c(-1,1), space="Lab", name="correlation") +
  theme_minimal() + theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1, size = 15),
                          axis.text.y = element_text(size=15),
                          axis.title.y = element_text(size=15),
                          axis.title.x = element_text(size=15),
                          panel.grid = element_blank()) + 
  xlab("Life stage") + ylab("Life stage")

ggsave(paste(path,"lsfigs/figS6.png",sep=""),
       dpi=300, width = 5, height=4, 
       scale = 1.7, units="in")


#***********************************************************************
##########*MAKE FIGS1 TO SHOW ASSAY CHARACTERISTICS AND GENE SET OUTPUTS
#***********************************************************************
#*read in DF with coding length info 
props$X <- gsub("location=", "", props$COORD)
props$X <- sapply(props$X, function(x) substr(x, 0, nchar(x)-3))
props <- props %>% arrange(GENE, desc(TOTAL_CODING_LENGTH))
props$ag <- props$GENE %in% ags
propls <- left_join(props, stage_labels_m3, by = 'GENE') 
propls50 <- left_join(props, stage_labels50, by = 'GENE') 

propls50$stage[propls50$GENE %in% spz_full$V1] <- 'sporozoite'
propls$stage <- factor(propls$stage, levels=c('sporozoite', 'ring', 'trophozoite', 'schizont', 
                                              'gametocyte', 'ookinete'))
propls50$stage <- factor(propls50$stage, levels=c('sporozoite', 'ring', 'trophozoite', 'schizont', 
                                                  'gametocyte', 'ookinete'))
propls50 %>% group_by(stage) %>% summarize(n())

csub <- filter(propls, !is.na(stage), !ag)
a <- pairwise.wilcox.test(csub$TOTAL_CODING_LENGTH, csub$stage,p.adjust.method = "BH") %>% tidy() %>% rowwise() %>% filter(p.value < 0.05)
a$y.position <- seq(max(log10(csub$TOTAL_CODING_LENGTH))*1.1, max(log10(csub$TOTAL_CODING_LENGTH))*1.25, length.out=nrow(a))
a$p.signif <- sapply(a$p.value, function(x) ifelse(x < 0.05, ifelse(x < 0.01, ifelse(x < 0.001, '***', '**'), '*'), ' '))

coding_length_m3 <- ggplot(data = filter(propls, !is.na(stage))) + 
  geom_point(aes(y=log10(TOTAL_CODING_LENGTH), x = stage, fill=stage), position = position_jitterdodge(), alpha=0.5) +  
  geom_boxplot(alpha=0.5, outlier.color = NA, aes(y=log10(TOTAL_CODING_LENGTH), x = stage, fill=stage))  + xlab("") + 
  ylab(expression(log10("Total Coding Length"))) + theme_bw() + 
  theme(axis.text = element_text(size=8),   
        axis.title.x=element_text(size=16), axis.title.y=element_text(size=12), 
        axis.text.x = element_text(size=12, angle = 45, hjust=1),
        plot.title = element_text(size=30), panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
        legend.text=element_text(size=15), legend.title=element_text(size=15), legend.position = "none") + 
  scale_fill_viridis(discrete=TRUE, labels=c('sporozoite', 'ring', 'trophozoite', 'schizont', 
                                             'gametocyte', 'ookinete'), name='Life stage') + 
  scale_color_viridis(discrete=TRUE, labels=c('sporozoite', 'ring', 'trophozoite', 'schizont', 
                                              'gametocyte', 'ookinete'), name='Life stage') + 
  stat_pvalue_manual(data.frame(a), label = "p.signif", size=5, bracket.size = 0.7)
# same plot for alt. gene sets 
csub <- filter(propls50, !is.na(stage), !ag)
a <- pairwise.wilcox.test(csub$TOTAL_CODING_LENGTH, csub$stage,p.adjust.method = "BH") %>% tidy() %>% rowwise() %>% filter(p.value < 0.05)
a$y.position <- seq(max(log10(csub$TOTAL_CODING_LENGTH))*1.1, max(log10(csub$TOTAL_CODING_LENGTH))*1.25, length.out=nrow(a))
a$p.signif <- sapply(a$p.value, function(x) ifelse(x < 0.05, ifelse(x < 0.01, ifelse(x < 0.001, '***', '**'), '*'), ' '))

coding_length_50 <- ggplot(data = filter(propls50, !is.na(stage))) + 
  geom_point(aes(y=log10(TOTAL_CODING_LENGTH), x = stage, fill=stage), position = position_jitterdodge(), alpha=0.5) +  
  geom_boxplot(alpha=0.5, outlier.color = NA, aes(y=log10(TOTAL_CODING_LENGTH), x = stage, fill=stage))  + xlab("") + 
  ylab(expression(log10("Total Coding Length"))) + theme_bw() + 
  theme(axis.text = element_text(size=8),   
        axis.title.x=element_text(size=16), axis.title.y=element_text(size=12), 
        axis.text.x = element_text(size=12, angle = 45, hjust=1),
        plot.title = element_text(size=30), panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
        legend.text=element_text(size=15), legend.title=element_text(size=15), legend.position = "none") + 
  scale_fill_viridis(discrete=TRUE, labels=c('sporozoite', 'ring', 'trophozoite', 'schizont', 
                                             'gametocyte', 'ookinete'), name='Life stage') + 
  scale_color_viridis(discrete=TRUE, labels=c('sporozoite', 'ring', 'trophozoite', 'schizont', 
                                              'gametocyte', 'ookinete'), name='Life stage') + 
  stat_pvalue_manual(data.frame(a), label = "p.signif", size=5, bracket.size = 0.7)

scaled = 8/5.5
ggarrange(ss2_counts, ss2_features, coding_length_m3,
          coding_length_50, font.label = list(size = 15),
          labels = c("A", "B", "C", "D"), ncol = 2, nrow = 2)
ggsave(filename=paste(path,"lsfigs/figS1.png",sep=""), dpi=300,
       width = 5.5, height=5.5, scale = scaled, units="in")

spz_full_non_ags <- data.frame(GENE=spz_full, stage=rep("sporozoite",length(spz_full))) %>% 
  filter(!GENE %in% ags)
stage_labels_50_spz_corrected <- stage_labels_50 %>% filter(stage!="sporozoite") %>% 
  dplyr::select(-X) %>% 
  rbind(spz_full_non_ags)
stage_labels_50_spz_corrected %>% group_by(stage) %>% summarize(n())
fn_out_gs50_corr <- paste(path,"stage_specific_spz_corrected.csv",sep="")
write.csv(stage_labels_50_spz_corrected, fn_out_gs50_corr, row.names = FALSE)

