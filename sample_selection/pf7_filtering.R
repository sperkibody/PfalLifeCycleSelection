library(dplyr)
setwd('~')
path <- '' # input path to data directory 
pfdf <- read.csv(paste(path,'Pf7_samples.txt',sep=''), sep = '\t')
# initial filtering steps 
pfdf %>% filter(Country=='Cambodia', QC.pass=='True', Sample.type=='gDNA', X..callable > 50) %>% group_by(Study) %>% summarize(n())
pfdf %>% filter(Country=='Democratic Republic of the Congo', QC.pass=='True', Sample.type=='gDNA', X..callable > 50) %>% group_by(Study) %>% summarize(n())
pfdf %>% filter(Country=='Ghana', QC.pass=='True', Sample.type=='gDNA', X..callable > 50) %>% group_by(Study) %>% summarize(n())
pfdf %>% filter(Country=='Tanzania', QC.pass=='True', Sample.type=='gDNA', X..callable > 50) %>% group_by(Study) %>% summarize(n())

# extract sample IDs for each. 
cambodia <- pfdf %>% filter(Country=='Cambodia', QC.pass=='True', Sample.type=='gDNA', X..callable > 50) %>% dplyr::select(Sample)
drc <- pfdf %>% filter(Country=='Democratic Republic of the Congo', QC.pass=='True', Sample.type=='gDNA', X..callable > 50) %>% dplyr::select(Sample)
ghana <- pfdf %>% filter(Country=='Ghana', QC.pass=='True', Sample.type=='gDNA', X..callable > 50) %>% dplyr::select(Sample) 
tanzania <- pfdf %>% filter(Country=='Tanzania', QC.pass=='True', Sample.type=='gDNA', X..callable > 50) %>% dplyr::select(Sample)


write.table(cambodia, file = paste(path,"cambodia.txt",sep=''), 
            quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(tanzania, file = paste(path,"tanzania.txt",sep=''), 
            quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(drc, file = paste(path,"drc.txt",sep=''), 
            quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(ghana, file = paste(path,"ghana.txt",sep=''), 
            quote = FALSE, row.names = FALSE, col.names = FALSE)
######################################################
#FILTER BY COI (Fws)
######################################################
coi_ls <- read.delim(paste(path, "Pf7_fws.txt", sep=''), sep = '\t')
pfdf <- left_join(pfdf, coi_ls, by="Sample") %>% filter(Fws > 0.95)

################################################################
#CHOOSE REPRESENTATIVE FROM HMMIBD CLUSTERS BASED ON COVERAGE
################################################################

cambodia <- pfdf %>% filter(Country=='Cambodia', QC.pass=='True', Sample.type=='gDNA', X..callable > 50) 
drc <- pfdf %>% filter(Country=='Democratic Republic of the Congo', QC.pass=='True', Sample.type=='gDNA', X..callable > 50)
ghana <- pfdf %>% filter(Country=='Ghana', QC.pass=='True', Sample.type=='gDNA', X..callable > 50) 
tanzania <- pfdf %>% filter(Country=='Tanzania', QC.pass=='True', Sample.type=='gDNA', X..callable > 50)

cghana <- read.csv(paste(path,'ghana.clusters.0.25.txt',sep=''), sep = '\t')
cghana$Sample <- row.names(cghana)
cghana <- cghana %>% dplyr::rename(cluster = x)
ghana <- left_join(ghana, cghana, by = 'Sample')

cdrc <- data.frame(read.csv(paste(path,'drc.clusters.0.25.txt',sep=''), sep = '\t'))
cdrc$Sample <- row.names(cdrc)
cdrc <- cdrc %>% dplyr::rename(cluster = x)
drc <- left_join(drc, cdrc, by = 'Sample')

ctan <- read.csv(paste(path,'tanzania.clusters.0.25.txt',sep=''), sep = '\t')
ctan$Sample <- row.names(ctan)
ctan <- ctan %>% dplyr::rename(cluster = x)
tanzania <- left_join(tanzania, ctan, by = 'Sample')

ccam <- read.csv(paste(path,'cambodia.clusters.0.25.txt',sep=''), sep = '\t')
ccam$Sample <- row.names(ccam)
ccam <- ccam %>% dplyr::rename(cluster = x)
cambodia <- left_join(cambodia, ccam, by = 'Sample')

nrow(ghana)
ghana$cluster_rep <- sapply(ghana$Sample, function(x) ifelse(ghana[ghana$Sample==x,'X..callable']== 
                                          max(ghana[ghana$cluster==ghana[ghana$Sample==x,'cluster'],'X..callable']), TRUE, FALSE))
ghana <- ghana %>% filter(cluster_rep == TRUE)
ghana %>% group_by(cluster) %>% summarize(n=n()) %>% filter(n > 1) %>% nrow()
nrow(ghana)

nrow(drc)
drc$cluster_rep <- sapply(drc$Sample, function(x) ifelse(drc[drc$Sample==x,'X..callable']== 
                                                               max(drc[drc$cluster==drc[drc$Sample==x,'cluster'],'X..callable']), TRUE, FALSE))
drc <- drc %>% filter(cluster_rep == TRUE)
drc %>% group_by(cluster) %>% summarize(n=n()) %>% filter(n > 1) %>% nrow()
nrow(drc)

nrow(cambodia)
cambodia$cluster_rep <- sapply(cambodia$Sample, function(x) ifelse(cambodia[cambodia$Sample==x,'X..callable']== 
                                                               max(cambodia[cambodia$cluster==cambodia[cambodia$Sample==x,'cluster'],'X..callable']), TRUE, FALSE))
cambodia <- cambodia %>% filter(cluster_rep == TRUE)
cambodia %>% group_by(cluster) %>% summarize(n=n()) %>% filter(n > 1) %>% nrow()
nrow(cambodia)

nrow(tanzania)
tanzania$cluster_rep <- sapply(tanzania$Sample, function(x) ifelse(tanzania[tanzania$Sample==x,'X..callable']== 
                                                               max(tanzania[tanzania$cluster==tanzania[tanzania$Sample==x,'cluster'],'X..callable']), TRUE, FALSE))
tanzania <- tanzania %>% filter(cluster_rep == TRUE)
tanzania %>% group_by(cluster) %>% summarize(n=n()) %>% filter(n > 1) %>% nrow()
nrow(tanzania)

ghana <- ghana %>% dplyr::select(Sample)
cambodia <- cambodia %>% dplyr::select(Sample)
drc <- drc %>% dplyr::select(Sample)
tanzania <- tanzania %>% dplyr::select(Sample)

write.table(cambodia, file = paste(path,"cambodia.txt",sep=''), 
            quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(tanzania, file = paste(path,"tanzania.txt",sep=''), 
            quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(drc, file = paste(path,"drc.txt",sep=''), 
            quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(ghana, file = paste(path,"ghana.txt",sep=''), 
            quote = FALSE, row.names = FALSE, col.names = FALSE)

# print out summary of sample number for all 
print(paste("Ghana N =", nrow(ghana),"| Tanzania N =", nrow(tanzania), "| DRC N =", nrow(drc), "Cambodia N=", nrow(cambodia)))
# obtain summary of sample collection dates 
print(summary(unlist(c(tanzania$Year, ghana$Year, drc$Year, cambodia$Year))))

