library(dplyr)
library(stringr)
library(vsn)
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(readxl)
library(data.table)
library(tximport)
library(org.Hs.eg.db)

## clin <- read.csv('clinical-211205.csv', na.strings='#N/A', row.names=1)

## rownames(clin) <- str_replace_all(rownames(clin), '-', '.')

## head(clin)

## clin$Batch <- ifelse(str_detect(rownames(clin), 'KH'), 'KH', NA)
## clin$Batch <- ifelse(str_detect(rownames(clin), 'MJ'), 'MJ', clin$Batch)
## clin$Batch <- ifelse(str_detect(rownames(clin), 'Mini.DP'), 'Mini-DP', clin$Batch)

## clin$Condition. <- ifelse(clin$Condition. %in% c('PTD','Severe PE','FGR','Control','Control PTD'), clin$Condition., 'FGR+HDP')

## clin$PrePregBMI <- clin$PrePregWt_Kg / clin$Height_Meters.HeightMeters^2

## clin$Smoking <- ifelse(str_detect(clin$Smoke, 'Yes'), 'Yes', NA)
## clin$Smoking <- ifelse(str_detect(clin$Smoke, 'Never'), 'No', clin$Smoking)

## clin$Race[!clin$Race %in% c('W', 'B')] <- 'other'
## clin$FDELTYPE[clin$FDELTYPE==0] <- NA

## clin$DiscWksGest <- arules::discretize(clin$WksGest, method='interval', breaks=4)
## levels(clin$DiscWksGest)
## hist(clin$WksGest)
## hist(sapply(clin$DiscWksGest, function(x) which(levels(clin$DiscWksGest)==x)))

## clin$DiscPrePregBMI <- arules::discretize(clin$PrePregBMI, method='interval', breaks=4)
## levels(clin$DiscPrePregBMI)
## hist(clin$PrePregBMI)
## hist(unlist(sapply(clin$DiscPrePregBMI, function(x) which(levels(clin$DiscPrePregBMI)==x)), use.names=FALSE))

## final.clin <- clin %>% dplyr::select(Condition., Batch, DiscWksGest, InfSex, Race,
##                                      DiscPrePregBMI, WksGest, PrePregBMI,
##                                      Smoking, FDELTYPE, Labor.initiation)

clin <- read_excel('Primary clinical variables clustered.xlsx', 1) %>% as.data.frame
colnames(clin) <- make.names(colnames(clin))
clin <- clin[-1,]
rownames(clin) <- make.names(clin$Study_ID.ID)
head(clin)

clin <- clin %>% mutate_if(function(x) mean(is.na(as.numeric(x))) < 0.5, as.numeric)

clin$Batch <- ifelse(str_detect(rownames(clin), 'KH'), 'KH', NA)
clin$Batch <- ifelse(str_detect(rownames(clin), 'MJ'), 'MJ', clin$Batch)
clin$Batch <- ifelse(str_detect(rownames(clin), 'Mini.DP'), 'Mini-DP', clin$Batch)

clin$Condition. <- ifelse(clin$Condition. %in% c('PTD','Severe PE','FGR','Control','Control PTD'), clin$Condition., 'FGR+HDP')

clin$Condition. <- make.names(clin$Condition.)

clin$PrePregBMI <- clin$PrePregWt_Kg / clin$Height_Meters.HeightMeters^2

clin$Smoking <- ifelse(str_detect(clin$Smoke, 'Yes'), 'Yes', NA)
clin$Smoking <- ifelse(str_detect(clin$Smoke, 'Never'), 'No', clin$Smoking)

final.clin <- clin %>% dplyr::select(Condition., Batch, WksGest, InfSex, Race, PrePregBMI,
                                     Smoking, FDELTYPE, Labor.initiation, Labor.no.labor)

final.clin$Race[!final.clin$Race %in% c('W', 'B')] <- 'other'
final.clin$FDELTYPE[final.clin$FDELTYPE==0] <- NA

final.clin <- final.clin %>% mutate_at(c('Race', 'FDELTYPE', 'Labor.initiation', 'Labor.no.labor'), factor)

tail(final.clin)

counts <- read.csv('mdp3_pmat.csv', row.names=1)
dim(counts)
head(counts)

idMap <- read_excel('Transcriptomic sample information.xlsx')
colnames(idMap) <- make.names(colnames(idMap))
idMap <- idMap[!is.na(idMap$Sample.ID),]

idMap$Sample.ID <- make.names(sapply(idMap$Sample.ID, function(x) sub(' .*', '', x)))
idMap$Study.ID <- make.names(idMap$Study.ID)
tail(idMap)

setnames(counts, idMap$Sample.ID, idMap$Study.ID)

txdb <- org.Hs.eg.db
k <- keys(txdb, keytype = "ENSEMBL")
df <- select(txdb, keys = k, keytype = "ENSEMBL", columns = c("SYMBOL", "CHR"))
tx2gene <- df[, 1:3]
tx2gene <- tx2gene[tx2gene$ENSEMBL %in% rownames(counts),]

## dups <- tx2gene %>% group_by(ENSEMBL) %>% summarize(n()) %>% filter(`n()` > 1) %>% dplyr::select(ENSEMBL)

tx2gene <- tx2gene %>% group_by(ENSEMBL) %>% filter(row_number() == 1, CHR!='Y', !SYMBOL %in% c('XIST', 'TSIX'))

## tx2gene[tx2gene$ENSEMBL %in% dups$ENSEMBL,]

t.counts <- t(counts) %>% as.data.frame

setnames(t.counts, tx2gene$ENSEMBL, tx2gene$SYMBOL)

counts <- t(t.counts) %>% as.data.frame

## counts.sum <- summarizeToGene(counts, tx2gene)

dim(counts[,colnames(counts) %in% rownames(final.clin)])

final.clin$Batch <- make.names(final.clin$Batch)
final.clin$Condition. <- make.names(final.clin$Condition.)

final.clin <- final.clin %>% mutate_if(function(x) mean(is.na(as.numeric(x))) < 0.5,
                                       as.numeric)
final.clin <- final.clin %>% mutate_if(function(x) length(unique(x)) > 7,
                                       scale)

final.clin <- final.clin %>% mutate_at(c('FDELTYPE', 'Labor.initiation', 'Race', 'InfSex', 'Labor.no.labor', 'Smoking', 'Condition.', 'Batch'), as.factor)

head(final.clin)
dim(final.clin)
dim(final.clin[rowSums(is.na(final.clin))==0,])
final.clin <- final.clin[rowSums(is.na(final.clin))==0,]

samp.ids <- intersect(rownames(final.clin), colnames(counts))

## outliers <- c('MJ.0357', 'MJ.0376', 'MJ.0505', 'MJ.0581', 'MJ.0868', 'MJ.0978', 'Mini.DP.189', 'Mini.DP.237')

## outliers <- c('MJ.0042', 'MJ.0357', 'MJ.0376', 'MJ.0472', 'MJ.0505', 'MJ.0581', 'MJ.0587', 'MJ.0766', 'MJ.0868', 'MJ.0978', 'Mini.DP.109', 'Mini.DP.115', 'Mini.DP.139', 'Mini.DP.187', 'Mini.DP.189', 'Mini.DP.232', 'Mini.DP.237')

outliers <- c('MJ.0868', 'Mini.DP.189', 'MJ.0581', 'Mini.DP.237', 'MJ.0505', 'MJ.0357',
              'MJ.0376', 'Mini.DP.139', 'Mini.DP.109', 'MJ.0978', 'MJ.0587', 'Mini.DP.232',
              'MJ.0472')

outliers <- c('MJ.0788', 'KH.1103', 'MJ.1151', 'MJ.0286', 'MJ.0687', 'Mini.DP.025',
              'MJ.0909', 'Mini.DP.104', 'MJ.1065', 'KH.1100', 'Mini.DP.181', 'Mini.DP.013',
              'MJ.0899') # , 'MJ.0859', 'MJ.0868')

samp.ids <- samp.ids[!samp.ids %in% outliers]

## samp.ids <- rev(samp.ids)

dim(counts[rownames(counts) %in% tx2gene$SYMBOL,samp.ids])

dds <- DESeqDataSetFromMatrix(counts[rownames(counts) %in% tx2gene$SYMBOL,samp.ids],
                              final.clin[samp.ids,],
                              design=~Condition. + Batch + WksGest + InfSex +
                                  Race + PrePregBMI + Smoking + FDELTYPE +
                                  Labor.initiation + Labor.no.labor)

keep <- rowSums(counts(dds)) >= 500 & rowMeans(counts(dds)==0) < 0.8
## keep <- rowMeans(counts(dds)) > 10
dds <- dds[keep,]

dim(dds)

## colMeans(counts(dds)==0)

library(sva)

combat.counts <- ComBat_seq(counts(dds),
                            batch=colData(dds)$Batch,
                            covar_mod=model.matrix(~ Condition. + WksGest + InfSex +
                                                       Race + PrePregBMI + Smoking +
                                                       FDELTYPE + Labor.initiation,
                                                   colData(dds)))

dim(dds)

write.csv(t(combat.counts), 'data/combat_rna_counts.csv')

## combat.counts <- t(read.csv('data/combat_rna_counts.csv', row.names=1))

## samp.ids <- colnames(combat.counts)

dds <- DESeqDataSetFromMatrix(combat.counts[,samp.ids],
                              final.clin[samp.ids,],
                              design=~Condition. + WksGest + InfSex + Race +
                                  PrePregBMI + Smoking + FDELTYPE +
                                  Labor.initiation)


## library("BiocParallel")
## register(MulticoreParam(14))

dds <- DESeq(dds)

plotDispEsts(dds)

vsd <- vst(dds, blind=FALSE)

meanSdPlot(assay(vsd))

write.csv(t(assay(vsd)), 'data/combat_vst_rna_expression.csv')


dds <- DESeqDataSetFromMatrix(combat.counts[,samp.ids],
                              final.clin[samp.ids,],
                              design=~Condition. + WksGest + InfSex + Race +
                                  PrePregBMI + Smoking + FDELTYPE +
                                  Labor.initiation + Labor.no.labor)


## library("BiocParallel")
## register(MulticoreParam(14))

dds <- DESeq(dds)

plotDispEsts(dds)

vsd <- vst(dds, blind=FALSE)

meanSdPlot(assay(vsd))


## scaled.sig <- scale(t(assay(vsd)[unique(all.de.genes),]))

## group.meta <- colData(dds)

## annot_df <- data.frame(Condition = as.factor(group.meta$Condition.),
##                        InfantSex = as.factor(group.meta$InfSex),
##                        GestationalAge = as.factor(group.meta$DiscWksGest), 
##                        row.names=rownames(group.meta))

## dist_mat <- vegan::vegdist(scaled.vsd, method='euclidean')
## dim(dist_mat)
## pheatmap(scaled.sig, annotation_col=annot_df)

## de.all.pheat <- pheatmap(t(scaled.sig), breaks = seq(-range, range, length.out = 101),
##                          annotation_col=annot_df, main="All differentially expressed genes",
##                          ## cluster_rows=sum(de$adj.P.Val < 0.05)>1,
##                          clustering_method='ward.D', clustering_distance_rows = "manhattan",
##                          clustering_distance_cols = "manhattan",
##                          labels_row=unique(all.de.genes))

## plot(de.all.pheat$tree_col)

## rect.hclust(de.all.pheat$tree_col, 3)

## clusts <- cutree(de.all.pheat$tree_col, 3)

## table(colData(dds)[names(clusts), 'Condition.'], clusts)
## chisq.test(table(colData(dds)[names(clusts), 'Condition.'], clusts))

## res.nbclust <- NbClust(scaled.sig, distance = "manhattan",
##                   min.nc = 2, max.nc = 15, 
##                   method = "ward.D", index ="silhouette")

## res.nbclust

## library(ClusterR)

## opt_gmm = Optimal_Clusters_GMM(scaled.sig, max_clusters = 20, criterion = "BIC",
##                                seed_mode = "random_subset",
##                                km_iter = 10, em_iter = 5, var_floor = 1e-5,
##                                plot_data = T)

## gmm <- GMM(scaled.sig, 10, seed_mode = "random_subset",
##            km_iter = 10, em_iter = 5, verbose = T, var_floor=1e-5)          

## # predict centroids, covariance matrix and weights
## gmm.clusts <- predict(gmm, scaled.sig)

## table(colData(dds)[names(clusts), 'Condition.'], gmm.clusts)
## chisq.test(table(colData(dds)[names(clusts), 'Condition.'], gmm.clusts))

## avg_sil <- function(k) {
##   km.res <- kmeans(scaled.sig, centers = k, nstart = 25)
##   ss <- silhouette(km.res$cluster, dist(scaled.sig))
##   mean(ss[, 3])
## }

## # Compute and plot wss for k = 2 to k = 15
## k.values <- 2:15

## # extract avg silhouette for 2-15 clusters
## avg_sil_values <- map_dbl(k.values, avg_sil)

## plot(k.values, avg_sil_values,
##        type = "b", pch = 19, frame = FALSE, 
##        xlab = "Number of clusters K",
##        ylab = "Average Silhouettes")

## km.clust <- kmeans(scaled.sig, 6)

## table(colData(dds)[names(clusts), 'Condition.'], km.clust$cluster)
## chisq.test(table(colData(dds)[names(clusts), 'Condition.'], km.clust$cluster))

## cor_mat <- cor(log2(assay(dds, 'mu')))
## dim(cor_mat)
## pheatmap(cor_mat, annotation_col=annot_df)

groups <- c('FGR+HDP', 'Control')

comparisons <- list(c('Control PTD', 'Control'),
                    c('FGR', 'Control'),
                    c('Severe PE', 'Control PTD'),
                    c('Severe PE', 'Control'),
                    c('PTD', 'Control PTD'),
                    c('PTD', 'Control'),
                    c('FGR+HDP', 'Control PTD'),
                    c('FGR+HDP', 'Control'))

library(RColorBrewer)
brewer.colors <- brewer.pal(n=12, "Paired")## )(10)

brewer.colors <- c(brewer.colors[1:2], brewer.colors[11], brewer.colors[5:6], brewer.colors[8], brewer.colors[3:4], brewer.colors[9:10])

annot.colors <- list()
annot.colors[['Diagnosis']] <- brewer.colors[1:6]
names(annot.colors[['Diagnosis']]) <- c('Control', 'Control PTD', 'PTD',
                                        'FGR', 'FGR+HDP', 'Severe PE')

diag.map <- c('Control', 'Control PTD', 'PTD',
              'FGR', 'FGR+HDP', 'Severe PE')

names(diag.map) <- c('Control', 'Control.PTD', 'PTD',
                     'FGR', 'FGR.HDP', 'Severe.PE')

## all.de.genes <- c()

for (groups in comparisons) {

    print(groups)
    
    groups <- make.names(groups)

    res <- results(dds, contrast=c('Condition.', groups), alpha=0.05)
    ## lfcShrink(dds, contrast=c('Condition.', groups), res=res, type='ashr')
    print(summary(res))

    ## resLFC <- lfcShrink(dds, contrast=c('Condition.', groups), res=res, type='ashr')
    ## print(summary(resLFC))
    
    de <- res[order(res$pvalue),] %>% as.data.frame
    de <- mutate(de, Gene=rownames(de), .before=1) %>% filter(baseMean > 0)
    ## de[!is.na(de$padj) & de$padj<0.05,]
    
    ## vsd <- vst(dds, blind=FALSE)

    group.meta <- colData(dds)[colData(dds)[,'Condition.'] %in% groups,]

    ## all.de.genes <- c(all.de.genes, de[!is.na(de$padj) & de$padj<0.05,'Gene'])

    if (nrow(de[!is.na(de$padj) & de$padj<0.05,])>1) {
        
        if (min(nrow(de[!is.na(de$padj) & de$padj<0.05,]), 500) == 500) {
            col.labels <- de[1:500,'Gene']
        } else {
            col.labels <- de[!is.na(de$padj) & de$padj<0.05,'Gene']
        }

        scaled.sig <- scale(t(assay(vsd))[rownames(group.meta), col.labels])
            
        
        range <- max(abs(scaled.sig))

        print(range)

        range <- min(max(abs(scaled.sig)),3)
        
        annot_df <- data.frame(Diagnosis = factor(diag.map[as.character(group.meta$Condition.)],
                                                  levels=diag.map[groups]),
                               `Infant Sex` = factor(group.meta$InfSex),
                               `Gestational Age` = clin[rownames(group.meta), 'WksGest'], 
                               row.names=rownames(group.meta))

        temp.annot.cols <- list()
        for (g in groups) {
            temp.annot.cols[['Diagnosis']][diag.map[g]] <- annot.colors[['Diagnosis']][diag.map[g]]
        }
       
        pdf(paste('plots-cond3/', str_replace(groups[1], ' ', '-'), '_',
                  str_replace(groups[2], ' ', '.'),
                  '_RNA_diffExp_heatmap.pdf',sep=''),
            width=10,
            height=3 + 0.12 * min(nrow(de[!is.na(de$padj) & de$padj<0.05,]), 500))

        pheatmap(t(scaled.sig), breaks = seq(-range, range, length.out = 101),
                 annotation_col=annot_df, main=paste(diag.map[groups[1]], 'vs.', diag.map[groups[2]]),
                 ## cluster_rows=sum(de$adj.P.Val < 0.05)>1,
                 clustering_method='ward.D', clustering_distance_rows = "euclidean",
                 clustering_distance_cols = "euclidean", annotation_colors=temp.annot.cols,
                 labels_row=col.labels, show_colnames=F)
        dev.off()

        png(paste('plots-cond3/', str_replace(groups[1], ' ', '-'), '_',
                  str_replace(groups[2], ' ', '.'),
                  '_RNA_diffExp_heatmap.png',sep=''),
            width=10,
            height=3 + 0.12 * min(nrow(de[!is.na(de$padj) & de$padj<0.05,]), 500),
            units='in', res=400)

        pheatmap(t(scaled.sig), breaks = seq(-range, range, length.out = 101),
                 annotation_col=annot_df, main=paste(diag.map[groups[1]], 'vs.', diag.map[groups[2]]),
                 ## cluster_rows=sum(de$adj.P.Val < 0.05)>1,
                 clustering_method='ward.D', clustering_distance_rows = "euclidean",
                 clustering_distance_cols = "euclidean", annotation_colors=temp.annot.cols,
                 labels_row=col.labels, show_colnames=F)
        dev.off()


        if (min(nrow(de[!is.na(de$padj) & de$padj<0.05,]), 25) == 25) {
            col.labels <- de[1:25,'Gene']
        } else {
            col.labels <- de[!is.na(de$padj) & de$padj<0.05,'Gene']
        }

        scaled.sig <- scale(t(assay(vsd))[rownames(group.meta), col.labels])
            
        
        range <- max(abs(scaled.sig))

        print(range)

        range <- min(max(abs(scaled.sig)),3)
        
        pdf(paste('plots-cond3/', str_replace(groups[1], ' ', '-'), '_',
                  str_replace(groups[2], ' ', '.'),
                  '_RNA_diffExp_top25_heatmap_unlabeled.pdf',sep=''),
            width=8,
            height=7)

        pheatmap(t(scaled.sig), breaks = seq(-range, range, length.out = 101),
                 annotation_col=annot_df, main=paste(diag.map[groups[1]], 'vs.', diag.map[groups[2]]),
                 ## cluster_rows=sum(de$adj.P.Val < 0.05)>1,
                 clustering_method='ward.D', clustering_distance_rows = "euclidean",
                 clustering_distance_cols = "euclidean", annotation_colors=temp.annot.cols,
                 labels_row=col.labels, show_colnames=F, show_rownames=F)
        dev.off()

        png(paste('plots-cond3/', str_replace(groups[1], ' ', '-'), '_',
                  str_replace(groups[2], ' ', '.'),
                  '_RNA_diffExp_top25_heatmap_unlabeled.png',sep=''),
            width=8,
            height=7,
            units='in', res=400)

        pheat <- pheatmap(t(scaled.sig), breaks = seq(-range, range, length.out = 101),
                 annotation_col=annot_df, main=paste(diag.map[groups[1]], 'vs.', diag.map[groups[2]]),
                 ## cluster_rows=sum(de$adj.P.Val < 0.05)>1,
                 clustering_method='ward.D', clustering_distance_rows = "euclidean",
                 clustering_distance_cols = "euclidean", annotation_colors=temp.annot.cols,
                 labels_row=col.labels, show_colnames=F, show_rownames=F)

        write.csv(data.frame(Labels=col.labels[pheat$tree_row[['order']]]),
                  paste('plots-cond3/', str_replace(groups[1], ' ', '-'), '_',
                        str_replace(groups[2], ' ', '.'),
                        '_RNA_diffExp_top25_heatmap_row_labels.csv',sep=''),
                  row.names=F)
        
        dev.off()

    }

    ## de[!is.na(de$padj) & de$padj<0.05,]

    write.csv(de, paste('DE-cond3/diffExpRNA_', str_replace(groups[1], ' ', '-'), '_',
                        str_replace(groups[2], ' ', '.'), '.csv', sep=''))
    write.csv(de[!is.na(de$padj) & de$padj<0.05,],
              paste('DE-cond3/diffExpRNA_', str_replace(groups[1], ' ', '-'), '_',
                    str_replace(groups[2], ' ', '.'), '_FDR05.csv', sep=''))
    write.csv(de[!is.na(de$padj) & de$padj<0.1,],
              paste('DE-cond3/diffExpRNA_', str_replace(groups[1], ' ', '-'), '_',
                    str_replace(groups[2], ' ', '.'), '_FDR1.csv', sep=''))
    
}

plotMA(res)

resultsNames(dds)

resLFC <- lfcShrink(dds, contrast= c('Condition.', 'FGR', 'Control'), type='ashr')

plotMA(resLFC)

summary(resLFC)
head(as.data.frame(res[order(res$pvalue),]), 45)
head(as.data.frame(resLFC[order(resLFC$pvalue),]), 45)

## vsd <- vst(dds, blind=FALSE)
meanSdPlot(assay(vsd))

pcaData <- plotPCA(vsd, 'Batch', returnData=TRUE)

summary(lm(pcaData$PC1 ~ Batch + Condition. + WksGest + InfSex + Race +
               PrePregBMI + Smoking + FDELTYPE + Labor.initiation, colData(dds)))
summary(lm(pcaData$PC2 ~ Batch + Condition. + WksGest + InfSex + Race +
               PrePregBMI + Smoking + FDELTYPE + Labor.initiation, colData(dds)))

hist(pcaData$PC1)
hist(pcaData$PC2)

boxplot(pcaData$PC1, range=2.22)
boxplot(pcaData$PC2, range=2.22)


box.pc1 <- boxplot.stats(pcaData$PC1, coef=2.22)
box.pc2 <- boxplot.stats(pcaData$PC2, coef=2.22)


ggplot(pcaData, aes(x=PC1, y=PC2, col=Batch)) +
    geom_point() +
    geom_vline(xintercept=mean(pcaData$PC1)-4*sd(pcaData$PC1)) +
    geom_vline(xintercept=mean(pcaData$PC1)+4*sd(pcaData$PC1)) +
    geom_hline(yintercept=mean(pcaData$PC2)-4*sd(pcaData$PC2)) +
    geom_hline(yintercept=mean(pcaData$PC2)+4*sd(pcaData$PC2)) +
    ggtitle('PC1 vs PC2 by Batch') +
    theme_bw()

ggplot(pcaData, aes(x=PC1, y=PC2, col=Condition.)) +
    geom_point() +
    geom_vline(xintercept=median(pcaData$PC1)-3*mad(pcaData$PC1, constant=1)) +
    geom_vline(xintercept=median(pcaData$PC1)+3*mad(pcaData$PC1, constant=1)) +
    geom_hline(yintercept=median(pcaData$PC2)-3*mad(pcaData$PC2, constant=1)) +
    geom_hline(yintercept=median(pcaData$PC2)+3*mad(pcaData$PC2, constant=1)) +
    ggtitle('PC1 vs PC2 by Batch') +
    theme_bw()


ggplot(pcaData, aes(x=PC1, y=PC2, col=Condition.)) +
    geom_point() +
    geom_vline(xintercept=box.pc1$stats[1]) +
    geom_vline(xintercept=box.pc1$stats[5]) +
    geom_hline(yintercept=box.pc2$stats[1]) +
    geom_hline(yintercept=box.pc2$stats[5]) +
    ggtitle('PC1 vs PC2 by Batch') +
    theme_bw()


pdf('uncorrected_iqrOutLines_pcaPlots_RNAseq.pdf', width=7, height=6.5)
ggplot(pcaData, aes(x=PC1, y=PC2, col=Condition.)) +
    geom_point() +
    geom_vline(xintercept=box.pc1$stats[1]) +
    geom_vline(xintercept=box.pc1$stats[5]) +
    geom_hline(yintercept=box.pc2$stats[1]) +
    geom_hline(yintercept=box.pc2$stats[5]) +
    ggtitle('PC1 vs PC2 by Batch') +
    theme_bw()
dev.off()

## outliers <- pcaData %>% filter(abs(PC1) > 3*sd(PC1) | abs(PC2) > 3*sd(PC2))
## rownames(outliers)

outliers <- pcaData %>% filter(PC1 < box.pc1$stats[1] | PC1 > box.pc1$stats[5] |
                               PC2 < box.pc2$stats[1] | PC2 > box.pc2$stats[5])
outliers
outliers <- rownames(outliers)

paste(rownames(outliers), collapse="', '")

## outliers <- c('MJ.0357', 'MJ.0376', 'MJ.0505', 'MJ.0581', 'MJ.0868', 'MJ.0978', 'Mini.DP.189', 'Mini.DP.237')

## outliers <- c('MJ.0042', 'MJ.0357', 'MJ.0376', 'MJ.0472', 'MJ.0505', 'MJ.0581', 'MJ.0587', 'MJ.0766', 'MJ.0868', 'MJ.0978', 'Mini.DP.109', 'Mini.DP.115', 'Mini.DP.139', 'Mini.DP.187', 'Mini.DP.189', 'Mini.DP.232', 'Mini.DP.237')

## outliers <- c('MJ.0042', 'MJ.0357', 'MJ.0376', 'MJ.0472', 'MJ.0505', 'MJ.0581', 'MJ.0587', 'MJ.0868', 'MJ.0978', 'Mini.DP.109', 'Mini.DP.139', 'Mini.DP.189', 'Mini.DP.232', 'Mini.DP.237')

library(EnvStats)
test <- rosnerTest(pcaData$PC2, k=length(outliers))
test$all.stats
test <- rosnerTest(pcaData$PC2, k=length(outliers))
outliers <- rownames(pcaData)[test$all.stats$Obs.Num[test$all.stats$Outlier]]

paste(outliers, collapse="', '")

outliers <- c('MJ.0868', 'Mini.DP.189', 'MJ.0581', 'Mini.DP.237', 'MJ.0505', 'MJ.0357', 'MJ.0376', 'Mini.DP.139', 'Mini.DP.109', 'MJ.0978', 'MJ.0587', 'Mini.DP.232', 'MJ.0472')

pdf('corrected_QC_pcaPlots_RNAseq.pdf', width=7, height=5.5)
plotPCA(vsd, 'Batch') +
    theme_bw() +
    guides(color=guide_legend(title="Batch"))

plotPCA(vsd, 'Condition.') +
    theme_bw() +
    guides(color=guide_legend(title="Condition"))

plotPCA(vsd, 'DiscWksGest') +
    theme_bw() +
    ## guides(color=guide_legend(title="Gestational Age"))
    labs(color="Gestational Age")

plotPCA(vsd, 'InfSex') +
    theme_bw()  +
    guides(color=guide_legend(title="Infant Sex"))

plotPCA(vsd, 'Race') +
    theme_bw() +
    guides(color=guide_legend(title="Race"))

plotPCA(vsd, 'DiscPrePregBMI') +
    theme_bw() +
    labs(color="Pre-Pregnancy BMI")

plotPCA(vsd, 'Smoking') +
    theme_bw() +
    guides(color=guide_legend(title="Smoking"))

plotPCA(vsd, 'FDELTYPE') +
    theme_bw()+
    guides(color=guide_legend(title="Delivery Type"))

plotPCA(vsd, 'Labor.initiation') +
    theme_bw() +
    guides(color=guide_legend(title="Labor Initiation"))

## plotPCA(vsd, 'zeroFrac') +
##     theme_bw() +
##     labs(color="Zero Fraction")

dev.off()


colData(vsd)$FLT1overPGF <- assay(vsd)['FLT1',]/assay(vsd)['PGF',]

library(rstatix)

FLT1overPGF.tests <- data.frame()

for (groups in comparisons) {
    dataset_group <- colData(vsd)[colData(vsd)$Condition. %in% make.names(groups),] %>% as.data.frame
    dataset_group$Condition. <- factor(dataset_group$Condition.)
    FLT1overPGF.tests <- rbind(FLT1overPGF.tests,
                               wilcox_test(dataset_group, FLT1overPGF ~ Condition.))
}

FLT1overPGF.tests <- FLT1overPGF.tests %>% adjust_pvalue(method = 'holm') %>%
    filter(p.adj < 0.05) %>%
    mutate(y.position = c(1.7, 1.8, 1.9, 2.0, 2.1, 2.2))

library(ggpubr)
ggviolin(as.data.frame(colData(vsd)), 'Condition.', 'FLT1overPGF',
         fill='Condition.', add = "boxplot") +
    stat_pvalue_manual(FLT1overPGF.tests)

ggsave('RNA-FLT1overPGF-dists.pdf',
       ggviolin(as.data.frame(colData(vsd)), 'Condition.', 'FLT1overPGF',
                fill='Condition.', ylab='FLT1 / PGF', add = "boxplot") +
       stat_pvalue_manual(FLT1overPGF.tests), height=6, width=6)

ggsave('RNA-FLT1overPGF-dists.png',
       ggviolin(as.data.frame(colData(vsd)), 'Condition.', 'FLT1overPGF',
                fill='Condition.', ylab='FLT1 / PGF', add = "boxplot") +
       stat_pvalue_manual(FLT1overPGF.tests), height=6, width=6)



res <- results(dds, contrast=c('InfSex', c('M', 'F')), alpha=0.05)
print(summary(res))

de <- res[order(res$pvalue),] %>% as.data.frame
de <- mutate(de, Gene=rownames(de), .before=1) %>% filter(baseMean > 0)

head(de, n=50)

de[!is.na(de$padj) & de$padj < 0.05,]

de[!is.na(de$padj) & de$padj < 0.05 & abs(de$log2FoldChange)>1,]


sex.de.genes <- tx2gene[tx2gene$SYMBOL %in% rownames(de[!is.na(de$padj) & de$padj < 0.1,]),]

sex.de.genes <- sex.de.genes %>% filter(!SYMBOL %in% c('ARHGAP11B', 'ERICH2', 'POLA2')) %>% as.data.frame

rownames(sex.de.genes) <- sex.de.genes$SYMBOL

de.x.genes <- (sex.de.genes[rownames(de[!is.na(de$padj) & de$padj < 0.01,]),] %>% filter(CHR=='X'))$SYMBOL

paste(de.x.genes, collapse="', '")

## de.x.genes <- c('XIST', 'TSIX', 'KDM6A', 'ZFX', 'KDM5C', 'EIF2S3', 'SMC1A', 'DDX3X', 'PUDP', 'EIF1AX', 'UBA1', 'CD99', 'VAMP7', 'STS', 'NUDT10', 'NAA10', 'CDK16', 'BRCC3', 'BCLAF3', 'DIPK2B', 'JPX', 'INE1', 'HSD17B10', 'MIR4767', 'CA5BP1', 'CHM', 'IL13RA2', 'TXLNG', 'ZRSR2', 'ASMTL', 'CA5B', 'ARMCX3', 'CXorf38', 'MBTPS2', 'LAS1L', 'USP11', 'FHL1', 'NEXMIF', 'HDX', 'SERTM2', 'SMS', 'SCML1', 'CD99P1', 'ARMCX6', 'RPS4X', 'LOC101928228', 'ZMAT1', 'PDK3', 'HDAC8', 'YIPF6', 'PHF8', 'NKAP', 'LINC00632', 'AMMECR1', 'ARMCX1', 'TMEM47', 'RPS6KA6', 'TRAPPC2', 'PPP1R3F', 'HMGB1P12', 'ZFP92', 'FOXP3', 'RAB40AL', 'IRS4', 'PAICSP7', 'LINC01284', 'KLHL13', 'GABRE', 'HCFC1', 'GSPT2', 'MIR6895', 'ARX', 'ZBED1', 'KLHL4', 'ALAS2', 'ZXDA', 'LOC100420869', 'YY2', 'FAM133A', 'AP1S2', 'MXRA5', 'CYTH1P1', 'MORC4', 'DHRSX', 'CCNB3', 'TFDP3', 'OFD1', 'FIRRE', 'TMEM187', 'BRWD3', 'RBM41', 'FUNDC2')
