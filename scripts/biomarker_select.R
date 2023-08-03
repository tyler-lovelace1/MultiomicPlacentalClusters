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

clustMarkersWilcox <- function(data, clusters, nMarkers=10, pval=0.05) {
    clusts <- as.numeric(sort(unique(clusters)))

    clust.counts <- table(clusters)

    means <- matrix(NA, nrow=ncol(data), ncol=length(clusts))
    rownames(means) <- colnames(data)
    colnames(means) <- paste('Expression', clusts, sep='.')

    otherMeans <- matrix(NA, nrow=ncol(data), ncol=length(clusts))
    rownames(otherMeans) <- colnames(data)
    colnames(otherMeans) <- paste('Expression', clusts, sep='.')

    log2fc <- matrix(NA, nrow=ncol(data), ncol=length(clusts))
    rownames(log2fc) <- colnames(data)
    colnames(log2fc) <- paste('log2FoldChange', clusts, sep='.')

    stats <- matrix(NA, nrow=ncol(data), ncol=length(clusts))
    rownames(stats) <- colnames(data)
    colnames(stats) <- paste('W', clusts, sep='.')
    
    pvals <- matrix(NA, nrow=ncol(data), ncol=length(clusts))
    rownames(pvals) <- colnames(data)
    colnames(pvals) <- paste('p.value', clusts)

    adj.pvals <- matrix(NA, nrow=ncol(data), ncol=length(clusts))
    rownames(adj.pvals) <- colnames(data)
    colnames(adj.pvals) <- paste('adj.p.value', clusts)
    
    for (c in clusts) {
        for (feat in 1:ncol(data)) {
            res <- wilcox.test(data[,feat] ~ factor(clusters==c, levels=c(TRUE, FALSE)),
                               alternative='two.sided')
            ## aucs[feat,c] <- auc(clusters==c, data[,feat], direction='<'), alternative='greater'
            means[feat,c] <- mean(data[clusters==c,feat])
            otherMeans[feat,c] <- mean(data[clusters!=c,feat])
            log2fc[feat,c] <- means[feat,c] - otherMeans[feat,c]
            stats[feat,c] <- res$statistic
            pvals[feat,c] <- res$p.value
        }
        adj.pvals[,c] <- p.adjust(pvals[,c], method='fdr')
    }

    marker.stats <- data.frame()
    signif.marker.stats <- data.frame()

    aucs <- stats

    colnames(aucs) <- paste('AUC', clusts, sep='.')

    for (c in clusts) {
        aucs[,c] <- aucs[,c] / (clust.counts[c] * (sum(clust.counts) - clust.counts[c]))
        aucs[,c] <- pmax(abs(aucs[,c]-0.5) + 0.5, aucs[,c])
    }

    for (c in clusts) {        
        o <- order(aucs[,c], decreasing=T)[1:nMarkers]
        marker.stats <- rbind(marker.stats, cbind(Feature=rownames(stats)[o],
                                                  Cluster=rep(c,nMarkers),
                                                  ExpressionInCluster=means[o,c],
                                                  ExpressionOutCluster=otherMeans[o,c],
                                                  log2FoldChange=log2fc[o,c],
                                                  W=stats[o,c],
                                                  AUC=aucs[o,c],
                                                  maxOtherAUC=apply(
                                                      aucs[o,-c], 1,
                                                      function(x) {
                                                          max(abs(x-0.5))+0.5
                                                      }),
                                                  p.value=pvals[o,c],
                                                  adj.p.value=adj.pvals[o,c]))

        o <- order(aucs[,c], decreasing=T)
        o <- o[adj.pvals[o,c] < pval]

        if (length(o) > nMarkers) {
            o <- o[1:nMarkers]
        }

        signif.marker.stats <- rbind(signif.marker.stats,
                                     cbind(Feature=rownames(stats)[o],
                                           Cluster=rep(c,length(o)),
                                           ExpressionInCluster=means[o,c],
                                           ExpressionOutCluster=otherMeans[o,c],
                                           log2FoldChange=log2fc[o,c],
                                           W=stats[o,c],
                                           AUC=aucs[o,c],
                                           maxOtherAUC=apply(
                                               aucs[o,-c], 1,
                                               function(x) {
                                                   max(abs(x-0.5))+0.5
                                               }),
                                           p.value=pvals[o,c],
                                           adj.p.value=adj.pvals[o,c]))
    }

    marker.stats <- marker.stats %>%
        mutate_at(c('ExpressionInCluster', 'ExpressionOutCluster', 'W', 'AUC',
                    'maxOtherAUC', 'p.value', 'adj.p.value', 'log2FoldChange'), as.numeric)

    signif.marker.stats <- signif.marker.stats %>%
        mutate_at(c('ExpressionInCluster', 'ExpressionOutCluster', 'W', 'AUC',
                    'maxOtherAUC', 'p.value', 'adj.p.value', 'log2FoldChange'), as.numeric)

    list(all=marker.stats, signif=signif.marker.stats)
}


rna <- read.csv('data/combat_vst_rna_expression.csv', row.names=1)
## colnames(rna) <- paste0('RNA_', colnames(rna))
head(rna)

mirna <- read.csv('data/combat_vst_mirna_expression.csv', row.names=1)
## colnames(mirna) <- paste0('miRNA_', colnames(mirna))
head(mirna)

rna.means <- apply(rna[samp.ids,], 2, mean)
rna.sds <- apply(rna[samp.ids,], 2, sd)
rna.kwp <- apply(rna[samp.ids,], 2, function(x) kruskal.test(x, clin$hvar.Clusters)$p.value)

ggplot(data.frame(Mean=rna.means, SD=rna.sds, KW=rna.kwp), aes(Mean, SD)) +
    geom_point() +
    geom_smooth() +
    theme_bw()

rna.means.p <- 1-ecdf(rna.means)(rna.means) + 1/ncol(rna)
rna.sds.p <- ecdf(rna.sds)(rna.sds)
## rna.kwstat.p <- 1-ecdf(rna.kwstat)(rna.kwstat)

rna.hkstat <- -2*rowSums(log(data.frame(Mean=rna.means.p, SD=rna.sds.p)))
names(rna.hkstat) <- colnames(rna)
rna.hk.p <- pchisq(rna.hkstat, 4, lower.tail=F)

rna.fisherstat <- -2*rowSums(log(data.frame(HK=rna.hk.p, KW=1-rna.kwp)))
names(rna.fisherstat) <- colnames(rna)

rna.fisherp <- pchisq(rna.fisherstat, 4, lower.tail=F)

housekeeping.rna <- names(sort(rna.fisherp[rna.means>med.rna & rna.kwp>0.5]))[1:10]

library(reshape2)

plotdf <- melt(data.frame(Clusters=clin$hvar.Clusters, rna[samp.ids,housekeeping.rna]), value.name='Expression', variable.name='RNA')

library(ggpubr)

ggplot(plotdf, aes(x=1, y=Expression, fill=Clusters)) +
    geom_violin() +
    stat_compare_means() +
    facet_wrap(~RNA, ncol=5) +
    scale_x_continuous(labels=NULL, breaks=NULL, name=NULL) +
    theme_bw()

ggsave('housekeeping_candidates_rna.pdf', width=11, height=7.5)
ggsave('housekeeping_candidates_rna.png', width=11, height=7.5, dpi=400)

mirna.means <- apply(mirna[samp.ids,], 2, mean)
mirna.sds <- apply(mirna[samp.ids,], 2, sd)
mirna.kwp <- apply(mirna[samp.ids,], 2, function(x) kruskal.test(x, clin$hvar.Clusters)$p.value)

ggplot(data.frame(Mean=mirna.means, SD=mirna.sds, KW=mirna.kwp), aes(Mean, SD)) +
    geom_point() +
    geom_smooth() +
    theme_bw()

mirna.means.p <- 1-ecdf(mirna.means)(mirna.means) + 1/ncol(mirna)
mirna.sds.p <- ecdf(mirna.sds)(mirna.sds)
## mirna.kwstat.p <- 1-ecdf(mirna.kwstat)(mirna.kwstat)

mirna.hkstat <- -2*rowSums(log(data.frame(Mean=mirna.means.p, SD=mirna.sds.p)))
names(mirna.hkstat) <- colnames(mirna)
mirna.hk.p <- pchisq(mirna.hkstat, 4, lower.tail=F)

mirna.fisherstat <- -2*rowSums(log(data.frame(HK=mirna.hk.p, KW=1-mirna.kwp)))
names(mirna.fisherstat) <- colnames(mirna)

mirna.fisherp <- pchisq(mirna.fisherstat, 4, lower.tail=F)

housekeeping.mirna <- names(sort(mirna.fisherp[mirna.means>med.mirna & mirna.kwp>0.5]))[1:10]

library(reshape2)

plotdf <- melt(data.frame(Clusters=clin$hvar.Clusters, mirna[samp.ids,housekeeping.mirna]), value.name='Expression', variable.name='miRNA')

library(ggpubr)

ggplot(plotdf, aes(x=1, y=Expression, fill=Clusters)) +
    geom_violin() +
    stat_compare_means() +
    facet_wrap(~miRNA, ncol=5) +
    scale_x_continuous(labels=NULL, breaks=NULL, name=NULL) +
    theme_bw()

ggsave('housekeeping_candidates_mirna.pdf', width=11, height=7.5)
ggsave('housekeeping_candidates_mirna.png', width=11, height=7.5, dpi=400)

med.rna <- quantile(as.matrix(rna), 0.75)
med.mirna <- quantile(as.matrix(mirna), 0.75)

write.csv(cbind(RNA=housekeeping.rna, miRNA=housekeeping.mirna), 'housekeeping_candidates.csv')

clin <- read_excel('Primary clinical variables doppler.xlsx', 1) %>% as.data.frame
colnames(clin) <- make.names(colnames(clin))
clin <- clin[-1,]
rownames(clin) <- make.names(clin$Study_ID.ID)
head(clin)

## clin[incomp.samp.ids,'pred.hvar.Clusters'] <- pred.final

## write.csv(clin, 'doppler_clinical.csv')

clin <- clin %>% mutate_if(function(x) mean(is.na(as.numeric(x))) < 0.5, as.numeric)

clin$Condition. <- ifelse(clin$Condition. %in% c('PTD','Severe PE','FGR','Control',
                                                 'Control PTD'), clin$Condition., 'FGR+HDP')

clin$PrePregBMI <- clin$PrePregWt_Kg / clin$Height_Meters.HeightMeters^2

clin$Smoking <- ifelse(str_detect(clin$Smoke, 'Yes'), 'Yes', NA)
clin$Smoking <- ifelse(str_detect(clin$Smoke, 'Never'), 'No', clin$Smoking)

clin$Race[!clin$Race %in% c('W', 'B')] <- 'other'
clin$FDELTYPE[clin$FDELTYPE==0] <- NA

clin <- clin %>% dplyr::select(c('Condition.', 'WksGest', 'InfSex', 'Race', 'PrePregBMI',
                                 'Smoking', 'FDELTYPE', 'Labor.initiation', 'hvar.Clusters'))

clin <- clin %>% mutate_at(c('Condition.', 'InfSex', 'Race', 'Smoking', 'FDELTYPE',
                             'Labor.initiation', 'hvar.Clusters'), factor)

#### Pathology data (slides)

pathology.slides <- read_excel('path features (from slides) for analysis.xlsx' , 1, na='n/a') %>% as.data.frame
colnames(pathology.slides) <- make.names(colnames(pathology.slides))
colnames(pathology.slides)[6] <- "Syncytial.knots"
rownames(pathology.slides) <- make.names(pathology.slides$Study.ID)
pathology.slides <- pathology.slides[,4:6]
pathology.slides <- pathology.slides %>% mutate_all(factor)

clin <- clin[!is.na(clin$hvar.Clusters),]

samp.ids <- rownames(clin)

med.rna <- quantile(as.matrix(rna), 0.75)
med.mirna <- quantile(as.matrix(mirna), 0.75)

mirna.markers <- clustMarkersWilcox(
    mirna[samp.ids[clin$Condition. %in% c('Severe PE', 'FGR+HDP')],],
    clin$hvar.Clusters[clin$Condition. %in% c('Severe PE', 'FGR+HDP')],
    nMarkers=ncol(mirna))

mirna.marker.list <- lapply(mirna.markers$signif %>% group_by(Cluster) %>%
                            filter(ExpressionInCluster > med.mirna |
                                   ExpressionOutCluster > med.mirna,
                                   log2FoldChange > 0) %>%
                            group_split(),
                            as.data.frame)


rna.markers <- clustMarkersWilcox(
    rna[samp.ids[clin$Condition. %in% c('Severe PE', 'FGR+HDP')],],
    clin$hvar.Clusters[clin$Condition. %in% c('Severe PE', 'FGR+HDP')],
    nMarkers=ncol(rna))

rna.marker.list <- lapply(rna.markers$signif %>% group_by(Cluster) %>%
                          filter(ExpressionInCluster > med.rna |
                                 ExpressionOutCluster > med.rna,
                                 log2FoldChange > 0) %>% group_split(),
                          as.data.frame)

placenta.markers <- read.csv('tissue_category_rna_placenta_Tissue.tsv', sep='\t')

placenta.genes <- placenta.markers$Gene[placenta.markers$RNA.blood.cell.distribution=='Not detected']

placenta.rna.marker.list <- lapply(rna.markers$signif %>% group_by(Cluster) %>%
                                   filter(Feature %in% paste0('RNA_',placenta.genes),
                                          ExpressionInCluster > med.rna |
                                          ExpressionOutCluster > med.rna,
                                          log2FoldChange > 0) %>% group_split(),
                                   as.data.frame)

top5.plac.rna.genes <- do.call(base::c, lapply(placenta.rna.marker.list,
                                               function(x) x$Feature[1:5]))

top5.plac.rna.genes <- na.omit(top5.plac.rna.genes)

top5.plac.rna.expr <- melt(data.frame(Marker=clin$hvar.Clusters,
                                      rna[samp.ids,top5.plac.rna.genes]),
                           variable.name='Gene', value.name='Expression')


top5.plac.rna.expr$Cluster <- c(rep(1, 5*271), rep(2, 5*271), rep(3, 5*271), rep(4, 3*271))
top5.plac.rna.expr$Marker <- ifelse(top5.plac.rna.expr$Marker==top5.plac.rna.expr$Cluster,
                                    'Yes', 'No')

top5.plac.rna.expr <- na.omit(top5.plac.rna.expr)
## top5.plac.rna.expr$Gene[top5.plac.rna.expr$Gene=='RNA_LAMC1.1'] <- 'RNA_LAMC1'
top5.plac.rna.expr$Gene <- factor(gsub('RNA_','',top5.plac.rna.expr$Gene),
                                  levels=c(gsub('RNA_','',na.omit(top5.plac.rna.genes))[-18], 'LAMC1.1'),
                                  labels=c(gsub('RNA_','',na.omit(top5.plac.rna.genes))[-18], 'LAMC1'))

p1 <- ggplot(top5.plac.rna.expr %>% filter(Cluster==1), aes(Gene, Expression, fill=Marker)) +
    ## geom_violin() +
    geom_boxplot() +
    ggtitle('Cluster 1 Markers') +
    theme_bw()

p2 <- ggplot(top5.plac.rna.expr %>% filter(Cluster==2), aes(Gene, Expression, fill=Marker)) +
    ## geom_violin() +
    geom_boxplot() +
    ggtitle('Cluster 2 Markers') +
    theme_bw()

p3 <- ggplot(top5.plac.rna.expr %>% filter(Cluster==3), aes(Gene, Expression, fill=Marker)) +
    ## geom_violin() +
    geom_boxplot() +
    ggtitle('Cluster 3 Markers') +
    theme_bw()

p4 <- ggplot(top5.plac.rna.expr %>% filter(Cluster==4), aes(Gene, Expression, fill=Marker)) +
    ## geom_violin() +
    geom_boxplot() +
    ggtitle('Cluster 4 Markers') +
    theme_bw()

library(cowplot)

plot_grid(p1,p2,p3,p4, nrow=2)
ggsave('marker_candidates/AnyPE_placenta_rna_markers_noBloodGene_topQuartileExpression.png', width=8, height=8, dpi=400)

for (clust in 1:4) {
    write.csv(head(placenta.rna.marker.list[[clust]], 20), paste0('marker_candidates/AnyPE_placenta_rna_markers_noBloodGene_topQuartileExpression_clust', clust, '_top20.csv'), row.names=FALSE)
}


top5.rna.genes <- do.call(base::c, lapply(rna.marker.list,
                                          function(x) x$Feature[1:5]))

top5.rna.expr <- melt(data.frame(Marker=clin$hvar.Clusters,
                                      rna[samp.ids,top5.rna.genes]),
                           variable.name='Gene', value.name='Expression')

top5.rna.expr$Cluster <- c(rep(1, 5*271), rep(2, 5*271), rep(3, 5*271), rep(4, 5*271))
top5.rna.expr$Marker <- ifelse(top5.rna.expr$Marker==top5.rna.expr$Cluster,
                                    'Yes', 'No')

top5.rna.expr$Gene <- factor(gsub('RNA_','',top5.rna.expr$Gene),
                                  levels=gsub('RNA_','',top5.rna.genes))

p1 <- ggplot(top5.rna.expr %>% filter(Cluster==1), aes(Gene, Expression, fill=Marker)) +
    ## geom_violin() +
    geom_boxplot() +
    ggtitle('Cluster 1 Markers') +
    theme_bw()

p2 <- ggplot(top5.rna.expr %>% filter(Cluster==2), aes(Gene, Expression, fill=Marker)) +
    ## geom_violin() +
    geom_boxplot() +
    ggtitle('Cluster 2 Markers') +
    theme_bw()

p3 <- ggplot(top5.rna.expr %>% filter(Cluster==3), aes(Gene, Expression, fill=Marker)) +
    ## geom_violin() +
    geom_boxplot() +
    ggtitle('Cluster 3 Markers') +
    theme_bw()

p4 <- ggplot(top5.rna.expr %>% filter(Cluster==4), aes(Gene, Expression, fill=Marker)) +
    ## geom_violin() +
    geom_boxplot() +
    ggtitle('Cluster 4 Markers') +
    theme_bw()

library(cowplot)

plot_grid(p1,p2,p3,p4, nrow=2)
ggsave('marker_candidates/AnyPE_all_rna_markers_topQuartileExpression.png', width=8, height=8, dpi=400)

for (clust in 1:4) {
    write.csv(head(rna.marker.list[[clust]], 20), paste0('marker_candidates/AnyPE_all_rna_markers_topQuartileExpression_clust', clust, '_top20.csv'), row.names=FALSE)
}


top5.mirna.genes <- do.call(base::c, lapply(mirna.marker.list,
                                          function(x) x$Feature[1:5]))

top5.mirna.expr <- melt(data.frame(Marker=clin$hvar.Clusters,
                                      mirna[samp.ids,top5.mirna.genes]),
                           variable.name='Gene', value.name='Expression')

top5.mirna.expr$Cluster <- c(rep(1, 5*271), rep(2, 5*271), rep(3, 5*271), rep(4, 5*271))
top5.mirna.expr$Marker <- ifelse(top5.mirna.expr$Marker==top5.mirna.expr$Cluster,
                                    'Yes', 'No')

top5.mirna.expr$Gene <- factor(gsub('[.]','-',gsub('miRNA_hsa[.]','',top5.mirna.expr$Gene)),
                               levels=gsub('[.]','-',gsub('miRNA_hsa[.]','',top5.mirna.genes)))

p1 <- ggplot(top5.mirna.expr %>% filter(Cluster==1), aes(Gene, Expression, fill=Marker)) +
    ## geom_violin() +
    geom_boxplot() +
    ggtitle('Cluster 1 Markers') +
    theme_bw()

p2 <- ggplot(top5.mirna.expr %>% filter(Cluster==2), aes(Gene, Expression, fill=Marker)) +
    ## geom_violin() +
    geom_boxplot() +
    ggtitle('Cluster 2 Markers') +
    theme_bw()

p3 <- ggplot(top5.mirna.expr %>% filter(Cluster==3), aes(Gene, Expression, fill=Marker)) +
    ## geom_violin() +
    geom_boxplot() +
    ggtitle('Cluster 3 Markers') +
    theme_bw()

p4 <- ggplot(top5.mirna.expr %>% filter(Cluster==4), aes(Gene, Expression, fill=Marker)) +
    ## geom_violin() +
    geom_boxplot() +
    ggtitle('Cluster 4 Markers') +
    theme_bw()

library(cowplot)

plot_grid(p1,p2,p3,p4, nrow=2)
ggsave('marker_candidates/AnyPE_mirna_markers_topQuartileExpression.png', width=8, height=8, dpi=400)

for (clust in 1:4) {
    write.csv(head(mirna.marker.list[[clust]], 20), paste0('marker_candidates/AnyPE_mirna_markers_topQuartileExpression_clust', clust, '_top20.csv'), row.names=FALSE)
}
