library(SNFtool)
library(stringr)
library(readxl)
library(dplyr)
library(ggplot2)
library(limma)
library(cowplot)
library(scales)
library(ggpubr)
library(pheatmap)
library(RColorBrewer)
library(ComplexHeatmap)

snf.cluster <- function(dataL, C, K, Sigma) {

    dataL.dist <- lapply(dataL, dist, diag=TRUE, upper=TRUE, method='euclidean')

    res <- list()
    
    for (k in K) {
        for (sigma in Sigma) {
            dataL.aff <- lapply(dataL.dist,
                                function(x) affinityMatrix(as.matrix(x),
                                                           K=k, sigma=sigma))

            W <- SNF(dataL.aff, K=k, t=100)
            
            for (c in C) {
                res[[paste0('K', k, '_Sigma', sigma, '_C', c)]] <- spectralClustering(W, c)
            }
        }
    }
    
    res
}


snf.instabs <- function(dataL, C, K, Sigma, snf.full, nSub, subFrac, seed) {
    set.seed(seed)

    n <- nrow(dataL[[1]])
    N <- floor(subFrac * n)

    res <- list()
    
    for (i in 1:nSub) {

        ss <- sample(1:n, N)

        dataLsub <- lapply(dataL, function(x, ii) x[ii,], ii=ss)

        snf.sub <- snf.cluster(dataLsub, C, K, Sigma)

        for (k in K) {
            for (sigma in Sigma) {       
                for (c in C) {
                    name <- paste0('K', k, '_Sigma', sigma, '_C', c)
                    res[[name]] <- c(res[[name]],
                                     aricode::AMI(snf.full[[name]][ss], snf.sub[[name]]))
                }
            }
        }
    }

    for (k in K) {
        for (sigma in Sigma) {       
            for (c in C) {
                name <- paste0('K', k, '_Sigma', sigma, '_C', c)
                res[[paste0(name, '_Mean')]] <- mean(res[[name]])
                res[[paste0(name, '_StDev')]] <- sd(res[[name]])
            }
        }
    }

    res
}


clustMarkersAUC <- function(data, clusters, nMarkers=10, pval=0.05) {
    require(pROC)

    clusts <- as.numeric(sort(unique(clusters)))

    aucs <- matrix(NA, nrow=ncol(data), ncol=length(clusts))
    rownames(aucs) <- colnames(data)
    colnames(aucs) <- paste('AUROC', clusts)
    
    pvals <- matrix(NA, nrow=ncol(data), ncol=length(clusts))
    rownames(pvals) <- colnames(data)
    colnames(pvals) <- paste('p.value', clusts)

    adj.pvals <- matrix(NA, nrow=ncol(data), ncol=length(clusts))
    rownames(adj.pvals) <- colnames(data)
    colnames(adj.pvals) <- paste('adj.p.value', clusts)
    
    for (c in clusts) {
        for (feat in 1:ncol(data)) {
            res <- roc.test(clusters==c, data[,feat], rep(1,length(clusters)), direction='<', alternative='greater')
            ## aucs[feat,c] <- auc(clusters==c, data[,feat], direction='<')
            aucs[feat,c] <- res$roc1$auc
            pvals[feat,c] <- res$p.value
        }
        adj.pvals[,c] <- p.adjust(pvals[,c], method='fdr')
    }

    marker.aucs <- data.frame()
    signif.marker.aucs <- data.frame()

    for (c in clusts) {
        o <- order(aucs[,c], decreasing=T)[1:nMarkers]
        marker.aucs <- rbind(marker.aucs, cbind(Feature=rownames(aucs)[o],
                                                Cluster=rep(c,nMarkers),
                                                aucs[o,],
                                                p.value=pvals[o,c],
                                                adj.p.value=adj.pvals[o,c]))

        o <- order(aucs[,c], decreasing=T)
        o <- o[adj.pvals[o,c] < pval]

        if (length(o) > nMarkers) {
            o <- o[1:nMarkers]
        }

        signif.marker.aucs <- rbind(signif.marker.aucs, cbind(Feature=rownames(aucs)[o],
                                                              Cluster=rep(c,length(o)),
                                                              aucs[o,],
                                                              p.value=pvals[o,c],
                                                              adj.p.value=adj.pvals[o,c]))
    }

    list(all=marker.aucs, signif=signif.marker.aucs)
}

clustMarkersWilcox <- function(data, clusters, nMarkers=10, pval=0.05) {
    clusts <- as.numeric(sort(unique(clusters)))

    clust.counts <- table(clusters)

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
                               alternative='greater')
            ## aucs[feat,c] <- auc(clusters==c, data[,feat], direction='<')
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
    }

    for (c in clusts) {        
        o <- order(stats[,c], decreasing=T)[1:nMarkers]
        marker.stats <- rbind(marker.stats, cbind(Feature=rownames(stats)[o],
                                                  Cluster=rep(c,nMarkers),
                                                  stats[o,],
                                                  aucs[o,],
                                                  p.value=pvals[o,c],
                                                  adj.p.value=adj.pvals[o,c]))

        o <- order(stats[,c], decreasing=T)
        o <- o[adj.pvals[o,c] < pval]

        if (length(o) > nMarkers) {
            o <- o[1:nMarkers]
        }

        signif.marker.stats <- rbind(signif.marker.stats, cbind(Feature=rownames(stats)[o],
                                                                Cluster=rep(c,length(o)),
                                                                stats[o,],
                                                                aucs[o,],
                                                                p.value=pvals[o,c],
                                                                adj.p.value=adj.pvals[o,c]))
    }

    list(all=marker.stats, signif=signif.marker.stats)
}

clustMarkers <- function(data, clusters, nMarkers=10, alpha=0.05, sign=c('positive', 'negative')) {

    if (sign != 'both') {
        sign <- ifelse(sign[1]=='negative', -1, 1)
    }
    clusters <- as.factor(clusters)

    markers <- data.frame()

    for (l in levels(clusters)) {

        comparison <- factor(ifelse(clusters==l, 'clust', 'other'),
                             levels=c('clust', 'other'))

        design <- model.matrix(~0 + comparison)

        colnames(design) <- c('clust', 'other')

        contrast <- makeContrasts(Diff = 'clust - other', levels=design)

        fit <- lmFit(t(data), design=design, maxit=1000)

        contrast_fit<-contrasts.fit(fit, contrast)
        
        ebays_fit<-eBayes(contrast_fit, robust=T)
        ## summary
        print(summary(decideTests(ebays_fit)))
        ## extract DE results
        de<-topTable(ebays_fit, n=ncol(data),
                     adjust.method="fdr", confint=TRUE)

        if (sign != 'both') {
            clust.markers <- de[sign * de$logFC > 0 & de$adj.P.Val < alpha,]
        } else {
            clust.markers <- de[de$adj.P.Val < alpha,]
        }
        
        if (nrow(clust.markers)==0) next

        clust.markers <- cbind(cluster=l, marker=rownames(clust.markers), clust.markers)

        markers <- rbind(markers, clust.markers)

    }

    de.res <- markers %>%
        group_by(cluster) %>%
        arrange(abs(logFC), .by_group=TRUE) %>%
        group_split()

    markers <- markers %>%
        group_by(marker) %>%
        arrange(abs(logFC), .by_group=TRUE) %>%
        filter(row_number()==1) %>%
        group_by(cluster) %>%
        arrange(abs(logFC), .by_group=TRUE) %>%
        slice_max(order_by=abs(logFC), n=nMarkers) %>%
        as.data.frame

    rownames(markers) <- markers$marker

    markers

    list(markers=markers, de.results=de.res)
}


clin <- read_excel('Primary clinical variables doppler.xlsx', 1) %>% as.data.frame
colnames(clin) <- make.names(colnames(clin))
clin <- clin[-1,]
rownames(clin) <- make.names(clin$Study_ID.ID)
head(clin)

clin <- clin %>% mutate_if(function(x) mean(is.na(as.numeric(x))) < 0.5, as.numeric)

## clin <- read.csv('clinical-211205.csv', na.strings='#N/A', row.names=1)

## rownames(clin) <- str_replace_all(rownames(clin), '-', '.')

## head(clin)

## clin$Batch <- ifelse(str_detect(rownames(clin), 'KH'), 'KH', NA)
## clin$Batch <- ifelse(str_detect(rownames(clin), 'MJ'), 'MJ', clin$Batch)
## clin$Batch <- ifelse(str_detect(rownames(clin), 'Mini.DP'), 'Mini-DP', clin$Batch)

clin$Condition. <- ifelse(clin$Condition. %in% c('PTD','Severe PE','FGR','Control','Control PTD'), clin$Condition., 'FGR+HDP')

clin$PrePregBMI <- clin$PrePregWt_Kg / clin$Height_Meters.HeightMeters^2

clin$Smoking <- ifelse(str_detect(clin$Smoke, 'Yes'), 'Yes', NA)
clin$Smoking <- ifelse(str_detect(clin$Smoke, 'former'), 'Former', clin$Smoking)
clin$Smoking <- ifelse(str_detect(clin$Smoke, 'Never'), 'No', clin$Smoking)

clin$Race[!clin$Race %in% c('W', 'B')] <- 'other'
clin$FDELTYPE[clin$FDELTYPE==0] <- NA

clin$MENDDIAB <- as.numeric(clin$MENDDIAB>0)
clin$MDEPRESS <- as.numeric(clin$MDEPRESS>0)


pathology <- read_excel('Path diagnoses for analysis.xlsx' , 1) %>% as.data.frame
colnames(pathology)
rownames(pathology) <- make.names(pathology$case)
pathology <- pathology[,11:15]
colnames(pathology) <- c("High.grade.MVM", "FVM",
                         "Acute.inflammation", "Chronic.inflammation", "VUE")

pathology <- pathology %>% mutate_at(c("High.grade.MVM", "FVM",
                                       "Acute.inflammation", "Chronic.inflammation",
                                       "VUE"),
                                     as.factor)
head(pathology)
dim(pathology)

clin[rownames(pathology),
     c("High.grade.MVM", "FVM", "Acute.inflammation", "Chronic.inflammation", "VUE")] <- pathology


## Multi-omics clustering

prot <- read.csv('data/combat_prot_data.csv', row.names=1)
colnames(prot) <- paste0('Prot_', colnames(prot))
head(prot)

metab <- read.csv('data/combat_metab_data.csv', row.names=1)
colnames(metab) <- paste0('Metab_', colnames(metab))
head(metab)

rna <- read.csv('data/combat_vst_rna_expression.csv', row.names=1)
colnames(rna) <- paste0('RNA_', colnames(rna))
head(rna)

mirna <- read.csv('data/combat_vst_mirna_expression.csv', row.names=1)
colnames(mirna) <- paste0('miRNA_', colnames(mirna))
head(mirna)


#### Highly variant features, top quartile

hvar.feats <- list()

prot.vars <- apply(prot, 2, stats::var)
hvar.feats[['prot']] <- colnames(prot)[which(prot.vars >= quantile(prot.vars, 0.75))]

metab.vars <- apply(metab, 2, stats::var)
hvar.feats[['metab']] <- colnames(metab)[which(metab.vars >= quantile(metab.vars, 0.75))]

rna.vars <- apply(rna, 2, stats::var)
hvar.feats[['rna']] <- colnames(rna)[which(rna.vars >= quantile(rna.vars, 0.75))]

mirna.vars <- apply(mirna, 2, stats::var)
hvar.feats[['mirna']] <- colnames(mirna)[which(mirna.vars >= quantile(mirna.vars, 0.75))]



rna.outliers <- c('MJ.0788', 'KH.1103', 'MJ.1151', 'MJ.0286', 'MJ.0687', 'Mini.DP.025',
                  'MJ.0909', 'Mini.DP.104', 'MJ.1065', 'KH.1100', 'Mini.DP.181',
                  'Mini.DP.013', 'MJ.0899', 'Mini.DP.045')

samp.ids <- intersect(rownames(prot), rownames(metab))

samp.ids <- intersect(samp.ids, rownames(rna))

samp.ids <- intersect(samp.ids, rownames(mirna))

samp.ids <- setdiff(samp.ids, rna.outliers)

data.views <- list(metab=metab[samp.ids,hvar.feats[['metab']]],
                   prot=prot[samp.ids, hvar.feats[['prot']]],
                   rna=rna[samp.ids, hvar.feats[['rna']]],
                   mirna=mirna[samp.ids, hvar.feats[['mirna']]])

clin <- clin[samp.ids,]

data.views <- lapply(data.views, standardNormalization)


Cs <- 4
Ks <- c(5:10, 15, 20, 25, 30)
Sigmas <- seq(0.3, 0.8, 0.05)

snf.res.full <- snf.cluster(data.views, Cs, Ks, Sigmas)

instabs <- snf.instabs(data.views, Cs, Ks, Sigmas, snf.res.full, 50, 0.8, 20220407)

instab.df <- data.frame()

for (k in Ks) {
    for (sigma in Sigmas) {       
        for (c in Cs) {
            name <- paste0('K', k, '_Sigma', sigma, '_C', c)
            print(name)

            instab.df <- rbind(instab.df,
                               data.frame(k, sigma, c,
                                          instabs[[paste0(name, '_Mean')]],
                                          instabs[[paste0(name, '_StDev')]]))
        }
    }
}

colnames(instab.df) <- c('K', 'Sigma', 'C', 'AMI', 'AMI_sd')


gg <- ggplot(instab.df %>% filter(C>1), aes(x=C, y=AMI)) +
    geom_line() +
    geom_point() +
    facet_grid(rows=vars(Sigma), cols=vars(K)) +
    theme_bw()
gg

ggsave('instability_grid_search_hvar2.pdf', gg, width=16, height=12, dpi=400)
ggsave('instability_grid_search_hvar2.png', gg, width=16, height=12, dpi=400)

instab.df %>% filter(C>1) %>% group_by(C) %>% arrange(desc(AMI), .by_group=TRUE) %>% as.data.frame

instab.df$AMI_SE <- instab.df$AMI_sd/sqrt(50)

instab.df %>% filter(C>1) %>% group_by(C) %>% arrange(desc(AMI)) %>% as.data.frame


data.views.dist <- lapply(data.views, dist, diag=TRUE, upper=TRUE, method='euclidean') ## euclidean

#### hvar top quartile
##   Best:                  C=6, Sigma=0.35, K=30
##   1se:                   C=6, Sigma=0.30, K=30
##   Best with different C: C=4, Sigma=0.30, K=20
##   Best Eigen Gap: 

lapply(data.views.dist, head)

## Best Stability
k <- 30
sigma <- 0.35
C <- 6

## Best Eigen Gap
k <- 20
sigma <- 0.30
C <- 4

## hist(data.views.dist[[4]])

data.views.aff <- lapply(data.views.dist, function(x) affinityMatrix(as.matrix(x), K=k, sigma=sigma))

W <- SNF(data.views.aff, K=k, t=50)

hist(log(2*W[W!=0.5]))

best.C <- estimateNumberOfClustersGivenGraph(W, 3:20)

best.C

labels <- spectralClustering(W, C)

temp <- labels

labels[temp==2] <- 3
labels[temp==3] <- 2
# labels[temp==2] <- 3

ind <- sort(labels, index.return = TRUE)
ind <- ind$ix

snf.heatmap <- function(W) {
    normalize <- function(X) X/rowSums(X)
    ind <- sort(as.vector(labels), index.return = TRUE)
    ind <- ind$ix
    heatmap.W <- W
    diag(heatmap.W) <- median(as.vector(W))
    heatmap.W <- normalize(heatmap.W)
    heatmap.W <- heatmap.W + t(heatmap.W)
    heatmap.W
}

heatmap(snf.heatmap(W)[ind, ind], scale = "none", Rowv = NA, Colv = NA)

pheatmap(log(snf.heatmap(W)[ind, ind]), cluster_rows=F, cluster_cols=F)
pheatmap(snf.heatmap(W)[ind, ind], breaks=c(0, 0.05), cluster_rows=F, cluster_cols=F)

print(calNMI(labels, sapply(clin$Condition., function(x) which(unique(clin$Condition.)==x))))

print(concordanceNetworkNMI(data.views.aff, C))

distinct.colors <- read.csv('../distinct_colors/colors.10.txt', header=F, sep='')
colnames(distinct.colors)[2:8] <- c('index', 'R', 'G', 'B', 'L', 'C', 'H')
distinct.colors$hex <- apply(distinct.colors, 1,
                             function(x) rgb(x['R'], x['G'], x['B'], maxColorValue=255))
distinct.colors <- distinct.colors[sample(10),]

clin.clusters <- read.csv('clustered_clinical_hvar.csv', row.names=1)

clin[,'Clusters'] <- clin.clusters[samp.ids,'hvar.Clusters']

print(ggplot(clin, aes(x=Clusters, fill=Condition.)) +
      geom_bar(position='fill') +
      ylab('Percentage') +
      scale_fill_manual(values=distinct.colors$hex[1:6]) +
      ## scale_fill_brewer(palette = 'Dark2') + 
      theme_bw())

print(ggplot(clin, aes(x=Clusters, fill=Condition.)) +
      geom_bar(position='stack') +
      ylab('Count') +
      scale_fill_manual(values=distinct.colors$hex[1:6]) +
      ## scale_fill_brewer(palette = 'Dark2') + 
      theme_bw())

clin$Prot_FLT1overPGF <- prot[samp.ids,'Prot_FLT1'] / prot[samp.ids,'Prot_PGF']

clin$RNA_FLT1overPGF <- rna[samp.ids,'RNA_FLT1'] / rna[samp.ids,'RNA_PGF']

ggplot(clin %>% mutate(Clusters=factor(Clusters)),
                    aes(x=Clusters, y=Prot_FLT1overPGF, fill=Clusters)) +
    geom_violin() +
    geom_boxplot(width=0.2) +
    xlab('Clusters') +
    ylab('Protein FLT1/PGF Ratio') +
    stat_compare_means() + 
    ## scale_fill_brewer(palette = 'Dark2') + 
    theme_bw()

ggplot(clin %>% mutate(Clusters=factor(Clusters)),
                    aes(x=Clusters, y=RNA_FLT1overPGF, fill=Clusters)) +
    geom_violin() +
    geom_boxplot(width=0.2) +
    xlab('Clusters') +
    ylab('RNA FLT1/PGF Ratio') +
    stat_compare_means() + 
    ## scale_fill_brewer(palette = 'Dark2') + 
    theme_bw()


table(clin$Clusters, clin$Condition.)
chisq.test(table(clin$Clusters, clin$Condition.))

## library(circlize)

## chordDiagram(W)

clin$Clusters <- factor(clin$Clusters)

pstack3 <- ggplot(clin, aes(x=Clusters, fill=Condition.)) +
    geom_bar(position='stack') +
    ylab('Count') +
    scale_fill_manual(values=distinct.colors$hex[1:6]) +
    ## scale_fill_brewer(palette = 'Dark2') + 
    theme_bw()

pfill3 <- ggplot(clin, aes(x=Clusters, fill=Condition.)) +
    geom_bar(position='fill') +
    ylab('Percentage') +
    scale_fill_manual(values=distinct.colors$hex[1:6]) +
    ## scale_fill_brewer(palette = 'Dark2') + 
    theme_bw()

samp.order <- samp.ids[ind]

annot.df <- data.frame(Clusters=factor(clin[samp.order,'Clusters']),
                       row.names=samp.order)

annot.colors <- list()
annot.colors[['Clusters']] <- distinct.colors$hex[7:10]
names(annot.colors[['Clusters']]) <- levels(clin$Clusters)

clustmap3 <- pheatmap(snf.heatmap(data.views.aff[['prot']])[ind, ind], cluster_rows=F, cluster_cols=F,
                      fontsize_row = 2, fontsize_col = 2, cellwidth=2, cellheight=2,
                      show_rownames=F, show_colnames=F, annotation_col=annot.df)


clustmap3 <- pheatmap(snf.heatmap(W)[ind, ind], cluster_rows=F, cluster_cols=F,
                      fontsize_row = 2, fontsize_col = 2, cellwidth=2, cellheight=2,
                      show_rownames=F, show_colnames=F, annotation_col=annot.df,
                      annotation_colors = annot.colors)


plot_grid(clustmap3$gtable, pstack3, pfill3, nrow=1, rel_widths=c(1,0.5,0.5))

ggsave('snf_hvar_clusts/cluster_comp_colorsV1.pdf', plot_grid(clustmap3$gtable, pstack3, pfill3, nrow=1, rel_widths=c(1,0.5,0.5)), width=18, height=8)

ggsave('snf_hvar_clusts/cluster_comp_colorsV1.png', plot_grid(clustmap3$gtable, pstack3, pfill3, nrow=1, rel_widths=c(1,0.5,0.5)), width=18, height=8)

### COLOR V2


clin$Diagnosis <- factor(clin$Condition., levels=c('Control', 'Control PTD', 'Severe PE',
                                                   'FGR', 'FGR+HDP', 'PTD'))


brewer.colors <- brewer.pal(n=12, "Paired")## )(10)

brewer.colors <- c(brewer.colors[1:2], brewer.colors[8], brewer.colors[5:6], brewer.colors[11], brewer.colors[3:4], brewer.colors[9:10])


clin$Clusters <- factor(labels)

pstack3 <- ggplot(clin, aes(x=Clusters, fill=Diagnosis)) +
    geom_bar(position='stack', color='black') +
    ylab('Count') +
    scale_fill_manual(values=brewer.colors[1:6]) +
    ## scale_fill_brewer(palette = 'Dark2') + 
    theme_bw()

pfill3 <- ggplot(clin, aes(x=Clusters, fill=Diagnosis)) +
    geom_bar(position='fill', color='black') +
    ylab('Percentage') +
    scale_fill_manual(values=brewer.colors[1:6]) +
    ## scale_fill_brewer(palette = 'Dark2') + 
    theme_bw()

samp.order <- samp.ids[ind]

annot.df <- data.frame(Clusters=factor(clin[samp.order,'Clusters']),
                       row.names=samp.order)

annot.colors <- list()
annot.colors[['Clusters']] <- brewer.colors[7:10]
names(annot.colors[['Clusters']]) <- levels(clin$Clusters)


clustmap3 <- pheatmap::pheatmap(snf.heatmap(W)[ind, ind],
                                breaks=seq(exp(-4.4), exp(-2.25), length.out=101), 
                                cluster_rows=F, cluster_cols=F,
                                fontsize_row = 2, fontsize_col = 2,
                                cellwidth=1.95, cellheight=1.95,
                                show_rownames=F, show_colnames=F,
                                annotation_col=annot.df,
                                annotation_colors = annot.colors)


plot_grid(clustmap3$gtable, pstack3, pfill3, nrow=1, rel_widths=c(1,0.5,0.5))

ggsave('snf_hvar_clusts/cluster_comp_finalColors_contrast.pdf', plot_grid(clustmap3$gtable, pstack3, pfill3, nrow=1, rel_widths=c(1,0.5,0.5)), width=18, height=8)

ggsave('snf_hvar_clusts/cluster_comp_finalColors_contrast.png', plot_grid(clustmap3$gtable, pstack3, pfill3, nrow=1, rel_widths=c(1,0.5,0.5)), width=18, height=8)


## Color V3

## iwanthue.colors <- c("#6ed154",
## "#8444c6",
## "#c5c351",
## "#d764af",
## "#62cbaa",
## "#c34944",
## "#6e83cb",
## "#bf8140",
## "#723164",
## "#587d38")

iwanthue.colors <- c("#63d048","#854dc9","#ccba39","#ce549e","#55c3a1",
                     "#b7534c","#6788ce","#c4692d","#8d4077","#799b57")


clin$Clusters <- factor(clin$Clusters)

pstack3 <- ggplot(clin, aes(x=Clusters, fill=Condition.)) +
    geom_bar(position='stack') +
    ylab('Count') +
    scale_fill_manual(values=iwanthue.colors[1:6]) +
    ## scale_fill_brewer(palette = 'Dark2') + 
    theme_bw()

pfill3 <- ggplot(clin, aes(x=Clusters, fill=Condition.)) +
    geom_bar(position='fill') +
    ylab('Percentage') +
    scale_fill_manual(values=iwanthue.colors[1:6]) +
    ## scale_fill_brewer(palette = 'Dark2') + 
    theme_bw()

samp.order <- samp.ids[ind]

annot.df <- data.frame(Clusters=factor(clin[samp.order,'Clusters']),
                       row.names=samp.order)

annot.colors <- list()
annot.colors[['Clusters']] <- iwanthue.colors[7:10]
names(annot.colors[['Clusters']]) <- levels(clin$Clusters)

## clustmap3 <- pheatmap(snf.heatmap(data.views.aff[['prot']])[ind, ind], cluster_rows=F, cluster_cols=F,
##                       fontsize_row = 2, fontsize_col = 2, cellwidth=2, cellheight=2,
##                       show_rownames=F, show_colnames=F, annotation_col=annot.df)


clustmap3 <- pheatmap(snf.heatmap(W)[ind, ind], cluster_rows=F, cluster_cols=F,
                      fontsize_row = 2, fontsize_col = 2, cellwidth=2, cellheight=2,
                      show_rownames=F, show_colnames=F, annotation_col=annot.df,
                      annotation_colors = annot.colors)


plot_grid(clustmap3$gtable, pstack3, pfill3, nrow=1, rel_widths=c(1,0.5,0.5))

ggsave('snf_hvar_clusts/cluster_comp_colorsV3.pdf', plot_grid(clustmap3$gtable, pstack3, pfill3, nrow=1, rel_widths=c(1,0.5,0.5)), width=18, height=8)

ggsave('snf_hvar_clusts/cluster_comp_colorsV3.png', plot_grid(clustmap3$gtable, pstack3, pfill3, nrow=1, rel_widths=c(1,0.5,0.5)), width=18, height=8)



prot.markers <- clustMarkersWilcox(prot[samp.ids,], clin$Clusters)

write.csv(prot.markers, 'snf_hvar_clusts/markers/wilcox_auc_prot_markers.csv')

## samp.order <- c()
## for (clust in 1:3) {
##     inclust.tree <- hclust(dist(scale(prot[samp.ids,])[clin$Clusters==clust, prot.markers$Feature[((clust-1)*10+1):(clust*10)]]), method='ward.D')
##     samp.order <- c(samp.order, inclust.tree$labels[inclust.tree$order])
## }

samp.order <- samp.ids[ind]

annot.df <- data.frame(Clusters=factor(clin[samp.order,'Clusters']),
                       row.names=samp.order)

col.gaps <- c(rev(which(annot.df$Clusters==1))[1],
              rev(which(annot.df$Clusters==2))[1],
              rev(which(annot.df$Clusters==3))[1])

pdf('snf_hvar_clusts/markers/auc_prot_markers.pdf', width=6, height=7)
pheatmap(t(prot[samp.order,prot.markers$signif$Feature]), breaks=seq(-3,3,length.out=101),
         cluster_rows=F, cluster_cols=F, scale='row', annotation_col=annot.df,
         clustering_method='ward.D', labels_row=prot.markers$Feature, show_colnames=FALSE,
         gaps_row=c(10,20,30), gaps_col=col.gaps)
dev.off()

png('snf_hvar_clusts/markers/auc_prot_markers.png', width=6, height=7, units='in', res=400)
pheatmap(t(prot[samp.order,prot.markers$Feature]), breaks=seq(-3,3,length.out=101), cluster_rows=F, cluster_cols=F, scale='row', annotation_col=annot.df, clustering_method='ward.D', labels_row=prot.markers$Feature, show_colnames=FALSE, gaps_row=c(10,20,30), gaps_col=col.gaps)
dev.off()


metab.markers <- clustMarkersAUC(metab[samp.ids,], clin$Clusters)

write.csv(metab.markers, 'snf_hvar_clusts/markers/auc_metab_markers.csv')

## samp.order <- c()
## for (clust in 1:3) {
##     inclust.tree <- hclust(dist(scale(metab[samp.ids,])[clin$Clusters==clust, metab.markers$Feature[((clust-1)*10+1):(clust*10)]]), method='ward.D')
##     samp.order <- c(samp.order, inclust.tree$labels[inclust.tree$order])
## }

samp.order <- samp.ids[ind]

annot.df <- data.frame(Clusters=factor(clin[samp.order,'Clusters']),
                       row.names=samp.order)

pdf('snf_hvar_clusts/markers/auc_metab_markers.pdf', width=8, height=7)
pheatmap(t(metab[samp.order,metab.markers$Feature]), breaks=seq(-3,3,length.out=101), cluster_rows=F, cluster_cols=F, scale='row', annotation_col=annot.df, clustering_method='ward.D', labels_row=metab.markers$Feature, show_colnames=FALSE, gaps_row=c(10,20,30), gaps_col=col.gaps)
dev.off()

png('snf_hvar_clusts/markers/auc_metab_markers.png', width=8, height=7, units='in', res=400)
pheatmap(t(metab[samp.order,metab.markers$Feature]), breaks=seq(-3,3,length.out=101), cluster_rows=F, cluster_cols=F, scale='row', annotation_col=annot.df, clustering_method='ward.D', labels_row=metab.markers$Feature, show_colnames=FALSE, gaps_row=c(10,20,30), gaps_col=col.gaps)
dev.off()


rna.markers <- clustMarkersAUC(rna[samp.ids,], clin$Clusters)

write.csv(rna.markers, 'snf_hvar_clusts/markers/auc_rna_markers.csv')

pdf('snf_hvar_clusts/markers/auc_rna_markers.pdf', width=6, height=7)
pheatmap(t(rna[samp.ids[ind],rna.markers$Feature]), breaks=seq(-3,3,length.out=101), cluster_rows=F, cluster_cols=F, scale='row', annotation_col=annot.df, clustering_method='ward.D', labels_row=rna.markers$Feature, show_colnames=FALSE, gaps_row=c(10,20,30), gaps_col=col.gaps)
dev.off()

png('snf_hvar_clusts/markers/auc_rna_markers.png', width=6, height=7, units='in', res=400)
pheatmap(t(rna[samp.ids[ind],rna.markers$Feature]), breaks=seq(-3,3,length.out=101), cluster_rows=F, cluster_cols=F, scale='row', annotation_col=annot.df, clustering_method='ward.D', labels_row=rna.markers$Feature, show_colnames=FALSE, gaps_row=c(10,20,30), gaps_col=col.gaps)
dev.off()


mirna.markers <- clustMarkersAUC(mirna[samp.ids,], clin$Clusters)

write.csv(mirna.markers, 'snf_hvar_clusts/markers/auc_mirna_markers.csv')

pdf('snf_hvar_clusts/markers/auc_mirna_markers.pdf', width=6, height=7)
pheatmap(t(mirna[samp.ids[ind],mirna.markers$Feature]), breaks=seq(-3,3,length.out=101), cluster_rows=F, cluster_cols=F, scale='row', annotation_col=annot.df, clustering_method='ward.D', labels_row=mirna.markers$Feature, show_colnames=FALSE, gaps_row=c(10,20,30), gaps_col=col.gaps)
dev.off()

png('snf_hvar_clusts/markers/auc_mirna_markers.png', width=6, height=7, units='in', res=400)
pheatmap(t(mirna[samp.ids[ind],mirna.markers$Feature]), breaks=seq(-3,3,length.out=101), cluster_rows=F, cluster_cols=F, scale='row', annotation_col=annot.df, clustering_method='ward.D', labels_row=mirna.markers$Feature, show_colnames=FALSE, gaps_row=c(10,20,30), gaps_col=col.gaps)
dev.off()



clustered.clin <- read.csv('clustered_clinical.csv', row.names=1)

table(de=clustered.clin$Clusters, hvar=clin$Clusters)

clustered.clin$hvar.Clusters <- clin$Clusters

clustered.clin$Smoking.detailed <- clin$Smoking

write.csv(clustered.clin, 'clustered_clinical_hvar.csv')


###### FULL AND SIGNIFICANT AUC BASED CLUSTER MARKERS ######

ind <- sort(as.vector(clin.clusters$hvar.Clusters), index.return = TRUE)
ind <- ind$ix

samp.order <- samp.ids[ind]

annot.df <- data.frame(Clusters=factor(clin.clusters[samp.order,'hvar.Clusters']),
                       row.names=samp.order)

col.gaps <- c(rev(which(annot.df$Clusters==1))[1],
              rev(which(annot.df$Clusters==2))[1],
              rev(which(annot.df$Clusters==3))[1])

### Proteins

prot.markers <- clustMarkersWilcox(prot[samp.ids,], clin.clusters$hvar.Clusters, nMarkers=ncol(prot))

## prot.markers$signif <- prot.markers$signif[prot.markers$signif$adj.p.value < 0.05,]

table(prot.markers$all$Cluster)

table(prot.markers$signif$Cluster)

write.csv(prot.markers$all, 'snf_hvar_clusts/markers/wilcox_auc_prot_all_markers.csv')

write.csv(prot.markers$signif, 'snf_hvar_clusts/markers/wilcox_auc_prot_signif_markers.csv')

row.gaps <- c(rev(which(prot.markers$signif$Cluster==1))[1],
              rev(which(prot.markers$signif$Cluster==2))[1],
              rev(which(prot.markers$signif$Cluster==3))[1])

pdf('snf_hvar_clusts/markers/wilcox_auc_prot_signif_markers.pdf', width=15, height=3 + 0.13*nrow(prot.markers$signif))
pheatmap(t(prot[samp.order,prot.markers$signif$Feature]), breaks=seq(-3,3,length.out=101), cluster_rows=F, cluster_cols=F, scale='row', annotation_col=annot.df, clustering_method='ward.D', labels_row=prot.markers$Feature, show_colnames=FALSE, gaps_row=row.gaps, gaps_col=col.gaps, annotation_colors=annot.colors)
dev.off()

for ( c in 1:4 ) {
    write.csv(prot.markers$all %>% filter(Cluster==c),
              paste0('snf_hvar_clusts/markers/clust', c,
                     '/clust', c, '_wilcox_auc_prot_all_markers.csv'))
    write.csv(prot.markers$signif %>% filter(Cluster==c),
              paste0('snf_hvar_clusts/markers/clust', c,
                     '/clust', c, '_wilcox_auc_prot_signif_markers.csv'))
}


prot.markers <- clustMarkersWilcox(prot[samp.ids,], clin.clusters$hvar.Clusters, nMarkers=10)

table(prot.markers$all$Cluster)

cumsum(table(prot.markers$all$Cluster))

row.gaps <- cumsum(table(prot.markers$all$Cluster))[-4]

annot.colors <- list()
annot.colors[['Clusters']] <- brewer.colors[7:10]
names(annot.colors[['Clusters']]) <- levels(clin$Clusters)

prot.meta <- read.tcsv('olink-prot-211205.csv', na.strings=c('#N/A', ''), sep=',')
prot.meta <- prot.meta %>% select(c('Panel', 'Assay', 'Uniprot.ID', 'OlinkID', 'Missing.Data.freq.', 'LOD'))
rownames(prot.meta) <- make.names(prot.meta$Assay, unique=T)

row.labels <- prot.markers$all$Feature

idx <- 1
for (strlist in str_split(prot.markers$all$Feature,'_')) {
    if (strlist[1] == 'Prot') {
        row.labels[idx] <- prot.meta[strlist[2],'Assay']
    }
    idx <- idx+1
}

write.csv(data.frame(Labels=row.labels),
          'snf_hvar_clusts/markers/wilcox_auc_prot_markers_row_labels.csv',
          row.names=F)

pdf('snf_hvar_clusts/markers/wilcox_auc_prot_markers.pdf',
    width=6, height=7)

pheatmap::pheatmap(t(prot[samp.order,prot.markers$all$Feature]),
                   breaks=seq(-3,3,length.out=101), cluster_rows=F,
                   cluster_cols=F, scale='row', show_rownames=FALSE, 
                   annotation_col=annot.df, clustering_method='ward.D',
                   labels_row=row.labels, show_colnames=FALSE, gaps_row=row.gaps,
                   gaps_col=col.gaps, annotation_colors=annot.colors)

dev.off()

png('snf_hvar_clusts/markers/wilcox_auc_prot_markers.png',
    width=6, height=7, units='in', res=400)

pheatmap::pheatmap(t(prot[samp.order,prot.markers$all$Feature]),
         breaks=seq(-3,3,length.out=101), cluster_rows=F, cluster_cols=F, scale='row',
         annotation_col=annot.df, clustering_method='ward.D', show_rownames=FALSE,
         labels_row=row.labels, show_colnames=FALSE, gaps_row=row.gaps,
         gaps_col=col.gaps, annotation_colors=annot.colors)

dev.off()


## png('snf_hvar_clusts/markers/auc_prot_signif_markers.png', width=15, height=3 + 0.13*nrow(prot.markers$signif), units='in', res=400)
## pheatmap(t(prot[samp.order,prot.markers$signif$Feature]), breaks=seq(-3,3,length.out=101), cluster_rows=F, cluster_cols=F, scale='row', annotation_col=annot.df, clustering_method='ward.D', labels_row=prot.markers$Feature, show_colnames=FALSE, gaps_row=row.gaps, gaps_col=col.gaps)
## dev.off()

### Metabolites

metab.markers <- clustMarkersWilcox(metab[samp.ids,], clin.clusters$hvar.Clusters, nMarkers=ncol(metab))

table(metab.markers$all$Cluster)

table(metab.markers$signif$Cluster)

write.csv(metab.markers$all, 'snf_hvar_clusts/markers/wilcox_auc_metab_all_markers.csv')

write.csv(metab.markers$signif, 'snf_hvar_clusts/markers/wilcox_auc_metab_signif_markers.csv')

row.gaps <- c(rev(which(metab.markers$signif$Cluster==1))[1],
              rev(which(metab.markers$signif$Cluster==2))[1],
              rev(which(metab.markers$signif$Cluster==3))[1])

pdf('snf_hvar_clusts/markers/wilcox_auc_metab_signif_markers.pdf', width=18, height=3 + 0.13*nrow(metab.markers$signif))
pheatmap(t(metab[samp.order,metab.markers$signif$Feature]), breaks=seq(-3,3,length.out=101), cluster_rows=F, cluster_cols=F, scale='row', annotation_col=annot.df, clustering_method='ward.D', labels_row=metab.markers$Feature, show_colnames=FALSE, gaps_row=row.gaps, gaps_col=col.gaps)
dev.off()


for ( c in 1:4 ) {
    write.csv(metab.markers$all %>% filter(Cluster==c),
              paste0('snf_hvar_clusts/markers/clust', c,
                     '/clust', c, '_wilcox_auc_metab_all_markers.csv'))
    write.csv(metab.markers$signif %>% filter(Cluster==c),
              paste0('snf_hvar_clusts/markers/clust', c,
                     '/clust', c, '_wilcox_auc_metab_signif_markers.csv'))
}



metab.markers <- clustMarkersWilcox(metab[samp.ids,], clin.clusters$hvar.Clusters, nMarkers=10)

table(metab.markers$all$Cluster)

cumsum(table(metab.markers$all$Cluster))

row.gaps <- cumsum(table(metab.markers$all$Cluster))[-4]

annot.colors <- list()
annot.colors[['Clusters']] <- brewer.colors[7:10]
names(annot.colors[['Clusters']]) <- levels(clin$Clusters)

metab.meta <- read_excel('UPIT-02-20PHML+ DATA TABLES OB.xlsx', 2) %>% as.data.frame
rownames(metab.meta) <- make.names(metab.meta$CHEMICAL_NAME)

row.labels <- metab.markers$all$Feature

idx <- 1
for (strlist in str_split(metab.markers$all$Feature,'_')) {
    if (strlist[1] == 'Metab') {
        row.labels[idx] <- metab.meta[strlist[2],'CHEMICAL_NAME']
    }
    idx <- idx+1
}

write.csv(data.frame(Labels=row.labels),
          'snf_hvar_clusts/markers/wilcox_auc_metab_markers_row_labels.csv',
          row.names=F)

pdf('snf_hvar_clusts/markers/wilcox_auc_metab_markers.pdf',
    width=6, height=7)

pheatmap::pheatmap(t(metab[samp.order,metab.markers$all$Feature]),
         breaks=seq(-3,3,length.out=101), cluster_rows=F, cluster_cols=F, scale='row',
         annotation_col=annot.df, clustering_method='ward.D', show_rownames=FALSE,
         labels_row=row.labels, show_colnames=FALSE, gaps_row=row.gaps,
         gaps_col=col.gaps, annotation_colors=annot.colors)

dev.off()

png('snf_hvar_clusts/markers/wilcox_auc_metab_markers.png',
    width=6, height=7, units='in', res=400)

pheatmap::pheatmap(t(metab[samp.order,metab.markers$all$Feature]),
         breaks=seq(-3,3,length.out=101), cluster_rows=F, cluster_cols=F, scale='row',
         annotation_col=annot.df, clustering_method='ward.D', show_rownames=FALSE,
         labels_row=row.labels, show_colnames=FALSE, gaps_row=row.gaps,
         gaps_col=col.gaps, annotation_colors=annot.colors)

dev.off()



## png('snf_hvar_clusts/markers/auc_metab_signif_markers.png', width=15, height=3 + 0.13*nrow(metab.markers$signif), units='in', res=400)
## pheatmap(t(metab[samp.order,metab.markers$signif$Feature]), breaks=seq(-3,3,length.out=101), cluster_rows=F, cluster_cols=F, scale='row', annotation_col=annot.df, clustering_method='ward.D', labels_row=metab.markers$Feature, show_colnames=FALSE, gaps_row=row.gaps, gaps_col=col.gaps)
## dev.off()

### RNA

rna.markers <- clustMarkersWilcox(rna[samp.ids,], clin.clusters$hvar.Clusters, nMarkers=ncol(rna))

table(rna.markers$all$Cluster)

table(rna.markers$signif$Cluster)

write.csv(rna.markers$all, 'snf_hvar_clusts/markers/wilcox_auc_rna_all_markers.csv')

write.csv(rna.markers$signif, 'snf_hvar_clusts/markers/wilcox_auc_rna_signif_markers.csv')

row.gaps <- c(rev(which(rna.markers$signif$Cluster==1))[1],
              rev(which(rna.markers$signif$Cluster==2))[1],
              rev(which(rna.markers$signif$Cluster==3))[1])

pdf('snf_hvar_clusts/markers/wilcox_auc_rna_signif_markers.pdf', width=15, height=3 + 0.13*nrow(rna.markers$signif))
pheatmap(t(rna[samp.order,rna.markers$signif$Feature]), breaks=seq(-3,3,length.out=101), cluster_rows=F, cluster_cols=F, scale='row', annotation_col=annot.df, clustering_method='ward.D', labels_row=rna.markers$Feature, show_colnames=FALSE, gaps_row=row.gaps, gaps_col=col.gaps)
dev.off()


for ( c in 1:4 ) {
    write.csv(rna.markers$all %>% filter(Cluster==c),
              paste0('snf_hvar_clusts/markers/clust', c,
                     '/clust', c, '_wilcox_auc_rna_all_markers.csv'))
    write.csv(rna.markers$signif %>% filter(Cluster==c),
              paste0('snf_hvar_clusts/markers/clust', c,
                     '/clust', c, '_wilcox_auc_rna_signif_markers.csv'))
}


rna.markers <- clustMarkersWilcox(rna[samp.ids,], clin.clusters$hvar.Clusters, nMarkers=10)

table(rna.markers$all$Cluster)

cumsum(table(rna.markers$all$Cluster))

row.gaps <- cumsum(table(rna.markers$all$Cluster))[-4]

annot.colors <- list()
annot.colors[['Clusters']] <- brewer.colors[7:10]
names(annot.colors[['Clusters']]) <- levels(clin$Clusters)

row.labels <- rna.markers$all$Feature

idx <- 1
for (strlist in str_split(rna.markers$all$Feature,'_')) {
    if (strlist[1] == 'RNA') {
        row.labels[idx] <- strlist[2]
    }
    idx <- idx+1
}

write.csv(data.frame(Labels=row.labels),
          'snf_hvar_clusts/markers/wilcox_auc_rna_markers_row_labels.csv',
          row.names=F)

pdf('snf_hvar_clusts/markers/wilcox_auc_rna_markers.pdf',
    width=6, height=7)

pheatmap::pheatmap(t(rna[samp.order,rna.markers$all$Feature]),
         breaks=seq(-3,3,length.out=101), cluster_rows=F, cluster_cols=F, scale='row',
         annotation_col=annot.df, clustering_method='ward.D', show_rownames=FALSE,
         labels_row=row.labels, show_colnames=FALSE, gaps_row=row.gaps,
         gaps_col=col.gaps, annotation_colors=annot.colors)

dev.off()

png('snf_hvar_clusts/markers/wilcox_auc_rna_markers.png',
    width=6, height=7, units='in', res=400)

pheatmap::pheatmap(t(rna[samp.order,rna.markers$all$Feature]),
         breaks=seq(-3,3,length.out=101), cluster_rows=F, cluster_cols=F, scale='row',
         annotation_col=annot.df, clustering_method='ward.D', show_rownames=FALSE, 
         labels_row=row.labels, show_colnames=FALSE, gaps_row=row.gaps,
         gaps_col=col.gaps, annotation_colors=annot.colors)

dev.off()


## png('snf_hvar_clusts/markers/wilcox_auc_rna_signif_markers.png', width=15, height=3 + 0.13*nrow(rna.markers$signif), units='in', res=400)
## pheatmap(t(rna[samp.order,rna.markers$signif$Feature]), breaks=seq(-3,3,length.out=101), cluster_rows=F, cluster_cols=F, scale='row', annotation_col=annot.df, clustering_method='ward.D', labels_row=rna.markers$Feature, show_colnames=FALSE, gaps_row=row.gaps, gaps_col=col.gaps)
## dev.off()

### miRNA

mirna.markers <- clustMarkersWilcox(mirna[samp.ids,], clin.clusters$hvar.Clusters, nMarkers=ncol(mirna))

table(mirna.markers$all$Cluster)

table(mirna.markers$signif$Cluster)

write.csv(mirna.markers$all, 'snf_hvar_clusts/markers/wilcox_auc_mirna_all_markers.csv')

write.csv(mirna.markers$signif, 'snf_hvar_clusts/markers/wilcox_auc_mirna_signif_markers.csv')

row.gaps <- c(rev(which(mirna.markers$signif$Cluster==1))[1],
              rev(which(mirna.markers$signif$Cluster==2))[1],
              rev(which(mirna.markers$signif$Cluster==3))[1])

pdf('snf_hvar_clusts/markers/wilcox_auc_mirna_signif_markers.pdf', width=15, height=3 + 0.13*nrow(mirna.markers$signif))
pheatmap(t(mirna[samp.order,mirna.markers$signif$Feature]), breaks=seq(-3,3,length.out=101), cluster_rows=F, cluster_cols=F, scale='row', annotation_col=annot.df, clustering_method='ward.D', labels_row=mirna.markers$Feature, show_colnames=FALSE, gaps_row=row.gaps, gaps_col=col.gaps)
dev.off()

for ( c in 1:4 ) {
    write.csv(mirna.markers$all %>% filter(Cluster==c),
              paste0('snf_hvar_clusts/markers/clust', c,
                     '/clust', c, '_wilcox_auc_mirna_all_markers.csv'))
    write.csv(mirna.markers$signif %>% filter(Cluster==c),
              paste0('snf_hvar_clusts/markers/clust', c,
                     '/clust', c, '_wilcox_auc_mirna_signif_markers.csv'))
}


mirna.markers <- clustMarkersWilcox(mirna[samp.ids,], clin.clusters$hvar.Clusters, nMarkers=10)

table(mirna.markers$all$Cluster)

cumsum(table(mirna.markers$all$Cluster))

row.gaps <- cumsum(table(mirna.markers$all$Cluster))[-4]

annot.colors <- list()
annot.colors[['Clusters']] <- brewer.colors[7:10]
names(annot.colors[['Clusters']]) <- levels(clin$Clusters)

row.labels <- mirna.markers$all$Feature

idx <- 1
for (strlist in str_split(mirna.markers$all$Feature,'_')) {
    if (strlist[1] == 'miRNA') {
        row.labels[idx] <- gsub('[.]', '-', gsub('hsa[.]', '', strlist[2]))
    }
    idx <- idx+1
}

write.csv(data.frame(Labels=row.labels),
          'snf_hvar_clusts/markers/wilcox_auc_mirna_markers_row_labels.csv',
          row.names=F)

pdf('snf_hvar_clusts/markers/wilcox_auc_mirna_markers.pdf',
    width=6, height=7)

pheatmap::pheatmap(t(mirna[samp.order,mirna.markers$all$Feature]),
         breaks=seq(-3,3,length.out=101), cluster_rows=F, cluster_cols=F, scale='row',
         annotation_col=annot.df, clustering_method='ward.D', show_rownames=FALSE,
         labels_row=row.labels, show_colnames=FALSE, gaps_row=row.gaps,
         gaps_col=col.gaps, annotation_colors=annot.colors)

dev.off()

png('snf_hvar_clusts/markers/wilcox_auc_mirna_markers.png',
    width=6, height=7, units='in', res=400)

pheatmap::pheatmap(t(mirna[samp.order,mirna.markers$all$Feature]),
         breaks=seq(-3,3,length.out=101), cluster_rows=F, cluster_cols=F, scale='row',
         annotation_col=annot.df, clustering_method='ward.D', show_rownames=FALSE,
         labels_row=row.labels, show_colnames=FALSE, gaps_row=row.gaps,
         gaps_col=col.gaps, annotation_colors=annot.colors)

dev.off()



## png('snf_hvar_clusts/markers/auc_mirna_signif_markers.png', width=15, height=3 + 0.13*nrow(mirna.markers$signif), units='in', res=400)
## pheatmap(t(mirna[samp.order,mirna.markers$signif$Feature]), breaks=seq(-3,3,length.out=101), cluster_rows=F, cluster_cols=F, scale='row', annotation_col=annot.df, clustering_method='ward.D', labels_row=mirna.markers$Feature, show_colnames=FALSE, gaps_row=row.gaps, gaps_col=col.gaps)
## dev.off()

############################################################

ind <- sort(as.vector(clin.clusters$Clusters), index.return = TRUE)
ind <- ind$ix

prot.markers <- clustMarkersAUC(prot[samp.ids,], clin.clusters$Clusters)

write.csv(prot.markers, 'snfClust_DE/markers/auc_prot_markers.csv')

samp.order <- samp.ids[ind]

annot.df <- data.frame(Clusters=factor(clin.clusters[samp.order,'Clusters']),
                       row.names=samp.order)

col.gaps <- c(rev(which(annot.df$Clusters==1))[1],
              rev(which(annot.df$Clusters==2))[1])

pdf('snfClust_DE/markers/auc_prot_markers.pdf', width=6, height=7)
pheatmap(t(prot[samp.order,prot.markers$Feature]), breaks=seq(-3,3,length.out=101), cluster_rows=F, cluster_cols=F, scale='row', annotation_col=annot.df, clustering_method='ward.D', labels_row=prot.markers$Feature, show_colnames=FALSE, gaps_row=c(10,20), gaps_col=col.gaps)
dev.off()

png('snfClust_DE/markers/auc_prot_markers.png', width=6, height=7, units='in', res=400)
pheatmap(t(prot[samp.order,prot.markers$Feature]), breaks=seq(-3,3,length.out=101), cluster_rows=F, cluster_cols=F, scale='row', annotation_col=annot.df, clustering_method='ward.D', labels_row=prot.markers$Feature, show_colnames=FALSE, gaps_row=c(10,20), gaps_col=col.gaps)
dev.off()


metab.markers <- clustMarkersAUC(metab[samp.ids,], clin.clusters$Clusters)

write.csv(metab.markers, 'snfClust_DE/markers/auc_metab_markers.csv')

## samp.order <- c()
## for (clust in 1:3) {
##     inclust.tree <- hclust(dist(scale(metab[samp.ids,])[clin$Clusters==clust, metab.markers$Feature[((clust-1)*10+1):(clust*10)]]), method='ward.D')
##     samp.order <- c(samp.order, inclust.tree$labels[inclust.tree$order])
## }

samp.order <- samp.ids[ind]

annot.df <- data.frame(Clusters=factor(clin.clusters[samp.order,'Clusters']),
                       row.names=samp.order)

pdf('snfClust_DE/markers/auc_metab_markers.pdf', width=8, height=7)
pheatmap(t(metab[samp.order,metab.markers$Feature]), breaks=seq(-3,3,length.out=101), cluster_rows=F, cluster_cols=F, scale='row', annotation_col=annot.df, clustering_method='ward.D', labels_row=metab.markers$Feature, show_colnames=FALSE, gaps_row=c(10,20), gaps_col=col.gaps)
dev.off()

png('snfClust_DE/markers/auc_metab_markers.png', width=8, height=7, units='in', res=400)
pheatmap(t(metab[samp.order,metab.markers$Feature]), breaks=seq(-3,3,length.out=101), cluster_rows=F, cluster_cols=F, scale='row', annotation_col=annot.df, clustering_method='ward.D', labels_row=metab.markers$Feature, show_colnames=FALSE, gaps_row=c(10,20), gaps_col=col.gaps)
dev.off()


rna.markers <- clustMarkersAUC(rna[samp.ids,], clin.clusters$Clusters)

write.csv(rna.markers, 'snfClust_DE/markers/auc_rna_markers.csv')

pdf('snfClust_DE/markers/auc_rna_markers.pdf', width=6, height=7)
pheatmap(t(rna[samp.ids[ind],rna.markers$Feature]), breaks=seq(-3,3,length.out=101), cluster_rows=F, cluster_cols=F, scale='row', annotation_col=annot.df, clustering_method='ward.D', labels_row=rna.markers$Feature, show_colnames=FALSE, gaps_row=c(10,20), gaps_col=col.gaps)
dev.off()

png('snfClust_DE/markers/auc_rna_markers.png', width=6, height=7, units='in', res=400)
pheatmap(t(rna[samp.ids[ind],rna.markers$Feature]), breaks=seq(-3,3,length.out=101), cluster_rows=F, cluster_cols=F, scale='row', annotation_col=annot.df, clustering_method='ward.D', labels_row=rna.markers$Feature, show_colnames=FALSE, gaps_row=c(10,20), gaps_col=col.gaps)
dev.off()


mirna.markers <- clustMarkersAUC(mirna[samp.ids,], clin.clusters$Clusters)

write.csv(mirna.markers, 'snfClust_DE/markers/auc_mirna_markers.csv')

pdf('snfClust_DE/markers/auc_mirna_markers.pdf', width=6, height=7)
pheatmap(t(mirna[samp.ids[ind],mirna.markers$Feature]), breaks=seq(-3,3,length.out=101), cluster_rows=F, cluster_cols=F, scale='row', annotation_col=annot.df, clustering_method='ward.D', labels_row=mirna.markers$Feature, show_colnames=FALSE, gaps_row=c(10,20), gaps_col=col.gaps)
dev.off()

png('snfClust_DE/markers/auc_mirna_markers.png', width=6, height=7, units='in', res=400)
pheatmap(t(mirna[samp.ids[ind],mirna.markers$Feature]), breaks=seq(-3,3,length.out=101), cluster_rows=F, cluster_cols=F, scale='row', annotation_col=annot.df, clustering_method='ward.D', labels_row=mirna.markers$Feature, show_colnames=FALSE, gaps_row=c(10,20), gaps_col=col.gaps)
dev.off()
