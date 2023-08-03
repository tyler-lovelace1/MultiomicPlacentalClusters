library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(ReactomePA)
library(enrichplot)
library(ggplot2)
library(SetRank)

comparisons <- list(c('Control PTD', 'Control'),
                    c('FGR', 'Control'),
                    c('Severe PE', 'Control PTD'),
                    c('Severe PE', 'Control'),
                    c('PTD', 'Control PTD'),
                    c('PTD', 'Control'),
                    c('FGR+HDP', 'Control PTD'),
                    c('FGR+HDP', 'Control'))

groups <- comparisons[[7]]

de <- read.csv(paste('DE-cond3/diffExpRNA_', make.names(groups[1]), '_',
                     make.names(groups[2]), '.csv', sep=''), row.names=1)

ids <- bitr(de$Gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")



pathnames <- AnnotationDbi::keys(reactome.db::reactome.db, 'PATHID')
org.keep <- grepl('HSA', pathnames)
pathnames <- pathnames[org.keep]

reactome.df <- AnnotationDbi::select(reactome.db::reactome.db,
                                     columns=c('PATHID', 'PATHNAME', 'ENTREZID'),
                                     keys=pathnames,
                                     keytype='PATHID')

colnames(reactome.df) <- c('termID', 'termName', 'geneID')

reactome.df <- reactome.df[reactome.df$geneID %in% ids$ENTREZID,]

reactome.df$termName <- gsub('Homo sapiens\r: ', '', reactome.df$termName)
reactome.df$description <- reactome.df$termName

reactome.df$dbName <- 'REACTOME'

dim(reactome.df)

reactome.df <- reactome.df %>% group_by(termID) %>% filter(n() >= 10)

dim(reactome.df)

options(mc.cores=8)

setCollection <- buildSetCollection(reactome.df)



set.seed(20220817)

for (groups in comparisons) {

    de <- read.csv(paste('DE-cond3/diffExpRNA_', make.names(groups[1]), '_',
                         make.names(groups[2]), '.csv', sep=''), row.names=1)

    de.stats <- as.vector(matrix(NA, nrow(ids)))
    names(de.stats) <- ids$ENTREZID

    for ( i in 1:nrow(ids) ) {
        de.stats[i] <- de[ids$SYMBOL[i], 'stat']
    }

    de.stats <- sort(de.stats, decreasing=T)


    gsea.out <- gsePathway(de.stats, pvalueCutoff=1,
                           pAdjustMethod='none', scoreType='pos')

    gsea.out <- setReadable(gsea.out, OrgDb = org.Hs.eg.db, keyType="ENTREZID")

    gsea.out@result$Description <- gsub('Homo sapiens\r: ', '', gsea.out$Description)

    ## gsea.out@result <- gsea.out@result[gsea.out$NES > 0,]

    sr.out <- setRankAnalysis(names(de.stats), setCollection)

    exportSingleResult(sr.out, names(de.stats),
                       setCollection, paste0('clust',c),
                       outputPath='reactome/setrank')

    pathways <- read.table(paste0('reactome/setrank/clust', c, '_pathways.txt'),
                           header=T, row.names=1, sep='\t')
    
    pathways$p.adjust <- p.adjust(pathways$correctedPValue, 'fdr', n=1438)

    head(pathways)
    dim(pathways)

    write.csv(pathways, paste0('reactome/setrank/', make.names(groups[1]), '_',
                               make.names(groups[2]), '_pathways.csv'))

    filtPathways <- pathways[pathways$pSetRank < 0.05 |
                             (pathways$p.adjust < 0.05),]

    head(filtPathways)
    dim(filtPathways)

    write.csv(filtPathways, paste0('reactome/setrank/', make.names(groups[1]), '_',
                                   make.names(groups[2]), '_filtered_pathways.csv'))

    gsea.out@result <- gsea.out@result[rownames(filtPathways),]

    gsea.out@result$p.adjust <- filtPathways$p.adjust

    gsea.out@result$Description <- str_trunc(gsea.out$Description, 35)

    
    p1 <- dotplot(gsea.out, showCategory=20) +
        ggtitle(paste0('Reactome GSEA: ', paste(groups, collapse=' vs. ')))
    p1

    ggsave(paste0('reactome/setrank/dotplot_', make.names(groups[1]), '_',
                  make.names(groups[2]), '_20.pdf'),
           p1, width=7, height=9)

    ggsave(paste0('reactome/setrank/dotplot_', make.names(groups[1]), '_',
                  make.names(groups[2]), '_20.png'),
           p1, width=7, height=9, dpi=400, bg='white')

    p1 <- dotplot(gsea.out, showCategory=10) +
        ggtitle(paste0('Reactome GSEA: ', paste(groups, collapse=' vs. ')))
    p1

    ggsave(paste0('reactome/setrank/dotplot_', make.names(groups[1]), '_',
                  make.names(groups[2]), '_10.pdf'),
           p1, width=7, height=7)

    ggsave(paste0('reactome/setrank/dotplot_', make.names(groups[1]), '_',
                  make.names(groups[2]), '_10.png'),
           p1, width=7, height=7, dpi=400, bg='white')

    gsea.out <- pairwise_termsim(gsea.out, showCategory=20)

    ## Sys.sleep(5)

    p2 <- emapplot(gsea.out, showCategory=20, max.overlaps=0, repel=T,
                   force=5, node_label='category', layout='fr',
                   min_edge=0.05) +
        ggtitle(paste0('Reactome GSEA: ', paste(groups, collapse=' vs. ')))
    p2
    
    ggsave(paste0('reactome/setrank/emapplot_', make.names(groups[1]), '_',
                  make.names(groups[2]), '.pdf'),
           p2, width=12, height=12)

    ggsave(paste0('reactome/setrank/emapplot_', make.names(groups[1]), '_',
                  make.names(groups[2]), '.png'),
           p2, width=12, height=12, dpi=400, bg='white')
}


set.seed(20220819)

for (c in 1:4) {

    print(paste('Cluster', c))

    auc.rna <- read.csv(paste0('snf_hvar_clusts/markers/clust', c,'/clust', c,'_auc_rna_all_markers.csv'))

    auc.rna <- auc.rna[,-1]

    auc.rna <- auc.rna %>% filter(Cluster==c)

    rownames(auc.rna) <- gsub('RNA_', '', auc.rna$Feature)
    
    auc.stats <- as.vector(matrix(NA, nrow(ids)))
    names(auc.stats) <- ids$ENTREZID

    for ( i in 1:nrow(ids) ) {
        ## auc.stats[i] <- -log10(auc.rna[ids$SYMBOL[i], 'p.value'])
        auc.stats[i] <- auc.rna[ids$SYMBOL[i], paste0('AUROC.', c)]
    }

    auc.stats <- sort(auc.stats, decreasing=T)


    gsea.out <- gsePathway(auc.stats, pvalueCutoff=1,
                           pAdjustMethod='none', scoreType='pos')

    gsea.out <- setReadable(gsea.out, OrgDb = org.Hs.eg.db, keyType="ENTREZID")

    gsea.out@result$Description <- gsub('Homo sapiens\r: ', '', gsea.out$Description)

    ## gsea.out@result <- gsea.out@result[gsea.out$NES > 0,]

    sr.out <- setRankAnalysis(names(auc.stats), setCollection)

    exportSingleResult(sr.out, names(auc.stats),
                       setCollection, paste0('clust',c),
                       outputPath='reactome/setrank')

    pathways <- read.table(paste0('reactome/setrank/clust', c, '_pathways.txt'),
                           header=T, row.names=1, sep='\t')
    
    pathways$p.adjust <- p.adjust(pathways$correctedPValue, 'fdr', n=1438)

    head(pathways)
    dim(pathways)

    write.csv(pathways, paste0('reactome/setrank/clust', c, '_pathways.csv'))

    filtPathways <- pathways[pathways$pSetRank < 0.05 |
                             (pathways$p.adjust < 0.05),]

    head(filtPathways)
    dim(filtPathways)

    write.csv(filtPathways, paste0('reactome/setrank/clust', c, '_filtered_pathways.csv'))

    gsea.out@result <- gsea.out@result[rownames(filtPathways),]

    gsea.out@result$p.adjust <- filtPathways$p.adjust

    gsea.out@result$Description <- str_trunc(gsea.out$Description, 35)

    gsea.out

    ## write.csv(gsea.out@result, paste0('reactome/gsea_clust', c, '.csv'))

    p1 <- dotplot(gsea.out, showCategory=20) +
        ggtitle(paste0('Reactome Enrichment: Cluster ', c))
    p1
    
    ggsave(paste0('reactome/setrank/dotplot_clust', c, '_20.pdf'),
           p1, width=7, height=9)

    ggsave(paste0('reactome/setrank/dotplot_clust', c, '_20.png'),
           p1, width=7, height=9, dpi=400, bg='white')

    
    p1 <- dotplot(gsea.out, showCategory=10) +
        ggtitle(paste0('Reactome Enrichment: Cluster ', c))
    p1
    
    ggsave(paste0('reactome/setrank/dotplot_clust', c, '_10.pdf'),
           p1, width=7, height=7)

    ggsave(paste0('reactome/setrank/dotplot_clust', c, '_10.png'),
           p1, width=7, height=7, dpi=400, bg='white')

    gsea.out <- pairwise_termsim(gsea.out, showCategory=min(20, nrow(gsea.out)))

    ## Sys.sleep(5)

    p2 <- emapplot(gsea.out, showCategory=min(20, nrow(gsea.out)),
                   max.overlaps=0, repel=T, force=5, layout='fr',
                   min_edge=0.05) +
        ggtitle(paste0('Reactome Enrichment: Cluster ', c))
    p2
    
    ggsave(paste0('reactome/setrank/emapplot_clust', c, '.pdf'),
           p2, width=12, height=12)

    ggsave(paste0('reactome/setrank/emapplot_clust', c, '.png'),
           p2, width=12, height=12, dpi=400, bg='white')
}


#### setrank test

pathnames <- AnnotationDbi::keys(reactome.db::reactome.db, 'PATHID')
org.keep <- grepl('HSA', pathnames)
pathnames <- pathnames[org.keep]

reactome.df <- AnnotationDbi::select(reactome.db::reactome.db,
                                     columns=c('PATHID', 'PATHNAME', 'ENTREZID'),
                                     keys=pathnames,
                                     keytype='PATHID')

colnames(reactome.df) <- c('termID', 'termName', 'geneID')

reactome.df <- reactome.df[reactome.df$geneID %in% names(auc.stats),]

reactome.df$termName <- gsub('Homo sapiens\r: ', '', reactome.df$termName)
reactome.df$description <- reactome.df$termName

reactome.df$dbName <- 'REACTOME'

dim(reactome.df)

reactome.df <- reactome.df %>% group_by(termID) %>% filter(n() >= 10)

options(mc.cores=16)

setCollection <- buildSetCollection(reactome.df)

sr.out <- setRankAnalysis(names(auc.stats), setCollection)

exportSingleResult(sr.out, names(auc.stats), setCollection, paste0('clust',c), outputPath='reactome/setrank')

pathways <- read.table(paste0('reactome/setrank/clust', c, '_pathways.txt'), header=T, row.names=1, sep='\t')

head(pathways)
dim(pathways)

write.csv(pathways, paste0('reactome/setrank/clust', c, '_pathways.csv'))

filtPathways <- pathways[pathways$pSetRank < 0.05 | (pathways$correctedPValue < 0.001 & pathways$adjustedPValue < 0.001),]

head(filtPathways)
dim(filtPathways)

write.csv(filtPathways, paste0('reactome/setrank/clust', c, '_filtered_pathways.csv'))

dotplot(gsea.out)

sr.out.df <- igraph::as_data_frame(sr.out, what='vertices')

sr.out.df <- sr.out.df[order(sr.out.df$pSetRank, sr.out.df$adjustedPValue, sr.out.df$correctedPValue),]

sr.out.df

sr.out.df[sr.out.df$pSetRank<0.05,]

dim(sr.out.df[sr.out.df$pSetRank<0.05 | (sr.out.df$adjustedPValue<0.001 & sr.out.df$correctedPValue<0.001),])

sr.out.edges <- igraph::as_data_frame(sr.out, what='edges')

disc.sources <- sr.out.edges$from[
