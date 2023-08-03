library(dplyr)
library(stringr)
library(fgsea)
library(limma)
library(clusterProfiler)
library(ggplot2)
library(RColorBrewer)
library(forcats)

comparisons <- list(c('Control PTD', 'Control'),
                    c('FGR', 'Control'),
                    c('Severe PE', 'Control PTD'),
                    c('Severe PE', 'Control'),
                    c('PTD', 'Control PTD'),
                    c('PTD', 'Control'),
                    c('FGR+HDP', 'Control PTD'),
                    c('FGR+HDP', 'Control'))

groups <- c('FGR+HDP', 'Control PTD')

de.res <- read.csv(paste('DE-cond3/diffExpMetab_', str_replace(groups[1], ' ', '-'), '_',
                         str_replace(groups[2], ' ', '.'), '.csv', sep=''), row.names=1)

pathways <- unique(de.res$SUB_PATHWAY)
pathways <- pathways[!is.na(pathways)]

metab.sets <- list()
metab.sets.df <- data.frame()

for (i in 1:length(pathways)) {
    ## idxs <- which(de.res$SUB_PATHWAY==pathways[i])
    ## if (length(idxs) >= 3) {
    ##     metab.sets[[pathways[i]]] <- idxs
    ## }
    ## if (de.res$SUPER_PATHWAY[which(de.res$SUB_PATHWAY==pathways[i])][1] == 'Lipid') {
    if (TRUE) {    
        metab.sets[[pathways[i]]] <- rownames(de.res)[which(de.res$SUB_PATHWAY==pathways[i])]
    

        metab.sets.df <- rbind(metab.sets.df,
                               data.frame(pathway=rep(pathways[i],
                                                      length(metab.sets[[pathways[i]]])),
                                          metabolite=metab.sets[[pathways[i]]]))
    }
}

sphingo.pathways <- c('Sphingolipid Synthesis', 'Dihydrosphingomyelins', 'Sphingosines', 'Sphingomyelins', 'Hexosylceramides (HCER)', 'Ceramides', 'Dihydroceramides', 'Lactosylceramides (LCER)', 'Ceramide PEs')

metab.sets[['Sphingolipids']] <- unlist(metab.sets[sphingo.pathways], use.names=F)

metab.sets.df <- rbind(metab.sets.df,
                       data.frame(pathway=rep('Sphingolipids',
                                              length(metab.sets[['Sphingolipids']])),
                                  metabolite=metab.sets[['Sphingolipids']]))

## pvals <- sapply(metab.sets, function(set) gset(S=set, r=stats))

## head(sort(p.adjust(pvals, 'fdr')),10)
 
## plotdf <- data.frame(Pathway=c(), Comparison=c(), NES=c(), padj=c(), Size=c())
## sig.mask <- rep(FALSE, 89)
set.seed(20220215)
for (groups in comparisons) {

    print(paste(groups, collapse=' vs '))

    de.res <- read.csv(paste('DE-cond3/diffExpMetab_', str_replace(groups[1], ' ', '-'), '_',
                             str_replace(groups[2], ' ', '.'), '.csv', sep=''), row.names=1)

    ## stats <- ifelse(de.res$logFC > 0,
    ##          ifelse(de.res$CI.L < 0, 0, de.res$CI.L),
    ##          ifelse(de.res$CI.R > 0, 0, -de.res$CI.R))
    stats <- de.res$t
    names(stats) <- rownames(de.res)

    ## msea.res <- fgsea(metab.sets, stats, minSize=3)

    msea.res <- GSEA(sort(stats, decreasing=T), minGSSize=3, TERM2GENE=metab.sets.df, pvalueCutoff=0.05)
    
    ## msea.res <- msea.res[order(msea.res$pval),]
    ## msea.res$pathway <- factor(msea.res$pathway)
    print(msea.res[,1:7], 10)

    ## sig.mask <- sig.mask | (msea.res$padj<0.1 & !is.na(msea.res$padj))
    ## temp <- data.frame(Pathway=msea.res$pathway[sig.mask],
    ##                    Comparison=rep(paste(groups, collapse=' vs '), sum(sig.mask)),
    ##                    NES=msea.res$NES[sig.mask],
    ##                    padj=msea.res$NES[sig.mask],
    ##                    Size=msea.res$NES[sig.mask])

    ## plotdf <- rbind(plotdf, temp)

    if (nrow(msea.res) != 0) {

        colors <- rev(brewer.pal(n = 7, name = "RdYlBu"))
        pdf(paste(paste(make.names(groups), collapse='.v.'), 'sphingolipid_group_pathway_results-cond3.pdf', sep='_'), width=6, height=6)
        gg <- ggplot(msea.res %>% as.data.frame,
                     aes(x=-log10(pvalue),
                         y=fct_reorder(Description, desc(pvalue)),
                         fill=NES)) +
            geom_bar(stat='identity') +
            scale_fill_gradientn(colors=colors, limits=c(-2.6, 2.6)) +
            scale_y_discrete(labels = function(x) str_wrap(x, width=16)) + 
            ylab('Pathway') + 
            theme_bw()

        print(gg)
        
        dev.off()

    }
}


nes.mat <- matrix(0, nrow=sum(sig.mask), ncol=length(comparisons)) %>% as.data.frame
rownames(nes.mat) <- msea.res$pathway[sig.mask]
comp.idx <- 1
for (groups in comparisons) {

    print(paste(groups, collapse=' vs '))

    de.res <- read.csv(paste('DE-cond2/diffExpMetab_', str_replace(groups[1], ' ', '-'), '_',
                             str_replace(groups[2], ' ', '.'), '.csv', sep=''), row.names=1)

    stats <- ifelse(de.res$logFC > 0,
             ifelse(de.res$CI.L < 0, de.res$logFC/10, de.res$CI.L),
             ifelse(de.res$CI.R > 0, -de.res$logFC/10, -de.res$CI.R))
    ## stats <- de.res$t
    names(stats) <- rownames(de.res)

    msea.res <- fgsea(metab.sets, stats, minSize=3)

    ## msea.res <- msea.res[order(msea.res$pval),]
    print(head(msea.res[order(msea.res$pval),1:7], 10))

    nes.mat[,comp.idx] <- msea.res$NES[sig.mask]

    colnames(nes.mat)[comp.idx] <- paste(groups, collapse=' vs ')

    ## sig.mask <- sig.mask | (msea.res$padj<0.1 & !is.na(msea.res$padj))
    ## temp <- data.frame(Pathway=msea.res$pathway[sig.mask],
    ##                    Comparison=rep(paste(groups, collapse=' vs '), sum(sig.mask)),
    ##                    NES=msea.res$NES[sig.mask],
    ##                    padj=msea.res$NES[sig.mask],
    ##                    Size=msea.res$NES[sig.mask])

    ## plotdf <- rbind(plotdf, temp)
    comp.idx <- comp.idx + 1
}
