library(dplyr)
library(ggplot2)
library(ggrepel)

comparisons <- list(c('Control PTD', 'Control'),
                    c('FGR', 'Control'),
                    c('Severe PE', 'Control PTD'),
                    c('Severe PE', 'Control'),
                    c('PTD', 'Control PTD'),
                    c('PTD', 'Control'),
                    c('FGR+HDP', 'Control PTD'),
                    c('FGR+HDP', 'Control'))


#### Proteins

for (groups in comparisons) {
    de.prot <- read.csv(paste0('DE-cond3/diffExpProt_',
                               gsub(' ', '-', groups[1]),
                               '_',
                               gsub(' ', '.', groups[2]),
                               '.csv'))

    volplot.df <- data.frame(log2FoldChange=de.prot$logFC,
                             pvalue=de.prot$P.Value,
                             Significant=factor(ifelse(de.prot$adj.P.Val < 0.05,
                                                       'Significant', 'Not Significant'),
                                               levels=c('Not Significant', 'Significant')),
                             Label=ifelse(1:nrow(de.prot) <= 5 & de.prot$adj.P.Val < 0.05, de.prot$Assay, NA))

    max.sig.idx <- length(de.metab$adj.P.Val[de.metab$adj.P.Val < 0.05])

    groups <- gsub('Control PTD', 'Control PT', groups)

    xrange <- max(volplot.df$log2FoldChange) - min(volplot.df$log2FoldChange)

    ggvol <- ggplot(volplot.df, aes(x=log2FoldChange, y=-log10(pvalue),
                                    col=Significant, label=Label)) +
        geom_point() +
        scale_color_manual(values=c('black', 'red')) +
        geom_vline(xintercept=c(-1, 1), col='darkgrey', lty=5) +
        ## geom_hline(yintercept=-log10(mean(volplot.df[max.sig.idx:(max.sig.idx+1),'pvalue'])),
        ##            col='darkgrey', lty=5) +
        labs(title=paste('Proteins', groups[1], 'vs', groups[2])) +
        guides(col=FALSE) +
        expand_limits(y=1.1 * max(-log10(volplot.df$pvalue)),
                      x=c(min(volplot.df$log2FoldChange) - 0.05*xrange,
                          max(volplot.df$log2FoldChange) + 0.05*xrange)) + 
        theme_bw()

    if (max.sig.idx > 0) {
        ggvol <- ggvol +
            geom_hline(yintercept=-log10(mean(volplot.df[max.sig.idx:(max.sig.idx+1),'pvalue'])),
                       col='darkgrey', lty=5)
    }

    ggvol <- ggvol + geom_text_repel(col='black', nudge_y = 0, force=5, 
                                     size=3.5, point.padding=0.25, ## box.padding=1,
                                     max.overlaps=Inf)

    groups <- gsub('Control PT', 'Control PTD', groups)

    ggsave(paste0('volcano_plots/diffExpProt_',
                  gsub(' ', '-', groups[1]),
                  '_',
                  gsub(' ', '.', groups[2]),
                  '_volPlot.png'), ggvol, width=5, height=5, dpi=400)

        ggsave(paste0('volcano_plots/diffExpProt_',
                  gsub(' ', '-', groups[1]),
                  '_',
                  gsub(' ', '.', groups[2]),
                  '_volPlot.pdf'), ggvol, width=5, height=5)
}


#### Metabolites

for (groups in comparisons) {
    de.metab <- read.csv(paste0('DE-cond3/diffExpMetab_',
                               gsub(' ', '-', groups[1]),
                               '_',
                               gsub(' ', '.', groups[2]),
                               '.csv'))

    volplot.df <- data.frame(log2FoldChange=de.metab$logFC,
                             pvalue=de.metab$P.Value,
                             Significant=factor(ifelse(de.metab$adj.P.Val < 0.05,
                                                       'Significant', 'Not Significant'),
                                               levels=c('Not Significant', 'Significant')),
                             Label=ifelse(1:nrow(de.metab) <= 5 & de.metab$adj.P.Val < 0.05,
                                          de.metab$CHEMICAL_NAME, NA))

    max.sig.idx <- length(de.metab$adj.P.Val[de.metab$adj.P.Val < 0.05])

    groups <- gsub('Control PTD', 'Control PT', groups)

    xrange <- max(volplot.df$log2FoldChange) - min(volplot.df$log2FoldChange)

    ggvol <- ggplot(volplot.df, aes(x=log2FoldChange, y=-log10(pvalue),
                                    col=Significant, label=Label)) +
        geom_point() +
        scale_color_manual(values=c('black', 'red')) +
        geom_vline(xintercept=c(-1, 1), col='darkgrey', lty=5) +
        ## geom_hline(yintercept=-log10(mean(volplot.df[max.sig.idx:(max.sig.idx+1),'pvalue'])),
        ##            col='darkgrey', lty=5) +
        labs(title=paste('Metabolites', groups[1], 'vs', groups[2])) +
        guides(col=FALSE) +
        expand_limits(y=1.1 * max(-log10(volplot.df$pvalue)),
                      x=c(min(volplot.df$log2FoldChange) - 0.05*xrange,
                          max(volplot.df$log2FoldChange) + 0.05*xrange)) + 
        theme_bw()

    if (max.sig.idx > 0) {
        ggvol <- ggvol +
            geom_hline(yintercept=-log10(mean(volplot.df[max.sig.idx:(max.sig.idx+1),'pvalue'])),
                       col='darkgrey', lty=5)
    }

    ggvol <- ggvol + geom_text_repel(col='black', nudge_y = 0, force=5, 
                                     size=3.5, point.padding=0.25, ## box.padding=1,
                                     max.overlaps=Inf)

    groups <- gsub('Control PT', 'Control PTD', groups)

    ggsave(paste0('volcano_plots/diffExpMetab_',
                  gsub(' ', '-', groups[1]),
                  '_',
                  gsub(' ', '.', groups[2]),
                  '_volPlot.png'), ggvol, width=5, height=5, dpi=400)

        ggsave(paste0('volcano_plots/diffExpMetab_',
                  gsub(' ', '-', groups[1]),
                  '_',
                  gsub(' ', '.', groups[2]),
                  '_volPlot.pdf'), ggvol, width=5, height=5)
}



#### miRNA

for (groups in comparisons) {
    de.mirna <- read.csv(paste0('DE-cond3/diffExpMiRNA_',
                               make.names(groups[1]),
                               '_',
                               make.names(groups[2]),
                               '.csv'))

    de.mirna

    volplot.df <- data.frame(log2FoldChange=de.mirna$log2FoldChange,
                             pvalue=de.mirna$pvalue,
                             Significant=factor(ifelse(!is.na(de.mirna$padj) &
                                                       de.mirna$padj < 0.05,
                                                       'Significant', 'Not Significant'),
                                               levels=c('Not Significant', 'Significant')),
                             Label=ifelse(1:nrow(de.mirna) <= 5 & de.mirna$padj < 0.05,
                                          de.mirna$miRNA, NA))

    max.sig.idx <- length(de.mirna$padj[!is.na(de.mirna$padj) & de.mirna$padj < 0.05])

    groups <- gsub('Control PTD', 'Control PT', groups)

    xrange <- max(volplot.df$log2FoldChange) - min(volplot.df$log2FoldChange)

    ggvol <- ggplot(volplot.df, aes(x=log2FoldChange, y=-log10(pvalue),
                                    col=Significant, label=Label)) +
        geom_point() +
        scale_color_manual(values=c('black', 'red')) +
        geom_vline(xintercept=c(-1, 1), col='darkgrey', lty=5) +
        ## geom_hline(yintercept=-log10(mean(volplot.df[max.sig.idx:(max.sig.idx+1),'pvalue'])),
        ##            col='darkgrey', lty=5) +
        labs(title=paste('miRNA', groups[1], 'vs', groups[2])) +
        guides(col=FALSE) +
        expand_limits(y=1.1 * max(-log10(volplot.df$pvalue)),
                      x=c(min(volplot.df$log2FoldChange) - 0.05*xrange,
                          max(volplot.df$log2FoldChange) + 0.05*xrange)) + 
        theme_bw()

    if (max.sig.idx > 0) {
        ggvol <- ggvol +
            geom_hline(yintercept=-log10(mean(volplot.df[max.sig.idx:(max.sig.idx+1),'pvalue'])),
                       col='darkgrey', lty=5)
    }

    ggvol <- ggvol + geom_text_repel(col='black', nudge_y = 0, force=5, 
                                     size=3.5, point.padding=0.25, ## box.padding=1,
                                     max.overlaps=Inf)

    groups <- gsub('Control PT', 'Control PTD', groups)

    ggsave(paste0('volcano_plots/diffExpMiRNA_',
                  make.names(groups[1]),
                  '_',
                  make.names(groups[2]),
                  '_volPlot.png'), ggvol, width=5, height=5, dpi=400)

        ggsave(paste0('volcano_plots/diffExpMiRNA_',
                  make.names(groups[1]),
                  '_',
                  make.names(groups[2]),
                  '_volPlot.pdf'), ggvol, width=5, height=5)
}



#### RNA

for (groups in comparisons) {
    de.rna <- read.csv(paste0('DE-cond3/diffExpRNA_',
                               make.names(groups[1]),
                               '_',
                               make.names(groups[2]),
                               '.csv'))

    de.rna

    volplot.df <- data.frame(log2FoldChange=de.rna$log2FoldChange,
                             pvalue=de.rna$pvalue,
                             Significant=factor(ifelse(!is.na(de.rna$padj) &
                                                       de.rna$padj < 0.05,
                                                       'Significant', 'Not Significant'),
                                               levels=c('Not Significant', 'Significant')),
                             Label=ifelse(1:nrow(de.rna) <= 5 & de.rna$padj < 0.05,
                                          de.rna$Gene, NA))

    max.sig.idx <- length(de.rna$padj[!is.na(de.rna$padj) & de.rna$padj < 0.05])

    groups <- gsub('Control PTD', 'Control PT', groups)

    xrange <- max(volplot.df$log2FoldChange) - min(volplot.df$log2FoldChange)

    ggvol <- ggplot(volplot.df, aes(x=log2FoldChange, y=-log10(pvalue),
                                    col=Significant, label=Label)) +
        geom_point() +
        scale_color_manual(values=c('black', 'red')) +
        geom_vline(xintercept=c(-1, 1), col='darkgrey', lty=5) +
        ## geom_hline(yintercept=-log10(mean(volplot.df[max.sig.idx:(max.sig.idx+1),'pvalue'])),
        ##            col='darkgrey', lty=5) +
        labs(title=paste('RNA', groups[1], 'vs', groups[2])) +
        guides(col=FALSE) +
        expand_limits(y=1.1 * max(-log10(volplot.df$pvalue)),
                      x=c(min(volplot.df$log2FoldChange) - 0.05*xrange,
                          max(volplot.df$log2FoldChange) + 0.05*xrange)) + 
        theme_bw()

    if (max.sig.idx > 0) {
        ggvol <- ggvol +
            geom_hline(yintercept=-log10(mean(volplot.df[max.sig.idx:(max.sig.idx+1),'pvalue'])),
                       col='darkgrey', lty=5)
    }

    ggvol <- ggvol + geom_text_repel(col='black', nudge_y = 0, force=5, 
                                     size=3.5, point.padding=0.25, ## box.padding=1,
                                     max.overlaps=Inf)

    groups <- gsub('Control PT', 'Control PTD', groups)

    ggsave(paste0('volcano_plots/diffExpRNA_',
                  make.names(groups[1]),
                  '_',
                  make.names(groups[2]),
                  '_volPlot.png'), ggvol, width=5, height=5, dpi=400)

        ggsave(paste0('volcano_plots/diffExpRNA_',
                  make.names(groups[1]),
                  '_',
                  make.names(groups[2]),
                  '_volPlot.pdf'), ggvol, width=5, height=5)
}
