library(dplyr)
library(stringr)
library(vsn)
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(readxl)
library(data.table)

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

## final.clin <- clin %>% dplyr::select(Condition., Batch, DiscWksGest, InfSex, Race, WksGest,
##                                      DiscPrePregBMI, PrePregBMI, Smoking, FDELTYPE,
##                                      Labor.initiation)

## final.clin <- clin %>% select(Condition., Batch, WksGest, InfSex, Race, PrePregBMI,
##                               Smoking, FDELTYPE, Labor.initiation)

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

counts <- read.csv('miniDP3_miRNA_rds.csv', row.names=1)
dim(counts)
head(counts)

idMap <- read_excel('Transcriptomic sample information.xlsx')
colnames(idMap) <- make.names(colnames(idMap))
idMap <- idMap[!is.na(idMap$Sample.ID),]

idMap$Sample.ID <- make.names(sapply(idMap$Sample.ID, function(x) sub(' .*', '', x)))
idMap$Study.ID <- make.names(idMap$Study.ID)
tail(idMap)

setnames(counts, idMap$Sample.ID, idMap$Study.ID)

dim(counts[,colnames(counts) %in% rownames(final.clin)])

final.clin$Batch <- make.names(final.clin$Batch)
final.clin$Condition. <- make.names(final.clin$Condition.)

## final.clin <- final.clin %>% mutate_if(function(x) mean(is.na(as.numeric(x))) < 0.5,
##                                        as.numeric)
## final.clin <- final.clin %>% mutate_if(function(x) length(unique(x)) > 7,
##                                        scale)

final.clin <- final.clin %>% mutate_at(c('Batch', 'Condition.', 'FDELTYPE', 'Labor.initiation', 'Labor.no.labor', 'Race', 'InfSex', 'Smoking'), factor)

head(final.clin)

dim(final.clin[rowSums(is.na(final.clin))==0,])
final.clin <- final.clin[rowSums(is.na(final.clin))==0,]

final.clin <- final.clin %>% mutate_at(c('WksGest', 'PrePregBMI'), function(x) scale(x)[,1])

samp.ids <- intersect(rownames(final.clin), colnames(counts))

dds <- DESeqDataSetFromMatrix(counts[,samp.ids], final.clin[samp.ids,],
                              design=~Condition. + WksGest + InfSex +
                                  Race + PrePregBMI + Smoking + FDELTYPE +
                                  Labor.initiation + Labor.no.labor)

## keep <- rowSums(counts(dds)) >= 50 & rowMeans(counts(dds)==0) < 0.95
keep <- rowMeans(counts(dds)) > 10
dds <- dds[keep,]

dim(dds)

library(sva)

combat.counts <- ComBat_seq(counts(dds),
                            batch=colData(dds)$Batch,
                            covar_mod=model.matrix(~ Condition. + WksGest + InfSex +
                                                       Race + PrePregBMI + Smoking +
                                                       FDELTYPE + Labor.initiation,
                                                   colData(dds)))

dim(dds)

write.csv(t(combat.counts), 'data/combat_mirna_counts.csv')

## combat.counts <- t(read.csv('data/combat_mirna_counts.csv', header=T, row.names=1))

## samp.ids <- colnames(combat.counts)

dds <- DESeqDataSetFromMatrix(combat.counts[,samp.ids],
                              final.clin[samp.ids,],
                              design=~Condition. + WksGest + InfSex + Race +
                                  PrePregBMI + Smoking + FDELTYPE +
                                  Labor.initiation)

dim(dds)

dds <- DESeq(dds, fitType='local')

plotDispEsts(dds)

vsd <- vst(dds, nsub=min(1000, nrow(dds)), blind=FALSE)

write.csv(t(assay(vsd)), 'data/combat_vst_mirna_expression.csv')

dds <- DESeqDataSetFromMatrix(combat.counts[,samp.ids],
                              final.clin[samp.ids,],
                              design=~Condition. + WksGest + InfSex + Race +
                                  PrePregBMI + Smoking + FDELTYPE +
                                  Labor.initiation + Labor.no.labor)

dim(dds)

dds <- DESeq(dds, fitType='local')

plotDispEsts(dds)

vsd <- vst(dds, nsub=min(1000, nrow(dds)), blind=FALSE)


meanSdPlot(assay(vsd))

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


for (groups in comparisons) {

    print(groups)
    groups <- make.names(groups)

    res <- results(dds, contrast=c('Condition.', groups), alpha=0.05, independentFiltering=TRUE)
    print(summary(res))

    ## resLFC <- lfcShrink(dds, contrast=c('Condition.', groups), res=res, type='ashr')
    ## print(summary(resLFC))
    
    de <- res[order(res$pvalue),] %>% as.data.frame
    de <- mutate(de, miRNA=rownames(de), .before=1) %>% filter(baseMean > 0)

    ## vsd <- vst(dds, blind=FALSE)

    group.meta <- colData(dds)[colData(dds)[,'Condition.'] %in% groups,]
        

    if (nrow(de[!is.na(de$padj) & de$padj<0.05,])>1) {


        if (min(nrow(de[!is.na(de$padj) & de$padj<0.05,]), 500) == 500) {
            col.labels <- de[1:500,'miRNA']
        } else {
            col.labels <- de[!is.na(de$padj) & de$padj<0.05,'miRNA']
        }

        scaled.sig <- scale(t(assay(vsd))[rownames(group.meta), col.labels])

        col.labels <- gsub('hsa[-]', '', col.labels)
        
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
                  '_miRNA_diffExp_heatmap.pdf',sep=''),
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
                  '_miRNA_diffExp_heatmap.png',sep=''),
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
            col.labels <- de[1:25,'miRNA']
        } else {
            col.labels <- de[!is.na(de$padj) & de$padj<0.05,'miRNA']
        }

        scaled.sig <- scale(t(assay(vsd))[rownames(group.meta), col.labels])

        col.labels <- gsub('hsa[-]', '', col.labels)
        
        range <- max(abs(scaled.sig))

        print(range)

        range <- min(max(abs(scaled.sig)),3)
        
        pdf(paste('plots-cond3/', str_replace(groups[1], ' ', '-'), '_',
                  str_replace(groups[2], ' ', '.'),
                  '_miRNA_diffExp_top25_heatmap_unlabeled.pdf',sep=''),
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
                  '_miRNA_diffExp_top25_heatmap_unlabeled.png',sep=''),
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
                        '_miRNA_diffExp_top25_heatmap_row_labels.csv',sep=''),
                  row.names=F)
        
        dev.off()

        
        ## pdf(paste('plots-cond3/', str_replace(groups[1], ' ', '-'), '_',
        ##           str_replace(groups[2], ' ', '.'),
        ##           '_miRNA_diffExp_heatmap.pdf',sep=''),
        ##     width=5 + 0.12 * nrow(group.meta),
        ##     height=3 + 0.12 * nrow(de[!is.na(de$padj) & de$padj<0.05,]))

        ## scaled.sig <- scale(t(assay(vsd))[rownames(group.meta),
        ##                                   de[!is.na(de$padj) & de$padj<0.05,'miRNA']])
        ## range <- min(max(abs(scaled.sig)),3)
        
        ## annot_df <- data.frame(Diagnosis = factor(diag.map[as.character(group.meta$Condition.)],
        ##                                           levels=diag.map[groups]),
        ##                        `Infant Sex` = factor(group.meta$InfSex),
        ##                        `Gestational Age` = group.meta$WksGest, 
        ##                        row.names=rownames(group.meta))

        ## temp.annot.cols <- list()
        ## for (g in groups) {
        ##     temp.annot.cols[['Diagnosis']][diag.map[g]] <- annot.colors[['Diagnosis']][diag.map[g]]
        ## }
        
        ## pheatmap(t(scaled.sig), breaks = seq(-range, range, length.out = 101),
        ##          annotation_col=annot_df, main=paste(diag.map[groups[2]], 'vs.', diag.map[groups[1]]),
        ##          ## cluster_rows=sum(de$adj.P.Val < 0.05)>1,
        ##          clustering_method='ward.D', clustering_distance_rows = "euclidean",
        ##          clustering_distance_cols = "euclidean", annotation_colors=temp.annot.cols,
        ##          labels_row=gsub('hsa[-]', '', de[!is.na(de$padj) & de$padj<0.05,'miRNA']))
        ## dev.off()


        ## png(paste('plots-cond3/', str_replace(groups[1], ' ', '-'), '_',
        ##           str_replace(groups[2], ' ', '.'),
        ##           '_miRNA_diffExp_heatmap.png',sep=''),
        ##     width=5 + 0.12 * nrow(group.meta),
        ##     height=3 + 0.12 * nrow(de[!is.na(de$padj) & de$padj<0.05,]),
        ##     units='in', res=400)

        ## pheatmap(t(scaled.sig), breaks = seq(-range, range, length.out = 101),
        ##          annotation_col=annot_df, main=paste(diag.map[groups[2]], 'vs.', diag.map[groups[1]]),
        ##          ## cluster_rows=sum(de$adj.P.Val < 0.05)>1,
        ##          clustering_method='ward.D', clustering_distance_rows = "euclidean",
        ##          clustering_distance_cols = "euclidean", annotation_colors=temp.annot.cols,
        ##          labels_row=gsub('hsa[-]', '', de[!is.na(de$padj) & de$padj<0.05,'miRNA']))
        ## dev.off()
    
    }

    ## de[!is.na(de$padj) & de$padj<0.05,]

    write.csv(de, paste('DE-cond3/diffExpMiRNA_', str_replace(groups[1], ' ', '-'), '_',
                        str_replace(groups[2], ' ', '.'), '.csv', sep=''))
    write.csv(de[!is.na(de$padj) & de$padj<0.05,],
              paste('DE-cond3/diffExpMiRNA_', str_replace(groups[1], ' ', '-'), '_',
                    str_replace(groups[2], ' ', '.'), '_FDR05.csv', sep=''))
    write.csv(de[!is.na(de$padj) & de$padj<0.1,],
              paste('DE-cond3/diffExpMiRNA_', str_replace(groups[1], ' ', '-'), '_',
                    str_replace(groups[2], ' ', '.'), '_FDR1.csv', sep=''))
    
}

dev.off()

plotMA(res)

resultsNames(dds)

resLFC <- lfcShrink(dds, contrast= c('Condition.', 'FGR', 'Control'), type='ashr')

plotMA(resLFC)

summary(resLFC)
head(as.data.frame(res[order(res$pvalue),]), 45)
head(as.data.frame(resLFC[order(resLFC$pvalue),]), 45)

vsd <- vst(dds, blind=FALSE)
meanSdPlot(assay(vsd))

pcaData <- plotPCA(vsd, 'Batch', returnData=TRUE)

kruskal.test(pcaData$PC1, pcaData$Batch)
kruskal.test(pcaData$PC2, pcaData$Batch)

pdf('corrected_pcaPlots_miRNA.pdf', width=7, height=5.5)
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
