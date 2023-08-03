library(Seurat)
library(InstaPrism)


raw.filt <- readRDS('../data/placenta_atlas_C_PE_annotated.rds')

dim(raw.filt)

table(raw.filt$Annotation)

raw.filt <- subset(raw.filt,
                   subset=Annotation!='B Cells' &
                       Annotation!='Erythroid Cells' &
                       Annotation!='Monocytes' &
                       Annotation!='T Cells' &
                       Annotation!='Extravillous Trophoblasts'
                   )

dim(raw.filt)

raw.filt$Annotation <- factor(raw.filt$Annotation)

table(raw.filt$Annotation)

## samp.names <- c('GSM5261695_C1', 'GSM5261696_C2', 'GSM5261699_P1', 'GSM5261700_P2')
## sample.meta <- data.frame(condition=c('C', 'C', 'P', 'P'), sample=samp.names)
## rownames(sample.meta) <- samp.names

## de.results <- pseudobulkDE(raw.filt, raw.filt$Sample, raw.filt$Annotation, sample.meta)

## de.results

## set.seed(230620)
## celltype.profiles <- pseudobulkProfiles(raw.filt, raw.filt$Annotation, 5)

## head(celltype.profiles)

## cell.counts <- as.matrix(raw.filt@assays$RNA@counts)
## colnames(cell.counts) <- raw.filt$Annotation

## write.table(cell.counts, 'cibersort/placenta_atlas_annotated_raw_counts.txt', sep='\t')


## write.table(celltype.profiles, 'cibersort/placenta_atlas_annotated_pseudobulk_replicates.txt', sep='\t')


## raw.filt <- raw.filt %>% subset(Annotation != "Low Quality Cells")

raw.filt$Condition <- factor(raw.filt$Condition, levels=c('C','P'),
                             labels=c('Control', 'PE')) 

DimPlot(raw.filt, label=T, repel=T, group.by='Annotation') + DimPlot(raw.filt, group.by='Condition')

cell.types <- ifelse(grepl('Cytotrophoblasts', as.character(raw.filt$Annotation)),
                     'Cytotrophoblasts',
              ifelse(grepl('Neutrophils', as.character(raw.filt$Annotation)),
                     'Granulocytes', as.character(raw.filt$Annotation)))

cell.states <- as.character(raw.filt$Annotation)

table(cell.types)
table(cell.states)

raw.filt@meta.data$`Cell.Type` <- cell.types

DimPlot(raw.filt, label=T, repel=T, group.by='Cell.Type')

ggsave('../out/scReference_deconv_celltypes_noEVT.png', width=10, height=9, dpi=400)
ggsave('../out/scReference_deconv_celltypes_noEVT.pdf', width=10, height=9)

DimPlot(raw.filt, label=T, group.by='Annotation') + DimPlot(raw.filt, group.by='Sample')


bulk <- t(read.csv('../data/combat_rna_counts.csv', row.names=1))

head(bulk)

require(snowfall)

ip.res.up <- InstaPrism_legacy(input_type='raw',
                     sc_Expr=as.matrix(raw.filt@assays$RNA@counts),
                     bulk_Expr=bulk,
                     cell.type.labels=cell.types,
                     cell.state.labels=cell.states,
                     return.Z.cs=FALSE, return.Z.ct=FALSE, n.core=1,
                     convergence.plot=TRUE, n.iter=100,
                     update=TRUE)

par(mfrow=c(4,2))
for (ct in unique(cell.types)) {
    plot(ip.res.up@Post.ini.ct@theta[ct,],
         ip.res.up@Post.updated.ct@theta[ct,],
         main=ct, xlab="Initial", ylab="Updated")
}


library(ggplot2)
library(reshape2)
library(stringr)

clin <- read_excel('../data/Primary clinical variables clustered.xlsx') %>% as.data.frame
colnames(clin) <- make.names(colnames(clin))
rownames(clin) <- make.names(clin$Study_ID.ID)
clin <- clin[-1,,drop=T]
head(clin)

clin <- clin %>% mutate_if(function(x) mean(is.na(as.numeric(x))) < 0.5, as.numeric)

clin$Condition. <- ifelse(clin$Condition. %in% c('PTD','Severe PE','FGR','Control','Control PTD'), clin$Condition., 'FGR+HDP')

## clin$Condition. <- make.names(clin$Condition.)

clin$PrePregBMI <- clin$PrePregWt_Kg / clin$Height_Meters.HeightMeters^2

clin$Smoking <- ifelse(str_detect(clin$Smoke, 'Yes'), 'Yes', NA)
clin$Smoking <- ifelse(str_detect(clin$Smoke, 'Never'), 'No', clin$Smoking)

clin$Race[!clin$Race %in% c('W', 'B')] <- 'other'
clin$FDELTYPE[clin$FDELTYPE==0] <- NA


## clin <- read.csv('clustered_clinical_hvar.csv', row.names=1)

## clin$Early.PE <- as.integer(clin$Condition. %in% c('Severe PE', 'FGR+HDP') & clin$WksGest < 34)

clin$Early.PE <- as.integer(clin$early.onset.PE=='1')

clin$Early.PTD <- as.integer(clin$WksGest < 32)

clin$Doppler.flow[clin$Doppler.flow==0] <- NA
clin$AEDF.REDF <- as.integer(clin$Doppler.flow==3 | clin$Doppler.flow==4)

clin$Nulliparity <- as.integer(clin$Para == 0)

clin$NICU.admission <- as.integer(clin$NICU_LOS > 0)

clin$Diabetes <- as.integer(clin$MENDDIAB != 0)

## clin <- clin %>% select(-c('High.grade.MVM','FVM',
##                            'Acute.inflammation','Chronic.inflammation'))

head(clin)

pathology.diag <- read_excel('../data/Path diagnoses for analysis.xlsx' , 1) %>% as.data.frame
colnames(pathology.diag)
rownames(pathology.diag) <- make.names(pathology.diag$case)
path.amsterdam <- pathology.diag %>% select(MVM, FVM, AI, CI)
path.amsterdam$MVM <- ifelse(path.amsterdam$MVM=='Y' & pathology.diag$`MVM only small placenta`==1, NA, path.amsterdam$MVM)

path.amsterdam['Mini.DP.282','MVM'] <- 'Y'
path.amsterdam[is.na(path.amsterdam)] <- 'N'
path.amsterdam <- path.amsterdam %>% mutate_all(factor)

head(path.amsterdam)
dim(path.amsterdam)

path.samp.ids <- intersect(rownames(clin), rownames(path.amsterdam))

clin[path.samp.ids, colnames(path.amsterdam)] <- path.amsterdam[path.samp.ids,]


pathology <- read_excel('../data/path features (from reports) for analysis.xlsx' , 1, na='n/a') %>% as.data.frame
colnames(pathology)
rownames(pathology) <- make.names(pathology$case)
pathology <- pathology[,c(-2:-1, -16)]

pathology[!(is.na(pathology[,1]) | pathology[,1] %in% c('<3', '3 to 5', '5 to 10')),1] <- '>=10'

pathology <- pathology %>% mutate_all(factor)
head(pathology)
dim(pathology)

path.samp.ids <- intersect(rownames(clin), rownames(pathology))

## clin$VUE <- NA

clin[path.samp.ids,
     colnames(pathology)] <- pathology[path.samp.ids,]

head(clin)

pathology.slides <- read_excel('../data/path features (from slides) for analysis.xlsx' , 1, na='n/a') %>% as.data.frame
colnames(pathology.slides)
rownames(pathology.slides) <- make.names(pathology.slides$`Study ID`)
pathology.slides <- pathology.slides[,4:6]

clin[intersect(rownames(clin), rownames(pathology.slides)),colnames(pathology.slides)] <- pathology.slides[intersect(rownames(clin),rownames(pathology.slides)),]

head(clin)

clin <- clin[!is.na(clin$hvar.Clusters),]

samp.ids <- intersect(rownames(clin), colnames(bulk))

dim(clin)

plotdf <- melt(ip.res.up@Post.updated.ct@theta,
               value.name='Proportion',
               varnames=c('Cell Type', 'Sample'))

cell.types <- c('Cytotrophoblasts', 'Granulocytes', 
                'Hofbauer Cells', 'NK Cells', 'Syncitiotrophoblasts',
                'Endothelial Cells', 'Fibroblasts')

clin.ct <- cbind(clin[samp.ids,], t(ip.res.up@Post.updated.ct@theta)[samp.ids,])

clin.ct <- clin.ct %>%
    mutate_at(rownames(ip.res.up@Post.updated.ct@theta), function(x) ifelse(x < 0.001, 0, x))


plotdf <- melt(clin.ct[,c('Study_ID.ID', 'hvar.Clusters', cell.types)],
               value.name='Proportion',
               variable.name='Cell Type', id.vars=c('Study_ID.ID', 'hvar.Clusters'))


ordered.means <- plotdf %>% group_by(`Cell Type`) %>% summarize_at('Proportion', mean) %>% arrange(Proportion)

plotdf[,'Cell Type'] <- factor(plotdf[,'Cell Type'],
                               levels=as.vector(ordered.means$`Cell Type`),
                               labels=c("Granulocytes", "NK Cells", "Endothelial Cells",
                                        "Hofbauer Cells", "Cytotrophoblasts", "Fibroblasts",
                                        'Syncytiotrophoblasts'))

ordered.means

ordered.samp.ids <- rownames(clin.ct)[order(clin.ct[,'Syncitiotrophoblasts'], decreasing=T)]

plotdf <- melt(clin.ct[ordered.samp.ids,c('Study_ID.ID', 'hvar.Clusters', cell.types)],
               value.name='Proportion',
               variable.name='Cell Type', id.vars=c('Study_ID.ID', 'hvar.Clusters'))

plotdf[,'Cell Type'] <- factor(plotdf[,'Cell Type'],
                               levels=rev(as.vector(ordered.means$`Cell Type`)),
                               labels=rev(c("Granulocytes", "NK Cells", "Endothelial Cells",
                                            "Hofbauer Cells", "Cytotrophoblasts",
                                            "Fibroblasts", 'Syncytiotrophoblasts')))

plotdf$Sample <- factor(make.names(plotdf$Study_ID.ID), levels=ordered.samp.ids)


pbar <- ggplot(plotdf,
               aes(x=Sample, y=Proportion, fill=`Cell Type`)) +
    geom_bar(stat='identity') +
    scale_x_discrete(labels = NULL, breaks = NULL) +
    scale_fill_brewer(palette='Set3') + 
    theme_bw()

pbar

ggsave('../out/instaprism_celltype_proportions_bar_noEVT.pdf', pbar, width=10, height=5)
ggsave('../out/instaprism_celltype_proportions_bar_noEVT.png', pbar, width=10, height=5, dpi=400)

pbar <- ggplot(plotdf %>% filter(hvar.Clusters==1),
               aes(x=Sample, y=Proportion, fill=`Cell Type`)) +
    geom_bar(stat='identity') +
    scale_x_discrete(labels = NULL, breaks = NULL) +
    scale_fill_brewer(palette='Set3') + 
    theme_bw()

pbar

ggsave('../out/instaprism_celltype_proportions_bar_noEVT_clust1.pdf', pbar, width=10, height=5)
ggsave('../out/instaprism_celltype_proportions_bar_noEVT_clust1.png', pbar, width=10, height=5, dpi=400)

pbar <- ggplot(plotdf %>% filter(hvar.Clusters==2),
               aes(x=Sample, y=Proportion, fill=`Cell Type`)) +
    geom_bar(stat='identity') +
    scale_x_discrete(labels = NULL, breaks = NULL) +
    scale_fill_brewer(palette='Set3') + 
    theme_bw()

pbar

ggsave('../out/instaprism_celltype_proportions_bar_noEVT_clust2.pdf', pbar, width=10, height=5)
ggsave('../out/instaprism_celltype_proportions_bar_noEVT_clust2.png', pbar, width=10, height=5, dpi=400)

pbar <- ggplot(plotdf %>% filter(hvar.Clusters==3),
               aes(x=Sample, y=Proportion, fill=`Cell Type`)) +
    geom_bar(stat='identity') +
    scale_x_discrete(labels = NULL, breaks = NULL) +
    scale_fill_brewer(palette='Set3') + 
    theme_bw()

pbar

ggsave('../out/instaprism_celltype_proportions_bar_noEVT_clust3.pdf', pbar, width=10, height=5)
ggsave('../out/instaprism_celltype_proportions_bar_noEVT_clust3.png', pbar, width=10, height=5, dpi=400)

pbar <- ggplot(plotdf %>% filter(hvar.Clusters==4),
               aes(x=Sample, y=Proportion, fill=`Cell Type`)) +
    geom_bar(stat='identity') +
    scale_x_discrete(labels = NULL, breaks = NULL) +
    scale_fill_brewer(palette='Set3') + 
    theme_bw()

pbar

ggsave('../out/instaprism_celltype_proportions_bar_noEVT_clust4.pdf', pbar, width=10, height=5)
ggsave('../out/instaprism_celltype_proportions_bar_noEVT_clust4.png', pbar, width=10, height=5, dpi=400)


plotdf <- melt(cbind(t(ip.res.up@Post.updated.ct@theta)[samp.ids,],
                     clin[samp.ids,c('hvar.Clusters', 'Condition.', 'Study_ID.ID')]),
               value.name='Proportion',
               variable.name='Cell Type',
               id.vars=c('hvar.Clusters', 'Condition.', 'Study_ID.ID'))

cell.types <- c('Cytotrophoblasts', 'Granulocytes', 'Extravillous Trophoblasts',
                'Hofbauer Cells', 'Neutrophils', 'NK Cells', 'Syncitiotrophoblasts',
                'Endothelial Cells', 'Fibroblasts')

plotdf <- melt(cbind(clin.ct[samp.ids,cell.types],
                     clin[samp.ids,c('hvar.Clusters', 'Condition.', 'Study_ID.ID')]),
               value.name='Proportion',
               variable.name='Cell Type',
               id.vars=c('hvar.Clusters', 'Condition.', 'Study_ID.ID'))


colnames(plotdf)[1:3] <- c('Clusters', 'Diagnosis', 'Sample')

plotdf$Clusters <- factor(plotdf$Clusters, levels=1:4)

plotdf$Diagnosis <- factor(plotdf$Diagnosis,
                           levels=c('Control', 'Control PTD', 'PTD',
                                    'FGR', 'FGR+HDP', 'Severe PE'),
                           labels=c('Control', 'Control PT', 'PTD',
                                    'FGR', 'FGR+HDP', 'Severe PE'))

plotdf$`Cell Type` <- factor(plotdf$`Cell Type`,
                             levels=rev(as.character(ordered.means$`Cell Type`)))

plotdf$`Cell Type` <- factor(plotdf$`Cell Type`,
                             levels=rev(names(ordered.means)))

plotdf$Proportion <- ifelse(plotdf$Proportion < 0.001, 0, plotdf$Proportion)

library(ggpubr)
library(RColorBrewer)

brewer.colors <- brewer.pal(n=12, "Paired")## )(10)

brewer.colors <- c(brewer.colors[1:2], brewer.colors[11], brewer.colors[5:6], brewer.colors[8], brewer.colors[3:4], brewer.colors[9:10])

library(rstatix)

cluster.tests <- rbind(clin.ct %>% kruskal_test(Cytotrophoblasts ~ hvar.Clusters),
                       clin.ct %>% kruskal_test(Granulocytes ~ hvar.Clusters),
                       clin.ct %>% kruskal_test(`Hofbauer Cells` ~ hvar.Clusters),
                       clin.ct %>% kruskal_test(`NK Cells` ~ hvar.Clusters),
                       clin.ct %>% kruskal_test(Syncitiotrophoblasts ~ hvar.Clusters),
                       clin.ct %>% kruskal_test(`Endothelial Cells` ~ hvar.Clusters),
                       clin.ct %>% kruskal_test(Fibroblasts ~ hvar.Clusters))
## condition.tests$padj <- p.adjust(condition.tests$p, 'holm')

cluster.tests <- cluster.tests %>% adjust_pvalue()
cluster.tests

write.csv(cluster.tests, '../out/cluster_instaprism_celltype_proportions_noEVT.csv')

ph.cluster.tests <- rbind(clin.ct %>% dunn_test(Cytotrophoblasts ~ hvar.Clusters),
                          clin.ct %>% dunn_test(Granulocytes ~ hvar.Clusters),
                          clin.ct %>% dunn_test(`Hofbauer Cells` ~ hvar.Clusters),
                          clin.ct %>% dunn_test(`NK Cells` ~ hvar.Clusters),
                          clin.ct %>% dunn_test(Syncitiotrophoblasts ~ hvar.Clusters),
                          clin.ct %>% dunn_test(`Endothelial Cells` ~ hvar.Clusters),
                          clin.ct %>% dunn_test(Fibroblasts ~ hvar.Clusters))

ph.cluster.tests
ph.cluster.tests[ph.cluster.tests$p.adj<0.05,] %>% as.data.frame

write.csv(ph.cluster.tests, '../out/cluster_instaprism_celltype_proportions_dunn_posthoc_wEVT.csv')

pdf('../out/instaprism_cluster_celltype_violin_noEVT.pdf', width=6, height=5)
cell.types <- ordered.means$`Cell Type`
for (ct in cell.types) {
    gg <- ggplot(plotdf %>% filter(`Cell Type`==ct),
                 aes(y=Proportion, fill=Clusters, x=Clusters)) +
        geom_violin() +
        geom_boxplot(width=0.1) +
        annotate(geom='text',
                 y=max(plotdf[which(plotdf$`Cell Type`==ct),'Proportion']) + 0.1 * max(plotdf[which(plotdf$`Cell Type`==ct),'Proportion']), x=2,
                 label=paste("Kruskal-Wallis, p =",
                             cluster.tests[which(cluster.tests$.y.==ct),'p.adj'])) + 
        facet_wrap(vars(`Cell Type`), scales='free', ncol=4) +
        scale_fill_manual(values=brewer.colors[7:10]) +
        scale_x_discrete(labels = NULL, breaks = NULL) +
        ## geom_pwc(method='dunn_test') + 
        theme_bw()
    print(gg)
}
dev.off()
gg


ggsave('../out/instaprism_cluster_celltype_proportions_violin.pdf', width=15, height=9)
ggsave('../out/instaprism_cluster_celltype_proportions_violin.png', width=15, height=9, dpi=400)

condition.tests <- rbind(clin.ct %>% kruskal_test(Cytotrophoblasts ~ Condition.),
                         clin.ct %>% kruskal_test(Granulocytes ~ Condition.),
                         clin.ct %>% kruskal_test(`Hofbauer Cells` ~ Condition.),
                         clin.ct %>% kruskal_test(`NK Cells` ~ Condition.),
                         clin.ct %>% kruskal_test(Syncitiotrophoblasts ~ Condition.),
                         clin.ct %>% kruskal_test(`Endothelial Cells` ~ Condition.),
                         clin.ct %>% kruskal_test(Fibroblasts ~ Condition.))

condition.tests <- condition.tests %>% adjust_pvalue()
condition.tests

write.csv(condition.tests, '../out/condition_instaprism_celltype_proportions_noEVT.csv')

ph.condition.tests <- rbind(clin.ct %>% dunn_test(Cytotrophoblasts ~ Condition.),
                          clin.ct %>% dunn_test(Granulocytes ~ Condition.),
                          clin.ct %>% dunn_test(`Hofbauer Cells` ~ Condition.),
                          clin.ct %>% dunn_test(`NK Cells` ~ Condition.),
                          clin.ct %>% dunn_test(Syncitiotrophoblasts ~ Condition.),
                          clin.ct %>% dunn_test(`Endothelial Cells` ~ Condition.),
                          clin.ct %>% dunn_test(Fibroblasts ~ Condition.))

ph.condition.tests
ph.condition.tests[ph.condition.tests$p.adj<0.05,] %>% as.data.frame

write.csv(ph.condition.tests, '../out/condition_instaprism_celltype_proportions_dunn_posthoc_noEVT.csv')

pdf('../out/instaprism_condition_celltype_violin_noEVT.pdf', width=6, height=5)
cell.types <- ordered.means$`Cell Type`
for (ct in cell.types) {
    gg <- ggplot(plotdf %>% filter(`Cell Type`==ct),
                 aes(y=Proportion, fill=Diagnosis, x=Diagnosis)) +
        geom_violin() +
        geom_boxplot(width=0.1) +
        annotate(geom='text',
                 y=max(plotdf[which(plotdf$`Cell Type`==ct),'Proportion']) + 0.1 * max(plotdf[which(plotdf$`Cell Type`==ct),'Proportion']), x=2,
                 label=paste("Kruskal-Wallis, p =",
                             condition.tests[which(condition.tests$.y.==ct),'p.adj'])) + 
        facet_wrap(vars(`Cell Type`), scales='free', ncol=4) +
        scale_fill_manual(values=brewer.colors[1:6]) +
        scale_x_discrete(labels = NULL, breaks = NULL) +
        ## geom_pwc(method='dunn_test') + 
        theme_bw()
    print(gg)
}
dev.off()
