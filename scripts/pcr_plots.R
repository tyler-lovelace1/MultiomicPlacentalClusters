library(dplyr)
library(reshape2)
library(readxl)
library(ggplot2)
library(outliers)

pcr.data <- read_excel('plasma pcr for Tyler.xlsx', sheet=1)

head(pcr.data)

pcr.data <- na.omit(pcr.data[!apply(sweep(pcr.data, 2, colnames(pcr.data), '=='), 1, any, na.rm=T),])

pcr.data <- pcr.data %>% mutate_at('Ave Fold', as.numeric)
## pcr.data <- pcr.data %>% mutate_at('Ave Fold', log2)

pcr.data$Primer <- sapply(strsplit(pcr.data$Primer, ' |[(]'), function(x) x[1])

genes <- unique(pcr.data$Primer)

marker.map <- c(ALPP=1, PAPPA=1, LGR5=2, DUSP9=2, HTRA4=3, FLT1=3, LYVE1=4, EDNRB=4)

pcr.data$Groups <- ifelse(pcr.data$cluster==marker.map[pcr.data$Primer],
                          'Cluster', 'Rest')

ggplot(pcr.data, aes(x=Groups, y=`Ave Fold`, fill=Groups)) +
    geom_boxplot() +
    facet_wrap(~factor(Primer, levels=genes), scales='free', ncol=2) +
    theme_classic()



grubbs.pvals <- c()
box.outlier.rows <- c()
z.outlier.rows <- c()
par(mfrow=c(4, 2))
for (gene in genes) {
    print(gene)
    expr.vec <- as.vector(pcr.data[pcr.data$Primer==gene & pcr.data$Groups=='Rest','Ave Fold'])$`Ave Fold`
    modZ <- (expr.vec-median(expr.vec))/mad(expr.vec)
    box.out <- boxplot(expr.vec, main=gene)
    print(paste('IQR Outliers:', paste(which(as.vector(pcr.data[,'Ave Fold'])$`Ave Fold` %in% box.out$out), collapse=', ')))
    print(paste('Robust Z-Score Outliers:', paste(which(as.vector(pcr.data[,'Ave Fold'])$`Ave Fold` %in% expr.vec[abs(modZ)>4]), collapse=', ')))
    box.outlier.rows <- c(box.outlier.rows, which(as.vector(pcr.data[,'Ave Fold'])$`Ave Fold` %in% box.out$out))
    z.outlier.rows <- c(z.outlier.rows, which(as.vector(pcr.data[,'Ave Fold'])$`Ave Fold` %in% expr.vec[abs(modZ)>4]))
    
    for (test in c('10')) {
        res1 <- grubbs.test(expr.vec, type=test)
        res2 <- grubbs.test(expr.vec, type=test, opposite=T)
        pvals <- c(res1$p.value)##, res2$p.value)
        names(pvals) <- c(paste(gene, test, sep='.'))## ,
                          ## paste(gene, stringi::stri_reverse(test), sep='.'))
        grubbs.pvals <- c(grubbs.pvals, pvals)
    }
}

grubbs.pvals[grubbs.pvals<0.05]

padj <- p.adjust(grubbs.pvals, 'holm')

padj[padj < 0.05]

grubbs.pvals <- c()
box.outlier.rows <- c()
z.outlier.rows <- c()
par(mfrow=c(4, 2))
for (gene in genes) {
    print(gene)
    expr.vec <- as.vector(pcr.data[pcr.data$Primer==gene & pcr.data$Groups=='Cluster','Ave Fold'])$`Ave Fold`
    modZ <- (expr.vec-median(expr.vec))/mad(expr.vec)
    box.out <- boxplot(expr.vec, main=gene)
    print(paste('IQR Outliers:', paste(which(as.vector(pcr.data[,'Ave Fold'])$`Ave Fold` %in% box.out$out), collapse=', ')))
    print(paste('Robust Z-Score Outliers:', paste(which(as.vector(pcr.data[,'Ave Fold'])$`Ave Fold` %in% expr.vec[abs(modZ)>4]), collapse=', ')))
    box.outlier.rows <- c(box.outlier.rows, which(as.vector(pcr.data[,'Ave Fold'])$`Ave Fold` %in% box.out$out))
    z.outlier.rows <- c(z.outlier.rows, which(as.vector(pcr.data[,'Ave Fold'])$`Ave Fold` %in% expr.vec[abs(modZ)>4]))
    
    for (test in c('10')) {
        res1 <- grubbs.test(expr.vec, type=test)
        res2 <- grubbs.test(expr.vec, type=test, opposite=T)
        pvals <- c(res1$p.value)##, res2$p.value)
        names(pvals) <- c(paste(gene, test, sep='.'))## ,
                          ## paste(gene, stringi::stri_reverse(test), sep='.'))
        grubbs.pvals <- c(grubbs.pvals, pvals)
    }
}

grubbs.pvals[grubbs.pvals<0.05]

padj <- p.adjust(grubbs.pvals, 'holm')

padj[padj < 0.05]

#### Only one outlier: Maximum DUSP9 value, index 63

pcr.data <- pcr.data[c(-63),]

rna <- read.csv('data/combat_vst_rna_expression.csv', row.names=1)

clin <- read.csv('doppler_clinical.csv', row.names=1)
colnames(clin) <- make.names(colnames(clin))
rownames(clin) <- make.names(clin$Study_ID.ID)
head(clin)

clin <- clin %>% mutate_if(function(x) mean(is.na(as.numeric(x))) < 0.5, as.numeric)

clin$Condition. <- ifelse(clin$Condition. %in% c('PTD','Severe PE','FGR','Control','Control PTD'), clin$Condition., 'FGR+HDP')

## clin$Condition. <- make.names(clin$Condition.)

clin$PrePregBMI <- clin$PrePregWt_Kg / clin$Height_Meters.HeightMeters^2

clin$Smoking <- ifelse(str_detect(clin$Smoke, 'Yes'), 'Yes', NA)
clin$Smoking <- ifelse(str_detect(clin$Smoke, 'Never'), 'No', clin$Smoking)

clin$Race[!clin$Race %in% c('W', 'B')] <- 'other'
clin$FDELTYPE[clin$FDELTYPE==0] <- NA

clin <- clin[!is.na(clin$hvar.Clusters),]

samp.ids <- rownames(clin)
samp.ids <- make.names(unique(pcr.data$Study_ID))

rna.markers <- melt(data.frame(
    rna[samp.ids,genes],
    cluster=clin[samp.ids,'hvar.Clusters'],
    Study_ID=clin[samp.ids, 'Study_ID.ID']),
    value.name='Placenta RNA Expression',
    variable.name='Primer', id.vars=c('cluster', 'Study_ID'))

head(rna.markers)

rownames(rna.markers) <- make.names(paste(rna.markers$Study_ID, rna.markers$Primer))

rna.markers[make.names(paste(pcr.data$Study_ID, pcr.data$Primer)),'Plasma PCR Expression'] <- pcr.data$`Ave Fold`

## rna.markers <- na.omit(rna.markers)

rna.markers$Groups <- ifelse(rna.markers$cluster==marker.map[rna.markers$Primer],
                             'Cluster', 'Rest')

## rna.markers$Groups <- pcr.data$Groups

rna.markers$cluster <- factor(rna.markers$cluster)



for (gene in genes) {
    print(gene)
    print(summary(MASS::rlm(x=rna.markers[rna.markers$Primer==gene,'Placenta RNA Expression'], y=rna.markers[rna.markers$Primer==gene,'Plasma PCR Expression'], method='MM')))
}

ggplot(rna.markers,
       aes(x=`Placenta RNA Expression`, y=`Plasma PCR Expression`, color=Groups)) +
    geom_point() +
    geom_smooth(method=MASS::rlm) +
    stat_cor() +
    facet_wrap(~factor(Primer, levels=genes), scales='free', ncol=2) +
    theme_classic()

## rna.markers$Groups <- ifelse(rna.markers$cluster==marker.map[rna.markers$Primer],
##                              'Cluster', 'Rest')


ggplot(rna.markers, aes(x=Groups, y=`Plasma PCR Expression`, fill=Groups)) +
    geom_boxplot() +
    stat_compare_means(method='wilcox') +
    facet_wrap(~Primer, scales='free', ncol=2) +
    theme_classic()

ggplot(rna.markers, aes(x=Groups, y=`Placenta RNA Expression`, fill=Groups)) +
    geom_boxplot() +
    stat_compare_means(method='wilcox') +
    facet_wrap(~Primer, scales='free', ncol=2) +
    theme_classic()

library(gridExtra)

## Cotinine and hydroxycotinine

pdf('plasma_placenta_marker_expression_matched_samples.pdf', width=6, height=4)
for (gene in genes) {
    single.marker <- rna.markers[rna.markers$Primer==gene,]
    clust <- single.marker$cluster[single.marker$Groups=='Cluster'][1]

    single.marker$Cluster <- factor(single.marker$Groups,
                                   levels=c('Cluster', 'Rest'),
                                   labels=c(paste('Cluster', clust), 'Other Clusters'))

    single.marker <- melt(single.marker, value.name='Expression',
                          id.vars=c('Primer', 'Groups', 'cluster', 'Study_ID', 'Cluster'),
                          variable.name='Source')

    single.marker <- na.omit(single.marker)

    panels <- unique(single.marker$Source)

    plot.list <- lapply(panels,
                        function (x) {
                            ggplot(single.marker %>% filter(Source==x),
                                   aes(x=Cluster, y=Expression, fill=Cluster)) +
                                geom_boxplot(width=0.5) +
                                facet_wrap(~Source, scales='free') +
                                theme_classic() + 
                                { if(x=='Placenta RNA Expression') scale_y_continuous(position = "left") else scale_y_continuous(position = "right") } +
                                theme(axis.title.x = element_blank(), legend.position='none')
                        })

    print(do.call("grid.arrange", c(plot.list, ncol = 2, bottom='Cluster', top=gene)))
}
dev.off()

gg <- ggplot(single.marker, aes(x=Groups, y=Expression, fill=Groups)) +
    geom_boxplot(width=0.5) +
    facet_wrap(~Source, ncol=2, scales='free') +
    scale_y_continuous(position = "right") +
    theme_classic() +
    ggtitle(gene)

print(gg)


par(mfrow=c(4,2))
for (gene in genes) {
    plot(rna[samp.ids,gene], pcr.cast[samp.ids,gene],
         main=gene, xlab='RNA-seq', ylab='PCR')
    print(gene)
    print(cor.test(rna[samp.ids,gene], pcr.cast[samp.ids,gene]))
}

for (gene in genes) {
    print(gene)
    print(
        MASS::cov.rob(rna.markers[rna.markers$Primer==gene,
                                  c('Placenta RNA Expression','Plasma PCR Expression')],
                  cor=T, method='mve', nsamp='best')$cor)
}



pcr.cast <- dcast(pcr.data, Study_ID ~ Primer, value.var='Ave Fold')
rownames(pcr.cast) <- make.names(pcr.cast$Study_ID)
pcr.cast <- pcr.cast[,-1]
pcr.cast
pcr.cast$cluster <- factor(clin[rownames(pcr.cast), 'hvar.Clusters'])

pcr.cast[pcr.cast > 100] <- median(pcr.cast[,'DUSP9'])

## pcr.cast <- na.omit(pcr.cast)

dim(pcr.cast)

library(nnet)

res <- multinom(cluster ~ ., pcr.cast %>% mutate_if(is.numeric, scale), decay=2)
res

library(lmtest)

waldtest(glm((cluster==1) ~ 1, pcr.cast, family=binomial()),
         glm((cluster==1) ~ ALPP + PAPPA, pcr.cast, family=binomial()))

waldtest(glm((cluster==2) ~ 1, na.omit(pcr.cast), family=binomial()),
         glm((cluster==2) ~ LGR5 + DUSP9, pcr.cast, family=binomial()))

waldtest(glm((cluster==3) ~ 1, pcr.cast, family=binomial()),
         glm((cluster==3) ~ HTRA4 + FLT1, pcr.cast, family=binomial()))

waldtest(glm((cluster==4) ~ 1, pcr.cast, family=binomial()),
         glm((cluster==4) ~ LYVE1 + EDNRB, pcr.cast, family=binomial()))


lrtest(glm((cluster==1) ~ 1, pcr.cast, family=binomial()),
       glm((cluster==1) ~ ALPP + PAPPA, pcr.cast, family=binomial()))

lrtest(glm((cluster==2) ~ 1, na.omit(pcr.cast), family=binomial()),
       glm((cluster==2) ~ LGR5 + DUSP9, pcr.cast, family=binomial()))

lrtest(glm((cluster==3) ~ 1, pcr.cast, family=binomial()),
       glm((cluster==3) ~ HTRA4 + FLT1, pcr.cast, family=binomial()))

lrtest(glm((cluster==4) ~ 1, pcr.cast, family=binomial()),
       glm((cluster==4) ~ LYVE1 + EDNRB, pcr.cast, family=binomial()))

