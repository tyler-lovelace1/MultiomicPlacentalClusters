library(tidyverse)
library(readxl)
source('cv_utils.R')

chisq.posthoc.test <- function(count.table, method='holm') {
    require(gtools)
    
    if( dim(count.table)[1] > 1 ) {
        combMat1 <- combinations(dim(count.table)[1], 2)
    } else {
        combMat1 <- combinations(dim(count.table)[1], 1)
    }
    
    combMat2 <- combinations(dim(count.table)[2], 2)

    res <- data.frame(matrix(NA, nrow(combMat1)*nrow(combMat2), 5))
    colnames(res) <- c(paste(names(dimnames(count.table)), 'Comparison', sep='.'),
                       'chisq', 'pval', 'padj')

    idx <- 1
    for (comb1 in 1:nrow(combMat1)) {
        for (comb2 in 1:nrow(combMat2)) {
            if (sum(count.table[combMat1[comb1,], combMat2[comb2,]]) > 0) {
                temp.res <- chisq.test(count.table[combMat1[comb1,], combMat2[comb2,]])
                res[idx,1:4] <- c(paste(dimnames(count.table)[[1]][combMat1[comb1,]],
                                        collapse='-'),
                                  paste(dimnames(count.table)[[2]][combMat2[comb2,]],
                                        collapse='/'),
                                  temp.res$statistic, temp.res$p.value)
            } else {
                res[idx,1:4] <- c(paste(dimnames(count.table)[[1]][combMat1[comb1,]],
                                        collapse='-'),
                                  paste(dimnames(count.table)[[2]][combMat2[comb2,]],
                                        collapse='/'),
                                  NA, NA)
            }
            ## print(temp.res)
            idx <- idx+1
        }
    }

    res$padj <- p.adjust(res$pval, method)
    res
}


fisher.posthoc.test <- function(count.table, method='holm') {
    require(gtools)
    
    if( dim(count.table)[1] > 1 ) {
        combMat1 <- combinations(dim(count.table)[1], 2)
    } else {
        combMat1 <- combinations(dim(count.table)[1], 1)
    }
    
    combMat2 <- combinations(dim(count.table)[2], 2)

    res <- data.frame(matrix(NA, nrow(combMat1)*nrow(combMat2), 5))
    colnames(res) <- c(paste(names(dimnames(count.table)), 'Comparison', sep='.'),
                       'odds.ratio', 'pval', 'padj')

    idx <- 1
    for (comb1 in 1:nrow(combMat1)) {
        for (comb2 in 1:nrow(combMat2)) {
            if (sum(count.table[combMat1[comb1,], combMat2[comb2,]]) > 0) {
                temp.res <- fisher.test(count.table[combMat1[comb1,], combMat2[comb2,]])
                res[idx,1:4] <- c(paste(dimnames(count.table)[[1]][combMat1[comb1,]],
                                        collapse='-'),
                                  paste(dimnames(count.table)[[2]][combMat2[comb2,]],
                                        collapse='/'),
                                  temp.res$estimate, temp.res$p.value)
            } else {
                res[idx,1:4] <- c(paste(dimnames(count.table)[[1]][combMat1[comb1,]],
                                        collapse='-'),
                                  paste(dimnames(count.table)[[2]][combMat2[comb2,]],
                                        collapse='/'),
                                  NA, NA)
            }
            ## print(temp.res)
            idx <- idx+1
        }
    }

    res$padj <- p.adjust(res$pval, method)
    res
}


## clin <- read_excel('Primary clinical variables doppler.xlsx', 1) %>% as.data.frame
## colnames(clin) <- make.names(colnames(clin))
## rownames(clin) <- make.names(clin$Study_ID.ID)
## clin <- clin[-1,,drop=T]
## head(clin)

clin <- read.csv('doppler_clinical.csv', row.names=1)
colnames(clin) <- make.names(colnames(clin))
rownames(clin) <- make.names(clin$Study_ID.ID)
## clin <- clin[,-1,drop=T]
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

pathology.diag <- read_excel('Path diagnoses for analysis.xlsx' , 1) %>% as.data.frame
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


pathology <- read_excel('path features (from reports) for analysis.xlsx' , 1, na='n/a') %>% as.data.frame
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

pathology.slides <- read_excel('path features (from slides) for analysis.xlsx' , 1, na='n/a') %>% as.data.frame
colnames(pathology.slides)
rownames(pathology.slides) <- make.names(pathology.slides$`Study ID`)
pathology.slides <- pathology.slides[,4:6]

clin[intersect(rownames(clin), rownames(pathology.slides)),colnames(pathology.slides)] <- pathology.slides[intersect(rownames(clin),rownames(pathology.slides)),]

head(clin)

## %>% filter(Condition. %in% c('Severe PE', 'FGR+HDP'))

cond.clust.table <- data.frame()

cond.clust.posthoc <- data.frame()

clin$hvar.Clusters <- factor(clin$hvar.Clusters, levels=1:4)

for (cond in unique(clin$Condition)) {
    print(cond)
    res <- chisq.test(table(clin[clin$Condition.==cond,c('Condition.', 'hvar.Clusters')]))
    print(res)

    cond.clust.table[cond,c('Condition', 'chisq.stat', 'pval', 'padj')] <-
        c(cond, res$statistic, res$p.value, NA)

    cond.clust.posthoc <- rbind(cond.clust.posthoc,
                                chisq.posthoc.test(table(clin[clin$Condition.==cond,c('Condition.', 'hvar.Clusters')])))
}

cond.clust.table

cond.clust.table$padj <- p.adjust(cond.clust.table$pval, 'holm')

cond.clust.table

cond.clust.posthoc[cond.clust.posthoc$padj < 0.05,]

write.csv(cond.clust.table, 'conditions_v_clusters_table.csv', quote=F)

write.csv(cond.clust.posthoc, 'conditions_v_clusters_chisq_subset_posthoc.csv', quote=F)

#### DE cluster summaries

clin.cont <- clin %>% select(Clusters, MaternalAge.Age, PrePregBMI, WksGest, Birthweight)

clin.cont.summary <- clin.cont %>%
    group_by(Clusters) %>%
    summarize_all(c(mean=mean, sd=sd)) %>%
    as.data.frame

clin.cont.summary <- clin.cont.summary[,-1]

clin.cont.summary <- clin.cont.summary[,c(seq(1, ncol(clin.cont.summary), by=4),
                                          seq(2, ncol(clin.cont.summary), by=4),
                                          seq(3, ncol(clin.cont.summary), by=4),
                                          seq(4, ncol(clin.cont.summary), by=4))] ## ,
                                          ## seq(5, ncol(clin.cont.summary), by=6),
                                          ## seq(6, ncol(clin.cont.summary), by=6))]
clin.cont.summary

clin.cont.table <- data.frame(matrix(NA, 4, (ncol(clin.cont)-1)))

colnames(clin.cont.table) <- colnames(clin.cont)[-1]
rownames(clin.cont.table) <- c(paste('Cluster ', 1:3, ' (n = ',
                                     as.vector(table(clin.cont$Clusters)),
                                     ')', sep=''), 'pval')

for (var in colnames(clin.cont)[-1]) {
    print(var)
    res <- kruskal.test(clin.cont[,var], clin.cont[,'Clusters'])
    print(res)

    clin.cont.table[1:3,var] <- paste(signif(clin.cont.summary[,paste(var, 'mean', sep='_')], 3),
                                   ' (',
                                   signif(clin.cont.summary[,paste(var, 'sd', sep='_')], 3),
                                   ')', sep='')
    clin.cont.table['pval',var] <- res$p.value
}

clin.cont.table

clin.cat <- clin %>% select(Clusters, Condition., Early.PE, Early.PTD, Smoking, Race,
                            Labor.initiation, FDELTYPE, InfSex, MVM, FVM, AI, CI,
                            `placental weight percentile`, `multipal infarcts`, DV,
                            `Syncytial knots - (Tony: considered with accelerated villous maturation)`,
                            RPH, `Accelerated villous maturity`,
                            `Distal villous hypoplasia focal/diffuse`,
                            `Segmental Avascular villi`) %>% mutate_all(factor)


clin.cat.table <- data.frame()

cat.table.rows <- c(paste('Cluster ', 1:3, ' (n = ',
                          as.vector(table(clin.cat$Clusters)),
                          ')', sep=''), 'pval')

for (var in colnames(clin.cat)[-1]) {

    print(var)

    count.table <- table(clin.cat[,c('Clusters', var)])

    clust.counts <- rowSums(count.table)

    print(count.table)
    
    if (ncol(count.table) >= 4) {
        res <- chisq.test(count.table)
    } else {
        res <- fisher.test(count.table)
    }
    print(res)

    if (ncol(count.table) == 2) {
        clin.cat.table[cat.table.rows,var] <- c(paste(signif(100 * count.table[,2] / clust.counts, 3), ' (', count.table[,2], '/', clust.counts, ')', sep=''), res$p.value)
    } else {
        clin.cat.table[cat.table.rows,var] <- c('','','',res$p.value)
        vals <- levels(clin.cat[,var])
        print(var)
        print(vals)
        for (i in 1:ncol(count.table)) {
            clin.cat.table[cat.table.rows, paste(var, vals[i], sep=': ')] <- c(paste(signif(100 * count.table[,i] / clust.counts, 3), ' (', count.table[,i], '/', clust.counts, ')', sep=''), NA)
        }
    }
}

t(clin.cat.table)

clin.table <- as.data.frame(t(cbind(clin.cont.table, clin.cat.table)))

clin.table$pval <- as.numeric(clin.table$pval)

clin.table[!is.na(clin.table$pval),
           'adj.pval'] <- p.adjust(clin.table$pval[!is.na(clin.table$pval)], method='fdr')

clin.table[clin.table$adj.pval < 0.05,]

clin.table

write.csv(clin.table, 'clinical_de_cluster_table.csv', quote=F)


#### Highly variable cluster summaries

clin.cont <- clin %>% select(hvar.Clusters, MaternalAge.Age, PrePregBMI, WksGest, Birthweight) %>% filter(!is.na(hvar.Clusters))

clin.cont.summary <- clin.cont %>%
    group_by(hvar.Clusters) %>%
    summarize_all(c(mean=mean, sd=sd)) %>%
    as.data.frame

clin.cont.summary <- clin.cont.summary[,-1]

clin.cont.summary <- clin.cont.summary[,c(seq(1, ncol(clin.cont.summary), by=4),
                                          seq(2, ncol(clin.cont.summary), by=4),
                                          seq(3, ncol(clin.cont.summary), by=4),
                                          seq(4, ncol(clin.cont.summary), by=4))] ## ,
                                          ## seq(5, ncol(clin.cont.summary), by=6),
                                          ## seq(6, ncol(clin.cont.summary), by=6))]
clin.cont.summary

library(FSA)

clin.cont.table <- data.frame(matrix(NA, 5, (ncol(clin.cont)-1)))

clin.cont.posthoc <- data.frame()

colnames(clin.cont.table) <- colnames(clin.cont)[-1]
rownames(clin.cont.table) <- c(paste('Cluster ', 1:4, ' (n = ',
                                     as.vector(table(clin.cont$hvar.Clusters)),
                                     ')', sep=''), 'pval')

for (var in colnames(clin.cont)[-1]) {
    print(var)
    res <- kruskal.test(clin.cont[,var], clin.cont[,'hvar.Clusters'])
    print(res)

    clin.cont.table[1:4,var] <- paste(signif(clin.cont.summary[,paste(var, 'mean', sep='_')], 3),
                                   ' (',
                                   signif(clin.cont.summary[,paste(var, 'sd', sep='_')], 3),
                                   ')', sep='')
    clin.cont.table['pval',var] <- res$p.value

    if (res$p.value < 0.05) {
        res.posthoc <- dunnTest(clin.cont[,var],
                                factor(clin.cont[,'hvar.Clusters']),
                                method='holm')
        print(res.posthoc)

        clin.cont.posthoc <- rbind(clin.cont.posthoc,
                                   cbind(Feature=rep(var, nrow(res.posthoc$res)),
                                         res.posthoc$res))
    }
}

clin.cont.table

clin.cont.posthoc

clin.cat <- clin %>% select(hvar.Clusters, Condition., Early.PE, Early.PTD, Smoking, Race,
                            Labor.no.labor, Labor.initiation, FDELTYPE, InfSex,
                            MVM, FVM, AI, CI, Diabetes, Nulliparity, AEDF.REDF, NICU.admission,
                            `placental weight percentile`, `multipal infarcts`, DV,
                            `Syncytial knots - (Tony: considered with accelerated villous maturation)`,
                            RPH, `Accelerated villous maturity`,
                            `Distal villous hypoplasia focal/diffuse`,
                            `Segmental Avascular villi`) %>%
    filter(!is.na(hvar.Clusters)) %>%
    mutate_all(factor)


clin.cat.table <- data.frame()
clin.cat.posthoc <- data.frame()

cat.table.rows <- c(paste('Cluster ', 1:4, ' (n = ',
                          as.vector(table(clin.cat$hvar.Clusters)),
                          ')', sep=''), 'pval')

for (var in colnames(clin.cat)[-1]) {

    print(var)

    count.table <- table(clin.cat[,c('hvar.Clusters', var)])

    clust.counts <- rowSums(count.table)

    print(count.table)

    res <- chisq.test(count.table)
    ## if (ncol(count.table) >= 4) {
    ##     res <- chisq.test(count.table)
    ## } else {
    ##     res <- fisher.test(count.table)
    ## }
    print(res)

    if (res$p.value < 0.05) {
        posthoc.res <- fisher.posthoc.test(count.table, method='holm')
        print(posthoc.res)
        ## posthoc.res <- rbind(colnames(posthoc.res), posthoc.res)
        colnames(posthoc.res) <- c('Feat1.Comparison', 'Feat2.Comparison',
                                   'odds.ratio', 'pval', 'padj')
        ## if (nrow(clin.cat.posthoc) > 0) {
        ##     posthoc.res <- posthoc.res[,-1:-2]
        ## }
        ## clin.cat.posthoc <- rbind(clin.cat.posthoc, cbind(Feature=rep(var,ncol(posthoc.res)), t(posthoc.res)))
        clin.cat.posthoc <- rbind(clin.cat.posthoc,
                                  cbind(Feature1=rep('hvar.Clusters',nrow(posthoc.res)),
                                        Feature2=rep(var,nrow(posthoc.res)),
                                        posthoc.res))
    }
    ## print(fisher.multcomp(t(count.table), p.method='holm'))

    if (ncol(count.table) == 2) {
        clin.cat.table[cat.table.rows,var] <- c(paste(signif(100 * count.table[,2] / clust.counts, 3), ' (', count.table[,2], '/', clust.counts, ')', sep=''), res$p.value)
    } else {
        clin.cat.table[cat.table.rows,var] <- c('','','','',res$p.value)
        vals <- levels(clin.cat[,var])
        print(var)
        print(vals)
        for (i in 1:ncol(count.table)) {
            clin.cat.table[cat.table.rows, paste(var, vals[i], sep=': ')] <- c(paste(signif(100 * count.table[,i] / clust.counts, 3), ' (', count.table[,i], '/', clust.counts, ')', sep=''), NA)
        }
    }
}

t(clin.cat.table)

clin.cat.posthoc

clin.table <- as.data.frame(t(cbind(clin.cont.table, clin.cat.table)))

clin.table$pval <- as.numeric(clin.table$pval)

clin.table[!is.na(clin.table$pval),
           'adj.pval'] <- p.adjust(clin.table$pval[!is.na(clin.table$pval)], method='fdr')

clin.table[clin.table$adj.pval < 0.05,]

sig.vars <- rownames(clin.table[clin.table$adj.pval < 0.05,])[str_detect(rownames(clin.table[clin.table$adj.pval < 0.05,]), 'NA[.]', negate=T)]

sig.vars

clin.table

write.csv(clin.table, 'clinical_hvar_cluster_table.csv', quote=F)

write.csv(clin.cont.posthoc, 'clinical_hvar_cluster_dunn_posthoc.csv', quote=F, row.names=F)

write.table(clin.cat.posthoc, 'clinical_hvar_cluster_fisher_subset_posthoc.csv', sep=',', quote=F, row.names=F)


#### Any PE Highly variable cluster summaries

clin.cont <- clin %>%
    filter(Condition. %in% c('Severe PE', 'FGR+HDP')) %>%
    select(hvar.Clusters, MaternalAge.Age, PrePregBMI, WksGest, Birthweight) %>%
    filter(!is.na(hvar.Clusters))

clin.cont.summary <- clin.cont %>%
    group_by(hvar.Clusters) %>%
    summarize_all(c(mean=mean, sd=sd)) %>%
    as.data.frame

clin.cont.summary <- clin.cont.summary[,-1]

clin.cont.summary <- clin.cont.summary[,c(seq(1, ncol(clin.cont.summary), by=4),
                                          seq(2, ncol(clin.cont.summary), by=4),
                                          seq(3, ncol(clin.cont.summary), by=4),
                                          seq(4, ncol(clin.cont.summary), by=4))] ## ,
                                          ## seq(5, ncol(clin.cont.summary), by=6),
                                          ## seq(6, ncol(clin.cont.summary), by=6))]
clin.cont.summary

clin.cont.posthoc <- data.frame()

clin.cont.table <- data.frame(matrix(NA, 5, (ncol(clin.cont)-1)))

colnames(clin.cont.table) <- colnames(clin.cont)[-1]
rownames(clin.cont.table) <- c(paste('Cluster ', 1:4, ' (n = ',
                                     as.vector(table(clin.cont$hvar.Clusters)),
                                     ')', sep=''), 'pval')

for (var in colnames(clin.cont)[-1]) {
    print(var)
    res <- kruskal.test(clin.cont[,var], clin.cont[,'hvar.Clusters'])
    print(res)

    clin.cont.table[1:4,var] <- paste(signif(clin.cont.summary[,paste(var, 'mean', sep='_')], 3),
                                   ' (',
                                   signif(clin.cont.summary[,paste(var, 'sd', sep='_')], 3),
                                   ')', sep='')
    clin.cont.table['pval',var] <- res$p.value

    if (res$p.value < 0.05) {
        res.posthoc <- dunnTest(clin.cont[,var],
                                factor(clin.cont[,'hvar.Clusters']),
                                method='holm')
        print(res.posthoc)

        clin.cont.posthoc <- rbind(clin.cont.posthoc,
                                   cbind(Feature=rep(var, nrow(res.posthoc$res)),
                                         res.posthoc$res))
    }

}

clin.cont.table

clin.cont.posthoc

clin.cat <- clin %>%
    filter(Condition. %in% c('Severe PE', 'FGR+HDP')) %>%
    select(hvar.Clusters, Condition., Early.PE, Early.PTD, Smoking, Race,
           Labor.no.labor, Labor.initiation, FDELTYPE, InfSex,
           MVM, FVM, AI, CI, Diabetes, Nulliparity, AEDF.REDF, NICU.admission,
           `placental weight percentile`, `multipal infarcts`, DV,
           `Syncytial knots - (Tony: considered with accelerated villous maturation)`,
           RPH, `Accelerated villous maturity`,
           `Distal villous hypoplasia focal/diffuse`,
           `Segmental Avascular villi`) %>%
    filter(!is.na(hvar.Clusters)) %>%
    mutate_all(factor)


clin.cat.table <- data.frame()
clin.cat.posthoc <- data.frame()

cat.table.rows <- c(paste('Cluster ', 1:4, ' (n = ',
                          as.vector(table(clin.cat$hvar.Clusters)),
                          ')', sep=''), 'pval')

for (var in colnames(clin.cat)[-1]) {

    print(var)

    count.table <- table(clin.cat[,c('hvar.Clusters', var)])

    clust.counts <- rowSums(count.table)

    print(count.table)

    res <- chisq.test(count.table)
    ## print(res)
    ## if (ncol(count.table) >= 4) {
    ##     res <- chisq.test(count.table)
    ## } else {
    ##     res <- fisher.test(count.table)
    ## }
    print(res)

    if (res$p.value < 0.05) {
        posthoc.res <- fisher.posthoc.test(count.table, method='holm')
        print(posthoc.res)
        ## posthoc.res <- rbind(colnames(posthoc.res), posthoc.res)
        colnames(posthoc.res) <- c('Feat1.Comparison', 'Feat2.Comparison',
                                   'odds.ratio', 'pval', 'padj')
        ## if (nrow(clin.cat.posthoc) > 0) {
        ##     posthoc.res <- posthoc.res[,-1:-2]
        ## }
        ## clin.cat.posthoc <- rbind(clin.cat.posthoc, cbind(Feature=rep(var,ncol(posthoc.res)), t(posthoc.res)))
        clin.cat.posthoc <- rbind(clin.cat.posthoc,
                                  cbind(Feature1=rep('hvar.Clusters',nrow(posthoc.res)),
                                        Feature2=rep(var,nrow(posthoc.res)),
                                        posthoc.res))
    }
    

    if (ncol(count.table) == 2) {
        clin.cat.table[cat.table.rows,var] <- c(paste(signif(100 * count.table[,2] / clust.counts, 3), ' (', count.table[,2], '/', clust.counts, ')', sep=''), res$p.value)
    } else {
        clin.cat.table[cat.table.rows,var] <- c('','','','',res$p.value)
        vals <- levels(clin.cat[,var])
        print(var)
        print(vals)
        for (i in 1:ncol(count.table)) {
            clin.cat.table[cat.table.rows, paste(var, vals[i], sep=': ')] <- c(paste(signif(100 * count.table[,i] / clust.counts, 3), ' (', count.table[,i], '/', clust.counts, ')', sep=''), NA)
        }
    }
}

t(clin.cat.table)

clin.table <- as.data.frame(t(cbind(clin.cont.table, clin.cat.table)))

clin.table$pval <- as.numeric(clin.table$pval)

clin.table[!is.na(clin.table$pval),
           'adj.pval'] <- p.adjust(clin.table$pval[!is.na(clin.table$pval)], method='fdr')

clin.table[clin.table$adj.pval < 0.05,]

clin.table

write.csv(clin.table, 'AnyPE_clinical_hvar_cluster_table.csv', quote=F)

write.csv(clin.cont.posthoc, 'AnyPE_clinical_hvar_cluster_dunn_posthoc.csv', quote=F, row.names=F)

write.table(clin.cat.posthoc, 'AnyPE_clinical_hvar_cluster_fisher_subset_posthoc.csv', sep=',', quote=F, row.names=F)


#### PTD Highly variable cluster summaries

clin.cont <- clin %>%
    filter(Condition. %in% c('PTD')) %>%
    select(hvar.Clusters, MaternalAge.Age, PrePregBMI, WksGest, Birthweight) %>%
    filter(!is.na(hvar.Clusters))

clin.cont.summary <- clin.cont %>%
    group_by(hvar.Clusters) %>%
    summarize_all(c(mean=mean, sd=sd)) %>%
    as.data.frame

clin.cont.summary <- clin.cont.summary[,-1]

clin.cont.summary <- clin.cont.summary[,c(seq(1, ncol(clin.cont.summary), by=4),
                                          seq(2, ncol(clin.cont.summary), by=4),
                                          seq(3, ncol(clin.cont.summary), by=4),
                                          seq(4, ncol(clin.cont.summary), by=4))] ## ,
                                          ## seq(5, ncol(clin.cont.summary), by=6),
                                          ## seq(6, ncol(clin.cont.summary), by=6))]
clin.cont.summary

clin.cont.posthoc <- data.frame()

clin.cont.table <- data.frame(matrix(NA, 5, (ncol(clin.cont)-1)))

colnames(clin.cont.table) <- colnames(clin.cont)[-1]
rownames(clin.cont.table) <- c(paste('Cluster ', 1:4, ' (n = ',
                                     as.vector(table(clin.cont$hvar.Clusters)),
                                     ')', sep=''), 'pval')

for (var in colnames(clin.cont)[-1]) {
    print(var)
    res <- kruskal.test(clin.cont[,var], clin.cont[,'hvar.Clusters'])
    print(res)

    clin.cont.table[1:4,var] <- paste(signif(clin.cont.summary[,paste(var, 'mean', sep='_')], 3),
                                   ' (',
                                   signif(clin.cont.summary[,paste(var, 'sd', sep='_')], 3),
                                   ')', sep='')
    clin.cont.table['pval',var] <- res$p.value

    if (res$p.value < 0.05) {
        res.posthoc <- dunnTest(clin.cont[,var],
                                factor(clin.cont[,'hvar.Clusters']),
                                method='holm')
        print(res.posthoc)

        clin.cont.posthoc <- rbind(clin.cont.posthoc,
                                   cbind(Feature=rep(var, nrow(res.posthoc$res)),
                                         res.posthoc$res))
    }

}

clin.cont.table

clin.cont.posthoc

clin.cat <- clin %>%
    filter(Condition. %in% c('PTD')) %>%
    select(hvar.Clusters, Early.PTD, Smoking, Race,
           Labor.no.labor, FDELTYPE, InfSex,
           MVM, FVM, AI, CI, Diabetes, Nulliparity, NICU.admission,
           `placental weight percentile`, DV,
           `Syncytial knots - (Tony: considered with accelerated villous maturation)`,
           RPH, `Accelerated villous maturity`,
           `Distal villous hypoplasia focal/diffuse`,
           `Segmental Avascular villi`) %>%
    filter(!is.na(hvar.Clusters)) %>%
    mutate_all(factor)


clin.cat.table <- data.frame()
clin.cat.posthoc <- data.frame()

cat.table.rows <- c(paste('Cluster ', 1:4, ' (n = ',
                          as.vector(table(clin.cat$hvar.Clusters)),
                          ')', sep=''), 'pval')

for (var in colnames(clin.cat)[-1]) {

    print(var)

    count.table <- table(clin.cat[,c('hvar.Clusters', var)])

    clust.counts <- rowSums(count.table)

    print(count.table)

    res <- chisq.test(count.table)
    ## print(res)
    ## if (ncol(count.table) >= 4) {
    ##     res <- chisq.test(count.table)
    ## } else {
    ##     res <- fisher.test(count.table)
    ## }
    print(res)

    if (res$p.value < 0.05) {
        posthoc.res <- fisher.posthoc.test(count.table, method='holm')
        print(posthoc.res)
        ## posthoc.res <- rbind(colnames(posthoc.res), posthoc.res)
        colnames(posthoc.res) <- c('Feat1.Comparison', 'Feat2.Comparison',
                                   'odds.ratio', 'pval', 'padj')
        ## if (nrow(clin.cat.posthoc) > 0) {
        ##     posthoc.res <- posthoc.res[,-1:-2]
        ## }
        ## clin.cat.posthoc <- rbind(clin.cat.posthoc, cbind(Feature=rep(var,ncol(posthoc.res)), t(posthoc.res)))
        clin.cat.posthoc <- rbind(clin.cat.posthoc,
                                  cbind(Feature1=rep('hvar.Clusters',nrow(posthoc.res)),
                                        Feature2=rep(var,nrow(posthoc.res)),
                                        posthoc.res))
    }
    

    if (ncol(count.table) == 2) {
        clin.cat.table[cat.table.rows,var] <- c(paste(signif(100 * count.table[,2] / clust.counts, 3), ' (', count.table[,2], '/', clust.counts, ')', sep=''), res$p.value)
    } else {
        clin.cat.table[cat.table.rows,var] <- c('','','','',res$p.value)
        vals <- levels(clin.cat[,var])
        print(var)
        print(vals)
        for (i in 1:ncol(count.table)) {
            clin.cat.table[cat.table.rows, paste(var, vals[i], sep=': ')] <- c(paste(signif(100 * count.table[,i] / clust.counts, 3), ' (', count.table[,i], '/', clust.counts, ')', sep=''), NA)
        }
    }
}

t(clin.cat.table)
clin.cat.posthoc

clin.table <- as.data.frame(t(cbind(clin.cont.table, clin.cat.table)))

clin.table$pval <- as.numeric(clin.table$pval)

clin.table[!is.na(clin.table$pval),
           'adj.pval'] <- p.adjust(clin.table$pval[!is.na(clin.table$pval)], method='fdr')

clin.table[clin.table$adj.pval < 0.05,]

clin.table

write.csv(clin.table, 'PTD_clinical_hvar_cluster_table.csv', quote=F)

write.csv(clin.cont.posthoc, 'PTD_clinical_hvar_cluster_dunn_posthoc.csv', quote=F, row.names=F)

write.table(clin.cat.posthoc, 'PTD_clinical_hvar_cluster_fisher_subset_posthoc.csv', sep=',', quote=F, row.names=F)

#### Condition summaries

clin.cont <- clin %>% select(Condition., MaternalAge.Age, PrePregBMI, WksGest, Birthweight) %>% filter(!is.na(Condition.)) %>% filter(!is.na(MaternalAge.Age) & !is.na(PrePregBMI))

clin.cont.summary <- clin.cont %>%
    group_by(Condition.) %>%
    summarize_all(c(mean=mean, sd=sd)) %>%
    as.data.frame

clin.cont.summary <- clin.cont.summary[,-1]

clin.cont.summary <- clin.cont.summary[,c(seq(1, ncol(clin.cont.summary), by=4),
                                          seq(2, ncol(clin.cont.summary), by=4),
                                          seq(3, ncol(clin.cont.summary), by=4),
                                          seq(4, ncol(clin.cont.summary), by=4))] ## ,
                                          ## seq(5, ncol(clin.cont.summary), by=6),
                                          ## seq(6, ncol(clin.cont.summary), by=6))]
clin.cont.summary

clin.cont.posthoc <- data.frame()

clin.cont.table <- data.frame(matrix(NA, length(unique(clin$Condition.))+1, (ncol(clin.cont)-1)))

colnames(clin.cont.table) <- colnames(clin.cont)[-1]
rownames(clin.cont.table) <- c(paste(sort(unique(clin$Condition.)), ' (n = ',
                                     as.vector(table(clin$Condition.)),
                                     ')', sep=''), 'pval')

for (var in colnames(clin.cont)[-1]) {
    print(var)
    res <- kruskal.test(clin.cont[,var], clin.cont[,'Condition.'])
    print(res)

    clin.cont.table[1:length(unique(clin$Condition.)),var] <- paste(signif(clin.cont.summary[,paste(var, 'mean', sep='_')], 3),
                                   ' (',
                                   signif(clin.cont.summary[,paste(var, 'sd', sep='_')], 3),
                                   ')', sep='')
    clin.cont.table['pval',var] <- res$p.value

    if (res$p.value < 0.05) {
        res.posthoc <- dunnTest(clin.cont[,var],
                                factor(clin.cont[,'Condition.']),
                                method='holm')
        print(res.posthoc)

        clin.cont.posthoc <- rbind(clin.cont.posthoc,
                                   cbind(Feature=rep(var, nrow(res.posthoc$res)),
                                         res.posthoc$res))
    }
}

clin.cont.table

clin.cont.posthoc

clin.cat <- clin %>% select(Condition., Early.PE, Early.PTD, Smoking, Race,
                            Labor.no.labor, Labor.initiation, FDELTYPE, InfSex,
                            MVM, FVM, AI, CI, Diabetes, Nulliparity, AEDF.REDF, NICU.admission,
                            `placental weight percentile`, `multipal infarcts`, DV,
                            `Syncytial knots - (Tony: considered with accelerated villous maturation)`,
                            RPH, `Accelerated villous maturity`,
                            `Distal villous hypoplasia focal/diffuse`,
                            `Segmental Avascular villi`) %>%
    filter(!is.na(Condition.)) %>%
    mutate_all(factor)


clin.cat.posthoc <- data.frame()
clin.cat.table <- data.frame()

cat.table.rows <- c(paste(sort(unique(clin$Condition.)), ' (n = ',
                          as.vector(table(clin.cat$Condition.)),
                          ')', sep=''), 'pval')

for (var in colnames(clin.cat)[-1]) {

    print(var)

    if (var == 'AEDF.REDF') {
        count.table <- table(clin.cat[,c('Condition.', var)] %>%
                             filter(Condition. %in% c('FGR', 'FGR+HDP', 'Severe PE')) %>%
                             mutate_all(factor))
    } else {
        count.table <- table(clin.cat[,c('Condition.', var)])
    }
    
    clust.counts <- rowSums(count.table)

    print(count.table)
    
    res <- chisq.test(count.table)
    ## } else {
    ##     res <- fisher.test(count.table)
    ## }
    print(res)

    count.table <- table(clin.cat[,c('Condition.', var)])

    clust.counts <- rowSums(count.table)

    if (res$p.value < 0.05) {
        posthoc.res <- fisher.posthoc.test(count.table, method='holm')
        print(posthoc.res)
        ## posthoc.res <- rbind(colnames(posthoc.res), posthoc.res)
        colnames(posthoc.res) <- c('Feat1.Comparison', 'Feat2.Comparison',
                                   'odds.ratio', 'pval', 'padj')
        ## if (nrow(clin.cat.posthoc) > 0) {
        ##     posthoc.res <- posthoc.res[,-1:-2]
        ## }
        ## clin.cat.posthoc <- rbind(clin.cat.posthoc, cbind(Feature=rep(var,ncol(posthoc.res)), t(posthoc.res)))
        clin.cat.posthoc <- rbind(clin.cat.posthoc,
                                  cbind(Feature1=rep('Condition.',nrow(posthoc.res)),
                                        Feature2=rep(var,nrow(posthoc.res)),
                                        posthoc.res))    }
    
    if (ncol(count.table) == 2) {
        clin.cat.table[cat.table.rows,var] <- c(paste(signif(100 * count.table[,2] / clust.counts, 3), ' (', count.table[,2], '/', clust.counts, ')', sep=''), res$p.value)
    } else {
        clin.cat.table[cat.table.rows,var] <- c(rep('', length(unique(clin$Condition.))),
                                                res$p.value)
        vals <- levels(clin.cat[,var])
        print(var)
        print(vals)
        for (i in 1:ncol(count.table)) {
            clin.cat.table[cat.table.rows, paste(var, vals[i], sep=': ')] <- c(paste(signif(100 * count.table[,i] / clust.counts, 3), ' (', count.table[,i], '/', clust.counts, ')', sep=''), NA)
        }
    }
}

t(clin.cat.table)

clin.table <- as.data.frame(t(cbind(clin.cont.table, clin.cat.table)))

clin.table$pval <- as.numeric(clin.table$pval)

clin.table[!is.na(clin.table$pval),
           'adj.pval'] <- p.adjust(clin.table$pval[!is.na(clin.table$pval)], method='fdr')

clin.table[clin.table$adj.pval < 0.05,]

clin.table

write.csv(clin.table, 'clinical_condition_table.csv', quote=F)

write.csv(clin.cont.posthoc, 'clinical_condition_dunn_posthoc.csv', quote=F, row.names=F)

write.table(clin.cat.posthoc, 'clinical_condition_fisher_subset_posthoc.csv', sep=',', quote=F, row.names=F)


#### Predicted highly variable clusters summaries

clin.cont <- clin %>% select(pred.hvar.Clusters, MaternalAge.Age, PrePregBMI, WksGest, Birthweight) %>% filter(!is.na(pred.hvar.Clusters))

clin.cont.summary <- clin.cont %>%
    group_by(pred.hvar.Clusters) %>%
    summarize_all(c(mean=mean, sd=sd)) %>%
    as.data.frame

clin.cont.summary <- clin.cont.summary[,-1]

clin.cont.summary <- clin.cont.summary[,c(seq(1, ncol(clin.cont.summary), by=4),
                                          seq(2, ncol(clin.cont.summary), by=4),
                                          seq(3, ncol(clin.cont.summary), by=4),
                                          seq(4, ncol(clin.cont.summary), by=4))] ## ,
                                          ## seq(5, ncol(clin.cont.summary), by=6),
                                          ## seq(6, ncol(clin.cont.summary), by=6))]
clin.cont.summary

library(FSA)

clin.cont.table <- data.frame(matrix(NA, 5, (ncol(clin.cont)-1)))

clin.cont.posthoc <- data.frame()

colnames(clin.cont.table) <- colnames(clin.cont)[-1]
rownames(clin.cont.table) <- c(paste('Cluster ', 1:4, ' (n = ',
                                     as.vector(table(clin.cont$pred.hvar.Clusters)),
                                     ')', sep=''), 'pval')

for (var in colnames(clin.cont)[-1]) {
    print(var)
    res <- kruskal.test(clin.cont[,var], clin.cont[,'pred.hvar.Clusters'])
    print(res)

    clin.cont.table[1:4,var] <- paste(signif(clin.cont.summary[,paste(var, 'mean', sep='_')], 3),
                                   ' (',
                                   signif(clin.cont.summary[,paste(var, 'sd', sep='_')], 3),
                                   ')', sep='')
    clin.cont.table['pval',var] <- res$p.value

    if (res$p.value < 0.05) {
        res.posthoc <- dunnTest(clin.cont[,var],
                                factor(clin.cont[,'pred.hvar.Clusters']),
                                method='holm')
        print(res.posthoc)

        clin.cont.posthoc <- rbind(clin.cont.posthoc,
                                   cbind(Feature=rep(var, nrow(res.posthoc$res)),
                                         res.posthoc$res))
    }
}

clin.cont.table

clin.cont.posthoc

clin.cat <- clin %>% select(pred.hvar.Clusters, Condition., Early.PE, Early.PTD, Smoking, Race,
                            Labor.no.labor, Labor.initiation, FDELTYPE, InfSex,
                            MVM, FVM, AI, CI, Diabetes, Nulliparity, AEDF.REDF, NICU.admission,
                            `placental weight percentile`, `multipal infarcts`, DV,
                            `Syncytial knots - (Tony: considered with accelerated villous maturation)`,
                            RPH, `Accelerated villous maturity`,
                            `Distal villous hypoplasia focal/diffuse`,
                            `Segmental Avascular villi`) %>%
    filter(!is.na(pred.hvar.Clusters)) %>%
    mutate_all(factor)


clin.cat.table <- data.frame()
clin.cat.posthoc <- data.frame()

cat.table.rows <- c(paste('Cluster ', 1:4, ' (n = ',
                          as.vector(table(clin.cat$pred.hvar.Clusters)),
                          ')', sep=''), 'pval')

for (var in colnames(clin.cat)[-1]) {

    print(var)

    count.table <- table(clin.cat[,c('pred.hvar.Clusters', var)])

    clust.counts <- rowSums(count.table)

    print(count.table)

    res <- chisq.test(count.table)
    ## if (ncol(count.table) >= 4) {
    ##     res <- chisq.test(count.table)
    ## } else {
    ##     res <- fisher.test(count.table)
    ## }
    print(res)

    if (res$p.value < 0.05) {
        posthoc.res <- chisq.posthoc.test(count.table, method='holm')
        print(posthoc.res)
        ## posthoc.res <- rbind(colnames(posthoc.res), posthoc.res)
        colnames(posthoc.res) <- c('Feat1.Comparison', 'Feat2.Comparison',
                                   'chisq', 'pval', 'padj')
        ## if (nrow(clin.cat.posthoc) > 0) {
        ##     posthoc.res <- posthoc.res[,-1:-2]
        ## }
        ## clin.cat.posthoc <- rbind(clin.cat.posthoc, cbind(Feature=rep(var,ncol(posthoc.res)), t(posthoc.res)))
        clin.cat.posthoc <- rbind(clin.cat.posthoc,
                                  cbind(Feature1=rep('pred.hvar.Clusters',nrow(posthoc.res)),
                                        Feature2=rep(var,nrow(posthoc.res)),
                                        posthoc.res))
    }
    ## print(fisher.multcomp(t(count.table), p.method='holm'))

    if (ncol(count.table) == 2) {
        clin.cat.table[cat.table.rows,var] <- c(paste(signif(100 * count.table[,2] / clust.counts, 3), ' (', count.table[,2], '/', clust.counts, ')', sep=''), res$p.value)
    } else {
        clin.cat.table[cat.table.rows,var] <- c('','','','',res$p.value)
        vals <- levels(clin.cat[,var])
        print(var)
        print(vals)
        for (i in 1:ncol(count.table)) {
            clin.cat.table[cat.table.rows, paste(var, vals[i], sep=': ')] <- c(paste(signif(100 * count.table[,i] / clust.counts, 3), ' (', count.table[,i], '/', clust.counts, ')', sep=''), NA)
        }
    }
}

t(clin.cat.table)

clin.cat.posthoc

clin.table <- as.data.frame(t(cbind(clin.cont.table, clin.cat.table)))

clin.table$pval <- as.numeric(clin.table$pval)

clin.table[!is.na(clin.table$pval),
           'adj.pval'] <- p.adjust(clin.table$pval[!is.na(clin.table$pval)], method='fdr')

clin.table[clin.table$adj.pval < 0.05,]

sig.vars <- rownames(clin.table[clin.table$adj.pval < 0.05,])[str_detect(rownames(clin.table[clin.table$adj.pval < 0.05,]), 'NA[.]', negate=T)]

sig.vars

clin.table

write.csv(clin.table, 'clinical_pred_hvar_cluster_table.csv', quote=F)

write.csv(clin.cont.posthoc, 'clinical_pred_hvar_cluster_dunn_posthoc.csv', quote=F, row.names=F)

write.table(clin.cat.posthoc, 'clinical_pred_hvar_cluster_chisq_subset_posthoc.csv', sep=',', quote=F, row.names=F)


#### Comparisons of clusters and predicted clusters

clin.table <- read.csv('clinical_hvar_cluster_table.csv', row.names=1)
clin.table

sig.vars <- rownames(clin.table[clin.table$adj.pval < 0.05,])[str_detect(rownames(clin.table[clin.table$adj.pval < 0.05,]), 'NA', negate=T)]

sig.vars

pred.v.orig <- data.frame(Feature=c(),
                          Cluster=c(),
                          p.val=c(),
                          p.adj=c())

for (var in sig.vars) {
    for (clust in 1:4) {
        if (is.numeric(clin[,var])) {
            res <- wilcox.test(na.omit(clin[clin$hvar.Clusters==clust,var]),
                               na.omit(clin[clin$pred.hvar.Clusters==clust,var]))
            print(res)

            pred.v.orig <- rbind(pred.v.orig,
                                 data.frame(Feature=var,
                                            Cluster=clust,
                                            p.val=res$p.value,
                                            p.adj=NA))
        } else {
            count <- table(c(clin[!is.na(clin$hvar.Clusters) &
                                  clin$hvar.Clusters==clust,var],
                             clin[!is.na(clin$pred.hvar.Clusters) &
                                  clin$pred.hvar.Clusters==clust,var]),
                           c(rep('Original', sum(clin$hvar.Clusters==clust, na.rm=T)),
                             rep('Predicted', sum(clin$pred.hvar.Clusters==clust, na.rm=T))))
            res <- chisq.test(count)
            print(res)

            pred.v.orig <- rbind(pred.v.orig,
                                 data.frame(Feature=var,
                                            Cluster=clust,
                                            p.val=res$p.value,
                                            p.adj=NA))
        }
    }
}

pred.v.orig$p.adj <- p.adjust(pred.v.orig$p.val, 'fdr')
pred.v.orig

write.csv(pred.v.orig, 'hvar_cluster_v_pred_hvar_cluster_clinical_comparison.csv', row.names=F)

#### Comparisons of histopathology significance Vuong's test
library(nonnest2)

histopath.comparisons <- data.frame(Feature=c('MVM', 'FVM', 'CI', 'AI'),
                                    z.stat=rep(NA,4),
                                    p.val=rep(NA,4),
                                    p.adj=rep(NA,4))

mod.clust <- glm(MVM ~ factor(hvar.Clusters), family=binomial(), clin[!is.na(clin$hvar.Clusters),])
summary(mod.clust)
mod.cond <- glm(MVM ~ factor(Condition.), family=binomial(), clin[!is.na(clin$hvar.Clusters),])
summary(mod.cond)


res <- vuongtest(mod.clust, mod.cond, adj='aic')

res

histopath.comparisons[1,c('z.stat', 'p.val')] <- c(res$LRTstat, res$p_LRT$A)

mod.clust <- glm(FVM ~ factor(hvar.Clusters), family=binomial(), clin[!is.na(clin$hvar.Clusters),])
summary(mod.clust)
mod.cond <- glm(FVM ~ factor(Condition.), family=binomial(), clin[!is.na(clin$hvar.Clusters),])
summary(mod.cond)

res <- vuongtest(mod.clust, mod.cond, adj='aic')
res

histopath.comparisons[2,c('z.stat', 'p.val')] <- c(res$LRTstat, res$p_LRT$A)


mod.clust <- glm(CI ~ factor(hvar.Clusters), family=binomial(), clin[!is.na(clin$hvar.Clusters),])
summary(mod.clust)
mod.cond <- glm(CI ~ factor(Condition.), family=binomial(), clin[!is.na(clin$hvar.Clusters),])
summary(mod.cond)

res <- vuongtest(mod.clust, mod.cond, adj='aic')
res

histopath.comparisons[3,c('z.stat', 'p.val')] <- c(res$LRTstat, res$p_LRT$A)


mod.clust <- glm(AI ~ factor(hvar.Clusters), family=binomial(), clin[!is.na(clin$hvar.Clusters),])
summary(mod.clust)
mod.cond <- glm(AI ~ factor(Condition.), family=binomial(), clin[!is.na(clin$hvar.Clusters),])
summary(mod.cond)

res <- vuongtest(mod.clust, mod.cond, adj='aic')
res

histopath.comparisons[4,c('z.stat', 'p.val')] <- c(res$LRTstat, res$p_LRT$A)

histopath.comparisons$p.adj <- p.adjust(histopath.comparisons$p.val, 'fdr')

histopath.comparisons

write.csv(histopath.comparisons, 'hitopathology_hvar_clusters_v_condition.csv')


library(clarkeTest)

histopath.comparisons <- data.frame(Feature=c('MVM', 'FVM', 'CI', 'AI'),
                                    stat=rep(NA,4),
                                    p.val=rep(NA,4),
                                    p.adj=rep(NA,4))

mod.clust <- glm(MVM ~ factor(hvar.Clusters), family=binomial(), clin[!is.na(clin$hvar.Clusters),])
summary(mod.clust)
mod.cond <- glm(MVM ~ factor(Condition.), family=binomial(), clin[!is.na(clin$hvar.Clusters),])
summary(mod.cond)

res <- clarke_test(mod.clust, mod.cond, level=0.05, digits=10)
res
histopath.comparisons[1,c('stat', 'p.val')] <- c(res$stat, pbinom(sum(res$loglik1 > res$loglik2), res$nobs, 0.5, lower.tail=F))

mod.clust <- glm(FVM ~ factor(hvar.Clusters), family=binomial(), clin[!is.na(clin$hvar.Clusters),])
summary(mod.clust)
mod.cond <- glm(FVM ~ factor(Condition.), family=binomial(), clin[!is.na(clin$hvar.Clusters),])
summary(mod.cond)

res <- clarke_test(mod.clust, mod.cond, level=0.05, digits=10)
res
histopath.comparisons[2,c('stat', 'p.val')] <- c(res$stat, pbinom(sum(res$loglik1 > res$loglik2), res$nobs, 0.5, lower.tail=F))


mod.clust <- glm(CI ~ factor(hvar.Clusters), family=binomial(), clin[!is.na(clin$hvar.Clusters),])
summary(mod.clust)
mod.cond <- glm(CI ~ factor(Condition.), family=binomial(), clin[!is.na(clin$hvar.Clusters),])
summary(mod.cond)

res <- clarke_test(mod.clust, mod.cond, level=0.05, digits=10)
res
histopath.comparisons[3,c('stat', 'p.val')] <- c(res$stat, pbinom(sum(res$loglik1 > res$loglik2), res$nobs, 0.5, lower.tail=F))


mod.clust <- glm(AI ~ factor(hvar.Clusters), family=binomial(), clin[!is.na(clin$hvar.Clusters),])
summary(mod.clust)
mod.cond <- glm(AI ~ factor(Condition.), family=binomial(), clin[!is.na(clin$hvar.Clusters),])
summary(mod.cond)

res <- clarke_test(mod.clust, mod.cond, level=0.05, digits=10)
res
histopath.comparisons[4,c('stat', 'p.val')] <- c(res$stat, pbinom(sum(res$loglik1 > res$loglik2), res$nobs, 0.5, lower.tail=F))

histopath.comparisons$p.adj <- p.adjust(histopath.comparisons$p.val, 'fdr')

histopath.comparisons



#### Clinical distributions by Clusters for Severe PE only


clin.cont <- clin%>% filter(Condition. %in% c('Severe PE')) %>% select(Clusters, MaternalAge.Age, PrePregBMI, WksGest, Birthweight)

clin.cont.summary <- clin.cont %>%
    group_by(Clusters) %>%
    summarize_all(c(mean=mean, sd=sd)) %>%
    as.data.frame

clin.cont.summary <- clin.cont.summary[,-1]

clin.cont.summary <- clin.cont.summary[,c(seq(1, ncol(clin.cont.summary), by=4),
                                          seq(2, ncol(clin.cont.summary), by=4),
                                          seq(3, ncol(clin.cont.summary), by=4),
                                          seq(4, ncol(clin.cont.summary), by=4))]## ,
                                          ## seq(5, ncol(clin.cont.summary), by=6),
                                          ## seq(6, ncol(clin.cont.summary), by=6))]
clin.cont.summary

clin.cont.table <- data.frame(matrix(NA, 4, (ncol(clin.cont)-1)))

colnames(clin.cont.table) <- colnames(clin.cont)[-1]
rownames(clin.cont.table) <- c(paste('Cluster ', 1:3, ' (n = ',
                                     as.vector(table(clin.cont$Clusters)),
                                     ')', sep=''), 'pval')

for (var in colnames(clin.cont)[-1]) {
    print(var)
    res <- kruskal.test(clin.cont[,var], clin.cont[,'Clusters'])
    print(res)

    clin.cont.table[1:3,var] <- paste(signif(clin.cont.summary[,paste(var, 'mean', sep='_')], 3),
                                   ' (',
                                   signif(clin.cont.summary[,paste(var, 'sd', sep='_')], 3),
                                   ')', sep='')
    clin.cont.table['pval',var] <- res$p.value
}

clin.cont.table

clin.cat <- clin %>%
    filter(Condition. %in% c('Severe PE')) %>%
    select(Clusters, Early.PE, Smoking, Race, Labor.initiation,
           FDELTYPE, InfSex, High.grade.MVM, FVM, Acute.inflammation,
           Chronic.inflammation, `placental weight percentile`,
           `multipal infarcts`, DV, 
           `Syncytial knots - (Tony: considered with accelerated villous maturation)`,
           RPH, `Accelerated villous maturity`,
           `Distal villous hypoplasia focal/diffuse`,
           `Segmental Avascular villi`) %>%
    mutate_all(factor)


clin.cat.table <- data.frame()

cat.table.rows <- c(paste('Cluster ', 1:3, ' (n = ',
                          as.vector(table(clin.cat$Clusters)),
                          ')', sep=''), 'pval')

for (var in colnames(clin.cat)[-1]) {

    print(var)

    count.table <- table(clin.cat[,c('Clusters', var)])

    clust.counts <- rowSums(count.table)

    print(count.table)
    
    if (ncol(count.table) >= 5) {
        res <- chisq.test(count.table)
    } else {
        res <- fisher.test(count.table)
    }
    print(res)

    if (ncol(count.table) == 2) {
        clin.cat.table[cat.table.rows,var] <- c(paste(signif(100 * count.table[,2] / clust.counts, 3), ' (', count.table[,2], '/', clust.counts, ')', sep=''), res$p.value)
    } else {
        clin.cat.table[cat.table.rows,var] <- c('','','',res$p.value)
        vals <- levels(clin.cat[,var])
        print(var)
        print(vals)
        for (i in 1:ncol(count.table)) {
            clin.cat.table[cat.table.rows, paste(var, vals[i], sep=': ')] <- c(paste(signif(100 * count.table[,i] / clust.counts, 3), ' (', count.table[,i], '/', clust.counts, ')', sep=''), NA)
        }
    }
}

clin.cat.table

clin.table <- as.data.frame(t(cbind(clin.cont.table, clin.cat.table)))

clin.table$pval <- as.numeric(clin.table$pval)

clin.table[!is.na(clin.table$pval),
           'adj.pval'] <- p.adjust(clin.table$pval[!is.na(clin.table$pval)], method='fdr')

clin.table[clin.table$adj.pval < 0.05,]

clin.table

write.csv(clin.table, 'SeverePE_clinical_cluster_table.csv', quote=F)


#### Clinical distributions by Clusters for Severe PE only


clin.cont <- clin%>% filter(Condition. %in% c('Severe PE', 'FGR+HDP')) %>% select(Clusters, MaternalAge.Age, PrePregBMI, WksGest, Birthweight)

clin.cont.summary <- clin.cont %>%
    group_by(Clusters) %>%
    summarize_all(c(mean=mean, sd=sd)) %>%
    as.data.frame

clin.cont.summary <- clin.cont.summary[,-1]

clin.cont.summary <- clin.cont.summary[,c(seq(1, ncol(clin.cont.summary), by=4),
                                          seq(2, ncol(clin.cont.summary), by=4),
                                          seq(3, ncol(clin.cont.summary), by=4),
                                          seq(4, ncol(clin.cont.summary), by=4))]## ,
                                          ## seq(5, ncol(clin.cont.summary), by=6),
                                          ## seq(6, ncol(clin.cont.summary), by=6))]
clin.cont.summary

clin.cont.table <- data.frame(matrix(NA, 4, (ncol(clin.cont)-1)))

colnames(clin.cont.table) <- colnames(clin.cont)[-1]
rownames(clin.cont.table) <- c(paste('Cluster ', 1:3, ' (n = ',
                                     as.vector(table(clin.cont$Clusters)),
                                     ')', sep=''), 'pval')

for (var in colnames(clin.cont)[-1]) {
    print(var)
    res <- kruskal.test(clin.cont[,var], clin.cont[,'Clusters'])
    print(res)

    clin.cont.table[1:3,var] <- paste(signif(clin.cont.summary[,paste(var, 'mean', sep='_')], 3),
                                   ' (',
                                   signif(clin.cont.summary[,paste(var, 'sd', sep='_')], 3),
                                   ')', sep='')
    clin.cont.table['pval',var] <- res$p.value
}

clin.cont.table

clin.cat <- clin %>%
    filter(Condition. %in% c('Severe PE', 'FGR+HDP')) %>%
    select(Clusters, Early.PE, Smoking, Race, Labor.initiation,
           FDELTYPE, InfSex, High.grade.MVM, FVM, Acute.inflammation,
           Chronic.inflammation, `placental weight percentile`,
           `multipal infarcts`, DV,
           `Syncytial knots - (Tony: considered with accelerated villous maturation)`,
           RPH, `Accelerated villous maturity`,
           `Distal villous hypoplasia focal/diffuse`,
           `Segmental Avascular villi`) %>%
    mutate_all(factor)


clin.cat.table <- data.frame()

cat.table.rows <- c(paste('Cluster ', 1:3, ' (n = ',
                          as.vector(table(clin.cat$Clusters)),
                          ')', sep=''), 'pval')

for (var in colnames(clin.cat)[-1]) {

    print(var)

    count.table <- table(clin.cat[,c('Clusters', var)])

    clust.counts <- rowSums(count.table)

    print(count.table)
    
    if (ncol(count.table) >= 5) {
        res <- chisq.test(count.table)
    } else {
        res <- fisher.test(count.table)
    }
    print(res)

    if (ncol(count.table) == 2) {
        clin.cat.table[cat.table.rows,var] <- c(paste(signif(100 * count.table[,2] / clust.counts, 3), ' (', count.table[,2], '/', clust.counts, ')', sep=''), res$p.value)
    } else {
        clin.cat.table[cat.table.rows,var] <- c('','','',res$p.value)
        vals <- levels(clin.cat[,var])
        print(var)
        print(vals)
        for (i in 1:ncol(count.table)) {
            clin.cat.table[cat.table.rows, paste(var, vals[i], sep=': ')] <- c(paste(signif(100 * count.table[,i] / clust.counts, 3), ' (', count.table[,i], '/', clust.counts, ')', sep=''), NA)
        }
    }
}

clin.cat.table

clin.table <- as.data.frame(t(cbind(clin.cont.table, clin.cat.table)))

clin.table$pval <- as.numeric(clin.table$pval)

clin.table[!is.na(clin.table$pval),
           'adj.pval'] <- p.adjust(clin.table$pval[!is.na(clin.table$pval)], method='fdr')

clin.table[clin.table$adj.pval < 0.05,]

clin.table

write.csv(clin.table, 'AnyPE_clinical_cluster_table.csv', quote=F)


#### Clinical distributions by Clusters for PTD only


clin.cont <- clin%>% filter(Condition. %in% c('PTD')) %>% select(Clusters, MaternalAge.Age, PrePregBMI, WksGest, Birthweight)

clin.cont.summary <- clin.cont %>%
    group_by(Clusters) %>%
    summarize_all(c(mean=mean, sd=sd)) %>%
    as.data.frame

clin.cont.summary <- clin.cont.summary[,-1]

clin.cont.summary <- clin.cont.summary[,c(seq(1, ncol(clin.cont.summary), by=4),
                                          seq(2, ncol(clin.cont.summary), by=4),
                                          seq(3, ncol(clin.cont.summary), by=4),
                                          seq(4, ncol(clin.cont.summary), by=4))]## ,
                                          ## seq(5, ncol(clin.cont.summary), by=6),
                                          ## seq(6, ncol(clin.cont.summary), by=6))]
clin.cont.summary

clin.cont.table <- data.frame(matrix(NA, 4, (ncol(clin.cont)-1)))

colnames(clin.cont.table) <- colnames(clin.cont)[-1]
rownames(clin.cont.table) <- c(paste('Cluster ', 1:3, ' (n = ',
                                     as.vector(table(clin.cont$Clusters)),
                                     ')', sep=''), 'pval')

for (var in colnames(clin.cont)[-1]) {
    print(var)
    res <- kruskal.test(clin.cont[,var], clin.cont[,'Clusters'])
    print(res)

    clin.cont.table[1:3,var] <- paste(signif(clin.cont.summary[,paste(var, 'mean', sep='_')], 3),
                                   ' (',
                                   signif(clin.cont.summary[,paste(var, 'sd', sep='_')], 3),
                                   ')', sep='')
    clin.cont.table['pval',var] <- res$p.value
}

clin.cont.table

clin.cat <- clin %>%
    filter(Condition. %in% c('PTD')) %>%
    select(Clusters, Smoking, Race,
           FDELTYPE, InfSex, High.grade.MVM, FVM, Acute.inflammation,
           Chronic.inflammation, `placental weight percentile`, DV,
           `Syncytial knots - (Tony: considered with accelerated villous maturation)`,
           RPH, `Accelerated villous maturity`,
           `Distal villous hypoplasia focal/diffuse`,
           `Segmental Avascular villi`) %>%
    mutate_all(factor)


clin.cat.table <- data.frame()

cat.table.rows <- c(paste('Cluster ', 1:3, ' (n = ',
                          as.vector(table(clin.cat$Clusters)),
                          ')', sep=''), 'pval')

for (var in colnames(clin.cat)[-1]) {

    print(var)

    count.table <- table(clin.cat[,c('Clusters', var)])

    clust.counts <- rowSums(count.table)

    print(count.table)
    
    if (ncol(count.table) >= 5) {
        res <- chisq.test(count.table)
    } else {
        res <- fisher.test(count.table)
    }
    print(res)

    if (ncol(count.table) == 2) {
        clin.cat.table[cat.table.rows,var] <- c(paste(signif(100 * count.table[,2] / clust.counts, 3), ' (', count.table[,2], '/', clust.counts, ')', sep=''), res$p.value)
    } else {
        clin.cat.table[cat.table.rows,var] <- c('','','',res$p.value)
        vals <- levels(clin.cat[,var])
        print(var)
        print(vals)
        for (i in 1:ncol(count.table)) {
            clin.cat.table[cat.table.rows, paste(var, vals[i], sep=': ')] <- c(paste(signif(100 * count.table[,i] / clust.counts, 3), ' (', count.table[,i], '/', clust.counts, ')', sep=''), NA)
        }
    }
}

clin.cat.table

clin.table <- as.data.frame(t(cbind(clin.cont.table, clin.cat.table)))

clin.table$pval <- as.numeric(clin.table$pval)

clin.table[!is.na(clin.table$pval),
           'adj.pval'] <- p.adjust(clin.table$pval[!is.na(clin.table$pval)], method='fdr')

clin.table[clin.table$adj.pval < 0.05,]

clin.table

write.csv(clin.table, 'PTD_clinical_cluster_table.csv', quote=F)



#### Distribution of Marker Genes / Proteins

prot <- read.csv('data/combat_prot_data.csv', row.names=1)
colnames(prot) <- paste0('Prot_', colnames(prot))
colnames(prot)

rna <- read.csv('data/combat_vst_rna_expression.csv', row.names=1)
colnames(rna) <- paste0('RNA_', colnames(rna))
colnames(rna)

rna.markers <- paste0('RNA_', c('FLT1', 'PGF', 'VEGFA', 'ENG',
                                'LEP', 'FSTL3', 'IGFBP1', 'IGF1', 'HIF1A'))

prot.markers <- paste0('Prot_', c('FLT1', 'PGF', 'VEGFA', 'LEP', 'FSTL3', 'IGFBP.1'))

samp.ids <- rownames(clin)

#### DE Clusters

markers <- cbind(Clusters=clin$Clusters, rna[samp.ids,] %>% select(rna.markers), prot[samp.ids,] %>% select(prot.markers))

markers$Prot_FLT1overPGF <- prot[samp.ids,'Prot_FLT1'] / prot[samp.ids,'Prot_PGF']

markers$RNA_FLT1overPGF <- rna[samp.ids,'RNA_FLT1'] / rna[samp.ids,'RNA_PGF']

head(markers)


markers.summary <- markers %>%
    group_by(Clusters) %>%
    summarize_all(c(mean=mean, sd=sd)) %>%
    as.data.frame

markers.summary <- markers.summary[,-1]

markers.table <- data.frame(matrix(NA, 4, (ncol(markers)-1)))

colnames(markers.table) <- colnames(markers)[-1]
rownames(markers.table) <- c(paste('Cluster ', 1:3, ' (n = ',
                                     as.vector(table(markers$Clusters)),
                                     ')', sep=''), 'pval')

for (var in colnames(markers)[-1]) {
    print(var)
    res <- kruskal.test(markers[,var], markers[,'Clusters'])
    print(res)

    markers.table[1:3,var] <- paste(signif(markers.summary[,paste(var, 'mean', sep='_')], 3),
                                   ' (',
                                   signif(markers.summary[,paste(var, 'sd', sep='_')], 3),
                                   ')', sep='')
    markers.table['pval',var] <- res$p.value
}

markers.table

markers.table <- as.data.frame(t(markers.table))

markers.table$pval <- as.numeric(markers.table$pval)

markers.table[!is.na(markers.table$pval),
           'adj.pval'] <- p.adjust(markers.table$pval[!is.na(markers.table$pval)], method='fdr')

markers.table[markers.table$adj.pval < 0.05,]

markers.table

write.csv(markers.table, 'markers_de_cluster_table.csv', quote=F)

library(ggplot2)
library(ggpubr)
library(rstatix)

expr.df <- data.frame()

for (var in colnames(markers)[-1]) {

    if (markers.table[var, 'adj.pval'] > 0.05) {
        next
    }

    expr <- data.frame()
    expr[rownames(markers),'Clusters'] <- factor(markers[,c('Clusters')])
    expr$marker <- var
    expr$Expression <- scale(markers[,var])[,1]
    expr$y.position <- 1.05 * max(expr$Expression)
    expr$adj.pval <- markers.table[var,'adj.pval']

    expr.df <- rbind(expr.df, expr)

    ## res <- compare_means(as.formula(paste(var, 'Clusters', sep='~')), markers, method='kruskal.test')
    ## res$p.adj <- markers.table[var,'adj.pval']
    
    ## gg <- ggplot(markers %>% mutate(Clusters=factor(Clusters)),
    ##              aes(x=Clusters, y=markers[,var], fill=Clusters)) +
    ##     geom_violin() +
    ##     geom_boxplot(width=0.2) +
    ##     xlab('Clusters') +
    ##     ylab(var) +
    ##     stat_pvalue_manual(res) + 
    ##     ## geom_text(x=2, y=1.15 * max(markers[,var]), label=paste('Kruskal, size=11) + 
    ##     theme_bw()

    ## print(gg)
}

stat.test <- expr.df %>%
    group_by(marker) %>%
    wilcox_test(Expression ~ Clusters) %>%
    add_xy_position(x = "marker", dodge = 0.8)

stat.test

pdf('signif_markers_de_cluster_violin.pdf', width=15, height=7)

gg <- ggplot(expr.df,
             aes(x=Clusters, y=Expression, fill=Clusters)) +
    geom_violin() +
    geom_boxplot(width=0.2) +
    ## stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.01) + 
    xlab('Clusters') +
    ylab('Z-scored Expression') +
    scale_y_continuous(expand = expansion(mult = c(0.01, 0.1))) + 
    ## ylim(-4,4) + 
    facet_wrap(vars(marker), ncol=8, scales='free_y') +
    theme_bw()

gg

dev.off()


png('signif_markers_de_cluster_violin.png', width=15, height=7, units='in', res=400)

gg <- ggplot(expr.df,
             aes(x=Clusters, y=Expression, fill=Clusters)) +
    geom_violin() +
    geom_boxplot(width=0.2) +
    ## stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.01) + 
    xlab('Clusters') +
    ylab('Z-scored Expression') +
    scale_y_continuous(expand = expansion(mult = c(0.01, 0.1))) + 
    ## ylim(-4,4) + 
    facet_wrap(vars(marker), ncol=8, scales='free_y') +
    theme_bw()

gg

dev.off()



#### Highly variant

markers <- cbind(Clusters=clin$hvar.Clusters, rna[samp.ids,] %>% select(rna.markers), prot[samp.ids,] %>% select(prot.markers))

markers$Prot_FLT1overPGF <- prot[samp.ids,'Prot_FLT1'] / prot[samp.ids,'Prot_PGF']

markers$RNA_FLT1overPGF <- rna[samp.ids,'RNA_FLT1'] / rna[samp.ids,'RNA_PGF']

head(markers)


markers.summary <- markers %>%
    group_by(Clusters) %>%
    summarize_all(c(mean=mean, sd=sd)) %>%
    as.data.frame

markers.summary <- markers.summary[,-1]

markers.table <- data.frame(matrix(NA, 5, (ncol(markers)-1)))

colnames(markers.table) <- colnames(markers)[-1]
rownames(markers.table) <- c(paste('Cluster ', 1:4, ' (n = ',
                                     as.vector(table(markers$Clusters)),
                                     ')', sep=''), 'pval')

for (var in colnames(markers)[-1]) {
    print(var)
    res <- kruskal.test(markers[,var], markers[,'Clusters'])
    print(res)

    markers.table[1:4,var] <- paste(signif(markers.summary[,paste(var, 'mean', sep='_')], 3),
                                   ' (',
                                   signif(markers.summary[,paste(var, 'sd', sep='_')], 3),
                                   ')', sep='')
    markers.table['pval',var] <- res$p.value
}

markers.table

markers.table <- as.data.frame(t(markers.table))

markers.table$pval <- as.numeric(markers.table$pval)

markers.table[!is.na(markers.table$pval),
           'adj.pval'] <- p.adjust(markers.table$pval[!is.na(markers.table$pval)], method='fdr')

markers.table[markers.table$adj.pval < 0.05,]

markers.table

write.csv(markers.table, 'markers_hvar_cluster_table.csv', quote=F)

library(ggplot2)
library(ggpubr)
library(rstatix)

expr.df <- data.frame()

for (var in colnames(markers)[-1]) {

    if (markers.table[var, 'adj.pval'] > 0.05) {
        next
    }

    expr <- data.frame()
    expr[rownames(markers),'Clusters'] <- factor(markers[,c('Clusters')])
    expr$marker <- var
    expr$Expression <- scale(markers[,var])[,1]
    expr$y.position <- 1.05 * max(expr$Expression)
    expr$adj.pval <- markers.table[var,'adj.pval']

    expr.df <- rbind(expr.df, expr)

    ## res <- compare_means(as.formula(paste(var, 'Clusters', sep='~')), markers, method='kruskal.test')
    ## res$p.adj <- markers.table[var,'adj.pval']
    
    ## gg <- ggplot(markers %>% mutate(Clusters=factor(Clusters)),
    ##              aes(x=Clusters, y=markers[,var], fill=Clusters)) +
    ##     geom_violin() +
    ##     geom_boxplot(width=0.2) +
    ##     xlab('Clusters') +
    ##     ylab(var) +
    ##     stat_pvalue_manual(res) + 
    ##     ## geom_text(x=2, y=1.15 * max(markers[,var]), label=paste('Kruskal, size=11) + 
    ##     theme_bw()

    ## print(gg)
}

stat.test <- expr.df %>%
    group_by(marker) %>%
    wilcox_test(Expression ~ Clusters) %>%
    add_xy_position(x = "marker", dodge = 0.8)

stat.test %>% as.data.frame

pdf('signif_markers_hvar_cluster_violin.pdf', width=15, height=7)

gg <- ggplot(expr.df,
             aes(x=Clusters, y=Expression, fill=Clusters)) +
    geom_violin() +
    geom_boxplot(width=0.2) +
    ## stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.01) + 
    xlab('Clusters') +
    ylab('Z-scored Expression') +
    scale_y_continuous(expand = expansion(mult = c(0.01, 0.1))) + 
    ## ylim(-4,4) + 
    facet_wrap(vars(marker), ncol=8, scales='free_y') +
    theme_bw()

gg

dev.off()


png('signif_markers_hvar_cluster_violin.png', width=15, height=7, units='in', res=400)

gg <- ggplot(expr.df,
             aes(x=Clusters, y=Expression, fill=Clusters)) +
    geom_violin() +
    geom_boxplot(width=0.2) +
    ## stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.01) + 
    xlab('Clusters') +
    ylab('Z-scored Expression') +
    scale_y_continuous(expand = expansion(mult = c(0.01, 0.1))) + 
    ## ylim(-4,4) + 
    facet_wrap(vars(marker), ncol=8, scales='free_y') +
    theme_bw()

gg

dev.off()



#### Distribution of Markers in Severe PE


markers <- cbind(Clusters=clin$Clusters, rna[samp.ids,] %>% select(rna.markers), prot[samp.ids,] %>% select(prot.markers))

markers$Prot_FLT1overPGF <- prot[samp.ids,'Prot_FLT1'] / prot[samp.ids,'Prot_PGF']

markers$RNA_FLT1overPGF <- rna[samp.ids,'RNA_FLT1'] / rna[samp.ids,'RNA_PGF']


markers <- markers %>% filter(clin$Condition. %in% c('Severe PE'))

markers.summary <- markers %>%
    group_by(Clusters) %>%
    summarize_all(c(mean=mean, sd=sd)) %>%
    as.data.frame

markers.summary <- markers.summary[,-1]

markers.table <- data.frame(matrix(NA, 4, (ncol(markers)-1)))

colnames(markers.table) <- colnames(markers)[-1]
rownames(markers.table) <- c(paste('Cluster ', 1:3, ' (n = ',
                                     as.vector(table(markers$Clusters)),
                                     ')', sep=''), 'pval')

for (var in colnames(markers)[-1]) {
    print(var)
    res <- kruskal.test(markers[,var], markers[,'Clusters'])
    print(res)

    markers.table[1:3,var] <- paste(signif(markers.summary[,paste(var, 'mean', sep='_')], 3),
                                   ' (',
                                   signif(markers.summary[,paste(var, 'sd', sep='_')], 3),
                                   ')', sep='')
    markers.table['pval',var] <- res$p.value
}

markers.table

markers.table <- as.data.frame(t(markers.table))

markers.table$pval <- as.numeric(markers.table$pval)

markers.table[!is.na(markers.table$pval),
           'adj.pval'] <- p.adjust(markers.table$pval[!is.na(markers.table$pval)], method='fdr')

markers.table[markers.table$adj.pval < 0.05,]

markers.table

write.csv(markers.table, 'SeverePE_markers_cluster_table.csv', quote=F)


expr.df <- data.frame()

for (var in colnames(markers)[-1]) {

    if (markers.table[var, 'adj.pval'] > 0.05) {
        next
    }


    expr <- data.frame()
    expr[rownames(markers),'Clusters'] <- factor(markers[,c('Clusters')])
    expr$marker <- var
    expr$Expression <- scale(markers[,var])[,1]

    expr.df <- rbind(expr.df, expr)

    ## res <- compare_means(as.formula(paste(var, 'Clusters', sep='~')), markers, method='kruskal.test')
    ## res$p.adj <- markers.table[var,'adj.pval']
    
    ## gg <- ggplot(markers %>% mutate(Clusters=factor(Clusters)),
    ##              aes(x=Clusters, y=markers[,var], fill=Clusters)) +
    ##     geom_violin() +
    ##     geom_boxplot(width=0.2) +
    ##     xlab('Clusters') +
    ##     ylab(var) +
    ##     stat_pvalue_manual(res) + 
    ##     ## geom_text(x=2, y=1.15 * max(markers[,var]), label=paste('Kruskal, size=11) + 
    ##     theme_bw()

    ## print(gg)
}

pdf('SeverePE_signif_markers_cluster_violin.pdf', width=15, height=10)

gg <- ggplot(expr.df,
             aes(x=Clusters, y=Expression, fill=Clusters)) +
    geom_violin() +
    geom_boxplot(width=0.2) +
    xlab('Clusters') +
    ylab('Z-scored Expression') +
    ## ylim(-4,4) + 
    facet_wrap(vars(marker), ncol=5, scales='free_y') +
    theme_bw()

gg

dev.off()


png('SeverePE_signif_markers_cluster_violin.png', width=15, height=10, res=400, units='in')

gg <- ggplot(expr.df,
             aes(x=Clusters, y=Expression, fill=Clusters)) +
    geom_violin() +
    geom_boxplot(width=0.2) +
    xlab('Clusters') +
    ylab('Z-scored Expression') +
    ## ylim(-4,4) + 
    facet_wrap(vars(marker), ncol=5, scales='free_y') +
    theme_bw()

gg

dev.off()




#### Distribution of Markers in Any PE


markers <- cbind(Clusters=clin$Clusters, rna[samp.ids,] %>% select(rna.markers), prot[samp.ids,] %>% select(prot.markers))

markers$Prot_FLT1overPGF <- prot[samp.ids,'Prot_FLT1'] / prot[samp.ids,'Prot_PGF']

markers$RNA_FLT1overPGF <- rna[samp.ids,'RNA_FLT1'] / rna[samp.ids,'RNA_PGF']


markers <- markers %>% filter(clin$Condition. %in% c('Severe PE', 'FGR+HDP'))

markers.summary <- markers %>%
    group_by(Clusters) %>%
    summarize_all(c(mean=mean, sd=sd)) %>%
    as.data.frame

markers.summary <- markers.summary[,-1]

markers.table <- data.frame(matrix(NA, 4, (ncol(markers)-1)))

colnames(markers.table) <- colnames(markers)[-1]
rownames(markers.table) <- c(paste('Cluster ', 1:3, ' (n = ',
                                     as.vector(table(markers$Clusters)),
                                     ')', sep=''), 'pval')

for (var in colnames(markers)[-1]) {
    print(var)
    res <- kruskal.test(markers[,var], markers[,'Clusters'])
    print(res)

    markers.table[1:3,var] <- paste(signif(markers.summary[,paste(var, 'mean', sep='_')], 3),
                                   ' (',
                                   signif(markers.summary[,paste(var, 'sd', sep='_')], 3),
                                   ')', sep='')
    markers.table['pval',var] <- res$p.value
}

markers.table

markers.table <- as.data.frame(t(markers.table))

markers.table$pval <- as.numeric(markers.table$pval)

markers.table[!is.na(markers.table$pval),
           'adj.pval'] <- p.adjust(markers.table$pval[!is.na(markers.table$pval)], method='fdr')

markers.table[markers.table$adj.pval < 0.05,]

markers.table

write.csv(markers.table, 'AnyPE_markers_cluster_table.csv', quote=F)


expr.df <- data.frame()

for (var in colnames(markers)[-1]) {

    if (markers.table[var, 'adj.pval'] > 0.05) {
        next
    }

    expr <- data.frame()
    expr[rownames(markers),'Clusters'] <- factor(markers[,c('Clusters')])
    expr$marker <- var
    expr$Expression <- scale(markers[,var])[,1]

    expr.df <- rbind(expr.df, expr)

    ## res <- compare_means(as.formula(paste(var, 'Clusters', sep='~')), markers, method='kruskal.test')
    ## res$p.adj <- markers.table[var,'adj.pval']
    
    ## gg <- ggplot(markers %>% mutate(Clusters=factor(Clusters)),
    ##              aes(x=Clusters, y=markers[,var], fill=Clusters)) +
    ##     geom_violin() +
    ##     geom_boxplot(width=0.2) +
    ##     xlab('Clusters') +
    ##     ylab(var) +
    ##     stat_pvalue_manual(res) + 
    ##     ## geom_text(x=2, y=1.15 * max(markers[,var]), label=paste('Kruskal, size=11) + 
    ##     theme_bw()

    ## print(gg)
}

pdf('AnyPE_signif_markers_cluster_violin.pdf', width=15, height=10)

gg <- ggplot(expr.df,
             aes(x=Clusters, y=Expression, fill=Clusters)) +
    geom_violin() +
    geom_boxplot(width=0.2) +
    xlab('Clusters') +
    ylab('Z-scored Expression') +
    ## ylim(-4,4) + 
    facet_wrap(vars(marker), ncol=5, scales='free_y') +
    theme_bw()

gg

dev.off()


png('AnyPE_signif_markers_cluster_violin.png', width=15, height=10, res=400, units='in')

gg <- ggplot(expr.df,
             aes(x=Clusters, y=Expression, fill=Clusters)) +
    geom_violin() +
    geom_boxplot(width=0.2) +
    xlab('Clusters') +
    ylab('Z-scored Expression') +
    ## ylim(-4,4) + 
    facet_wrap(vars(marker), ncol=5, scales='free_y') +
    theme_bw()

gg

dev.off()




#### Distribution of Markers in PTD

markers <- cbind(Clusters=clin$Clusters, rna[samp.ids,] %>% select(rna.markers), prot[samp.ids,] %>% select(prot.markers))

markers$Prot_FLT1overPGF <- prot[samp.ids,'Prot_FLT1'] / prot[samp.ids,'Prot_PGF']

markers$RNA_FLT1overPGF <- rna[samp.ids,'RNA_FLT1'] / rna[samp.ids,'RNA_PGF']


markers <- markers %>% filter(clin$Condition. %in% c('PTD'))

markers.summary <- markers %>%
    group_by(Clusters) %>%
    summarize_all(c(mean=mean, sd=sd)) %>%
    as.data.frame

markers.summary <- markers.summary[,-1]

markers.table <- data.frame(matrix(NA, 4, (ncol(markers)-1)))

colnames(markers.table) <- colnames(markers)[-1]
rownames(markers.table) <- c(paste('Cluster ', 1:3, ' (n = ',
                                     as.vector(table(markers$Clusters)),
                                     ')', sep=''), 'pval')

for (var in colnames(markers)[-1]) {
    print(var)
    res <- kruskal.test(markers[,var], markers[,'Clusters'])
    print(res)

    markers.table[1:3,var] <- paste(signif(markers.summary[,paste(var, 'mean', sep='_')], 3),
                                   ' (',
                                   signif(markers.summary[,paste(var, 'sd', sep='_')], 3),
                                   ')', sep='')
    markers.table['pval',var] <- res$p.value
}

markers.table

markers.table <- as.data.frame(t(markers.table))

markers.table$pval <- as.numeric(markers.table$pval)

markers.table[!is.na(markers.table$pval),
           'adj.pval'] <- p.adjust(markers.table$pval[!is.na(markers.table$pval)], method='fdr')

markers.table[markers.table$adj.pval < 0.05,]

markers.table

write.csv(markers.table, 'PTD_markers_cluster_table.csv', quote=F)


expr.df <- data.frame()

for (var in colnames(markers)[-1]) {

    if (markers.table[var, 'adj.pval'] > 0.05) {
        next
    }

    expr <- data.frame()
    expr[rownames(markers),'Clusters'] <- factor(markers[,c('Clusters')])
    expr$marker <- var
    expr$Expression <- scale(markers[,var])[,1]

    expr.df <- rbind(expr.df, expr)

    ## res <- compare_means(as.formula(paste(var, 'Clusters', sep='~')), markers, method='kruskal.test')
    ## res$p.adj <- markers.table[var,'adj.pval']
    
    ## gg <- ggplot(markers %>% mutate(Clusters=factor(Clusters)),
    ##              aes(x=Clusters, y=markers[,var], fill=Clusters)) +
    ##     geom_violin() +
    ##     geom_boxplot(width=0.2) +
    ##     xlab('Clusters') +
    ##     ylab(var) +
    ##     stat_pvalue_manual(res) + 
    ##     ## geom_text(x=2, y=1.15 * max(markers[,var]), label=paste('Kruskal, size=11) + 
    ##     theme_bw()

    ## print(gg)
}

pdf('PTD_signif_markers_cluster_violin.pdf', width=15, height=10)

gg <- ggplot(expr.df,
             aes(x=Clusters, y=Expression, fill=Clusters)) +
    geom_violin() +
    geom_boxplot(width=0.2) +
    xlab('Clusters') +
    ylab('Z-scored Expression') +
    ## ylim(-4,4) + 
    facet_wrap(vars(marker), ncol=4, scales='free_y') +
    theme_bw()

gg

dev.off()


png('PTD_signif_markers_cluster_violin.png', width=15, height=10, res=400, units='in')

gg <- ggplot(expr.df,
             aes(x=Clusters, y=Expression, fill=Clusters)) +
    geom_violin() +
    geom_boxplot(width=0.2) +
    xlab('Clusters') +
    ylab('Z-scored Expression') +
    ## ylim(-4,4) + 
    facet_wrap(vars(marker), ncol=4, scales='free_y') +
    theme_bw()

gg

dev.off()



#### Differential Expression with Cluster 1 as reference

library(DESeq2)
library(vsn)

rna.counts <- read.csv('data/combat_rna_counts.csv', row.names=1)
head(rna.counts)

samp.ids <- rownames(clin)

clin <- clin %>% mutate_at(c('hvar.Clusters', 'InfSex', 'Race', 'Smoking', 'FDELTYPE',
                             'Labor.initiation'), factor)

clin <- clin %>% mutate_at(c('WksGest', 'PrePregBMI'), function(x) scale(x)[,1])

dds <- DESeqDataSetFromMatrix(t(rna.counts[samp.ids,]), clin,
                              design=~hvar.Clusters) ## + WksGest + InfSex + Race +
                                  ## PrePregBMI + Smoking + FDELTYPE + Labor.initiation)

dds <- DESeq(dds)

plotDispEsts(dds)

vsd <- vst(dds, blind=FALSE)

meanSdPlot(assay(vsd))

## 2 vs 1 Differential Expression

res <- results(dds, contrast=c('hvar.Clusters', '2', '1'), alpha=0.05, independentFiltering=TRUE)
print(summary(res))

## resLFC <- lfcShrink(dds, contrast=c('Condition.', groups), res=res, type='ashr')
## print(summary(resLFC))

de <- res[order(res$pvalue),] %>% as.data.frame
de <- mutate(de, RNA=rownames(de), .before=1) %>% filter(baseMean > 0)

## vsd <- vst(dds, blind=FALSE)

group.meta <- colData(dds)[colData(dds)[,'hvar.Clusters'] %in% c('1', '2'),]


if (nrow(de[!is.na(de$padj) & de$padj<0.05,])>1) {
    pdf(paste('snf_hvar_clusts/DE/plots/cluster2_diffExp_heatmap.pdf',sep=''),
        width=5 + 0.12 * nrow(group.meta),
        height=3 + 0.12 * 500)

    scaled.sig <- scale(t(assay(vsd))[rownames(group.meta), 1:500])
                                      ## de[!is.na(de$padj) & de$padj<0.05,'RNA']])
    range <- min(max(abs(scaled.sig)),3)
        
    annot_df <- data.frame(Condition = as.factor(group.meta$Condition.),
                           Clusters = as.factor(group.meta$hvar.Clusters),
                           GestationalAge = group.meta$WksGest, 
                           row.names=rownames(group.meta))
    
    pheatmap(t(scaled.sig), breaks = seq(-range, range, length.out = 101),
             annotation_col=annot_df, main=paste('Cluster 2', 'vs.', 'Cluster 1'),
             ## cluster_rows=sum(de$adj.P.Val < 0.05)>1,
             clustering_method='ward.D', clustering_distance_rows = "euclidean",
             clustering_distance_cols = "euclidean",
             labels_row=de[1:500,'RNA'])
    dev.off()
}

write.csv(de, 'snf_hvar_clusts/DE/diffExpRNA_1_2.csv')
write.csv(de[!is.na(de$padj) & de$padj<0.05,], 'snf_hvar_clusts/DE/diffExpRNA_1_2_FDR05.csv')
write.csv(de[!is.na(de$padj) & de$padj<0.1,], 'snf_hvar_clusts/DE/diffExpRNA_1_2_FDR1.csv')


library(enrichplot)
library(clusterProfiler)
library(org.Hs.eg.db)

geneList <- de$stat
keys <- bitr(de$RNA,
             fromType="SYMBOL",
             toType="ENTREZID",
             OrgDb="org.Hs.eg.db") %>%
    group_by(SYMBOL) %>%
    filter(row_number()==1)
names(geneList) <- unlist(as.vector(keys[,2]))
geneList <- sort(geneList, decreasing=TRUE)

set.seed(20220526)

gsea.go <- gseGO(geneList     = geneList,
                 OrgDb        = org.Hs.eg.db,
                 ont          = "BP",
                 minGSSize    = 50,
                 maxGSSize    = 500,
                 pvalueCutoff = 0.05,
                 verbose      = FALSE)

gsea.go <- simplify(gsea.go, cutoff=0.7, by='p.adjust', select_fun=min)

gsea.go

xx <- pairwise_termsim(gsea.go, showCategory=nrow(gsea.go))
p1 <- emapplot(xx, showCategory=nrow(gsea.go), max.overlaps=0)

p1

gsea.go.up <- gsea.go
gsea.go.up@result <- gsea.go.up@result[gsea.go$NES>0,]
gsea.go.up@geneSets <- gsea.go.up@geneSets[gsea.go$NES>0]

gsea.go.up

xx <- pairwise_termsim(gsea.go.up, showCategory=nrow(gsea.go.up))
p2 <- emapplot(xx, showCategory=nrow(gsea.go.up), max.overlaps=0, repel=TRUE, force=10, group_category=TRUE, node_label='group', nCluster=floor(2*sqrt(nrow(gsea.go.up))), nWords=5, clusterFunction=cluster::pam)
p2

gsea.go.down <- gsea.go
gsea.go.down@result <- gsea.go.down@result[gsea.go$NES<0,]
gsea.go.down@geneSets <- gsea.go.down@geneSets[gsea.go$NES<0]

gsea.go.down

xx <- pairwise_termsim(gsea.go.down, showCategory=nrow(gsea.go.down))
p3 <- emapplot(xx, showCategory=nrow(gsea.go.down), max.overlaps=0, repel=TRUE, force=10)

p3

gsea.go.df <- as.data.frame(gsea.go)

head(gsea.go.df[gsea.go.df$NES>0,1:10], 20)

head(gsea.go.df[gsea.go.df$NES<0,1:10], 20)

library(cowplot)
plot_grid(p2, p3, nrow=1, labels=c('Up', 'Down'))

ggsave('snf_hvar_clusts/DE/plots/cluster2_go_gsea.pdf',
       plot_grid(p2, p3, nrow=1, labels=c('Upregulated', 'Downregulated')),
       width=20, height=10, dpi=400)
ggsave('snf_hvar_clusts/DE/plots/cluster2_go_gsea.png',
       plot_grid(p2, p3, nrow=1, labels=c('Upregulated', 'Downregulated')),
       width=20, height=10, dpi=400)

write.csv(gsea.go.df, 'snf_hvar_clusts/DE/cluster2_go_gsea.csv')

## 3 vs 1 Differential Expression

res <- results(dds, contrast=c('hvar.Clusters', '3', '1'), alpha=0.05, independentFiltering=TRUE)
print(summary(res))

## resLFC <- lfcShrink(dds, contrast=c('Condition.', groups), res=res, type='ashr')
## print(summary(resLFC))

de <- res[order(res$pvalue),] %>% as.data.frame
de <- mutate(de, RNA=rownames(de), .before=1) %>% filter(baseMean > 0)

## vsd <- vst(dds, blind=FALSE)

group.meta <- colData(dds)[colData(dds)[,'hvar.Clusters'] %in% c('1', '3'),]


if (nrow(de[!is.na(de$padj) & de$padj<0.05,])>1) {
    pdf(paste('snf_hvar_clusts/DE/plots/cluster3_diffExp_heatmap.pdf',sep=''),
        width=5 + 0.12 * nrow(group.meta),
        height=3 + 0.12 * 500)

    scaled.sig <- scale(t(assay(vsd))[rownames(group.meta), 1:500])
                                      ## de[!is.na(de$padj) & de$padj<0.05,'RNA']])
    range <- min(max(abs(scaled.sig)),3)
        
    annot_df <- data.frame(Condition = as.factor(group.meta$Condition.),
                           Clusters = as.factor(group.meta$hvar.Clusters),
                           GestationalAge = group.meta$WksGest, 
                           row.names=rownames(group.meta))
    
    pheatmap(t(scaled.sig), breaks = seq(-range, range, length.out = 101),
             annotation_col=annot_df, main=paste('Cluster 3', 'vs.', 'Cluster 1'),
             ## cluster_rows=sum(de$adj.P.Val < 0.05)>1,
             clustering_method='ward.D', clustering_distance_rows = "euclidean",
             clustering_distance_cols = "euclidean",
             labels_row=de[1:500,'RNA'])
    dev.off()
}

write.csv(de, 'snf_hvar_clusts/DE/diffExpRNA_1_3.csv')
write.csv(de[!is.na(de$padj) & de$padj<0.05,], 'snf_hvar_clusts/DE/diffExpRNA_1_3_FDR05.csv')
write.csv(de[!is.na(de$padj) & de$padj<0.1,], 'snf_hvar_clusts/DE/diffExpRNA_1_3_FDR1.csv')


geneList <- de$stat
keys <- bitr(de$RNA,
             fromType="SYMBOL",
             toType="ENTREZID",
             OrgDb="org.Hs.eg.db") %>%
    group_by(SYMBOL) %>%
    filter(row_number()==1)
names(geneList) <- unlist(as.vector(keys[,2]))
geneList <- sort(geneList, decreasing=TRUE)

set.seed(20220526)

gsea.go <- gseGO(geneList     = geneList,
                 OrgDb        = org.Hs.eg.db,
                 ont          = "BP",
                 minGSSize    = 50,
                 maxGSSize    = 500,
                 pvalueCutoff = 0.05,
                 verbose      = FALSE)

gsea.go <- simplify(gsea.go, cutoff=0.7, by='p.adjust', select_fun=min)

gsea.go

xx <- pairwise_termsim(gsea.go, showCategory=nrow(gsea.go))
p1 <- emapplot(xx, showCategory=nrow(gsea.go), max.overlaps=0)

p1

gsea.go.up <- gsea.go
gsea.go.up@result <- gsea.go.up@result[gsea.go$NES>0,]
gsea.go.up@geneSets <- gsea.go.up@geneSets[gsea.go$NES>0]

gsea.go.up

xx <- pairwise_termsim(gsea.go.up, showCategory=nrow(gsea.go.up))
p2 <- emapplot(xx, showCategory=nrow(gsea.go.up), max.overlaps=0, repel=TRUE, force=10, group_category=TRUE, node_label='group', nCluster=floor(2*sqrt(nrow(gsea.go.up))), nWords=5, clusterFunction=cluster::pam)

p2

gsea.go.down <- gsea.go
gsea.go.down@result <- gsea.go.down@result[gsea.go$NES<0,]
gsea.go.down@geneSets <- gsea.go.down@geneSets[gsea.go$NES<0]

gsea.go.down

xx <- pairwise_termsim(gsea.go.down, showCategory=nrow(gsea.go.down))
p3 <- emapplot(xx, showCategory=nrow(gsea.go.down), max.overlaps=0, repel=TRUE, force=10)

p3

gsea.go.df <- as.data.frame(gsea.go)

head(gsea.go.df[gsea.go.df$NES>0,1:10], 20)

head(gsea.go.df[gsea.go.df$NES<0,1:10], 20)

library(cowplot)
plot_grid(p2, p3, nrow=1, labels=c('Up', 'Down'))

ggsave('snf_hvar_clusts/DE/plots/cluster3_go_gsea.pdf',
       plot_grid(p2, p3, nrow=1, labels=c('Upregulated', 'Downregulated')),
       width=20, height=10, dpi=400)
ggsave('snf_hvar_clusts/DE/plots/cluster3_go_gsea.png',
       plot_grid(p2, p3, nrow=1, labels=c('Upregulated', 'Downregulated')),
       width=20, height=10, dpi=400)

write.csv(gsea.go.df, 'snf_hvar_clusts/DE/cluster3_go_gsea.csv')


## 4 vs 1 Differential Expression

res <- results(dds, contrast=c('hvar.Clusters', '4', '1'), alpha=0.05, independentFiltering=TRUE)
print(summary(res))

## resLFC <- lfcShrink(dds, contrast=c('Condition.', groups), res=res, type='ashr')
## print(summary(resLFC))

de <- res[order(res$pvalue),] %>% as.data.frame
de <- mutate(de, RNA=rownames(de), .before=1) %>% filter(baseMean > 0)

## vsd <- vst(dds, blind=FALSE)

group.meta <- colData(dds)[colData(dds)[,'hvar.Clusters'] %in% c('1', '4'),]


if (nrow(de[!is.na(de$padj) & de$padj<0.05,])>1) {
    pdf(paste('snf_hvar_clusts/DE/plots/cluster4_diffExp_heatmap.pdf',sep=''),
        width=5 + 0.12 * nrow(group.meta),
        height=4 + 0.12 * 500)

    scaled.sig <- scale(t(assay(vsd))[rownames(group.meta), 1:500])
                                      ## de[!is.na(de$padj) & de$padj<0.05,'RNA']])
    range <- min(max(abs(scaled.sig)),4)
        
    annot_df <- data.frame(Condition = as.factor(group.meta$Condition.),
                           Clusters = as.factor(group.meta$hvar.Clusters),
                           GestationalAge = group.meta$WksGest, 
                           row.names=rownames(group.meta))
    
    pheatmap(t(scaled.sig), breaks = seq(-range, range, length.out = 101),
             annotation_col=annot_df, main=paste('Cluster 4', 'vs.', 'Cluster 1'),
             ## cluster_rows=sum(de$adj.P.Val < 0.05)>1,
             clustering_method='ward.D', clustering_distance_rows = "euclidean",
             clustering_distance_cols = "euclidean",
             labels_row=de[1:500,'RNA'])
    dev.off()
}

write.csv(de, 'snf_hvar_clusts/DE/diffExpRNA_1_4.csv')
write.csv(de[!is.na(de$padj) & de$padj<0.05,], 'snf_hvar_clusts/DE/diffExpRNA_1_4_FDR05.csv')
write.csv(de[!is.na(de$padj) & de$padj<0.1,], 'snf_hvar_clusts/DE/diffExpRNA_1_4_FDR1.csv')


geneList <- de$stat
keys <- bitr(de$RNA,
             fromType="SYMBOL",
             toType="ENTREZID",
             OrgDb="org.Hs.eg.db") %>%
    group_by(SYMBOL) %>%
    filter(row_number()==1)
names(geneList) <- unlist(as.vector(keys[,2]))
geneList <- sort(geneList, decreasing=TRUE)

set.seed(20220526)

gsea.go <- gseGO(geneList     = geneList,
                 OrgDb        = org.Hs.eg.db,
                 ont          = "BP",
                 minGSSize    = 50,
                 maxGSSize    = 500,
                 pvalueCutoff = 0.05,
                 verbose      = FALSE)

gsea.go <- simplify(gsea.go, cutoff=0.7, by='p.adjust', select_fun=min)

gsea.go

xx <- pairwise_termsim(gsea.go, showCategory=nrow(gsea.go))
p1 <- emapplot(xx, showCategory=nrow(gsea.go), max.overlaps=0)

p1

gsea.go.up <- gsea.go
gsea.go.up@result <- gsea.go.up@result[gsea.go$NES>0,]
gsea.go.up@geneSets <- gsea.go.up@geneSets[gsea.go$NES>0]

gsea.go.up

xx <- pairwise_termsim(gsea.go.up, showCategory=nrow(gsea.go.up))
p2 <- emapplot(xx, showCategory=nrow(gsea.go.up), max.overlaps=0, repel=TRUE, force=10, group_category=TRUE, node_label='group', nCluster=floor(2*sqrt(nrow(gsea.go.up))), nWords=5, clusterFunction=cluster::pam)

p2

gsea.go.down <- gsea.go
gsea.go.down@result <- gsea.go.down@result[gsea.go$NES<0,]
gsea.go.down@geneSets <- gsea.go.down@geneSets[gsea.go$NES<0]

gsea.go.down

xx <- pairwise_termsim(gsea.go.down, showCategory=nrow(gsea.go.down))
p3 <- emapplot(xx, showCategory=nrow(gsea.go.down), max.overlaps=0, repel=TRUE, force=10)

p3

gsea.go.df <- as.data.frame(gsea.go)

head(gsea.go.df[gsea.go.df$NES>0,1:10], 20)

head(gsea.go.df[gsea.go.df$NES<0,1:10], 20)

library(cowplot)
plot_grid(p2, p3, nrow=1, labels=c('Up', 'Down'))

ggsave('snf_hvar_clusts/DE/plots/cluster4_go_gsea.pdf',
       plot_grid(p2, p4, nrow=1, labels=c('Upregulated', 'Downregulated')),
       width=20, height=10, dpi=400)
ggsave('snf_hvar_clusts/DE/plots/cluster4_go_gsea.png',
       plot_grid(p2, p4, nrow=1, labels=c('Upregulated', 'Downregulated')),
       width=20, height=10, dpi=400)

write.csv(gsea.go.df, 'snf_hvar_clusts/DE/cluster4_go_gsea.csv')



#### Conditional DE by cluster


dds <- DESeqDataSetFromMatrix(t(rna.counts[samp.ids,]), clin,
                              design=~hvar.Clusters + WksGest + InfSex + Race +
                                  PrePregBMI + Smoking + FDELTYPE + Labor.initiation)

dds <- DESeq(dds)

plotDispEsts(dds)

vsd <- vst(dds, blind=FALSE)

meanSdPlot(assay(vsd))

## 2 vs 1 Differential Expression

res <- results(dds, contrast=c('hvar.Clusters', '2', '1'), alpha=0.05, independentFiltering=TRUE)
print(summary(res))

## resLFC <- lfcShrink(dds, contrast=c('Condition.', groups), res=res, type='ashr')
## print(summary(resLFC))

de <- res[order(res$pvalue),] %>% as.data.frame
de <- mutate(de, RNA=rownames(de), .before=1) %>% filter(baseMean > 0)

## vsd <- vst(dds, blind=FALSE)

group.meta <- colData(dds)[colData(dds)[,'hvar.Clusters'] %in% c('1', '2'),]


if (nrow(de[!is.na(de$padj) & de$padj<0.05,])>1) {
    pdf(paste('snf_hvar_clusts/DE-cond2/plots/cluster2_diffExp_heatmap.pdf',sep=''),
        width=5 + 0.12 * nrow(group.meta),
        height=3 + 0.12 * 500)

    scaled.sig <- scale(t(assay(vsd))[rownames(group.meta), 1:500])
                                      ## de[!is.na(de$padj) & de$padj<0.05,'RNA']])
    range <- min(max(abs(scaled.sig)),3)
        
    annot_df <- data.frame(Condition = as.factor(group.meta$Condition.),
                           Clusters = as.factor(group.meta$hvar.Clusters),
                           GestationalAge = group.meta$WksGest, 
                           row.names=rownames(group.meta))
    
    pheatmap(t(scaled.sig), breaks = seq(-range, range, length.out = 101),
             annotation_col=annot_df, main=paste('Cluster 2', 'vs.', 'Cluster 1'),
             ## cluster_rows=sum(de$adj.P.Val < 0.05)>1,
             clustering_method='ward.D', clustering_distance_rows = "euclidean",
             clustering_distance_cols = "euclidean",
             labels_row=de[1:500,'RNA'])
    dev.off()
}

write.csv(de, 'snf_hvar_clusts/DE-cond2/diffExpRNA_1_2.csv')
write.csv(de[!is.na(de$padj) & de$padj<0.05,], 'snf_hvar_clusts/DE-cond2/diffExpRNA_1_2_FDR05.csv')
write.csv(de[!is.na(de$padj) & de$padj<0.1,], 'snf_hvar_clusts/DE-cond2/diffExpRNA_1_2_FDR1.csv')


library(enrichplot)
library(clusterProfiler)
library(org.Hs.eg.db)

geneList <- de$stat
keys <- bitr(de$RNA,
             fromType="SYMBOL",
             toType="ENTREZID",
             OrgDb="org.Hs.eg.db") %>%
    group_by(SYMBOL) %>%
    filter(row_number()==1)
names(geneList) <- unlist(as.vector(keys[,2]))
geneList <- sort(geneList, decreasing=TRUE)

set.seed(20220526)

gsea.go <- gseGO(geneList     = geneList,
                 OrgDb        = org.Hs.eg.db,
                 ont          = "BP",
                 minGSSize    = 50,
                 maxGSSize    = 500,
                 pvalueCutoff = 0.05,
                 verbose      = FALSE)

gsea.go <- simplify(gsea.go, cutoff=0.7, by='p.adjust', select_fun=min)

gsea.go

xx <- pairwise_termsim(gsea.go, showCategory=nrow(gsea.go))
p1 <- emapplot(xx, showCategory=nrow(gsea.go), max.overlaps=0)

p1

gsea.go.up <- gsea.go
gsea.go.up@result <- gsea.go.up@result[gsea.go$NES>0,]
gsea.go.up@geneSets <- gsea.go.up@geneSets[gsea.go$NES>0]

gsea.go.up

xx <- pairwise_termsim(gsea.go.up, showCategory=nrow(gsea.go.up))
p2 <- emapplot(xx, showCategory=nrow(gsea.go.up), max.overlaps=0, repel=TRUE, force=10, group_category=TRUE, node_label='group', nCluster=floor(2*sqrt(nrow(gsea.go.up))), nWords=5, clusterFunction=cluster::pam)
p2

gsea.go.down <- gsea.go
gsea.go.down@result <- gsea.go.down@result[gsea.go$NES<0,]
gsea.go.down@geneSets <- gsea.go.down@geneSets[gsea.go$NES<0]

gsea.go.down

xx <- pairwise_termsim(gsea.go.down, showCategory=nrow(gsea.go.down))
p3 <- emapplot(xx, showCategory=nrow(gsea.go.down), max.overlaps=0, repel=TRUE, force=10)

p3

gsea.go.df <- as.data.frame(gsea.go)

head(gsea.go.df[gsea.go.df$NES>0,1:10], 20)

head(gsea.go.df[gsea.go.df$NES<0,1:10], 20)

library(cowplot)
plot_grid(p2, p3, nrow=1, labels=c('Up', 'Down'))

ggsave('snf_hvar_clusts/DE-cond2/plots/cluster2_go_gsea.pdf',
       plot_grid(p2, p3, nrow=1, labels=c('Upregulated', 'Downregulated')),
       width=20, height=10, dpi=400)
ggsave('snf_hvar_clusts/DE-cond2/plots/cluster2_go_gsea.png',
       plot_grid(p2, p3, nrow=1, labels=c('Upregulated', 'Downregulated')),
       width=20, height=10, dpi=400)

write.csv(gsea.go.df, 'snf_hvar_clusts/DE-cond2/cluster2_go_gsea.csv')

## 3 vs 1 Differential Expression

res <- results(dds, contrast=c('hvar.Clusters', '3', '1'), alpha=0.05, independentFiltering=TRUE)
print(summary(res))

## resLFC <- lfcShrink(dds, contrast=c('Condition.', groups), res=res, type='ashr')
## print(summary(resLFC))

de <- res[order(res$pvalue),] %>% as.data.frame
de <- mutate(de, RNA=rownames(de), .before=1) %>% filter(baseMean > 0)

## vsd <- vst(dds, blind=FALSE)

group.meta <- colData(dds)[colData(dds)[,'hvar.Clusters'] %in% c('1', '3'),]


if (nrow(de[!is.na(de$padj) & de$padj<0.05,])>1) {
    pdf(paste('snf_hvar_clusts/DE-cond2/plots/cluster3_diffExp_heatmap.pdf',sep=''),
        width=5 + 0.12 * nrow(group.meta),
        height=3 + 0.12 * 500)

    scaled.sig <- scale(t(assay(vsd))[rownames(group.meta), 1:500])
                                      ## de[!is.na(de$padj) & de$padj<0.05,'RNA']])
    range <- min(max(abs(scaled.sig)),3)
        
    annot_df <- data.frame(Condition = as.factor(group.meta$Condition.),
                           Clusters = as.factor(group.meta$hvar.Clusters),
                           GestationalAge = group.meta$WksGest, 
                           row.names=rownames(group.meta))
    
    pheatmap(t(scaled.sig), breaks = seq(-range, range, length.out = 101),
             annotation_col=annot_df, main=paste('Cluster 3', 'vs.', 'Cluster 1'),
             ## cluster_rows=sum(de$adj.P.Val < 0.05)>1,
             clustering_method='ward.D', clustering_distance_rows = "euclidean",
             clustering_distance_cols = "euclidean",
             labels_row=de[1:500,'RNA'])
    dev.off()
}

write.csv(de, 'snf_hvar_clusts/DE-cond2/diffExpRNA_1_3.csv')
write.csv(de[!is.na(de$padj) & de$padj<0.05,], 'snf_hvar_clusts/DE-cond2/diffExpRNA_1_3_FDR05.csv')
write.csv(de[!is.na(de$padj) & de$padj<0.1,], 'snf_hvar_clusts/DE-cond2/diffExpRNA_1_3_FDR1.csv')


geneList <- de$stat
keys <- bitr(de$RNA,
             fromType="SYMBOL",
             toType="ENTREZID",
             OrgDb="org.Hs.eg.db") %>%
    group_by(SYMBOL) %>%
    filter(row_number()==1)
names(geneList) <- unlist(as.vector(keys[,2]))
geneList <- sort(geneList, decreasing=TRUE)

set.seed(20220526)

gsea.go <- gseGO(geneList     = geneList,
                 OrgDb        = org.Hs.eg.db,
                 ont          = "BP",
                 minGSSize    = 50,
                 maxGSSize    = 500,
                 pvalueCutoff = 0.05,
                 verbose      = FALSE)

gsea.go <- simplify(gsea.go, cutoff=0.7, by='p.adjust', select_fun=min)

gsea.go

xx <- pairwise_termsim(gsea.go, showCategory=nrow(gsea.go))
p1 <- emapplot(xx, showCategory=nrow(gsea.go), max.overlaps=0)

p1

gsea.go.up <- gsea.go
gsea.go.up@result <- gsea.go.up@result[gsea.go$NES>0,]
gsea.go.up@geneSets <- gsea.go.up@geneSets[gsea.go$NES>0]

gsea.go.up

xx <- pairwise_termsim(gsea.go.up, showCategory=nrow(gsea.go.up))
p2 <- emapplot(xx, showCategory=nrow(gsea.go.up), max.overlaps=0, repel=TRUE, force=10, group_category=TRUE, node_label='group', nCluster=floor(2*sqrt(nrow(gsea.go.up))), nWords=5, clusterFunction=cluster::pam)

p2

gsea.go.down <- gsea.go
gsea.go.down@result <- gsea.go.down@result[gsea.go$NES<0,]
gsea.go.down@geneSets <- gsea.go.down@geneSets[gsea.go$NES<0]

gsea.go.down

xx <- pairwise_termsim(gsea.go.down, showCategory=nrow(gsea.go.down))
p3 <- emapplot(xx, showCategory=nrow(gsea.go.down), max.overlaps=0, repel=TRUE, force=10)

p3

gsea.go.df <- as.data.frame(gsea.go)

head(gsea.go.df[gsea.go.df$NES>0,1:10], 20)

head(gsea.go.df[gsea.go.df$NES<0,1:10], 20)

library(cowplot)
plot_grid(p2, p3, nrow=1, labels=c('Up', 'Down'))

ggsave('snf_hvar_clusts/DE-cond2/plots/cluster3_go_gsea.pdf',
       plot_grid(p2, p3, nrow=1, labels=c('Upregulated', 'Downregulated')),
       width=20, height=10, dpi=400)
ggsave('snf_hvar_clusts/DE-cond2/plots/cluster3_go_gsea.png',
       plot_grid(p2, p3, nrow=1, labels=c('Upregulated', 'Downregulated')),
       width=20, height=10, dpi=400)

write.csv(gsea.go.df, 'snf_hvar_clusts/DE-cond2/cluster3_go_gsea.csv')


## 4 vs 1 Differential Expression

res <- results(dds, contrast=c('hvar.Clusters', '4', '1'), alpha=0.05, independentFiltering=TRUE)
print(summary(res))

## resLFC <- lfcShrink(dds, contrast=c('Condition.', groups), res=res, type='ashr')
## print(summary(resLFC))

de <- res[order(res$pvalue),] %>% as.data.frame
de <- mutate(de, RNA=rownames(de), .before=1) %>% filter(baseMean > 0)

## vsd <- vst(dds, blind=FALSE)

group.meta <- colData(dds)[colData(dds)[,'hvar.Clusters'] %in% c('1', '4'),]


if (nrow(de[!is.na(de$padj) & de$padj<0.05,])>1) {
    pdf(paste('snf_hvar_clusts/DE-cond2/plots/cluster4_diffExp_heatmap.pdf',sep=''),
        width=5 + 0.12 * nrow(group.meta),
        height=4 + 0.12 * 500)

    scaled.sig <- scale(t(assay(vsd))[rownames(group.meta), 1:500])
                                      ## de[!is.na(de$padj) & de$padj<0.05,'RNA']])
    range <- min(max(abs(scaled.sig)),4)
        
    annot_df <- data.frame(Condition = as.factor(group.meta$Condition.),
                           Clusters = as.factor(group.meta$hvar.Clusters),
                           GestationalAge = group.meta$WksGest, 
                           row.names=rownames(group.meta))
    
    pheatmap(t(scaled.sig), breaks = seq(-range, range, length.out = 101),
             annotation_col=annot_df, main=paste('Cluster 4', 'vs.', 'Cluster 1'),
             ## cluster_rows=sum(de$adj.P.Val < 0.05)>1,
             clustering_method='ward.D', clustering_distance_rows = "euclidean",
             clustering_distance_cols = "euclidean",
             labels_row=de[1:500,'RNA'])
    dev.off()
}

write.csv(de, 'snf_hvar_clusts/DE-cond2/diffExpRNA_1_4.csv')
write.csv(de[!is.na(de$padj) & de$padj<0.05,], 'snf_hvar_clusts/DE-cond2/diffExpRNA_1_4_FDR05.csv')
write.csv(de[!is.na(de$padj) & de$padj<0.1,], 'snf_hvar_clusts/DE-cond2/diffExpRNA_1_4_FDR1.csv')


geneList <- de$stat
keys <- bitr(de$RNA,
             fromType="SYMBOL",
             toType="ENTREZID",
             OrgDb="org.Hs.eg.db") %>%
    group_by(SYMBOL) %>%
    filter(row_number()==1)
names(geneList) <- unlist(as.vector(keys[,2]))
geneList <- sort(geneList, decreasing=TRUE)

set.seed(20220526)

gsea.go <- gseGO(geneList     = geneList,
                 OrgDb        = org.Hs.eg.db,
                 ont          = "BP",
                 minGSSize    = 50,
                 maxGSSize    = 500,
                 pvalueCutoff = 0.05,
                 verbose      = FALSE)

gsea.go <- simplify(gsea.go, cutoff=0.7, by='p.adjust', select_fun=min)

gsea.go

xx <- pairwise_termsim(gsea.go, showCategory=nrow(gsea.go))
p1 <- emapplot(xx, showCategory=nrow(gsea.go), max.overlaps=0)

p1

gsea.go.up <- gsea.go
gsea.go.up@result <- gsea.go.up@result[gsea.go$NES>0,]
gsea.go.up@geneSets <- gsea.go.up@geneSets[gsea.go$NES>0]

gsea.go.up

xx <- pairwise_termsim(gsea.go.up, showCategory=nrow(gsea.go.up))
p2 <- emapplot(xx, showCategory=nrow(gsea.go.up), max.overlaps=0, repel=TRUE, force=10, group_category=TRUE, node_label='group', nCluster=floor(2*sqrt(nrow(gsea.go.up))), nWords=5, clusterFunction=cluster::pam)

p2

gsea.go.down <- gsea.go
gsea.go.down@result <- gsea.go.down@result[gsea.go$NES<0,]
gsea.go.down@geneSets <- gsea.go.down@geneSets[gsea.go$NES<0]

gsea.go.down

xx <- pairwise_termsim(gsea.go.down, showCategory=nrow(gsea.go.down))
p3 <- emapplot(xx, showCategory=nrow(gsea.go.down), max.overlaps=0, repel=TRUE, force=10)

p3

gsea.go.df <- as.data.frame(gsea.go)

head(gsea.go.df[gsea.go.df$NES>0,1:10], 20)

head(gsea.go.df[gsea.go.df$NES<0,1:10], 20)

library(cowplot)
plot_grid(p2, p3, nrow=1, labels=c('Up', 'Down'))

ggsave('snf_hvar_clusts/DE-cond2/plots/cluster4_go_gsea.pdf',
       plot_grid(p2, p4, nrow=1, labels=c('Upregulated', 'Downregulated')),
       width=20, height=10, dpi=400)
ggsave('snf_hvar_clusts/DE-cond2/plots/cluster4_go_gsea.png',
       plot_grid(p2, p4, nrow=1, labels=c('Upregulated', 'Downregulated')),
       width=20, height=10, dpi=400)

write.csv(gsea.go.df, 'snf_hvar_clusts/DE-cond2/cluster4_go_gsea.csv')



aucVector <- function(data, clusters) {
    require(pROC)
    aucs <- rep(NA, ncol(data))
    names(aucs) <- colnames(data)

    for (feat in 1:ncol(data)) {
        aucs[feat] <- auc(clusters, data[,feat], direction='<')
    }

    aucs
}


#### AUC based GO:BP enrichment

## Cluster 1
geneList <- aucVector(rna[rownames(clin),], clin$hvar.Clusters==4)
keys <- bitr(de$RNA,
             fromType="SYMBOL",
             toType="ENTREZID",
             OrgDb="org.Hs.eg.db") %>%
    group_by(SYMBOL) %>%
    filter(row_number()==1)
names(geneList) <- unlist(as.vector(keys[,2]))
geneList <- sort(geneList, decreasing=TRUE)

set.seed(20220526)

gsea.go <- gseGO(geneList     = scale(geneList)[,1],
                 OrgDb        = org.Hs.eg.db,
                 ont          = "BP",
                 minGSSize    = 50,
                 maxGSSize    = 500,
                 pvalueCutoff = 0.05,
                 verbose      = FALSE)

gsea.go <- simplify(gsea.go, cutoff=0.7, by='p.adjust', select_fun=min)

gsea.go

xx <- pairwise_termsim(gsea.go, showCategory=nrow(gsea.go))
p1 <- emapplot(xx, showCategory=nrow(gsea.go), max.overlaps=0)

p1

gsea.go.up <- gsea.go
gsea.go.up@result <- gsea.go.up@result[gsea.go$NES>0,]
gsea.go.up@geneSets <- gsea.go.up@geneSets[gsea.go$NES>0]

gsea.go.up

xx <- pairwise_termsim(gsea.go.up, showCategory=nrow(gsea.go.up))
p2 <- emapplot(xx, showCategory=nrow(gsea.go.up), max.overlaps=0, repel=TRUE, force=10, group_category=TRUE, node_label='group', nCluster=floor(2*sqrt(nrow(gsea.go.up))), nWords=5, clusterFunction=cluster::pam)

p2

gsea.go.down <- gsea.go
gsea.go.down@result <- gsea.go.down@result[gsea.go$NES<0,]
gsea.go.down@geneSets <- gsea.go.down@geneSets[gsea.go$NES<0]

gsea.go.down

xx <- pairwise_termsim(gsea.go.down, showCategory=nrow(gsea.go.down))
p3 <- emapplot(xx, showCategory=nrow(gsea.go.down), max.overlaps=0, repel=TRUE, force=10)

p3

gsea.go.df <- as.data.frame(gsea.go)

head(gsea.go.df[gsea.go.df$NES>0,1:10], 20)

head(gsea.go.df[gsea.go.df$NES<0,1:10], 20)

library(cowplot)
plot_grid(p2, p3, nrow=1, labels=c('Up', 'Down'))

ggsave('snf_hvar_clusts/DE-cond2/plots/cluster4_go_gsea.pdf',
       plot_grid(p2, p4, nrow=1, labels=c('Upregulated', 'Downregulated')),
       width=20, height=10, dpi=400)
ggsave('snf_hvar_clusts/DE-cond2/plots/cluster4_go_gsea.png',
       plot_grid(p2, p4, nrow=1, labels=c('Upregulated', 'Downregulated')),
       width=20, height=10, dpi=400)

write.csv(gsea.go.df, 'snf_hvar_clusts/DE-cond2/cluster4_go_gsea.csv')





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

data.views <- list()
data.views[['prot']] <- prot
data.views[['metab']] <- metab
data.views[['rna']] <- rna
data.views[['mirna']] <- mirna

concat <- clin %>% dplyr::select(hvar.Clusters, WksGest, InfSex, Race, PrePregBMI, Smoking,
                          FDELTYPE, Labor.initiation,
                          `Syncytial knots - (Tony: considered with accelerated villous maturation)`,
                          `Accelerated villous maturity`,
                          `Distal villous hypoplasia focal/diffuse`)

concat <- concat %>% mutate_at(c('hvar.Clusters', 'InfSex', 'Race', 'Smoking', 'FDELTYPE',
                                 'Labor.initiation',
                                 'Syncytial knots - (Tony: considered with accelerated villous maturation)',
                                 'Accelerated villous maturity',
                                 'Distal villous hypoplasia focal/diffuse'), factor)

for (modal in c('prot', 'metab', 'rna', 'mirna')) {
    markers <- read.csv(paste('snf_hvar_clusts/markers/auc', modal, 'markers.csv', sep='_'))
    feats <- unique(markers$Feature)
    print(feats)

    concat <- cbind(concat, data.views[[modal]][rownames(concat),feats])
}

## library(rCausalMGM)

colnames(concat) <- make.names(colnames(concat))

dim(concat)
ig.path <- mgmPath(concat, lambdas=0.8*(0.01)^((0:99)/99), rank=F, verbose=T)

## plot(log10(ig.path$lambdas), ig.path$AIC)
## points(log10(ig.path$lambdas), ig.path$BIC, col='red')

ig.bic <- ig.path$graphs[[which.min(ig.path$BIC)]]
ig.bic
ig.bic$markov.blankets[['hvar.Clusters']]


## ig.aic <- ig.path$graphs[[which.min(ig.path$AIC)]]
## ig.aic
## ig.aic$markov.blankets[['hvar.Clusters']]


g <- fciMax(concat, initialGraph=ig.bic, alpha=0.05, fdr=T, rank=F, verbose=T)

g$edges

saveGraph(g, 'snf_hvar_clusts/graph/snfClust_markers_mgmfcimax_fdr05.txt')
saveGraph(g, 'snf_hvar_clusts/graph/snfClust_markers_mgmfcimax_fdr05.sif')


g.boot <- bootstrap(concat,
                    algorithm='mgmfcimax',
                    lambda=ig.bic$lambda, alpha=0.05,
                    numBoots=100, verbose=T)

g.boot
g.boot$edges
g.boot$stabilities


library(stringr)

g.stabs <- data.frame()

count <- 1
for (edge in str_split(g$edges, ' ')) {
    
    for (i in 1:nrow(g.boot$stabilities)) {
        if ((edge[1] == g.boot$stabilities[i,1] & edge[3] == g.boot$stabilities[i,3]) |
            (edge[3] == g.boot$stabilities[i,1] & edge[1] == g.boot$stabilities[i,3])) {
            ## print(edge)
            ## print(g.boot$stabilities[i,])
            g.stabs <- rbind(g.stabs, g.boot$stabilities[i,])

            g.stabs[count,2] <- edge[2]

            if (edge[1] > edge[3]) {
                g.stabs[count,] <- g.stabs[count,c(3,2,1,4,6,5,7,9,8,10,11)]
            }
            
            break
        }
    }

    count <- count + 1
}

g.stabs

write.csv(g.stabs,
          'snf_hvar_clusts/graph/snfClust_marker_mgmfcimax_fdr05.csv',
          row.names=FALSE, quote=FALSE)




comparisons <- list(c('Control-PTD', 'Control'),
                    c('FGR', 'Control'),
                    c('Severe-PE', 'Control.PTD'),
                    c('Severe-PE', 'Control'),
                    c('PTD', 'Control.PTD'),
                    c('PTD', 'Control'),
                    c('FGR+HDP', 'Control.PTD'),
                    c('FGR+HDP', 'Control'))

de.feats <- list()
de.feats[['metab']] <- c()
de.feats[['prot']] <- c()
de.feats[['rna']] <- c()
de.feats[['mirna']] <- c()


for (groups in comparisons) {
    de.metab <- read.csv(paste('DE-cond3/diffExpMetab_', groups[1], '_', groups[2], '_FDR05.csv', sep=''))

    de.prot <- read.csv(paste('DE-cond3/diffExpProt_', groups[1], '_', groups[2], '_FDR05.csv', sep=''))

    ## print(de.prot)

    ## print(rownames(de.prot))
    
    de.rna <- read.csv(paste('DE-cond3/diffExpRNA_', make.names(groups[1]), '_', make.names(groups[2]), '_FDR05.csv', sep=''))

    de.mirna <- read.csv(paste('DE-cond3/diffExpMiRNA_', make.names(groups[1]), '_', make.names(groups[2]), '_FDR05.csv', sep=''))

    de.feats[['metab']] <- c(de.feats[['metab']],
                             paste0('Metab_', make.names(de.metab$CHEMICAL_NAME)))
    de.feats[['prot']] <- c(de.feats[['prot']],
                            paste0('Prot_', make.names(de.prot$Assay)))
    de.feats[['rna']] <- c(de.feats[['rna']],
                           paste0('RNA_', make.names(de.rna$Gene)))
    de.feats[['mirna']] <- c(de.feats[['mirna']],
                             paste0('miRNA_', make.names(de.mirna$miRNA)))
    
    ## de.feats <- c(de.feats, make.names(de.metab$CHEMICAL_NAME), make.names(rownames(de.prot)))
}



de.feats[['metab']] <- setdiff(unique(de.feats[['metab']]), 'Metab_')

length(de.feats[['metab']])

de.feats[['prot']] <- setdiff(unique(de.feats[['prot']]), 'Prot_')

length(de.feats[['prot']])

de.feats[['rna']] <- setdiff(unique(de.feats[['rna']]), 'RNA_')

length(de.feats[['rna']])

de.feats[['mirna']] <- setdiff(unique(de.feats[['mirna']]), 'miRNA_')

length(de.feats[['mirna']])



library(ggvenn)

de.feats <- list()
de.feats[['metab']] <- list()
de.feats[['prot']] <- list()
de.feats[['rna']] <- list()
de.feats[['mirna']] <- list()

venn.diag.tests <- data.frame(Condition=c(NA),
                              Dataset=c(NA),
                              Total=c(NA),
                              Control=c(NA),
                              `Control PTD`=c(NA),
                              Common=c(NA),
                              pval=c(NA),
                              padj=c(NA))

for (groups in comparisons[3:8]) {

    groups[2] <- gsub(' ', '.', groups[2])
    groups[1] <- gsub(' ', '-', groups[1])
    
    de.metab <- read.csv(paste('DE-cond3/diffExpMetab_', groups[1], '_', groups[2], '_FDR05.csv', sep=''))

    de.prot <- read.csv(paste('DE-cond3/diffExpProt_', groups[1], '_', groups[2], '_FDR05.csv', sep=''))

    print(groups)
    groups <- gsub('[-]', ' ', groups)
    groups <- gsub('[.]', ' ', groups)
    print(groups)
    ## print(de.prot)

    ## print(rownames(de.prot))
    
    de.rna <- read.csv(paste('DE-cond3/diffExpRNA_', make.names(groups[1]), '_', make.names(groups[2]), '_FDR05.csv', sep=''))

    de.mirna <- read.csv(paste('DE-cond3/diffExpMiRNA_', make.names(groups[1]), '_', make.names(groups[2]), '_FDR05.csv', sep=''))

    de.feats[['metab']][[groups[2]]] <- de.metab$CHEMICAL_NAME
    de.feats[['prot']][[groups[2]]] <- de.prot$Assay
    de.feats[['rna']][[groups[2]]] <- de.rna$Gene
    de.feats[['mirna']][[groups[2]]] <- de.mirna$miRNA
    
    ## de.feats <- c(de.feats, make.names(de.metab$CHEMICAL_NAME), make.names(rownames(de.prot)))

    if (groups[2] == 'Control') {

        common <- length(intersect(de.feats[['metab']][['Control']],
                                   de.feats[['metab']][['Control PTD']]))
        de.control <- length(de.feats[['metab']][['Control']])
        de.controlPTD <- length(de.feats[['metab']][['Control PTD']])
        total <- ncol(metab)

        pval <- phyper(common-1, de.control, total-de.control, de.controlPTD, lower.tail=F)

        test.res <- c(Condition=groups[1], Dataset='Metabolite', Total=total,
                      Control=de.control, `Control PTD`=de.controlPTD,
                      Common=common, pval=pval, padj=NA)

        venn.diag.tests <- rbind(venn.diag.tests, test.res)
        
        venn <- ggvenn(de.feats[['metab']],
                       fill_color=rev(brewer.colors[1:2])) +
            ggtitle(paste0('DE Metabolites: ', groups[1], ' vs. Controls'))
        ggsave(paste0('common-DE-analytes/common-DE-metabs-',
                      make.names(groups[1]), '.pdf'), venn,
               width=5, height=5)
        ggsave(paste0('common-DE-analytes/common-DE-metabs-',
                      make.names(groups[1]), '.png'), venn,
               width=5, height=5, dpi=400, bg='white')

        common <- length(intersect(de.feats[['prot']][['Control']],
                                   de.feats[['prot']][['Control PTD']]))
        de.control <- length(de.feats[['prot']][['Control']])
        de.controlPTD <- length(de.feats[['prot']][['Control PTD']])
        total <- ncol(prot)

        pval <- phyper(common-1, de.control, total-de.control, de.controlPTD, lower.tail=F)

        test.res <- c(Condition=groups[1], Dataset='Protein', Total=total,
                      Control=de.control, `Control PTD`=de.controlPTD,
                      Common=common, pval=pval, padj=NA)

        venn.diag.tests <- rbind(venn.diag.tests, test.res)
        
        venn <- ggvenn(de.feats[['prot']],
                       fill_color=rev(brewer.colors[1:2])) +
            ggtitle(paste0('DE Proteins: ', groups[1], ' vs. Controls'))
        ggsave(paste0('common-DE-analytes/common-DE-prots-',
                      make.names(groups[1]), '.pdf'), venn,
               width=5, height=5)
        ggsave(paste0('common-DE-analytes/common-DE-prots-',
                      make.names(groups[1]), '.png'), venn,
               width=5, height=5, dpi=400, bg='white')

        common <- length(intersect(de.feats[['rna']][['Control']],
                                   de.feats[['rna']][['Control PTD']]))
        de.control <- length(de.feats[['rna']][['Control']])
        de.controlPTD <- length(de.feats[['rna']][['Control PTD']])
        total <- ncol(rna)

        pval <- phyper(common-1, de.control, total-de.control, de.controlPTD, lower.tail=F)

        test.res <- c(Condition=groups[1], Dataset='RNA', Total=total,
                      Control=de.control, `Control PTD`=de.controlPTD,
                      Common=common, pval=pval, padj=NA)

        venn.diag.tests <- rbind(venn.diag.tests, test.res)
        

        venn <- ggvenn(de.feats[['rna']],
                       fill_color=rev(brewer.colors[1:2])) +
            ggtitle(paste0('DE Genes: ', groups[1], ' vs. Controls'))
        ggsave(paste0('common-DE-analytes/common-DE-RNA-',
                      make.names(groups[1]), '.pdf'), venn,
               width=5, height=5)
        ggsave(paste0('common-DE-analytes/common-DE-RNA-',
                      make.names(groups[1]), '.png'), venn,
               width=5, height=5, dpi=400, bg='white')

        common <- length(intersect(de.feats[['mirna']][['Control']],
                                   de.feats[['mirna']][['Control PTD']]))
        de.control <- length(de.feats[['mirna']][['Control']])
        de.controlPTD <- length(de.feats[['mirna']][['Control PTD']])
        total <- ncol(mirna)

        pval <- phyper(common-1, de.control, total-de.control, de.controlPTD, lower.tail=F)

        test.res <- c(Condition=groups[1], Dataset='miRNA', Total=total,
                      Control=de.control, `Control PTD`=de.controlPTD,
                      Common=common, pval=pval, padj=NA)

        venn.diag.tests <- rbind(venn.diag.tests, test.res)
        
        
        venn <- ggvenn(de.feats[['mirna']],
                       fill_color=rev(brewer.colors[1:2])) +
            ggtitle(paste0('DE miRNA: ', groups[1], ' vs. Controls'))
        ggsave(paste0('common-DE-analytes/common-DE-miRNA-',
                      make.names(groups[1]), '.pdf'), venn,
               width=5, height=5)
        ggsave(paste0('common-DE-analytes/common-DE-miRNA-',
                      make.names(groups[1]), '.png'), venn,
               width=5, height=5, dpi=400, bg='white')
        
    }
}

venn.diag.tests <- venn.diag.tests[-1,]

venn.diag.tests$padj <- p.adjust(venn.diag.tests$pval, method='holm')

write.csv(venn.diag.tests, 'common-DE-analytes.csv', row.names=F)
