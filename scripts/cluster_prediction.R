library(limma)
library(DESeq2)
library(SNFtool)
library(dplyr)
library(stringr)
library(readxl)
library(glmnet)
library(pROC)
library(caret)
source('cv_utils.R')

tanimoto <- function(a, b) {
    intersection <- length(intersect(a, b))
    num <- length(a) + length(b) - 2 * intersection
    denom <- length(a) + length(b) - intersection
    return ( 1 - num / denom)
}

jaccard <- function(a, b) {
    intersection <- length(intersect(a, b))
    union <- length(a) + length(b) - intersection
    return ( intersection / union)
}


permLabels <- function(true, new.labels) {
    require(gtools)
    require(caret)
    
    K <- length(unique(true))
    if (K > 8) {
        stop("Too many classes for brute force permutation")
    }
    permMat <- permutations(K, K, v=sort(unique(true)))

    bestPerm <- 1
    bestDiagSum <- 0
    for (i in 1:nrow(permMat)) {
        temp.labels <- new.labels
        for (j in 1:ncol(permMat)) {
            temp.labels[new.labels==j] <- permMat[i,j]
        }
        confusion <- table(true, temp.labels)
        ## print(confusion)
        if (sum(diag(confusion)) > bestDiagSum) {
            bestDiagSum <- sum(diag(confusion))
            bestPerm <- i
        }
    }
    
    for (j in 1:ncol(permMat)) {
        temp.labels[new.labels==j] <- permMat[bestPerm,j]
    }

    ## print(confusionMatrix(factor(true, levels=sort(unique(true))),
    ##                       factor(temp.labels, levels=sort(unique(true)))))
    
    temp.labels
}


#### Clinical data

clin <- read_excel('Primary clinical variables doppler.xlsx', 1) %>% as.data.frame
colnames(clin) <- make.names(colnames(clin))
clin <- clin[-1,]
rownames(clin) <- make.names(clin$Study_ID.ID)
head(clin)

clin[incomp.samp.ids,'pred.hvar.Clusters'] <- pred.final

write.csv(clin, 'doppler_clinical.csv')

clin <- clin %>% mutate_if(function(x) mean(is.na(as.numeric(x))) < 0.5, as.numeric)

clin$Condition. <- ifelse(clin$Condition. %in% c('PTD','Severe PE','FGR','Control',
                                                 'Control PTD'), clin$Condition., 'FGR+HDP')

clin$PrePregBMI <- clin$PrePregWt_Kg / clin$Height_Meters.HeightMeters^2

clin$Smoking <- ifelse(str_detect(clin$Smoke, 'Yes'), 'Yes', NA)
clin$Smoking <- ifelse(str_detect(clin$Smoke, 'Never'), 'No', clin$Smoking)

clin$Race[!clin$Race %in% c('W', 'B')] <- 'other'
clin$FDELTYPE[clin$FDELTYPE==0] <- NA

clin <- clin %>% select(c('Condition.', 'WksGest', 'InfSex', 'Race', 'PrePregBMI',
                          'Smoking', 'FDELTYPE', 'Labor.initiation'))

clin <- clin %>% mutate_at(c('Condition.', 'InfSex', 'Race', 'Smoking', 'FDELTYPE',
                             'Labor.initiation'), factor)

#### Pathology data (slides)

pathology.slides <- read_excel('path features (from slides) for analysis.xlsx' , 1, na='n/a') %>% as.data.frame
colnames(pathology.slides) <- make.names(colnames(pathology.slides))
colnames(pathology.slides)[6] <- "Syncytial.knots"
rownames(pathology.slides) <- make.names(pathology.slides$Study.ID)
pathology.slides <- pathology.slides[,4:6]
pathology.slides <- pathology.slides %>% mutate_all(factor)

#### Multi-omics data

prot <- read.csv('data/combat_prot_data.csv', row.names=1)
colnames(prot) <- paste0('Prot_', colnames(prot))
## head(prot)

metab <- read.csv('data/combat_metab_data.csv', row.names=1)
colnames(metab) <- paste0('Metab_', colnames(metab))
## head(metab)

rna <- read.csv('data/combat_vst_rna_expression.csv', row.names=1)
colnames(rna) <- paste0('RNA_', colnames(rna))
## head(rna)

rna.counts <- read.csv('data/combat_rna_counts.csv', row.names=1)
colnames(rna.counts) <- paste0('RNA_', colnames(rna.counts))
## head(rna)

mirna <- read.csv('data/combat_vst_mirna_expression.csv', row.names=1)
colnames(mirna) <- paste0('miRNA_', colnames(mirna))
## head(mirna)

mirna.counts <- read.csv('data/combat_mirna_counts.csv', row.names=1)
colnames(mirna.counts) <- paste0('miRNA_', colnames(mirna.counts))


full.samp.ids <- unique(c(rownames(prot), rownames(metab), rownames(rna), rownames(mirna)))
clin <- clin[full.samp.ids,]
pathology.slides <- pathology.slides[full.samp.ids,]

samp.ids <- intersect(rownames(prot), rownames(metab))
samp.ids <- intersect(samp.ids, rownames(rna))
samp.ids <- intersect(samp.ids, rownames(mirna))
samp.ids <- intersect(samp.ids, rownames(pathology.slides))

replicates <- 10
K <- 10
C <- 4
neighbors <- 20
bandwidth <- 0.3

set.seed(20220714)

foldids <- list()


labels <- array(NA, c(length(samp.ids), replicates, K))
dimnames(labels)[[1]] <- samp.ids

concat.train <- list()
concat.test <- list()

clin.train <- list()
clin.test <- list()

dataL.train <- list()
dataL.test <- list()

alphas <- seq(0.5, 1, 0.05)
lambdas <- exp(seq(log(0.45), log(0.045), length.out=50))

prob.mat <- array(NA, c(length(samp.ids), C, length(lambdas), length(alphas), replicates))
pred.mat <- array(NA, c(length(samp.ids), length(lambdas), length(alphas), replicates))

score.auc <- array(NA, c(replicates, K, length(lambdas), length(alphas)))
score.bacc <- array(NA, c(replicates, K, length(lambdas), length(alphas)))

for (rep in 1:replicates) {

    print(paste0("Replicate ", rep))
    
    foldids[[rep]] <- stratifiedFolds(clin$Condition., K)

    names(foldids[[rep]]) <- full.samp.ids

    foldid <- foldids[[rep]]

    hvar.prots <- list()
    hvar.metabs <- list()
    hvar.rnas <- list()
    hvar.mirnas <- list()

    for (k in 1:K) {
        prot.vars <- apply(prot[rownames(prot)[foldid[rownames(prot)]!=k],], 2, stats::var)
        metab.vars <- apply(metab[rownames(metab)[foldid[rownames(metab)]!=k],],2,stats::var)
        rna.vars <- apply(rna[rownames(rna)[foldid[rownames(rna)]!=k],], 2, stats::var)
        mirna.vars <- apply(mirna[rownames(mirna)[foldid[rownames(mirna)]!=k],],2,stats::var)

        hvar.prots[[k]] <- colnames(prot)[which(prot.vars >= quantile(prot.vars, 0.75))]
        hvar.metabs[[k]] <- colnames(metab)[which(metab.vars >= quantile(metab.vars, 0.75))]
        hvar.rnas[[k]] <- colnames(rna)[which(rna.vars >= quantile(rna.vars, 0.75))]
        hvar.mirnas[[k]] <- colnames(mirna)[which(mirna.vars >= quantile(mirna.vars, 0.75))]
    }

    jac.prot <- matrix(NA, K, K)
    jac.metab <- matrix(NA, K, K)
    jac.mirna <- matrix(NA, K, K)
    jac.rna <- matrix(NA, K, K)

    for (ki in 1:K) {
        for (kj in 1:K) {
            jac.prot[ki,kj] <- jaccard(hvar.prots[[ki]], hvar.prots[[kj]])
            jac.metab[ki,kj] <- jaccard(hvar.metabs[[ki]], hvar.metabs[[kj]])
            jac.mirna[ki,kj] <- jaccard(hvar.mirnas[[ki]], hvar.mirnas[[kj]])
            jac.rna[ki,kj] <- jaccard(hvar.rnas[[ki]], hvar.rnas[[kj]])
        }
    }

    ## jac.prot
    print(mean(jac.prot[upper.tri(jac.prot)]))
    ## jac.metab
    print(mean(jac.metab[upper.tri(jac.metab)]))
    ## jac.rna
    print(mean(jac.rna[upper.tri(jac.rna)]))
    ## jac.mirna
    print(mean(jac.mirna[upper.tri(jac.mirna)]))

    concat.train[[rep]] <- list()
    concat.test[[rep]] <- list()

    clin.train[[rep]] <- list()
    clin.test[[rep]] <- list()

    dataL.train[[rep]] <- list()
    dataL.test[[rep]] <- list()

    for (k in 1:K) {

        dataL.train[[rep]][[k]] <- list()
        dataL.test[[rep]][[k]] <- list()
        
        dataL.train[[rep]][[k]][['prot']] <- prot[samp.ids[foldid[samp.ids]!=k],
                                                  hvar.prots[[k]]]
        dataL.train[[rep]][[k]][['metab']] <- metab[samp.ids[foldid[samp.ids]!=k],
                                                    hvar.metabs[[k]]]
        dataL.train[[rep]][[k]][['rna']] <- rna[samp.ids[foldid[samp.ids]!=k],
                                                hvar.rnas[[k]]]
        dataL.train[[rep]][[k]][['mirna']] <- mirna[samp.ids[foldid[samp.ids]!=k],
                                                     hvar.mirnas[[k]]]

        clin.train[[rep]][[k]] <- cbind(clin[samp.ids[foldid[samp.ids]!=k],-1],
                                        pathology.slides[samp.ids[foldid[samp.ids]!=k],])
        
        dataL.test[[rep]][[k]][['prot']] <- prot[samp.ids[foldid[samp.ids]==k],
                                                 hvar.prots[[k]]]
        dataL.test[[rep]][[k]][['metab']] <- metab[samp.ids[foldid[samp.ids]==k],
                                                   hvar.metabs[[k]]]
        dataL.test[[rep]][[k]][['rna']] <- rna[samp.ids[foldid[samp.ids]==k],
                                               hvar.rnas[[k]]]
        dataL.test[[rep]][[k]][['mirna']] <- mirna[samp.ids[foldid[samp.ids]==k],
                                                   hvar.mirnas[[k]]]

        clin.test[[rep]][[k]] <- cbind(clin[samp.ids[foldid[samp.ids]==k],-1],
                                       pathology.slides[samp.ids[foldid[samp.ids]==k],])

        for (idx in 1:4) {
            means <- apply(dataL.train[[rep]][[k]][[idx]], 2, mean)
            sds <- apply(dataL.train[[rep]][[k]][[idx]], 2, sd)

            for (feat in 1:ncol(dataL.train[[rep]][[k]][[idx]])) {
                dataL.train[[rep]][[k]][[idx]][,feat] <- (dataL.train[[rep]][[k]][[idx]][,feat] -
                                                          means[feat]) / sds[feat]
                dataL.test[[rep]][[k]][[idx]][,feat] <- (dataL.test[[rep]][[k]][[idx]][,feat] -
                                                         means[feat]) / sds[feat]
            }
        }

        concat.train[[rep]][[k]] <- cbind(dataL.train[[rep]][[k]][['prot']],
                                          dataL.train[[rep]][[k]][['metab']],
                                          dataL.train[[rep]][[k]][['rna']],
                                          dataL.train[[rep]][[k]][['mirna']],
                                          clin.train[[rep]][[k]])

        concat.test[[rep]][[k]] <- cbind(dataL.test[[rep]][[k]][['prot']],
                                         dataL.test[[rep]][[k]][['metab']],
                                         dataL.test[[rep]][[k]][['rna']],
                                         dataL.test[[rep]][[k]][['mirna']],
                                         clin.test[[rep]][[k]])

        snf.out <- SNFcluster(dataL.train[[rep]][[k]], dataL.test[[rep]][[k]],
                              C, neighbors, bandwidth)

        labels[samp.ids[foldid[samp.ids]!=k],rep,k] <- snf.out[['train.labels']]
        labels[samp.ids[foldid[samp.ids]==k],rep,k] <- snf.out[['test.labels']]
    }


    ami.labels <- matrix(NA, K, K)

    for (ki in 1:K) {
        for (kj in 1:K) {
            ami.labels[ki,kj] <- aricode::AMI(labels[,rep,ki], labels[,rep,kj])
        }
    }

    ## ami.labels
    print(paste0("Replicate ", rep, " Cluster Label AMI:"))
    print(mean(ami.labels[upper.tri(ami.labels)]))


    for (k in 1:K) {
        print(paste0("Fold ", k))
        idx <- which(foldid[samp.ids]==k)
        ## class.wts <- nrow(concat.train[[rep]][[k]]) /
        ##     (C * table(labels[-idx,rep,k]))
        ## case.wts <- class.wts[labels[-idx,rep,k]]
        case.wts <- rep(1, sum(foldid[samp.ids]!=k))

        for (a.idx in 1:length(alphas)) {
            
            res <- glmnet(model.matrix(~.,concat.train[[rep]][[k]])[,-1],
                          labels[-idx,rep,k],
                          lambda=lambdas/alphas[a.idx], alpha=alphas[a.idx],
                          family='multinomial', weights=as.vector(case.wts),
                          type.multinomial = 'grouped')

            prob.mat[idx,,,a.idx,rep] <- predict(res,
                                                 newx=model.matrix(~.,concat.test[[rep]][[k]])[,-1],
                                                 type='response')

            pred.mat[idx,,a.idx,rep] <- predict(res,
                                                newx=model.matrix(~.,concat.test[[rep]][[k]])[,-1],
                                                type='class')

            for (l.idx in 1:length(lambdas)) {

                probs <- prob.mat[idx,,l.idx,a.idx,rep]
                ## rownames(probs) <- samp.ids[idx]
                colnames(probs) <- c(1,2,3,4)
                
                roc.res <- multiclass.roc(labels[idx,rep,k], probs)

                print(roc.res$auc)

                score.auc[rep,k,l.idx,a.idx] <- as.numeric(roc.res$auc)

                cmat <- confusionMatrix(factor(pred.mat[idx,l.idx,a.idx,rep], levels=1:4),
                                        factor(labels[idx,rep,k], levels=1:4))

                print(cmat$table)

                score.bacc[rep,k,l.idx,a.idx] <- mean(cmat$byClass[,'Balanced Accuracy'], na.rm=TRUE)
            }

        }
    }
}

## hist(labels[,1])
## hist(labels[,2])
## hist(labels[,3])
## hist(labels[,4])
## hist(labels[,5])
## hist(labels[,6])
## hist(labels[,7])
## hist(labels[,8])
## hist(labels[,9])
## hist(labels[,10])


## prob.mat <- matrix(NA, nrow(labels), 3)
## rownames(prob.mat) <- rownames(labels)


cv.mean.bacc <- apply(score.bacc, c(3,4), median)
cv.sd.bacc <- apply(score.bacc, c(3,4), sd)

which(cv.mean.bacc == max(cv.mean.bacc), arr.ind=T)

which(cv.mean.bacc==min(cv.mean.bacc[which(cv.mean.bacc > min(max(cv.mean.bacc)-cv.sd.bacc[which(cv.mean.bacc == max(cv.mean.bacc), arr.ind=T)]/sqrt(replicates*K)), arr.ind=T)]), arr.ind=T)

## ungrouped: l.idx = 48, a.idx = 5, AUC = 0.9620043, BACC = 0.8770507
## grouped: l.idx = 46, a.idx = 7, AUC = 0.9656175, BACC = 0.8817613

cv.mean.auc <- apply(score.auc, c(3,4), mean)
cv.sd.auc <- apply(score.auc, c(3,4), sd)

which(cv.mean.auc == max(cv.mean.auc), arr.ind=T)

which(cv.mean.auc==min(cv.mean.auc[which(cv.mean.auc > min(max(cv.mean.auc)-cv.sd.auc[which(cv.mean.auc == max(cv.mean.auc), arr.ind=T)]/sqrt(replicates*K)), arr.ind=T)]), arr.ind=T)

## ungrouped: l.idx = 45, a.idx = 6, AUC = 0.9608601, BACC = 0.8758518
## grouped: l.idx = 45, a.idx = 6, AUC = 0.9652327, BACC = 0.8775397

aucMat <- matrix(NA, replicates*K, 5)
colnames(aucMat) <- 1:5

## combMat <- gtools::combinations(4, 2)

## for (comb in 1:nrow(combMat)) {
##     colnames(aucMat)[comb] <- paste(combMat[comb,], collapse='/')
## }

enet.feats <- list()

for (rep in 1:replicates) {
    foldid <- foldids[[rep]]
    enet.feats[[rep]] <- list()
    for (k in 1:K) {

        labels[samp.ids, rep, k] <- permLabels(clin.clusters[samp.ids, 'hvar.Clusters'], labels[samp.ids, rep, k])
        
        class.wts <- nrow(concat.train[[rep]][[k]]) /
            (4 * table(labels[samp.ids[foldid[samp.ids]!=k],rep,k]))
        case.wts <- class.wts[labels[samp.ids[foldid[samp.ids]!=k],rep,k]]
        ## case.wts <- (case.wts + rep(1, nrow(concat.train[[k]]))) / 2
        ## case.wts <- rep(1, nrow(concat.train[[rep]][[k]]))

        res <- glmnet(model.matrix(~.,concat.train[[rep]][[k]])[,-1],
                      labels[samp.ids[foldid[samp.ids]!=k],rep,k],
                      lambda=lambdas[l.idx]/alphas[a.idx], alpha=alphas[a.idx],
                      family='multinomial', weights=as.vector(case.wts),
                      type.multinomial = 'grouped')
        
        beta <- coef(res)
        print("Cluster 1")
        print(beta[[1]]@Dimnames[[1]][as.matrix(beta[[1]])!=0][-1])
        print(beta[[1]][as.matrix(beta[[1]])!=0][-1])
        print("Cluster 2")
        print(beta[[2]]@Dimnames[[1]][as.matrix(beta[[2]])!=0][-1])
        print(beta[[2]][as.matrix(beta[[2]])!=0][-1])
        print("Cluster 3")
        print(beta[[3]]@Dimnames[[1]][as.matrix(beta[[3]])!=0][-1])
        print(beta[[3]][as.matrix(beta[[3]])!=0][-1])
        print("Cluster 4")
        print(beta[[4]]@Dimnames[[1]][as.matrix(beta[[4]])!=0][-1])
        print(beta[[4]][as.matrix(beta[[4]])!=0][-1])

        enet.feats[[rep]][[k]] <- unique(c(beta[[1]]@Dimnames[[1]][as.matrix(beta[[1]])!=0][-1],
                                           beta[[2]]@Dimnames[[1]][as.matrix(beta[[2]])!=0][-1],
                                           beta[[3]]@Dimnames[[1]][as.matrix(beta[[3]])!=0][-1],
                                           beta[[4]]@Dimnames[[1]][as.matrix(beta[[4]])!=0][-1]))
        
        probs <- predict(res, newx=model.matrix(~.,concat.test[[rep]][[k]])[,-1], type='response')[,,1]

        roc.res <- multiclass.roc(factor(labels[samp.ids[foldid[samp.ids]==k],rep,k], levels=1:4), probs)

        ## combMat <- combinations(4, 2)

        ## for (comb in 1:nrow(combMat)) {
        ##     if (!is.null(roc.res$rocs[[paste(combMat[comb,], collapse='/')]])) {
        ##         aucMat[20*(rep-1) + 2 * (k-1) + 1:2, paste(combMat[comb,], collapse='/')] <- sapply(roc.res$rocs[[paste(combMat[comb,], collapse='/')]], auc)
        ##     }
        ## }

        for (clust in 1:4) {
            ovr.labels <- factor(labels[samp.ids[foldid[samp.ids]==k],rep,k]==clust, levels=c(TRUE, FALSE))
            if (!all(ovr.labels==FALSE)) {
                roc.ovr <- roc(ovr.labels, probs[,clust])
                print(roc.ovr$auc)
                aucMat[10*(rep-1) + k, clust] <- roc.ovr$auc
            } else {
                aucMat[10*(rep-1) + k, clust] <- NA
            }
        }

        print(roc.res$auc)

        aucMat[10*(rep-1) + k, 5] <- roc.res$auc

        ## print(BrierScore(model.matrix(~0+.,as.data.frame(factor(labels[samp.ids[foldid[samp.ids]==k],k]))), probs))

        print(confusionMatrix(factor(predict(res,
                                             newx=model.matrix(~.,concat.test[[rep]][[k]])[,-1],
                                             type='class')[,1], levels=1:4),
                              factor(labels[samp.ids[foldid[samp.ids]==k],rep,k], levels=1:4)))
        
        print("")
        print("")

        ## prob.mat[samp.ids[foldid[samp.ids]==k],] <- probs
    }
}

auc.plot.df <- data.frame(row.names=1:(5 * K * replicates))
colnames(aucMat) <- c(paste('Cluster', 1:4), 'Multiclass')
auc.plot.df$Comparison <- rep(colnames(aucMat), each=(K * replicates))
auc.plot.df$AUC <- c(aucMat)

## idx <- 1
## for ( rep in 1:replicates ) {
##     for (k in 1:K) {
##         auc.plot.df[(2 * 6 * K * replicates)+idx,] <- c('Multiclass', mean(aucMat[20*(rep-1) + 2 * (k-1) + 1:2,], na.rm=T))
##         idx <- idx+1
##     }
## }

auc.plot.df$Comparison <- factor(auc.plot.df$Comparison)
auc.plot.df$AUC <- as.numeric(auc.plot.df$AUC)

auc.summaries <- auc.plot.df %>% group_by(Comparison) %>%
    summarize_all(c(mean=mean, sd=sd,
                    n=function(x, ...) sum(!is.na(x)),
                    lCIquant=function(x, ...) quantile(x, probs=c(0.25), ...),
                    median=function(x, ...) quantile(x, probs=c(0.75), ...),
                    hCIquant=function(x, ...) quantile(x, probs=c(0.95), ...),
                    lCIgauss=function(x, ...) mean(x, ...) + qnorm(0.25) * sd(x, ...),
                    hCIgauss=function(x, ...) mean(x, ...) + qnorm(0.75) * sd(x, ...)),
                  na.rm=T)

auc.summaries

gg <- ggplot(auc.plot.df, aes(Comparison, AUC, fill=Comparison)) +
    geom_violin() +
    geom_hline(yintercept=0.5, lty=5, col='red') + 
    ## geom_boxplot(width=0.1) +
    scale_fill_manual(values=c(brewer.colors[7:10], '#fed9a6')) +
    theme_bw()

gg

ggsave('enet_prediction_aucs_ovr_wBaseline_v2.png', gg, width=5, height=5, dpi=400)
ggsave('enet_prediction_aucs_ovr_wBaseline_v2.pdf', gg, width=5, height=5)


gg <- ggplot(auc.plot.df, aes(Comparison, AUC, fill=Comparison)) +
    geom_violin() +
    ## geom_hline(yintercept=0.5, lty=5, col='red') + 
    ## geom_boxplot(width=0.1) +
    scale_fill_manual(values=c(brewer.colors[7:10], '#fed9a6')) +
    theme_bw()

gg

ggsave('enet_prediction_aucs_ovr_woBaseline.png', gg, width=5, height=5, dpi=400)
ggsave('enet_prediction_aucs_ovr_woBaseline.pdf', gg, width=5, height=5)


## enet.feat.union <- c()
## for (k in 1:K) {
##     ener.feat.union <- c(enet.feat.union, enet.feats[[k]])
## }
enet.feat.union <- unique(unlist(enet.feats))
enet.feat.union

enet.feat.counts <- matrix(0, length(enet.feat.union), 1)
rownames(enet.feat.counts) <- enet.feat.union

for (rep in 1:replicates) {
    for (k in 1:K) {
        enet.feat.counts[enet.feats[[rep]][[k]],] <- enet.feat.counts[enet.feats[[rep]][[k]],] + 1
    }
}

enet.feat.counts <- enet.feat.counts[order(enet.feat.counts, decreasing=T),]

enet.feat.counts[enet.feat.counts >= floor(0.5 * replicates * K)]

feat.bar.plot <- data.frame(Count=enet.feat.counts[enet.feat.counts >= floor(0.5 * replicates * K)])

feat.bar.plot$Feature <- rownames(feat.bar.plot)

metab.meta <- read_excel('UPIT-02-20PHML+ DATA TABLES OB.xlsx', 2) %>% as.data.frame
rownames(metab.meta) <- make.names(metab.meta$CHEMICAL_NAME)
head(metab.meta)

feat.bar.plot$Dataset <- "Clinical"

idx <- 1
for (strlist in str_split(rownames(feat.bar.plot),'_')) {
    if (length(strlist) > 1) {
        if (strlist[1] == 'Metab') {
            feat.bar.plot[idx, 'Dataset'] <- 'Metabolome'
            feat.bar.plot[idx, 'Feature'] <- metab.meta[strlist[2],'CHEMICAL_NAME']
        } else if (strlist[1] == 'RNA') {
            feat.bar.plot[idx, 'Dataset'] <- 'Transcriptome'
            feat.bar.plot[idx, 'Feature'] <- strlist[2]
        } else if (strlist[1] == 'Prot') {
            feat.bar.plot[idx, 'Dataset'] <- 'Proteome'
            feat.bar.plot[idx, 'Feature'] <- strlist[2]
        } else {
            feat.bar.plot[idx, 'Dataset'] <- strlist[1]
            feat.bar.plot[idx, 'Feature'] <- gsub('[.]', '-', gsub('hsa[.]', '', strlist[2]))
        }
    }
    idx <- idx + 1
}

feat.bar.plot$Feature[1] <- 'Gestational Age'

feat.bar.plot$Dataset <- factor(feat.bar.plot$Dataset,
                                levels=c('Clinical', 'Transcriptome', 'miRNA',
                                         'Proteome', 'Metabolome'))

featplot <- ggplot(feat.bar.plot,
       aes(x=Count,
           y=factor(Feature, levels=rev(Feature)),
           fill=Dataset)) +
    geom_bar(stat='identity', color='black') +
    scale_fill_manual(values=c('#F1E2CC', '#CBD5E8', '#FFFFCC', '#FDDAEC', '#E6F5C9')) +
    ylab("Feature") +
    theme_bw()
featplot

ggsave('10rep_10fold_hvar_enet_feats_bacc_grouped.png', featplot, width=7, height=10, dpi=400)
ggsave('10rep_10fold_hvar_enet_feats_bacc_grouped.pdf', featplot, width=7, height=10, dpi=400)


#### Final models

prot.vars.full <- apply(prot[samp.ids,], 2, stats::var)
metab.vars.full <- apply(metab[samp.ids,], 2, stats::var)
rna.vars.full <- apply(rna[samp.ids,], 2, stats::var)
mirna.vars.full <- apply(mirna[samp.ids,], 2, stats::var)

hvar.prots.full <- colnames(prot)[which(prot.vars.full >= quantile(prot.vars.full, 0.75))]
hvar.metabs.full <- colnames(metab)[which(metab.vars.full>=quantile(metab.vars.full, 0.75))]
hvar.rnas.full <- colnames(rna)[which(rna.vars.full >= quantile(rna.vars.full, 0.75))]
hvar.mirnas.full <- colnames(mirna)[which(mirna.vars.full>=quantile(mirna.vars.full, 0.75))]


concat <- cbind(prot[samp.ids,hvar.prots.full],
                metab[samp.ids,hvar.metabs.full],
                mirna[samp.ids,hvar.mirnas.full],
                rna[samp.ids,hvar.rnas.full],
                model.matrix(~., clin[samp.ids,] %>% select(-Condition.))[,-1],
                pathology.slides[samp.ids,] %>% mutate_all(as.numeric) - 1)

class.wts <- nrow(concat) /
    (4 * table(clin.clusters[samp.ids,'hvar.Clusters']))

case.wts <- (class.wts[clin.clusters[samp.ids,'hvar.Clusters']])
case.wts <- rep(1,nrow(concat))

res.final <- glmnet(as.matrix(concat),
                    factor(clin.clusters[samp.ids,'hvar.Clusters']),
                    lambda=lambdas[l.idx]/alphas[a.idx], alpha=alphas[a.idx],
                    family='multinomial', weights=as.vector(case.wts),
                    type.multinomial = 'grouped')

res.final

beta.final <- coef(res.final)
print("Cluster 1")
print(beta.final[[1]]@Dimnames[[1]][as.matrix(beta.final[[1]])!=0][-1])
print(beta.final[[1]][as.matrix(beta.final[[1]])!=0][-1])
print("Cluster 2")
print(beta.final[[2]]@Dimnames[[1]][as.matrix(beta.final[[2]])!=0][-1])
print(beta.final[[2]][as.matrix(beta.final[[2]])!=0][-1])
print("Cluster 3")
print(beta.final[[3]]@Dimnames[[1]][as.matrix(beta.final[[3]])!=0][-1])
print(beta.final[[3]][as.matrix(beta.final[[3]])!=0][-1])
print("Cluster 4")
print(beta.final[[4]]@Dimnames[[1]][as.matrix(beta.final[[4]])!=0][-1])
print(beta.final[[4]][as.matrix(beta.final[[4]])!=0][-1])

incomp.samp.ids <- setdiff(full.samp.ids, samp.ids)
concat.means <- colMeans(concat)

incomp <- matrix(rep(concat.means, length(incomp.samp.ids)), nrow=length(incomp.samp.ids), ncol=ncol(concat), byrow=TRUE)
## incomp <- matrix(NA, nrow=length(incomp.samp.ids), ncol=ncol(concat), byrow=TRUE)

rownames(incomp) <- incomp.samp.ids
colnames(incomp) <- colnames(concat)

incomp[intersect(rownames(prot),incomp.samp.ids),
       intersect(colnames(prot),colnames(concat))] <- as.matrix(prot[intersect(rownames(prot),
                                                                               incomp.samp.ids),
                                                                     intersect(colnames(prot),
                                                                               colnames(concat))])
incomp[intersect(rownames(metab),incomp.samp.ids),
       intersect(colnames(metab),colnames(concat))] <- as.matrix(metab[intersect(rownames(metab),
                                                                                 incomp.samp.ids),
                                                                       intersect(colnames(metab),
                                                                                 colnames(concat))])
incomp[intersect(rownames(rna),incomp.samp.ids),
       intersect(colnames(rna),colnames(concat))] <- as.matrix(rna[intersect(rownames(rna),
                                                                             incomp.samp.ids),
                                                                   intersect(colnames(rna),
                                                                             colnames(concat))])
incomp[intersect(rownames(mirna),incomp.samp.ids),
       intersect(colnames(mirna),colnames(concat))] <- as.matrix(mirna[intersect(rownames(mirna),
                                                                                 incomp.samp.ids),
                                                                       intersect(colnames(mirna),
                                                                                 colnames(concat))])

clin.mod.matrix <- model.matrix(~., clin[full.samp.ids,] %>% select(-Condition.))[,-1]

incomp[intersect(rownames(clin.mod.matrix),incomp.samp.ids),
       intersect(colnames(clin.mod.matrix),colnames(concat))] <- as.matrix(clin.mod.matrix[intersect(rownames(clin.mod.matrix),
                                                                                                     incomp.samp.ids),
                                                                                           intersect(colnames(clin.mod.matrix),
                                                                                                     colnames(concat))])

path.mat <- pathology.slides[full.samp.ids,] %>% mutate_all(as.numeric) - 1

incomp[intersect(rownames(path.mat),incomp.samp.ids),
       intersect(colnames(path.mat),colnames(concat))] <- as.matrix(path.mat[intersect(rownames(path.mat),
                                                                                       incomp.samp.ids),
                                                                             intersect(colnames(path.mat),
                                                                                       colnames(concat))])



pred.final <- predict(res.final, newx=incomp, type='class')
pred.final
table(pred.final)
table(clin[incomp.samp.ids, 'Condition.'], pred.final)

prob.final <- predict(res.final, newx=incomp, type='response')
## prob.final
par(mfrow=c(2,2))
for (c in 1:4) {
    hist(prob.final[,c,])
}


## modality specific models

res.prot <- glmnet(as.matrix(concat[,c('WksGest', names(enet.feat.counts[enet.feat.counts >= floor(0.5 * replicates * K)])[grepl('Prot_', names(enet.feat.counts[enet.feat.counts >= floor(0.5 * replicates * K)]))])]),
                   factor(clin.clusters[samp.ids,'hvar.Clusters']),
                   lambda=0.5, alpha=0,
                   family='multinomial', weights=as.vector(case.wts),
                   type.multinomial = 'grouped')

coef(res.prot)

pred.prot <- predict(res.prot,
                     newx=na.omit(incomp[intersect(incomp.samp.ids, rownames(prot)),c('WksGest', names(enet.feat.counts[enet.feat.counts >= floor(0.5 * replicates * K)])[grepl('Prot_', names(enet.feat.counts[enet.feat.counts >= floor(0.5 * replicates * K)]))])]),
                     type='class')
pred.prot
table(pred.prot)
table(clin[intersect(incomp.samp.ids, rownames(prot)), 'Condition.'], pred.prot)

prob.prot <- predict(res.prot,
                     newx=na.omit(incomp[intersect(incomp.samp.ids, rownames(prot)),c('WksGest', names(enet.feat.counts[enet.feat.counts >= floor(0.5 * replicates * K)])[grepl('Prot_', names(enet.feat.counts[enet.feat.counts >= floor(0.5 * replicates * K)]))])]),
                     type='response')
## prob.prot
par(mfrow=c(2,2))
for (c in 1:4) {
    hist(prob.prot[,c,])
}


res.metab <- glmnet(as.matrix(concat[,c('WksGest', names(enet.feat.counts[enet.feat.counts >= floor(0.5 * replicates * K)])[grepl('Metab_', names(enet.feat.counts[enet.feat.counts >= floor(0.5 * replicates * K)]))])]),
                   factor(clin.clusters[samp.ids,'hvar.Clusters']),
                   lambda=0.5, alpha=0,
                   family='multinomial', weights=as.vector(case.wts),
                   type.multinomial = 'grouped')

coef(res.metab)

pred.metab <- predict(res.metab,
                     newx=na.omit(incomp[intersect(incomp.samp.ids, rownames(metab)),c('WksGest', names(enet.feat.counts[enet.feat.counts >= floor(0.5 * replicates * K)])[grepl('Metab_', names(enet.feat.counts[enet.feat.counts >= floor(0.5 * replicates * K)]))])]),
                     type='class')
pred.metab
table(pred.metab)
table(clin[intersect(incomp.samp.ids, rownames(metab)), 'Condition.'], pred.metab)

prob.metab <- predict(res.metab,
                     newx=na.omit(incomp[intersect(incomp.samp.ids, rownames(metab)),c('WksGest', names(enet.feat.counts[enet.feat.counts >= floor(0.5 * replicates * K)])[grepl('Metab_', names(enet.feat.counts[enet.feat.counts >= floor(0.5 * replicates * K)]))])]),
                     type='response')
## prob.metab
par(mfrow=c(2,2))
for (c in 1:4) {
    hist(prob.metab[,c,])
}


res.rna <- glmnet(as.matrix(concat[,c('WksGest', names(enet.feat.counts[enet.feat.counts >= floor(0.5 * replicates * K)])[grepl('RNA_', names(enet.feat.counts[enet.feat.counts >= floor(0.5 * replicates * K)]))])]),
                   factor(clin.clusters[samp.ids,'hvar.Clusters']),
                   lambda=0.5, alpha=0,
                   family='multinomial', weights=as.vector(case.wts),
                   type.multinomial = 'grouped')

coef(res.rna)

pred.rna <- predict(res.rna,
                     newx=na.omit(incomp[intersect(incomp.samp.ids, rownames(rna)),c('WksGest', names(enet.feat.counts[enet.feat.counts >= floor(0.5 * replicates * K)])[grepl('RNA_', names(enet.feat.counts[enet.feat.counts >= floor(0.5 * replicates * K)]))])]),
                     type='class')
pred.rna
table(pred.rna)
table(clin[intersect(incomp.samp.ids, rownames(rna)), 'Condition.'], pred.rna)

prob.rna <- predict(res.rna,
                     newx=na.omit(incomp[intersect(incomp.samp.ids, rownames(rna)),c('WksGest', names(enet.feat.counts[enet.feat.counts >= floor(0.5 * replicates * K)])[grepl('RNA_', names(enet.feat.counts[enet.feat.counts >= floor(0.5 * replicates * K)]))])]),
                     type='response')
## prob.rna
par(mfrow=c(2,2))
for (c in 1:4) {
    hist(prob.rna[,c,])
}


res.mirna <- glmnet(as.matrix(concat[,c('WksGest', names(enet.feat.counts[enet.feat.counts >= floor(0.5 * replicates * K)])[grepl('miRNA_', names(enet.feat.counts[enet.feat.counts >= floor(0.5 * replicates * K)]))])]),
                   factor(clin.clusters[samp.ids,'hvar.Clusters']),
                   lambda=0.5, alpha=0,
                   family='multinomial', weights=as.vector(case.wts),
                   type.multinomial = 'grouped')

coef(res.mirna)

pred.mirna <- predict(res.mirna,
                     newx=na.omit(incomp[intersect(incomp.samp.ids, rownames(mirna)),c('WksGest', names(enet.feat.counts[enet.feat.counts >= floor(0.5 * replicates * K)])[grepl('miRNA_', names(enet.feat.counts[enet.feat.counts >= floor(0.5 * replicates * K)]))])]),
                     type='class')
pred.mirna
table(pred.mirna)
table(clin[intersect(incomp.samp.ids, rownames(mirna)), 'Condition.'], pred.mirna)

prob.mirna <- predict(res.mirna,
                     newx=na.omit(incomp[intersect(incomp.samp.ids, rownames(mirna)),c('WksGest', names(enet.feat.counts[enet.feat.counts >= floor(0.5 * replicates * K)])[grepl('miRNA_', names(enet.feat.counts[enet.feat.counts >= floor(0.5 * replicates * K)]))])]),
                     type='response')
## prob.mirna
par(mfrow=c(2,2))
for (c in 1:4) {
    hist(prob.mirna[,c,])
}



## modality specific elastic net models
library(glmnetUtils)

res.prot <- cv.glmnet(as.matrix(cbind(prot[samp.ids,hvar.prots.full],
                                      WksGest=clin[samp.ids,] %>% select(WksGest))),
                      factor(clin.clusters[samp.ids,'hvar.Clusters']),
                      alpha=1,
                      family='multinomial', weights=as.vector(case.wts),
                      type.multinomial = 'grouped')

plot(res.prot)

coef(res.prot, s=res.prot$lambda.1se)

pred.prot <- predict(res.prot,
                     newx=na.omit(incomp[intersect(incomp.samp.ids, rownames(prot)),c(hvar.prots.full, 'WksGest')]),
                     type='class', s=res.prot$lambda.1se)
## pred.prot
table(pred.prot)
table(clin[intersect(incomp.samp.ids, rownames(prot)), 'Condition.'], pred.prot)

prob.prot <- predict(res.prot,
                     newx=na.omit(incomp[intersect(incomp.samp.ids, rownames(prot)), c(hvar.prots.full, 'WksGest')]),
                     type='response', s=res.prot$lambda.1se)
## prob.prot
par(mfrow=c(2,2))
for (c in 1:4) {
    hist(prob.prot[,c,])
}


res.metab <- cv.glmnet(as.matrix(cbind(metab[samp.ids,hvar.metabs.full],
                                      WksGest=clin[samp.ids,] %>% select(WksGest))),
                      factor(clin.clusters[samp.ids,'hvar.Clusters']),
                      alpha=1,
                      family='multinomial', weights=as.vector(case.wts),
                      type.multinomial = 'grouped')

plot(res.metab)

coef(res.metab, s=res.metab$lambda.1se)

pred.metab <- predict(res.metab,
                     newx=na.omit(incomp[intersect(incomp.samp.ids, rownames(metab)),c(hvar.metabs.full, 'WksGest')]),
                     type='class', s=res.metab$lambda.1se)
## pred.metab
table(pred.metab)
table(clin[intersect(incomp.samp.ids, rownames(metab)), 'Condition.'], pred.metab)

prob.metab <- predict(res.metab,
                     newx=na.omit(incomp[intersect(incomp.samp.ids, rownames(metab)), c(hvar.metabs.full, 'WksGest')]),
                     type='response', s=res.metab$lambda.1se)
## prob.metab
par(mfrow=c(2,2))
for (c in 1:4) {
    hist(prob.metab[,c,])
}


res.rna <- cv.glmnet(as.matrix(cbind(rna[samp.ids,hvar.rnas.full],
                                      WksGest=clin[samp.ids,] %>% select(WksGest))),
                      factor(clin.clusters[samp.ids,'hvar.Clusters']),
                      alpha=1,
                      family='multinomial', weights=as.vector(case.wts),
                      type.multinomial = 'grouped')

plot(res.rna)

coef(res.rna, s=res.rna$lambda.1se)

pred.rna <- predict(res.rna,
                     newx=na.omit(incomp[intersect(incomp.samp.ids, rownames(rna)),c(hvar.rnas.full, 'WksGest')]),
                     type='class', s=res.rna$lambda.1se)
## pred.rna
table(pred.rna)
table(clin[intersect(incomp.samp.ids, rownames(rna)), 'Condition.'], pred.rna)

prob.rna <- predict(res.rna,
                     newx=na.omit(incomp[intersect(incomp.samp.ids, rownames(rna)), c(hvar.rnas.full, 'WksGest')]),
                     type='response', s=res.rna$lambda.1se)
## prob.rna
par(mfrow=c(2,2))
for (c in 1:4) {
    hist(prob.rna[,c,])
}


res.mirna <- cv.glmnet(as.matrix(cbind(mirna[samp.ids,hvar.mirnas.full],
                                      WksGest=clin[samp.ids,] %>% select(WksGest))),
                      factor(clin.clusters[samp.ids,'hvar.Clusters']),
                      alpha=1,
                      family='multinomial', weights=as.vector(case.wts),
                      type.multinomial = 'grouped')

plot(res.mirna)

coef(res.mirna, s=res.mirna$lambda.1se)

pred.mirna <- predict(res.mirna,
                     newx=na.omit(incomp[intersect(incomp.samp.ids, rownames(mirna)),c(hvar.mirnas.full, 'WksGest')]),
                     type='class', s=res.mirna$lambda.1se)
## pred.mirna
table(pred.mirna)
table(clin[intersect(incomp.samp.ids, rownames(mirna)), 'Condition.'], pred.mirna)

prob.mirna <- predict(res.mirna,
                     newx=na.omit(incomp[intersect(incomp.samp.ids, rownames(mirna)), c(hvar.mirnas.full, 'WksGest')]),
                     type='response', s=res.mirna$lambda.1se)
## prob.mirna
par(mfrow=c(2,2))
for (c in 1:4) {
    hist(prob.mirna[,c,])
}


incomp.pred.mat <- matrix(NA, nrow=length(incomp.samp.ids), ncol=5)
rownames(incomp.pred.mat) <- incomp.samp.ids
colnames(incomp.pred.mat) <- c('prot', 'metab', 'rna', 'mirna', 'final')

incomp.pred.mat[intersect(incomp.samp.ids, rownames(prot)), 'prot'] <- pred.prot
incomp.pred.mat[intersect(incomp.samp.ids, rownames(metab)), 'metab'] <- pred.metab
incomp.pred.mat[intersect(incomp.samp.ids, rownames(rna)), 'rna'] <- pred.rna
incomp.pred.mat[intersect(incomp.samp.ids, rownames(mirna)), 'mirna'] <- pred.mirna
incomp.pred.mat[incomp.samp.ids, 'final'] <- pred.final

incomp.pred.mat

ami.modality.mat <- matrix(0, 5, 5)
rownames(ami.modality.mat) <- c('prot', 'metab', 'rna', 'mirna', 'final')
colnames(ami.modality.mat) <- c('prot', 'metab', 'rna', 'mirna', 'final')

for (i in 1:4) {
    for (j in (i+1):5) {
        samps <- rownames(na.omit(incomp.pred.mat[,c(i,j)]))
        print(incomp.pred.mat[samps,])
        ami.modality.mat[i,j] <- aricode::AMI(incomp.pred.mat[samps,i],
                                              incomp.pred.mat[samps,j])
    }
}

ami.modality.mat



concat <- cbind(prot[samp.ids,],
                metab[samp.ids,],
                mirna[samp.ids,],
                rna[samp.ids,],
                clin[samp.ids,],
                pathology.slides[samp.ids,])

library(rCausalMGM)

clin.clusters <- read.csv('clustered_clinical_hvar.csv', row.names=1)

clin.clusters <- clin.clusters %>% mutate_at('hvar.Clusters', factor)

concat <- cbind(prot[samp.ids,],
                metab[samp.ids,],
                mirna[samp.ids,],
                rna[samp.ids,],
                hvar.Clusters=clin.clusters[samp.ids,'hvar.Clusters'],
                model.matrix(~0+hvar.Clusters, clin.clusters[samp.ids,]), 
                clin[samp.ids,],
                pathology.slides[samp.ids,])

concat <- concat %>% mutate_at(c('hvar.Clusters', colnames(pathology.slides),
                                 'Race', 'InfSex', 'Smoking',
                                 'Labor.initiation', 'FDELTYPE',
                                 paste0('hvar.Clusters', 1:4)), factor)

## , colnames(pathology.slides),
## 'Race', 'InfSex', 'PrePregBMI', 'Smoking',
## 'Labor.initiation', 'FDELTYPE'

ig.path <- mgmPath(concat[,c(names(enet.feat.counts[enet.feat.counts >=
                                                    floor(0.5*replicates*K)]),
                             'hvar.Clusters')],
                   lambdas=1*0.005^((0:199)/199), verbose=T, rank=F)

## nparams <- (ig.path$BIC + 2*ig.path$loglik) / log(270)
## ig.path$AIC <- -2 * ig.path$loglik + 2 * nparams

score.range <- max(ig.path$BIC) - min(-2*ig.path$loglik)
plot(log10(ig.path$lambda), ig.path$BIC, ylim = c(min(-2*ig.path$loglik) - 0.1*score.range,
                                                  max(ig.path$BIC) + 0.1*score.range))
points(log10(ig.path$lambda), ig.path$AIC, col='green')
points(log10(ig.path$lambda), -2*ig.path$loglik, col='blue')

ig.path$graphs[[which.min(ig.path$BIC)]]$markov.blankets[['hvar.Clusters']]
ig.path$graphs[[which.min(ig.path$AIC)]]$markov.blankets[['hvar.Clusters']]

ig.aic <- ig.path$graphs[[which.min(ig.path$AIC)]]
ig.aic

## ig.aicc <- ig.path$graphs[[which.min(ig.path$AICc)]]
## ig.aicc

ig.bic <- ig.path$graphs[[which.min(ig.path$BIC)]]
ig.bic

## ig.steps <- steps(concat[,c(names(enet.feat.counts[enet.feat.counts >=
##                                                    floor(0.5*replicates*K)]),
##                             'Clusters')], nLambda=50,
##                   numSub=50, g=0.1, verbose=T)

## ig <- mgm(concat[,c(names(enet.feat.counts[enet.feat.counts >=
##                                            floor(0.5*replicates*K)]),
##                     'Clusters')],
##           lambda=ig.steps$lambda[1])

g <- fciMax(concat[,c(names(enet.feat.counts[enet.feat.counts >=
                                             floor(0.5*replicates*K)]),
                      'hvar.Clusters')],
            initialGraph=ig.bic, 
            alpha=0.05, fdr=F, verbose=T, rank=F)

print(g)
print(g$edges)
print(g$markov.blankets$hvar.Clusters)

plot.graph(g, nodeAttr=list(shape='none', cex=0.9), edgeAttr=list(lwd=5))
plot.graph(g, unique(c('hvar.Clusters',
                       g$markov.blankets[['hvar.Clusters1']])),
           nodeAttr=list(shape='none', cex=0.9),
           edgeAttr=list(lwd=3))


for (clust in seq_len(4)) {

    ig.path <- mgmPath(concat[,c(names(enet.feat.counts[enet.feat.counts >=
                                                        floor(0.5*replicates*K)]),
                                 paste0('hvar.Clusters', clust))],
                       lambdas=1*0.005^((0:199)/199), verbose=T, rank=F)

    ## nparams <- (ig.path$BIC + 2*ig.path$loglik) / log(270)
    ## ig.path$AIC <- -2 * ig.path$loglik + 2 * nparams

    score.range <- max(ig.path$BIC) - min(-2*ig.path$loglik)
    print(plot(log10(ig.path$lambda), ig.path$BIC,
               ylim = c(min(-2*ig.path$loglik) - 0.1*score.range,
                        max(ig.path$BIC) + 0.1*score.range)))
    points(log10(ig.path$lambda), ig.path$AIC, col='green')
    points(log10(ig.path$lambda), -2*ig.path$loglik, col='blue')

    ig.path$graphs[[which.min(ig.path$BIC)]]$markov.blankets[[paste0('hvar.Clusters', clust)]]
    ig.path$graphs[[which.min(ig.path$AIC)]]$markov.blankets[[paste0('hvar.Clusters', clust)]]

    ig.aic <- ig.path$graphs[[which.min(ig.path$AIC)]]
    ig.aic

    ## ig.aicc <- ig.path$graphs[[which.min(ig.path$AICc)]]
    ## ig.aicc

    ig.bic <- ig.path$graphs[[which.min(ig.path$BIC)]]
    ig.bic

    ## ig.steps <- steps(concat[,c(names(enet.feat.counts[enet.feat.counts >=
    ##                                                    floor(0.5*replicates*K)]),
    ##                             'Clusters')], nLambda=50,
    ##                   numSub=50, g=0.1, verbose=T)

    ## ig <- mgm(concat[,c(names(enet.feat.counts[enet.feat.counts >=
    ##                                            floor(0.5*replicates*K)]),
    ##                     'Clusters')],
    ##           lambda=ig.steps$lambda[1])

    g <- fciMax(concat[,c(names(enet.feat.counts[enet.feat.counts >=
                                                 floor(0.5*replicates*K)]),
                          paste0('hvar.Clusters', clust))],
                initialGraph=ig.bic, 
                alpha=0.1, fdr=T, verbose=T, rank=F)


    saveGraph(g, paste0('snf_hvar_clusts/graph/hvar_snfClust',
                        clust, '_enetFeats_mgmBICfcimax_fdr1.txt'))
    
    saveGraph(g, paste0('snf_hvar_clusts/graph/hvar_snfClust',
                        clust, '_enetFeats_mgmBICfcimax_fdr1.sif'))


    g.boot <- bootstrap(concat[,c(names(enet.feat.counts[enet.feat.counts >=
                                                         floor(0.5*replicates*K)]),
                                  paste0('hvar.Clusters', clust))],
                        ensemble='highest',
                        algorithm='mgmfcimax',
                        lambda=ig.bic$lambda, alpha=0.1,
                        numBoots=100, verbose=T)

    g.table <- graphTable(g, g.boot$stabilities)

    write.csv(g.table,
              paste0('snf_hvar_clusts/graph/hvar_snfClust',
                     clust, '_enetFeats_mgmBICfcimax_fdr1.csv'),
              row.names=FALSE, quote=FALSE)

    
}



library(nnet)
f <- as.formula(paste0('hvar.Clusters ~ 1 + ',
                       paste(g$markov.blankets$hvar.Clusters, collapse=' + ')))

res <- multinom(f, concat)
print(res)
print(paste0('BIC: ', res$deviance + log(270) * 3 * (length(g$markov.blankets$hvar.Clusters) + 1)))

f2 <- hvar.Clusters ~ 1 + WksGest + RNA_MYO1A + RNA_H4C1 + RNA_EHD1 + RNA_BIN2 + Prot_FSTL3 + Prot_ACE2

res2 <- multinom(f2, concat)
print(res2)
print(paste0('BIC: ', res2$deviance + log(270) * 3 * (length(g$markov.blankets$hvar.Clusters) + 1)))

f3 <- as.formula(paste0('hvar.Clusters ~ 1 + ', paste(c(possPar, setdiff(g$markov.blankets$hvar.Clusters, c(possPar, adjSet)), adjSet), collapse=' + ')))
f3
res3 <- multinom(f3, concat)
print(res3)
print(paste0('BIC: ', res3$deviance + log(270) * 3 * (length(g$markov.blankets$hvar.Clusters) + 1)))
res3.sum <- summary(res3)


f4 <- as.formula(paste0('hvar.Clusters ~ 1 + ', paste(c(possPar, adjSet), collapse=' + ')))
f4
res4 <- multinom(f4, concat)
print(res4)
print(paste0('BIC: ', res4$deviance + log(270) * 3 * (length(c(possPar, adjSet)) + 1)))
res4.sum <- summary(res4)


getPossParents <- function(g, target) {
    possParents <- c()
    for (e in g$edges) {
        e <- strsplit(e, ' ')[[1]]
        if (e[3] == target & e[2] != '<->') {
            possParents <- c(possParents, e[1])
        } else if (e[1] == target & e[2] == 'o-o') {
            possParents <- c(possParents, e[3])
        }
    }
    return(possParents)
}

getConfounded <- function(g, target) {
    confounded <- c()
    for (e in g$edges) {
        e <- strsplit(e, ' ')[[1]]
        if (any(target == e) & e[2] == '<->') {
            confounded <- c(confounded, e[e != target & e != '<->'])
        } 
    }
    return(confounded)
}

getAdjustmentSet <- function(g, target) {
    visited <- as.vector(matrix(FALSE, nrow=length(g$nodes)))
    names(visited) <- g$nodes
    visited[target] <- TRUE
    adjSet <- c()
    possParents <- getPossParents(g, target)
    
    confounded <- getConfounded(g, target)
    while (length(possParents) > 0) {
        print(paste0('Possible parents of: ', possParents[1], ':'))
        print(getPossParents(g, possParents[1]))
        adjSet <- c(adjSet, getPossParents(g, possParents[1]))
        newConfounded <- getConfounded(g, possParents[1])
        newConfounded <- newConfounded[!visited[newConfounded]]
        visited[newConfounded] <- TRUE
        print(paste0('Confounders of: ', possParents[1], ':'))
        print(newConfounded)
        confounded <- c(confounded, newConfounded)
        possParents <- possParents[-1]
        while (length(confounded) > 0) {
            ## print(confounded)
            print(paste0('Possible parents of: ', confounded[1], ':'))
            print(getPossParents(g, confounded[1]))
            adjSet <- c(adjSet, getPossParents(g, confounded[1]))
            newConfounded <- getConfounded(g, confounded[1])
            newConfounded <- newConfounded[!visited[newConfounded]]
            visited[newConfounded] <- TRUE
            print(paste0('Confounders of: ', confounded[1], ':'))
            print(newConfounded)
            confounded <- c(confounded, newConfounded)
            confounded <- confounded[-1]
        }
    }
    return(setdiff(unique(adjSet), getPossParents(g, target)))
}

possPar <- getPossParents(g, 'hvar.Clusters')
possPar
getConfounded(g, 'hvar.Clusters')
adjSet <- getAdjustmentSet(g, 'hvar.Clusters')
adjSet

pheatmap::pheatmap(cor(concat[,names(enet.feat.counts[enet.feat.counts >= floor(0.5*replicates*K)])]) - diag(nrow=47), clustering_method='ward.D')

pheatmap::pheatmap(cor(concat[,g$markov.blankets$hvar.Clusters]) - diag(nrow=length(g$markov.blankets$hvar.Clusters)), clustering_method='ward.D')

saveGraph(g, 'snf_hvar_clusts/graph/hvar_snfClust_enetFeats_mgmBICfcimax_fdr05.txt')
saveGraph(g, 'snf_hvar_clusts/graph/hvar_snfClust_enetFeats_mgmBICfcimax_fdr05.sif')


g.boot <- bootstrap(concat[,c(names(enet.feat.counts[enet.feat.counts >=
                                                     floor(0.5*replicates*K)]),
                              'hvar.Clusters')],
                    ensemble='highest',
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

            if (edge[3] == g.boot$stabilities[i,1] & edge[1] == g.boot$stabilities[i,3]) {
                g.stabs[count,] <- g.stabs[count,c(3,2,1,4,6,5,7,9,8,10,11)]
            }
            
            break
        }
    }

    count <- count + 1
}

g.stabs[,4:11] <- pmax(sapply(g.stabs[,4:11], as.numeric), 0)

g.stabs

write.csv(g.stabs,
          'snf_hvar_clusts/graph/hvar_snfClust_enetFeats_mgmBICfcimax_fdr05.csv',
          row.names=FALSE, quote=FALSE)

saveRDS(score.auc, '10rep_10fold_auc_v2.rds')
saveRDS(score.bacc, '10rep_10fold_bacc_v2.rds')
saveRDS(enet.feat.counts, '10rep_10fold_feat_counts_v2.rds')


## miRNA log concentrations

for (clust in 1:4) {
    mirna.markers <- read.csv(paste0('snf_hvar_clusts/markers/clust',
                                     clust, '/clust', clust,
                                     '_auc_mirna_signif_markers.csv'),
                              row.names=1)

    mirna.markers <- mirna.markers[1:10,]

    mirna.means <- colMeans(mirna[samp.ids[clin.clusters[samp.ids,'hvar.Clusters']==clust],
                                  mirna.markers$Feature])

    
    mirna.means.table <- data.frame(miRNA=gsub('[.]', '-', gsub('miRNA_', '', mirna.markers$Feature)),
                                    ClusterMean=mirna.means)

    write.csv(mirna.means.table,
              paste0('snf_hvar_clusts/markers/clust',
                     clust, '/clust', clust, '_mirna_comir_table.csv'), row.names=FALSE)
}


## Marker Features heatmap

library(pheatmap)
library(ggplot2)

prot.markers <- c('FLT1', 'PGF', 'LEP', 'FSTL3', 'IGFBP-1', 'VEGFA', 'GDF-15')

rna.markers <- c('FLT1', 'PGF', 'ENG', 'LEP', 'FSTL3', 'IGFBP1', 'VEGFA', 'GDF15',
                 'TPBG', 'LGALS13', 'BTG2', 'CLDN1', 'LGALS14', 'ZNF554', 'SIGLEC6',
                 'HTRA1', 'ACVRL1', 'EGFL7')

mirna.markers <- c('hsa-miR-210-3p', 'hsa-miR-193b-5p', 'hsa-miR-24-3p',
                   'hsa-miR-1301-3p', 'hsa-miR-223-3p', 'hsa-miR-224-5p')

metab.markers <- c('ethanolamine', 'glycerophosphorylcholine (GPC)', 'aspartate',
                   'taurine', 'lysine', 'glucose', 'myo-inositol', 'glutamate',
                   'glycine', 'creatine')

prot.names <- paste0('Prot_', make.names(prot.markers))

metab.names <- paste0('Metab_', make.names(metab.markers))

rna.names <- paste0('RNA_', rna.markers)

mirna.names <- paste0('miRNA_', make.names(mirna.markers))

mirna.markers <- c('miR-210-3p', 'miR-193b-5p', 'miR-24-3p',
                   'miR-1301-3p', 'miR-223-3p', 'miR-224-5p')

marker.table <- cbind(rna[samp.ids,rna.names],
                      mirna[samp.ids,mirna.names],
                      prot[samp.ids,prot.names],
                      metab[samp.ids,metab.names])

scaled.table <- scale(marker.table)

samp.order <- order(clin.clusters[samp.ids,'hvar.Clusters'])

annot.df <- data.frame(Clusters=factor(clin.clusters[samp.ids,'hvar.Clusters'],
                                       levels=c(1,2,3,4)),
                       row.names=samp.ids)

col.breaks <- cumsum(table(clin.clusters[samp.ids,'hvar.Clusters']))[-4]
row.breaks <- cumsum(c(length(rna.names), length(mirna.names), length(prot.names)))

samp.clust.map <- pheatmap(t(scaled.table[samp.order,]), breaks = seq(-4, 4, length.out=101),
         cluster_rows = F, cluster_cols=F, gaps_col=col.breaks, gaps_row=row.breaks,
         annotation_col=annot.df, show_colnames=F,
         labels_row=c(rna.markers, mirna.markers, prot.markers, metab.markers))

ggsave('snf_hvar_clusts/known_markers_hvar_clusts_heatmap.pdf', samp.clust.map$gtable, width=6, height=8)
ggsave('snf_hvar_clusts/known_markers_hvar_clusts_heatmap.png', samp.clust.map$gtable, width=6, height=8, dpi=400)


samp.order <- order(factor(clin.clusters[samp.ids,'Condition.'],
                           levels=c('Control', 'Control PTD', 'Severe PE',
                                    'FGR', 'FGR+HDP', 'PTD')))

annot.df <- data.frame(Condition=factor(clin.clusters[samp.ids,'Condition.'],
                                        levels=c('Control', 'Control PTD', 'Severe PE',
                                                 'FGR', 'FGR+HDP', 'PTD')),
                       row.names=samp.ids)

col.breaks <- cumsum(table(clin.clusters[samp.ids,'Condition.']))[-6]
row.breaks <- cumsum(c(length(prot.names), length(rna.names), length(mirna.names)))

samp.cond.map <- pheatmap(t(scaled.table[samp.order,]), breaks = seq(-4, 4, length.out=101),
         cluster_rows = F, cluster_cols=F, gaps_col=col.breaks, gaps_row=row.breaks,
         annotation_col=annot.df, show_colnames=F,
         labels_row=c(prot.markers, rna.markers, mirna.markers, metab.markers))

ggsave('snf_hvar_clusts/known_markers_condition_heatmap.pdf', samp.cond.map$gtable, width=7, height=8)
ggsave('snf_hvar_clusts/known_markers_condition_heatmap.png', samp.cond.map$gtable, width=7, height=8, dpi=400)



summary.table <- scaled.table %>%
    as.data.frame %>%
    group_by(clin.clusters[samp.ids,'hvar.Clusters']) %>%
    summarize_all(mean) %>%
    as.data.frame

summary.table <- summary.table[,-1]

rownames(summary.table) <- c(1,2,3,4)

annot.df <- data.frame(Clusters=factor(c(1,2,3,4),
                                       levels=c(1,2,3,4)))

## row.breaks <- cumsum(c(length(prot.names), length(rna.names), length(mirna.names)))

ordered.names <- c(rna.names[order(summary.table[3,rna.names], decreasing=T)],
                   mirna.names[order(summary.table[3,mirna.names], decreasing=T)],
                   prot.names[order(summary.table[3,prot.names], decreasing=T)],
                   metab.names[order(summary.table[3,metab.names], decreasing=T)])

ordered.labels <- c(rna.markers[order(summary.table[3,rna.names], decreasing=T)],
                    mirna.markers[order(summary.table[3,mirna.names], decreasing=T)],
                    prot.markers[order(summary.table[3,prot.names], decreasing=T)],
                    metab.markers[order(summary.table[3,metab.names], decreasing=T)])

library(RColorBrewer)
brewer.colors <- brewer.pal(n=12, "Paired")## )(10)

brewer.colors <- c(brewer.colors[1:2], brewer.colors[5:6], brewer.colors[11], brewer.colors[8], brewer.colors[3:4], brewer.colors[9:10])

annot.colors <- list()
annot.colors[['Clusters']] <- brewer.colors[7:10]
names(annot.colors[['Clusters']]) <- levels(annot.df$Clusters)


max.val <- max(abs(summary.table))

mean.clust.map <- pheatmap(t(summary.table[,ordered.names]),
                           breaks = seq(-max.val, max.val, length.out=101),
         cluster_rows = F, cluster_cols=F, gaps_row=row.breaks, # scale='row',
         annotation_col=annot.df, show_colnames=F,
         labels_row=ordered.labels, annotation_colors=annot.colors)

ggsave('snf_hvar_clusts/known_markers_hvar_clusts_mean_heatmap.pdf', mean.clust.map$gtable, width=6, height=8)
ggsave('snf_hvar_clusts/known_markers_hvar_clusts_mean_heatmap.png', mean.clust.map$gtable, width=6, height=8, dpi=400)


ordered.names <- colnames(summary.table)[order(summary.table[3,], decreasing=T)]
ordered.labels <- c(rna.markers,
                    mirna.markers,
                    prot.markers,
                    metab.markers)[order(summary.table[3,], decreasing=T)]

col.annot.df <- data.frame(
    Modality=factor(c(rep("Transcriptome", length(rna.markers)),
                      rep("miRNA", length(mirna.markers)),
                      rep("Proteome", length(prot.markers)),
                      rep("Metabolome", length(metab.markers)))[order(summary.table[3,],
                                                                      decreasing=T)],
                    levels=c("Transcriptome", "miRNA", "Proteome", "Metabolome")))

rownames(col.annot.df) <- ordered.names

library(RColorBrewer)
brewer.colors <- brewer.pal(n=12, "Paired")## )(10)

brewer.colors <- c(brewer.colors[1:2], brewer.colors[5:6], brewer.colors[11], brewer.colors[8], brewer.colors[3:4], brewer.colors[9:10])

annot.colors <- list()
annot.colors[['Clusters']] <- brewer.colors[7:10]
names(annot.colors[['Clusters']]) <- levels(annot.df$Clusters)
annot.colors[["Modality"]] <- c('#CBD5E8', '#FFFFCC', '#FDDAEC', '#E6F5C9')
names(annot.colors[['Modality']]) <- levels(col.annot.df$Modality)


max.val <- max(abs(summary.table))

mean.clust.map <- pheatmap(t(summary.table[,ordered.names]),
                           breaks = seq(-max.val, max.val, length.out=101),
         cluster_rows = F, cluster_cols=F, # gaps_row=row.breaks, # scale='row',
         annotation_col=annot.df, show_colnames=F, annotation_row=col.annot.df,
         labels_row=ordered.labels, annotation_colors=annot.colors)

ggsave('snf_hvar_clusts/known_markers_hvar_clusts_mean_heatmap_ungrouped.pdf', mean.clust.map$gtable, width=6, height=8)
ggsave('snf_hvar_clusts/known_markers_hvar_clusts_mean_heatmap_ungrouped.png', mean.clust.map$gtable, width=6, height=8, dpi=400)



markers.summary <- scaled.table %>% as.data.frame %>% 
    group_by(clin.clusters[samp.ids,'hvar.Clusters']) %>%
    summarize_all(c(mean=mean, sd=sd)) %>%
    as.data.frame

markers.summary <- markers.summary[,-1]

rownames(markers.summary) <- c(1,2,3,4)

markers.summary

markers.means <- markers.summary[,paste0(c(rna.names, mirna.names,
                                           prot.names, metab.names),
                                         '_mean')]
colnames(markers.means) <- c(rna.names, mirna.names, prot.names, metab.names)

markers.sds <- markers.summary[,paste0(c(rna.names, mirna.names,
                                           prot.names, metab.names),
                                         '_sd')]
colnames(markers.sds) <- c(rna.names, mirna.names, prot.names, metab.names)


marker.posthoc <- data.frame()

marker.summary.table <- data.frame(matrix(NA, 5, ncol(markers.means)))

colnames(marker.summary.table) <- ordered.names
rownames(marker.summary.table) <- c(paste('Cluster ', 1:4, ' (n = ',
                                     as.vector(table(clin.clusters[samp.ids,
                                                                   'hvar.Clusters'])),
                                     ')', sep=''), 'pval')

library(FSA)

for (var in ordered.names) {
    print(var)
    res <- kruskal.test(scaled.table[,var], clin.clusters[samp.ids,'hvar.Clusters'])
    print(res)

    marker.summary.table[1:4,var] <- paste(signif(markers.means[,var], 3),
                                             ' (', signif(markers.sds[,var], 3),
                                             ')', sep='')
    marker.summary.table['pval',var] <- res$p.value

    if (res$p.value < 0.05) {
        res.posthoc <- dunnTest(scaled.table[,var],
                                factor(clin.clusters[samp.ids,'hvar.Clusters']),
                                method='holm')
        print(res.posthoc)

        marker.posthoc <- rbind(marker.posthoc,
                                cbind(Feature=rep(var, nrow(res.posthoc$res)),
                                      res.posthoc$res))
    }

}

marker.summary.table <- as.data.frame(t(marker.summary.table))

marker.summary.table[,'adj.pval'] <- p.adjust(marker.summary.table$pval, method='fdr')

marker.summary.table

marker.posthoc

write.csv(marker.summary.table,
          'snf_hvar_clusts/known_markers_hvar_cluster_table.csv', quote=F)

write.csv(marker.posthoc,
          'snf_hvar_clusts/known_markers_hvar_cluster_dunn_posthoc.csv',
          quote=F, row.names=F)


#### Any PE version for hvar clusters


summary.table <- cbind(clin.clusters[samp.ids, c('Condition.', 'hvar.Clusters')],
                       scaled.table) %>%
    as.data.frame %>%
    filter(Condition. %in% c('Severe PE', 'FGR+HDP')) %>%
    group_by(hvar.Clusters) %>%
    summarize_all(mean) %>%
    as.data.frame

summary.table <- summary.table[,-(1:2)]

rownames(summary.table) <- c(1,2,3,4)

annot.df <- data.frame(Clusters=factor(c(1,2,3,4),
                                       levels=c(1,2,3,4)))

## row.breaks <- cumsum(c(length(prot.names), length(rna.names), length(mirna.names)))

## ordered.names <- c(rna.names[order(summary.table[3,rna.names], decreasing=T)],
##                    mirna.names[order(summary.table[3,mirna.names], decreasing=T)],
##                    prot.names[order(summary.table[3,prot.names], decreasing=T)],
##                    metab.names[order(summary.table[3,metab.names], decreasing=T)])

## ordered.labels <- c(rna.markers[order(summary.table[3,rna.names], decreasing=T)],
##                     mirna.markers[order(summary.table[3,mirna.names], decreasing=T)],
##                     prot.markers[order(summary.table[3,prot.names], decreasing=T)],
##                     metab.markers[order(summary.table[3,metab.names], decreasing=T)])


max.val <- max(abs(summary.table))

mean.clust.map <- pheatmap(t(summary.table[,ordered.names]),
                           breaks = seq(-max.val, max.val, length.out=101),
         cluster_rows = F, cluster_cols=F, gaps_row=row.breaks, # scale='row',
         annotation_col=annot.df, show_colnames=F,
         labels_row=ordered.labels, annotation_colors=annot.colors)

ggsave('snf_hvar_clusts/AnyPE_known_markers_hvar_clusts_mean_heatmap.pdf', mean.clust.map$gtable, width=6, height=8)
ggsave('snf_hvar_clusts/AnyPE_known_markers_hvar_clusts_mean_heatmap.png', mean.clust.map$gtable, width=6, height=8, dpi=400)


## ordered.names <- colnames(summary.table)[order(summary.table[3,], decreasing=T)]
## ordered.labels <- c(rna.markers,
##                     mirna.markers,
##                     prot.markers,
##                     metab.markers)[order(summary.table[3,], decreasing=T)]

## col.annot.df <- data.frame(
##     Modality=factor(c(rep("Transcriptome", length(rna.markers)),
##                       rep("miRNA", length(mirna.markers)),
##                       rep("Proteome", length(prot.markers)),
##                       rep("Metabolome", length(metab.markers)))[order(summary.table[3,],
##                                                                       decreasing=T)],
##                     levels=c("Transcriptome", "miRNA", "Proteome", "Metabolome")))

## rownames(col.annot.df) <- ordered.names

library(RColorBrewer)
brewer.colors <- brewer.pal(n=12, "Paired")## )(10)

brewer.colors <- c(brewer.colors[1:2], brewer.colors[5:6], brewer.colors[11], brewer.colors[8], brewer.colors[3:4], brewer.colors[9:10])

annot.colors <- list()
annot.colors[['Clusters']] <- brewer.colors[7:10]
names(annot.colors[['Clusters']]) <- levels(annot.df$Clusters)
annot.colors[["Modality"]] <- c('#CBD5E8', '#FFFFCC', '#FDDAEC', '#E6F5C9')
names(annot.colors[['Modality']]) <- levels(col.annot.df$Modality)


max.val <- max(abs(summary.table))

mean.clust.map <- pheatmap(t(summary.table[,ordered.names]),
                           breaks = seq(-max.val, max.val, length.out=101),
         cluster_rows = F, cluster_cols=F, # gaps_row=row.breaks, # scale='row',
         annotation_col=annot.df, show_colnames=F, annotation_row=col.annot.df,
         labels_row=ordered.labels, annotation_colors=annot.colors)

ggsave('snf_hvar_clusts/AnyPE_known_markers_hvar_clusts_mean_heatmap_ungrouped.pdf', mean.clust.map$gtable, width=6, height=8)
ggsave('snf_hvar_clusts/AnyPE_known_markers_hvar_clusts_mean_heatmap_ungrouped.png', mean.clust.map$gtable, width=6, height=8, dpi=400)



markers.summary <- cbind(clin.clusters[samp.ids, c('Condition.', 'hvar.Clusters')],
                         scaled.table) %>%
    as.data.frame %>%
    filter(Condition. %in% c('Severe PE', 'FGR+HDP')) %>%
    group_by(hvar.Clusters) %>%
    select(-c('Condition.')) %>%
    summarize_all(c(mean=mean, sd=sd)) %>%
    as.data.frame

markers.summary <- markers.summary[,-1]

rownames(markers.summary) <- c(1,2,3,4)

markers.summary

markers.means <- markers.summary[,paste0(c(rna.names, mirna.names,
                                           prot.names, metab.names),
                                         '_mean')]
colnames(markers.means) <- c(rna.names, mirna.names, prot.names, metab.names)

markers.sds <- markers.summary[,paste0(c(rna.names, mirna.names,
                                           prot.names, metab.names),
                                         '_sd')]
colnames(markers.sds) <- c(rna.names, mirna.names, prot.names, metab.names)


marker.posthoc <- data.frame()

marker.summary.table <- data.frame(matrix(NA, 5, ncol(markers.means)))

colnames(marker.summary.table) <- ordered.names
rownames(marker.summary.table) <- c(paste('Cluster ', 1:4, ' (n = ',
                                          as.vector(table(clin.clusters[samp.ids[
                                              clin.clusters[samp.ids, 'Condition.'] %in%
                                              c('Severe PE', 'FGR+HDP')],
                                              'hvar.Clusters'])),
                                     ')', sep=''), 'pval')

library(FSA)

filt.scaled.table <- cbind(clin.clusters[samp.ids, c('Condition.', 'hvar.Clusters')],
                           scaled.table) %>%
    as.data.frame %>%
    filter(Condition. %in% c('Severe PE', 'FGR+HDP'))

for (var in ordered.names) {
    print(var)
    res <- kruskal.test(filt.scaled.table[,var], filt.scaled.table[,'hvar.Clusters'])
    print(res)

    marker.summary.table[1:4,var] <- paste(signif(markers.means[,var], 3),
                                             ' (', signif(markers.sds[,var], 3),
                                             ')', sep='')
    marker.summary.table['pval',var] <- res$p.value

    if (res$p.value < 0.05) {
        res.posthoc <- dunnTest(filt.scaled.table[,var],
                                factor(filt.scaled.table[,'hvar.Clusters']),
                                method='holm')
        print(res.posthoc)

        marker.posthoc <- rbind(marker.posthoc,
                                cbind(Feature=rep(var, nrow(res.posthoc$res)),
                                      res.posthoc$res))
    }

}

marker.summary.table <- as.data.frame(t(marker.summary.table))

marker.summary.table[,'adj.pval'] <- p.adjust(marker.summary.table$pval, method='fdr')

marker.summary.table

marker.summary.table[marker.summary.table$adj.pval < 0.05,]

marker.posthoc

write.csv(marker.summary.table,
          'snf_hvar_clusts/AnyPE_known_markers_hvar_cluster_table.csv', quote=F)

write.csv(marker.posthoc,
          'snf_hvar_clusts/AnyPE_known_markers_hvar_cluster_dunn_posthoc.csv',
          quote=F, row.names=F)


#### PTD only version for hvar clusters


summary.table <- cbind(clin.clusters[samp.ids, c('Condition.', 'hvar.Clusters')],
                       scaled.table) %>%
    as.data.frame %>%
    filter(Condition. %in% c('PTD')) %>%
    group_by(hvar.Clusters) %>%
    summarize_all(mean) %>%
    as.data.frame

summary.table <- summary.table[,-(1:2)]

rownames(summary.table) <- c(1,2,3,4)

annot.df <- data.frame(Clusters=factor(c(1,2,3,4),
                                       levels=c(1,2,3,4)))

## row.breaks <- cumsum(c(length(prot.names), length(rna.names), length(mirna.names)))

## ordered.names <- c(rna.names[order(summary.table[3,rna.names], decreasing=T)],
##                    mirna.names[order(summary.table[3,mirna.names], decreasing=T)],
##                    prot.names[order(summary.table[3,prot.names], decreasing=T)],
##                    metab.names[order(summary.table[3,metab.names], decreasing=T)])

## ordered.labels <- c(rna.markers[order(summary.table[3,rna.names], decreasing=T)],
##                     mirna.markers[order(summary.table[3,mirna.names], decreasing=T)],
##                     prot.markers[order(summary.table[3,prot.names], decreasing=T)],
##                     metab.markers[order(summary.table[3,metab.names], decreasing=T)])


max.val <- max(abs(summary.table))

mean.clust.map <- pheatmap(t(summary.table[,ordered.names]),
                           breaks = seq(-max.val, max.val, length.out=101),
         cluster_rows = F, cluster_cols=F, gaps_row=row.breaks, # scale='row',
         annotation_col=annot.df, show_colnames=F,
         labels_row=ordered.labels, annotation_colors=annot.colors)

ggsave('snf_hvar_clusts/PTD_known_markers_hvar_clusts_mean_heatmap.pdf', mean.clust.map$gtable, width=6, height=8)
ggsave('snf_hvar_clusts/PTD_known_markers_hvar_clusts_mean_heatmap.png', mean.clust.map$gtable, width=6, height=8, dpi=400)


markers.summary <- cbind(clin.clusters[samp.ids, c('Condition.', 'hvar.Clusters')],
                         scaled.table) %>%
    as.data.frame %>%
    filter(Condition. %in% c('PTD')) %>%
    group_by(hvar.Clusters) %>%
    select(-c('Condition.')) %>%
    summarize_all(c(mean=mean, sd=sd)) %>%
    as.data.frame

markers.summary <- markers.summary[,-1]

rownames(markers.summary) <- c(1,2,3,4)

markers.summary

markers.means <- markers.summary[,paste0(c(rna.names, mirna.names,
                                           prot.names, metab.names),
                                         '_mean')]
colnames(markers.means) <- c(rna.names, mirna.names, prot.names, metab.names)

markers.sds <- markers.summary[,paste0(c(rna.names, mirna.names,
                                           prot.names, metab.names),
                                         '_sd')]
colnames(markers.sds) <- c(rna.names, mirna.names, prot.names, metab.names)


marker.posthoc <- data.frame()

marker.summary.table <- data.frame(matrix(NA, 5, ncol(markers.means)))

colnames(marker.summary.table) <- ordered.names
rownames(marker.summary.table) <- c(paste('Cluster ', 1:4, ' (n = ',
                                          as.vector(table(clin.clusters[samp.ids[
                                              clin.clusters[samp.ids, 'Condition.'] %in%
                                              c('PTD')],
                                              'hvar.Clusters'])),
                                     ')', sep=''), 'pval')

library(FSA)

filt.scaled.table <- cbind(clin.clusters[samp.ids, c('Condition.', 'hvar.Clusters')],
                           scaled.table) %>%
    as.data.frame %>%
    filter(Condition. %in% c('PTD'))

for (var in ordered.names) {
    print(var)
    res <- kruskal.test(filt.scaled.table[,var], filt.scaled.table[,'hvar.Clusters'])
    print(res)

    marker.summary.table[1:4,var] <- paste(signif(markers.means[,var], 3),
                                             ' (', signif(markers.sds[,var], 3),
                                             ')', sep='')
    marker.summary.table['pval',var] <- res$p.value

    if (res$p.value < 0.05) {
        res.posthoc <- dunnTest(filt.scaled.table[,var],
                                factor(filt.scaled.table[,'hvar.Clusters']),
                                method='holm')
        print(res.posthoc)

        marker.posthoc <- rbind(marker.posthoc,
                                cbind(Feature=rep(var, nrow(res.posthoc$res)),
                                      res.posthoc$res))
    }

}

marker.summary.table <- as.data.frame(t(marker.summary.table))

marker.summary.table[,'adj.pval'] <- p.adjust(marker.summary.table$pval, method='fdr')

marker.summary.table

marker.summary.table[marker.summary.table$adj.pval < 0.05,]

marker.posthoc

write.csv(marker.summary.table,
          'snf_hvar_clusts/PTD_known_markers_hvar_cluster_table.csv', quote=F)

write.csv(marker.posthoc,
          'snf_hvar_clusts/PTD_known_markers_hvar_cluster_dunn_posthoc.csv',
          quote=F, row.names=F)

####

summary.table <- scaled.table %>% as.data.frame %>% group_by(clin.clusters[samp.ids,'Condition.']) %>% summarize_all(mean) %>% as.data.frame

summary.table <- summary.table[,-1]

rownames(summary.table) <- c('Control', 'Control PTD', 'FGR', 'FGR+HDP', 'PTD', 'Severe PE')

summary.table <- summary.table[c('Control', 'Control PTD', 'PTD',
                                 'FGR', 'FGR+HDP', 'Severe PE'),]

annot.df <- data.frame(Diagnosis=factor(c('Control', 'Control PTD', 'PTD',
                                          'FGR', 'FGR+HDP', 'Severe PE'),
                                       levels=c('Control', 'Control PTD', 'PTD',
                                                'FGR', 'FGR+HDP', 'Severe PE')),
                       row.names=c('Control', 'Control PTD', 'PTD',
                                   'FGR', 'FGR+HDP', 'Severe PE'))

row.breaks <- cumsum(c(length(rna.names), length(mirna.names), length(prot.names)))

ordered.names <- c(rna.names[order(summary.table['FGR+HDP',rna.names], decreasing=T)],
                   mirna.names[order(summary.table['FGR+HDP',mirna.names], decreasing=T)],
                   prot.names[order(summary.table['FGR+HDP',prot.names], decreasing=T)],
                   metab.names[order(summary.table['FGR+HDP',metab.names], decreasing=T)])

ordered.labels <- c(rna.markers[order(summary.table['FGR+HDP',rna.names], decreasing=T)],
                    mirna.markers[order(summary.table['FGR+HDP',mirna.names], decreasing=T)],
                    prot.markers[order(summary.table['FGR+HDP',prot.names], decreasing=T)],
                    metab.markers[order(summary.table['FGR+HDP',metab.names], decreasing=T)])

brewer.colors <- brewer.pal(n=12, "Paired")## )(10)

brewer.colors <- c(brewer.colors[1:2], brewer.colors[11], brewer.colors[5:6], brewer.colors[8], brewer.colors[3:4], brewer.colors[9:10])

annot.colors <- list()
annot.colors[['Diagnosis']] <- brewer.colors[1:6]
names(annot.colors[['Diagnosis']]) <- levels(annot.df$Diagnosis)

max.val <- max(abs(summary.table))

mean.cond.map <- pheatmap(t(summary.table[,ordered.names]),
                          breaks = seq(-max.val, max.val, length.out=101),
                          cluster_rows = F, cluster_cols=F, gaps_row=row.breaks,
                          ## scale='row',
                          annotation_col=annot.df, show_colnames=F,
                          labels_row=ordered.labels, annotation_colors=annot.colors)

ggsave('snf_hvar_clusts/known_markers_condition_mean_heatmap_FGR+HDPOrder.pdf', mean.cond.map$gtable, width=7, height=8)
ggsave('snf_hvar_clusts/known_markers_condition_mean_heatmap_FGR+HDPOrder.png', mean.cond.map$gtable, width=7, height=8, dpi=400)


markers.summary <- scaled.table %>% as.data.frame %>% 
    group_by(clin.clusters[samp.ids,'Condition.']) %>%
    summarize_all(c(mean=mean, sd=sd)) %>%
    as.data.frame

markers.summary <- markers.summary[,-1]

rownames(markers.summary) <- c('Control', 'Control PTD', 'FGR',
                               'FGR+HDP', 'PTD', 'Severe PE')

markers.summary <- markers.summary[c('Control', 'Control PTD', 'PTD',
                                     'FGR', 'FGR+HDP', 'Severe PE'),]


markers.summary

markers.means <- markers.summary[,paste0(c(rna.names, mirna.names,
                                           prot.names, metab.names),
                                         '_mean')]
colnames(markers.means) <- c(rna.names, mirna.names, prot.names, metab.names)

markers.sds <- markers.summary[,paste0(c(rna.names, mirna.names,
                                           prot.names, metab.names),
                                         '_sd')]
colnames(markers.sds) <- c(rna.names, mirna.names, prot.names, metab.names)


marker.posthoc <- data.frame()

marker.summary.table <- data.frame(matrix(NA, 7, ncol(markers.means)))

colnames(marker.summary.table) <- ordered.names
rownames(marker.summary.table) <- c(paste(c('Control', 'Control PTD', 'PTD',
                                            'FGR', 'FGR+HDP', 'Severe PE'),
                                          ' (n = ',
                                          as.vector(
                                              table(
                                                  clin.clusters[samp.ids,'Condition.'])[
                                                  c('Control', 'Control PTD', 'PTD',
                                                    'FGR', 'FGR+HDP', 'Severe PE')]),
                                          ')', sep=''), 'pval')

library(FSA)

for (var in ordered.names) {
    print(var)
    res <- kruskal.test(scaled.table[,var], clin.clusters[samp.ids,'Condition.'])
    print(res)

    marker.summary.table[1:6,var] <- paste(signif(markers.means[,var], 3),
                                             ' (', signif(markers.sds[,var], 3),
                                             ')', sep='')
    marker.summary.table['pval',var] <- res$p.value

    if (res$p.value < 0.05) {
        res.posthoc <- dunnTest(scaled.table[,var],
                                factor(clin.clusters[samp.ids,'Condition.']),
                                method='holm')
        print(res.posthoc)

        marker.posthoc <- rbind(marker.posthoc,
                                cbind(Feature=rep(var, nrow(res.posthoc$res)),
                                      res.posthoc$res))
    }

}

marker.summary.table <- as.data.frame(t(marker.summary.table))

marker.summary.table[,'adj.pval'] <- p.adjust(marker.summary.table$pval, method='fdr')

marker.summary.table

marker.summary.table[marker.summary.table$adj.pval > 0.05,]

marker.posthoc

write.csv(marker.summary.table,
          'snf_hvar_clusts/known_markers_condition_table.csv', quote=F)

write.csv(marker.posthoc,
          'snf_hvar_clusts/known_markers_condition_dunn_posthoc.csv',
          quote=F, row.names=F)



#### Cluster 2/3 feature correlation

prot <- read.csv('data/combat_prot_data.csv', row.names=1)
colnames(prot) <- paste0('Prot_', colnames(prot))

metab <- read.csv('data/combat_metab_data.csv', row.names=1)
colnames(metab) <- paste0('Metab_', colnames(metab))

rna <- read.csv('data/combat_vst_rna_expression.csv', row.names=1)
colnames(rna) <- paste0('RNA_', colnames(rna))

mirna <- read.csv('data/combat_vst_mirna_expression.csv', row.names=1)
colnames(mirna) <- paste0('miRNA_', colnames(mirna))

library(ggpubr)
library(cowplot)

prot.summary <- cbind(scale(prot[samp.ids,]), Clusters=clin[samp.ids,'hvar.Clusters']) %>%
    as.data.frame %>%
    group_by(Clusters) %>%
    summarize_all(mean)

prot.summary <- t(prot.summary[,-1]) %>% as.data.frame
colnames(prot.summary) <- paste0('Cluster', 1:4)

p1 <- ggplot(prot.summary, aes(Cluster1, Cluster4)) +
    geom_point(color='darkgrey') +
    stat_cor(method='spearman', cor.coef.name='R', label.y.npc = "bottom") +
    stat_smooth(method = "lm", color='black', lty=5, se=FALSE) + 
    labs(title='Protein Expression',
         x='Mean Expression in Cluster 1',
         y='Mean Expression in Cluster 4') +
    theme_bw()


rna.summary <- cbind(scale(rna[samp.ids,]), Clusters=clin[samp.ids,'hvar.Clusters']) %>%
    as.data.frame %>%
    group_by(Clusters) %>%
    summarize_all(mean)

rna.summary <- t(rna.summary[,-1]) %>% as.data.frame
colnames(rna.summary) <- paste0('Cluster', 1:4)

p2 <- ggplot(rna.summary, aes(Cluster1, Cluster4)) +
    geom_point(color='darkgrey') +
    stat_cor(method='spearman', cor.coef.name='R', label.y.npc = "bottom") +
    stat_smooth(method = "lm", color='black', lty=5, se=FALSE) + 
    labs(title='Gene Expression',
         x='Mean Expression in Cluster 1',
         y='Mean Expression in Cluster 4') +
    theme_bw()


metab.summary <- cbind(scale(metab[samp.ids,]), Clusters=clin[samp.ids,'hvar.Clusters']) %>%
    as.data.frame %>%
    group_by(Clusters) %>%
    summarize_all(mean)

metab.summary <- t(metab.summary[,-1]) %>% as.data.frame
colnames(metab.summary) <- paste0('Cluster', 1:4)

p3 <- ggplot(metab.summary, aes(Cluster1, Cluster4)) +
    geom_point(color='darkgrey') +
    stat_cor(method='spearman', cor.coef.name='R', label.y.npc = "bottom") +
    stat_smooth(method = "lm", color='black', lty=5, se=FALSE) + 
    labs(title='Metabolite Expression',
         x='Mean Expression in Cluster 1',
         y='Mean Expression in Cluster 4') +
    theme_bw()


mirna.summary <- cbind(scale(mirna[samp.ids,]), Clusters=clin[samp.ids,'hvar.Clusters']) %>%
    as.data.frame %>%
    group_by(Clusters) %>%
    summarize_all(mean)

mirna.summary <- t(mirna.summary[,-1]) %>% as.data.frame
colnames(mirna.summary) <- paste0('Cluster', 1:4)

p4 <- ggplot(mirna.summary, aes(Cluster1, Cluster4)) +
    geom_point(color='darkgrey') +
    stat_cor(method='spearman', cor.coef.name='R', label.y.npc = "bottom") +
    stat_smooth(method = "lm", color='black', lty=5, se=FALSE) + 
    labs(title='miRNA Expression',
         x='Mean Expression in Cluster 1',
         y='Mean Expression in Cluster 4') +
    theme_bw()


plot_grid(p1, p2, p3, p4, nrow=2, ncol=2)

ggsave('cluster_1v4_mean_expression.png', plot_grid(p1, p2, p3, p4, nrow=2, ncol=2), height=9, width=9)

ggsave('cluster_1v4_mean_expression.pdf', plot_grid(p1, p2, p3, p4, nrow=2, ncol=2), height=9, width=9)


prot.summary <- cbind(scale(prot[samp.ids,]), Clusters=clin[samp.ids,'hvar.Clusters']) %>%
    as.data.frame %>%
    group_by(Clusters) %>%
    summarize_all(mean)

prot.summary <- t(prot.summary[,-1]) %>% as.data.frame
colnames(prot.summary) <- paste0('Cluster', 1:4)

p1 <- ggplot(prot.summary, aes(Cluster2, Cluster3)) +
    geom_point(color='darkgrey') +
    stat_cor(method='spearman', cor.coef.name='R', label.y.npc = "bottom") +
    stat_smooth(method = "lm", color='black', lty=5, se=FALSE) + 
    labs(title='Protein Expression',
         x='Mean Expression in Cluster 2',
         y='Mean Expression in Cluster 3') +
    theme_bw()


rna.summary <- cbind(scale(rna[samp.ids,]), Clusters=clin[samp.ids,'hvar.Clusters']) %>%
    as.data.frame %>%
    group_by(Clusters) %>%
    summarize_all(mean)

rna.summary <- t(rna.summary[,-1]) %>% as.data.frame
colnames(rna.summary) <- paste0('Cluster', 1:4)

p2 <- ggplot(rna.summary, aes(Cluster2, Cluster3)) +
    geom_point(color='darkgrey') +
    stat_cor(method='spearman', cor.coef.name='R', label.y.npc = "bottom") +
    stat_smooth(method = "lm", color='black', lty=5, se=FALSE) + 
    labs(title='Gene Expression',
         x='Mean Expression in Cluster 2',
         y='Mean Expression in Cluster 3') +
    theme_bw()


metab.summary <- cbind(scale(metab[samp.ids,]), Clusters=clin[samp.ids,'hvar.Clusters']) %>%
    as.data.frame %>%
    group_by(Clusters) %>%
    summarize_all(mean)

metab.summary <- t(metab.summary[,-1]) %>% as.data.frame
colnames(metab.summary) <- paste0('Cluster', 1:4)

p3 <- ggplot(metab.summary, aes(Cluster2, Cluster3)) +
    geom_point(color='darkgrey') +
    stat_cor(method='spearman', cor.coef.name='R', label.y.npc = "bottom") +
    stat_smooth(method = "lm", color='black', lty=5, se=FALSE) + 
    labs(title='Metabolite Expression',
         x='Mean Expression in Cluster 2',
         y='Mean Expression in Cluster 3') +
    theme_bw()


mirna.summary <- cbind(scale(mirna[samp.ids,]), Clusters=clin[samp.ids,'hvar.Clusters']) %>%
    as.data.frame %>%
    group_by(Clusters) %>%
    summarize_all(mean)

mirna.summary <- t(mirna.summary[,-1]) %>% as.data.frame
colnames(mirna.summary) <- paste0('Cluster', 1:4)

p4 <- ggplot(mirna.summary, aes(Cluster2, Cluster3)) +
    geom_point(color='darkgrey') +
    stat_cor(method='spearman', cor.coef.name='R', label.y.npc = "bottom") +
    stat_smooth(method = "lm", color='black', lty=5, se=FALSE) + 
    labs(title='miRNA Expression',
         x='Mean Expression in Cluster 2',
         y='Mean Expression in Cluster 3') +
    theme_bw()


plot_grid(p1, p2, p3, p4, nrow=2, ncol=2)

ggsave('cluster_2v3_mean_expression.png', plot_grid(p1, p2, p3, p4, nrow=2, ncol=2), height=9, width=9)

ggsave('cluster_2v3_mean_expression.pdf', plot_grid(p1, p2, p3, p4, nrow=2, ncol=2), height=9, width=9)


alphas <- seq(0.5, 1, 0.05)
lambdas <- exp(seq(log(0.45), log(0.045), length.out=50))
## types <- c('ungrouped')

prob.mat <- array(NA, c(length(samp.ids), C, length(lambdas), length(alphas)))
pred.mat <- array(NA, c(length(samp.ids), length(lambdas), length(alphas)))

score.auc <- array(NA, c(K, length(lambdas), length(alphas)))
score.bacc <- array(NA, c(K, length(lambdas), length(alphas)))

for (t.idx in 1:length(types)) {
    for (a.idx in 1:length(alphas)) {
        for (k in 1:K) {
            idx <- which(foldid[samp.ids]==k)
            
            class.wts <- nrow(concat.train[[k]]) /
                (C * table(labels[-idx,k]))
            case.wts <- class.wts[labels[-idx,k]]

            res <- glmnet(model.matrix(~.,concat.train[[k]])[,-1],
                          labels[-idx,k],
                          lambda=lambdas/alphas[a.idx], alpha=alphas[a.idx],
                          family='multinomial', weights=as.vector(case.wts),
                          type.multinomial = types[t.idx])

            prob.mat[idx,,,a.idx,t.idx] <- predict(res,
                                                   newx=model.matrix(~.,concat.test[[k]])[,-1],
                                                   type='response')

            pred.mat[idx,,a.idx,t.idx] <- predict(res,
                                                  newx=model.matrix(~.,concat.test[[k]])[,-1],
                                                  type='class')

            for (l.idx in 1:length(lambdas)) {

                probs <- prob.mat[idx,,l.idx,a.idx,t.idx]
                ##rownames(probs) <- samp.ids[idx]
                colnames(probs) <- c(1,2,3)
                
                roc.res <- multiclass.roc(labels[idx,k], probs)

                print(roc.res$auc)

                score.auc[k,l.idx,a.idx,t.idx] <- as.numeric(roc.res$auc)

                cmat <- confusionMatrix(factor(pred.mat[idx,l.idx,a.idx,t.idx]),
                                        factor(labels[idx,k]))

                ## print(cmat)

                score.bacc[k,l.idx,a.idx,t.idx] <- mean(cmat$byClass[,'Balanced Accuracy'])
            }

        }
    }
}

cv.mean.bacc <- apply(score.bacc, c(3,4), mean)
cv.sd.bacc <- apply(score.bacc, c(3,4), sd)

which(cv.mean.bacc == max(cv.mean.bacc), arr.ind=T)

which(cv.mean.bacc==min(cv.mean.bacc[which(cv.mean.bacc > min(max(cv.mean.bacc)-cv.sd.bacc[which(cv.mean.bacc == max(cv.mean.bacc), arr.ind=T)]/sqrt(10)), arr.ind=T)]), arr.ind=T)


cv.mean.auc <- apply(score.auc, c(3,4), mean)
cv.sd.auc <- apply(score.auc, c(3,4), sd)

which(cv.mean.auc == max(cv.mean.auc), arr.ind=T)

which(cv.mean.auc==min(cv.mean.auc[which(cv.mean.auc > min(max(cv.mean.auc)-cv.sd.auc[which(cv.mean.auc == max(cv.mean.auc), arr.ind=T)]/sqrt(10)), arr.ind=T)]), arr.ind=T)

