library(dplyr)
library(stringr)
library(limma)
library(sva)
library(ggplot2)
library(pheatmap)
library(readxl)
source('utils.R')

## clin <- read.csv('clinical-211205.csv', na.strings='#N/A', row.names=1)

## rownames(clin) <- str_replace_all(rownames(clin), '-', '.')

## head(clin)

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

final.clin <- clin %>% select(Condition., Batch, WksGest, InfSex, Race, PrePregBMI,
                              Smoking, FDELTYPE, Labor.initiation, Labor.no.labor)

final.clin$Race[!final.clin$Race %in% c('W', 'B')] <- 'other'
final.clin$FDELTYPE[final.clin$FDELTYPE==0] <- NA

final.clin %>% mutate_at(c('Race', 'FDELTYPE', 'Labor.initiation', 'Labor.no.labor'), factor)

head(final.clin)


pathology <- read_excel('Path diagnoses for analysis.xlsx' , 1) %>% as.data.frame
rownames(pathology) <- make.names(pathology$case)
pathology <- pathology[,11:15]
colnames(pathology) <- c("High.grade.MVM", "FVM",
                         "Acute.inflammation", "Chronic.inflammation", "VUE")

pathology$VUE <- ifelse(pathology$VUE>0, 1, 0)



pathology <- pathology %>% mutate_at(c("High.grade.MVM", "FVM",
                                       "Acute.inflammation", "Chronic.inflammation",
                                       "VUE"),
                                     as.factor)
head(pathology)
dim(pathology)

final.clin[rownames(pathology),
           c("High.grade.MVM", "FVM", "Acute.inflammation", "Chronic.inflammation", "VUE")] <- pathology

final.clin$Control.VUE <- NA

final.clin$Control.VUE <- ifelse(final.clin$VUE==1, "VUE", "No VUE")
final.clin$Control.VUE <- ifelse(final.clin$VUE==0 & final.clin$Condition. %in% c('Control', 'Control PTD'), "Control No VUE", final.clin$Control.VUE)

final.clin$Control.VUE <- factor(final.clin$Control.VUE)

head(final.clin)

#### IMPORT AND FILTER PROTEIN DATA ####


prot <- read.tcsv('olink-prot-211205.csv', na.strings=c('#N/A', ''), sep=',')

prot.meta <- prot %>% select(c('Panel', 'Assay', 'Uniprot.ID', 'OlinkID', 'Missing.Data.freq.', 'LOD'))

samp.list <- intersect(rownames(final.clin), colnames(prot))

prot.exp <- prot %>% select(all_of(samp.list))

prot.quality <- colSums(tail(prot.exp,10)=='Warning')==0 # samples that passed all quality checks

prot.meta <- prot.meta[2:(nrow(prot.meta)-10),]

prot.meta$LOD <- as.numeric(prot.meta$LOD)

prot.exp <- prot.exp[2:(nrow(prot.exp)-10),]

rownames(prot.meta) <- prot.meta$OlinkID
rownames(prot.exp) <- prot.meta$OlinkID

prot.exp <- prot.exp %>% mutate_all(as.numeric)

samp.list <- samp.list[prot.quality]

prot.exp <- prot.exp[,samp.list]

head(prot.exp)

cond.totals <- table(factor(final.clin[samp.list,'Condition.']))

cond.totals

prots2remove <- c()

for (prot in rownames(prot.exp)) {
    
    temp <- table(ifelse(prot.exp[prot,] < prot.meta[prot,'LOD'], 'below', 'above'), factor(final.clin[samp.list,'Condition.']))

    if (nrow(temp) == 1) {
        if (rownames(temp) == 'below') {
            prots2remove <- c(prots2remove, prot)
        }
    } else {
        if(all(temp['below',] / cond.totals > 0.5)) {
            prots2remove <- c(prots2remove, prot)
        }
    }
}

prots2remove

prot.exp.clean <- t(prot.exp) %>% as.data.frame %>% select(-prots2remove)

prots2keep <- colnames(prot.exp.clean)

prots2keep

prot.exp.clean <- as.data.frame(t(avereps(t(prot.exp.clean), ID=prot.meta[colnames(prot.exp.clean),'Assay'])))

prot.data <- cbind(final.clin[samp.list,], prot.exp.clean[samp.list,])

pcs <- prcomp(t(scale(prot.data[,-c(1:16)])))

summary(pcs)

pc.plot <- pcs$rotation[,1:2]

pc.plot <- cbind(pc.plot, prot.data[,1:9])

ggplot(pc.plot, aes(x=PC1, y=PC2, col=Condition.)) +
    geom_point() +
    geom_vline(xintercept=mean(pc.plot$PC1)-4*sd(pc.plot$PC1)) +
    geom_vline(xintercept=mean(pc.plot$PC1)+4*sd(pc.plot$PC1)) +
    geom_hline(yintercept=mean(pc.plot$PC2)-4*sd(pc.plot$PC2)) +
    geom_hline(yintercept=mean(pc.plot$PC2)+4*sd(pc.plot$PC2)) +
    ggtitle('PC1 vs PC2 by Batch') +
    theme_bw()

outliers <- pc.plot %>% filter(abs(PC1) > 4*sd(PC1) | abs(PC2) > 4*sd(PC2))
rownames(outliers)

outlier.samps <- c('MJ.0435', 'MJ.0500', 'MJ.0622', 'MJ.1001', 'Mini.DP.053', 'Mini.DP.129')

prot.data.no.out <- prot.data[!rownames(prot.data) %in% rownames(outliers),]

write.csv(prot.data.no.out[,-(1:16)], 'data/raw_prot_data.csv')

pcs <- prcomp(t(scale(prot.data.no.out[,-c(1:16)])))

summary(pcs)

plot(pcs$sdev[1:50])

pc.plot <- pcs$rotation[,1:10]

pc.plot <- cbind(pc.plot, prot.data.no.out[,1:16])

summary(lm(WksGest ~ 1 + PC1 + PC2 + PC3 + PC4 + PC5 +
               PC6 + PC7 + PC8 + PC9 + PC10, pc.plot))

## Call:
## lm(formula = WksGest ~ 1 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + 
##     PC7 + PC8 + PC9 + PC10, data = pc.plot)

## Residuals:
##    Min     1Q Median     3Q    Max 
## -8.878 -1.036  0.112  1.399  6.055 

## Coefficients:
##              Estimate Std. Error t value Pr(>|t|)    
## (Intercept)  36.87900    0.11481 321.226  < 2e-16 ***
## PC1          -0.18868    2.09818  -0.090    0.928    
## PC2           8.68898    2.09818   4.141 4.41e-05 ***
## PC3          -0.02634    2.09818  -0.013    0.990    
## PC4           2.14691    2.09818   1.023    0.307    
## PC5         -21.24627    2.09818 -10.126  < 2e-16 ***
## PC6         -25.01389    2.09818 -11.922  < 2e-16 ***
## PC7          -3.03343    2.09818  -1.446    0.149    
## PC8          -2.28333    2.09818  -1.088    0.277    
## PC9         -23.86675    2.09818 -11.375  < 2e-16 ***
## PC10         -3.23130    2.09818  -1.540    0.125    
## ---
## Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## Residual standard error: 2.098 on 323 degrees of freedom
## Multiple R-squared:  0.552,	Adjusted R-squared:  0.5381 
## F-statistic: 39.79 on 10 and 323 DF,  p-value: < 2.2e-16

summary(glm(Batch=='Mini-DP' ~ 1 + PC1 + PC2 + PC3 + PC4 + PC5 +
                PC6 + PC7 + PC8 + PC9 + PC10, pc.plot, family='binomial'))

## Call:
## glm(formula = Batch == "Mini-DP" ~ 1 + PC1 + PC2 + PC3 + PC4 + 
##     PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family = "binomial", 
##     data = pc.plot)

## Deviance Residuals: 
##      Min        1Q    Median        3Q       Max  
## -3.12632  -0.03422   0.06403   0.21218   2.94199  

## Coefficients:
##             Estimate Std. Error z value Pr(>|z|)    
## (Intercept)   2.8537     0.4124   6.919 4.55e-12 ***
## PC1         -48.3938     6.9253  -6.988 2.79e-12 ***
## PC2          47.9125     8.5099   5.630 1.80e-08 ***
## PC3         -20.5880     5.4386  -3.786 0.000153 ***
## PC4          26.3042     6.0991   4.313 1.61e-05 ***
## PC5         -16.9117     5.0260  -3.365 0.000766 ***
## PC6           1.9396     4.6946   0.413 0.679491    
## PC7          -6.0008     5.5366  -1.084 0.278439    
## PC8          30.0616     5.8475   5.141 2.73e-07 ***
## PC9          -2.0211     4.1213  -0.490 0.623849    
## PC10          1.8959     5.0624   0.375 0.708027    
## ---
## Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## (Dispersion parameter for binomial family taken to be 1)

##     Null deviance: 387.25  on 333  degrees of freedom
## Residual deviance: 115.71  on 323  degrees of freedom
## AIC: 137.71

## Number of Fisher Scoring iterations: 7


pdf('uncorrected_prot_noOut_PCA_plots.pdf', width=6, height=5)
ggplot(pc.plot, aes(x=PC1, y=PC2, col=Batch)) +
    geom_point() +
    ggtitle('PC1 vs PC2 by Batch') +
    theme_bw()

ggplot(pc.plot, aes(x=PC1, y=PC2, col=Condition.)) +
    geom_point() +
    ggtitle('PC1 vs PC2 by Condition') +
    theme_bw()

ggplot(pc.plot, aes(x=PC1, y=PC2, col=WksGest)) +
    geom_point() +
    ggtitle('PC1 vs PC2 by Gestational Age') +
    theme_bw()

ggplot(pc.plot, aes(x=PC1, y=PC2, col=InfSex)) +
    geom_point() +
    ggtitle('PC1 vs PC2 by Infant Sex') +
    theme_bw()

ggplot(pc.plot, aes(x=PC1, y=PC2, col=factor(Labor.no.labor))) +
    geom_point() +
    ggtitle('PC1 vs PC2 by Labor') +
    theme_bw()

dev.off()


prot.data.no.out$Race[!prot.data.no.out$Race %in% c('W', 'B')] <- 'other'
prot.data.no.out$FDELTYPE[prot.data.no.out$FDELTYPE==0] <- NA

mod <- model.matrix(~Condition. + WksGest + InfSex + Race + PrePregBMI + Smoking + factor(FDELTYPE) + factor(Labor.initiation), prot.data.no.out)

## prot.data.corrected <- removeBatchEffect(t(prot.data[,-c(1:9)]), batch=prot.data$Batch, design=mod)

prot.data.corrected <- ComBat(t(prot.data.no.out[rownames(mod),-c(1:16)]), batch=prot.data.no.out[rownames(mod),'Batch'], mod=mod)

prot.data.combat <- cbind(prot.data.no.out[rownames(mod),1:16], t(prot.data.corrected))

write.csv(prot.data.combat[,17:ncol(prot.data.combat)], 'data/combat_prot_data.csv')

## prot.data.combat <- read.csv('data/combat_prot_data.csv', row.names=1)

## prot.data.combat <- cbind(final.clin[rownames(prot.data.combat),], prot.data.combat)


pcs <- prcomp(t(scale(t(prot.data.corrected))))

summary(pcs)

plot(pcs$sdev[1:50])

pc.plot <- pcs$rotation[,1:10]

pc.plot <- cbind(pc.plot, prot.data.combat[,1:16])


summary(lm(WksGest ~ 1 + PC1 + PC2 + PC3 + PC4 + PC5 +
               PC6 + PC7 + PC8 + PC9 + PC10 + Condition. +
               InfSex + Race + PrePregBMI + Smoking + factor(FDELTYPE) +
               factor(Labor.initiation),
           pc.plot))

## Call:
## lm(formula = WksGest ~ 1 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + 
##     PC7 + PC8 + PC9 + PC10 + Condition. + InfSex + Race + PrePregBMI + 
##     Smoking + factor(FDELTYPE), data = pc.plot)

## Residuals:
##     Min      1Q  Median      3Q     Max 
## -7.6716 -0.9459  0.1058  1.1144  3.9222 

## Coefficients:
##                        Estimate Std. Error t value Pr(>|t|)    
## (Intercept)            38.73200    0.58092  66.673  < 2e-16 ***
## PC1                     6.82458    1.88606   3.618 0.000349 ***
## PC2                    -2.95926    1.91932  -1.542 0.124200    
## PC3                    -2.47854    1.79918  -1.378 0.169385    
## PC4                     5.36059    1.85639   2.888 0.004171 ** 
## PC5                    10.16410    1.87568   5.419 1.26e-07 ***
## PC6                   -16.44244    2.21070  -7.438 1.15e-12 ***
## PC7                    -1.97151    1.88908  -1.044 0.297518    
## PC8                   -16.68198    2.02295  -8.246 5.62e-15 ***
## PC9                    -8.33267    1.88788  -4.414 1.43e-05 ***
## PC10                   -1.36092    1.82458  -0.746 0.456339    
## Condition.Control PTD  -0.75443    0.56809  -1.328 0.185209    
## Condition.FGR          -1.39457    0.38628  -3.610 0.000360 ***
## Condition.FGR+HDP      -3.83264    0.50431  -7.600 4.06e-13 ***
## Condition.PTD          -3.20037    0.33596  -9.526  < 2e-16 ***
## Condition.Severe PE    -2.26658    0.33551  -6.756 7.70e-11 ***
## InfSexM                 0.02624    0.20703   0.127 0.899242    
## Raceother               0.61543    0.42192   1.459 0.145742    
## RaceW                   0.10811    0.27782   0.389 0.697454    
## PrePregBMI             -0.01189    0.01629  -0.730 0.466193    
## SmokingYes              0.19947    0.22443   0.889 0.374863    
## factor(FDELTYPE)2      -0.16325    0.22605  -0.722 0.470750    
## ---
## Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## Residual standard error: 1.758 on 292 degrees of freedom
## Multiple R-squared:  0.701,	Adjusted R-squared:  0.6795 
## F-statistic: 32.61 on 21 and 292 DF,  p-value: < 2.2e-16

summary(glm(Batch=='Mini-DP' ~ 1 + PC1 + PC2 + PC3 + PC4 + PC5 +
                PC6 + PC7 + PC8 + PC9 + PC10 + Condition. + WksGest +
                InfSex + Race + PrePregBMI + Smoking + factor(FDELTYPE) +
                factor(Labor.initiation),
            pc.plot, family='binomial'))

## Call:
## glm(formula = Batch == "Mini-DP" ~ 1 + PC1 + PC2 + PC3 + PC4 + 
##     PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + Condition. + WksGest + 
##     InfSex + Race + PrePregBMI + Smoking + factor(FDELTYPE), 
##     family = "binomial", data = pc.plot)

## Deviance Residuals: 
##     Min       1Q   Median       3Q      Max  
## -2.6617  -0.3613   0.4115   0.6329   2.2568  

## Coefficients:
##                        Estimate Std. Error z value Pr(>|z|)    
## (Intercept)           -3.796151   3.556137  -1.067 0.285749    
## PC1                   -1.569439   2.898733  -0.541 0.588216    
## PC2                   -0.775615   3.033885  -0.256 0.798220    
## PC3                    1.891502   2.705843   0.699 0.484525    
## PC4                   -1.598881   2.876863  -0.556 0.578366    
## PC5                    1.957220   3.027927   0.646 0.518027    
## PC6                   -0.096095   3.873541  -0.025 0.980208    
## PC7                    1.037867   2.964915   0.350 0.726302    
## PC8                    1.208084   3.425300   0.353 0.724318    
## PC9                   -3.502921   3.049935  -1.149 0.250753    
## PC10                   0.091388   2.961405   0.031 0.975382    
## Condition.Control PTD -0.747852   0.757333  -0.987 0.323407    
## Condition.FGR         -1.236887   0.570136  -2.169 0.030048 *  
## Condition.FGR+HDP     -2.027711   0.825176  -2.457 0.013998 *  
## Condition.PTD          0.139961   0.625925   0.224 0.823063    
## Condition.Severe PE    0.452684   0.594453   0.762 0.446350    
## WksGest                0.129589   0.087992   1.473 0.140823    
## InfSexM               -0.043614   0.329533  -0.132 0.894708    
## Raceother              1.771328   0.715969   2.474 0.013360 *  
## RaceW                  1.234185   0.409779   3.012 0.002597 ** 
## PrePregBMI             0.009831   0.026236   0.375 0.707864    
## SmokingYes            -0.596797   0.329620  -1.811 0.070209 .  
## factor(FDELTYPE)2     -1.239814   0.342988  -3.615 0.000301 ***
## ---
## Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## (Dispersion parameter for binomial family taken to be 1)

##     Null deviance: 360.63  on 313  degrees of freedom
## Residual deviance: 265.25  on 291  degrees of freedom
## AIC: 311.25

## Number of Fisher Scoring iterations: 5


pdf('corrected_prot_noOut_PCA_plots.pdf', width=6, height=5)
gg <- ggplot(pc.plot, aes(x=PC1, y=PC2, col=Batch)) +
    geom_point() +
    ggtitle('PC1 vs PC2 by Batch') +
    theme_bw()
gg

ggplot(pc.plot, aes(x=PC1, y=PC2, col=Condition.)) +
    geom_point() +
    ggtitle('PC1 vs PC2 by Condition') +
    theme_bw()

ggplot(pc.plot, aes(x=PC1, y=PC2, col=WksGest)) +
    geom_point() +
    ggtitle('PC1 vs PC2 by Gestational Age') +
    theme_bw()

ggplot(pc.plot, aes(x=PC1, y=PC2, col=InfSex)) +
    geom_point() +
    ggtitle('PC1 vs PC2 by Infant Sex') +
    theme_bw()
dev.off()

groups <- c('Control PTD', 'Control')

run_limma<-function(groups) {
    ## extract data by group
    dataset_group<-prot.data.combat ## [prot.data.combat$Condition. %in% groups,]
    print(table(dataset_group$Condition.))
    ## design a model - control will always be on top as 0

    #### COND GROUP 1 ####
    ## design<-model.matrix(~0 + as.factor(dataset_group$Condition.) + dataset_group$WksGest + as.factor(dataset_group$InfSex))
    ## colnames(design)<-c("Control", "Control.PTD", "FGR", "FGR.HDP", "PTD", "Severe.PE", "GestationalAge", "InfantSex")

    #### COND GROUP 2 ####
    ## design<-model.matrix(~0 + as.factor(dataset_group$Condition.) + dataset_group$WksGest + as.factor(dataset_group$InfSex) + as.factor(dataset_group$Race) + dataset_group$PrePregBMI + as.factor(dataset_group$Smoking) + as.factor(dataset_group$FDELTYPE) + as.factor(dataset_group$Labor.initiation))
    ## ## colnames(design)<-c("case", "control", "GestationalAge", "InfantSex", "RaceOther", "RaceW", "PrePregBMI", "Smoking", "DeliveryType")
    ## colnames(design)<-c("Control", "Control.PTD", "FGR", "FGR.HDP", "PTD", "Severe.PE", "GestationalAge", "InfantSex", "RaceOther", "RaceW", "PrePregBMI", "Smoking", "DeliveryType", "LaborInitiation")

    #### COND GROUP 3 ####
    design<-model.matrix(~0 + as.factor(dataset_group$Condition.) + dataset_group$WksGest + as.factor(dataset_group$InfSex) + as.factor(dataset_group$Race) + dataset_group$PrePregBMI + as.factor(dataset_group$Smoking) + as.factor(dataset_group$FDELTYPE) + as.factor(dataset_group$Labor.initiation) + as.factor(dataset_group$Labor.no.labor))
    ## colnames(design)<-c("case", "control", "GestationalAge", "InfantSex", "RaceOther", "RaceW", "PrePregBMI", "Smoking", "DeliveryType")
    colnames(design)<-c("Control", "Control.PTD", "FGR", "FGR.HDP", "PTD", "Severe.PE", "GestationalAge", "InfantSex", "RaceOther", "RaceW", "PrePregBMI", "Smoking", "DeliveryType", "LaborInitiation", "Labor")

#### VUE all ####
    ## design<-model.matrix(~0 + VUE + WksGest + factor(InfSex) + factor(Race) + PrePregBMI + factor(Smoking) + factor(FDELTYPE) + factor(Labor.initiation), dataset_group)
    ## colnames(design)<-c("VUE0", "VUE1", "GestationalAge", "InfantSex", "RaceOther", "RaceW", "PrePregBMI", "Smoking", "DeliveryType", "LaborInitiation")
    ## make contrast - what to compare
    contrast<- makeContrasts(Diff = paste(str_replace(groups, '[+ ]', '.'), collapse=" - "), levels=design)
    ## contrast<- makeContrasts(Diff = "VUE1 - VUE0", levels=design)
    print(contrast)
    ## apply linear model to each protein
    ## Robust regression provides an alternative to least squares regression that works with less restrictive assumptions. Specifically, it provides much better regression coefficient estimates when outliers are present in the data
    fit<-lmFit(t(dataset_group[,17:ncol(dataset_group)]), design=design, maxit=1000)
    ## apply contrast
    contrast_fit<-contrasts.fit(fit, contrast)
    ## apply empirical Bayes smoothing to the SE
    ebays_fit<-eBayes(contrast_fit, robust=T)
    ## summary
    print(summary(decideTests(ebays_fit)))
    ## extract DE results
    DE_results<-topTable(ebays_fit, n=ncol(dataset_group), adjust.method="fdr", confint=TRUE)

    ## annot_df <- data.frame(Condition = as.factor(dataset_group$Condition.),
    ##                        Infant.Sex = as.factor(dataset_group$InfSex),
    ##                        Gestational.Age = dataset_group$WksGest, 
    ##                        row.names=rownames(dataset_group))

    ## if (sum(DE_results$adj.P.Val < 0.05)!=0) {
    ##     scaled.sig <- scale(dataset_group[rownames(DE_results)[DE_results$adj.P.Val < 0.05]])
    ##     range <- min(max(abs(scaled.sig)),4.5)
    ##     pheatmap(t(scaled.sig), breaks = seq(-range, range, length.out = 101),
    ##              annotation_col=annot_df, main=paste(groups[2], 'vs.', groups[1]),
    ##              clustering_method='ward', 
    ##              labels_row=prot.meta[rownames(DE_results)[DE_results$adj.P.Val < 0.05],'Assay'])
    ## }
    
    return(DE_results)
}


comparisons <- list(c('Control PTD', 'Control'),
                    c('FGR', 'Control'),
                    c('Severe PE', 'Control PTD'),
                    c('Severe PE', 'Control'),
                    c('PTD', 'Control PTD'),
                    c('PTD', 'Control'),
                    c('FGR+HDP', 'Control PTD'),
                    c('FGR+HDP', 'Control'))

## pdf('lmFit_robust_prot_diffExp_heatmaps.pdf', width=12, height=10)
## for (groups in comparisons) {
##     print(groups)
##     de <- run_limma(groups)
##     paste('diffExpProt_', str_replace(groups[1], ' ', '-'), '_',
##           str_replace(groups[2], ' ', '.'), '.csv', sep='')
##     de <- cbind(Protein=prot.meta[rownames(de),'Assay'], de)
##     write.csv(de, paste('diffExpProt_', str_replace(groups[1], ' ', '-'), '_',
##                         str_replace(groups[2], ' ', '.'), '.csv', sep=''))
## }
## dev.off()

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
    de <- run_limma(groups)
    print(paste('DE-cond3/diffExpProt_', str_replace(groups[1], ' ', '-'), '_',
                str_replace(groups[2], ' ', '.'), '.csv', sep=''))

    if (sum(de$adj.P.Val < 0.05)!=0) {

        if (min(nrow(de[!is.na(de$adj.P.Val) & de$adj.P.Val<0.05,]), 500) == 500) {
            col.labels <- rownames(de)[1:500]
            col.idxs <- rownames(de)[1:500]
        } else {
            col.labels <- rownames(de)[de$adj.P.Val < 0.05]
            col.idxs <- rownames(de)[de$adj.P.Val < 0.05]
        }

        scaled.sig <- scale(prot.data.combat[prot.data.combat$Condition. %in% make.names(groups), col.idxs])    
        
        range <- max(abs(scaled.sig))

        print(range)

        range <- min(max(abs(scaled.sig)),3)
        
        annot_df <- data.frame(Diagnosis = factor(diag.map[prot.data.combat[prot.data.combat$Condition.
                                                                             %in% make.names(groups), 'Condition.']], levels=groups),
                               `Infant Sex` = as.factor(prot.data.combat[prot.data.combat$Condition.
                                                                        %in% make.names(groups), 'InfSex']),
                               `Gestational Age` = prot.data.combat[prot.data.combat$Condition.
                                                                     %in% make.names(groups), 'WksGest'], 
                               row.names=rownames(prot.data.combat[prot.data.combat$Condition.
                                                                    %in% make.names(groups),]))

        temp.annot.cols <- list()
        for (g in groups) {
            temp.annot.cols[['Diagnosis']][g] <- annot.colors[['Diagnosis']][g]
        }
       
        pdf(paste('plots-cond3/', str_replace(groups[1], ' ', '-'), '_',
                  str_replace(groups[2], ' ', '.'),
                  '_prot_diffExp_heatmap.pdf',sep=''),
            width=10,
            height=3 + 0.12 * min(nrow(de[!is.na(de$adj.P.Val) & de$adj.P.Val<0.05,]), 500))

        pheatmap(t(scaled.sig), breaks = seq(-range, range, length.out = 101),
                 annotation_col=annot_df, main=paste(groups[1], 'vs.', groups[2]),
                 ## cluster_rows=sum(de$adj.P.Val < 0.05)>1,
                 clustering_method='ward.D', clustering_distance_rows = "euclidean",
                 clustering_distance_cols = "euclidean", annotation_colors=temp.annot.cols,
                 labels_row=col.labels, show_colnames=F)
        dev.off()

        png(paste('plots-cond3/', str_replace(groups[1], ' ', '-'), '_',
                  str_replace(groups[2], ' ', '.'),
                  '_prot_diffExp_heatmap.png',sep=''),
            width=10,
            height=3 + 0.12 * min(nrow(de[!is.na(de$adj.P.Val) & de$adj.P.Val<0.05,]), 500),
            units='in', res=400)

        pheatmap(t(scaled.sig), breaks = seq(-range, range, length.out = 101),
                 annotation_col=annot_df, main=paste(groups[1], 'vs.', groups[2]),
                 ## cluster_rows=sum(de$adj.P.Val < 0.05)>1,
                 clustering_method='ward.D', clustering_distance_rows = "euclidean",
                 clustering_distance_cols = "euclidean", annotation_colors=temp.annot.cols,
                 labels_row=col.labels, show_colnames=F)
        dev.off()


        if (min(nrow(de[!is.na(de$adj.P.Val) & de$adj.P.Val<0.05,]), 25) == 25) {
            col.labels <- rownames(de)[1:25]
            col.idxs <- rownames(de)[1:25]
        } else {
            col.labels <- rownames(de)[de$adj.P.Val < 0.05]
            col.idxs <- rownames(de)[de$adj.P.Val < 0.05]
        }

        scaled.sig <- scale(prot.data.combat[prot.data.combat$Condition. %in% make.names(groups), col.idxs])
            
        
        range <- max(abs(scaled.sig))

        print(range)

        range <- min(max(abs(scaled.sig)),3)
        
        pdf(paste('plots-cond3/', str_replace(groups[1], ' ', '-'), '_',
                  str_replace(groups[2], ' ', '.'),
                  '_prot_diffExp_top25_heatmap_unlabeled.pdf',sep=''),
            width=8,
            height=7)

        pheatmap(t(scaled.sig), breaks = seq(-range, range, length.out = 101),
                 annotation_col=annot_df, main=paste(groups[1], 'vs.', groups[2]),
                 ## cluster_rows=sum(de$adj.P.Val < 0.05)>1,
                 clustering_method='ward.D', clustering_distance_rows = "euclidean",
                 clustering_distance_cols = "euclidean", annotation_colors=temp.annot.cols,
                 labels_row=col.labels, show_colnames=F, show_rownames=F)
        
        dev.off()

        png(paste('plots-cond3/', str_replace(groups[1], ' ', '-'), '_',
                  str_replace(groups[2], ' ', '.'),
                  '_prot_diffExp_top25_heatmap_unlabeled.png',sep=''),
            width=8,
            height=7,
            units='in', res=400)

        pheat <- pheatmap(t(scaled.sig), breaks = seq(-range, range, length.out = 101),
                          annotation_col=annot_df, main=paste(groups[1], 'vs.', groups[2]),
                          ## cluster_rows=sum(de$adj.P.Val < 0.05)>1,
                          clustering_method='ward.D', clustering_distance_rows = "euclidean",
                          clustering_distance_cols = "euclidean",
                          annotation_colors=temp.annot.cols,
                          labels_row=col.labels, show_colnames=F, show_rownames=F)

        write.csv(data.frame(Labels=colnames(scaled.sig)[pheat$tree_row[['order']]]),
                  paste('plots-cond3/', str_replace(groups[1], ' ', '-'), '_',
                        str_replace(groups[2], ' ', '.'),
                        '_prot_diffExp_top25_heatmap_row_labels.csv',sep=''),
                  row.names=F)
        
        dev.off()
        
        ## pdf(paste('plots-cond3/', str_replace(groups[1], ' ', '-'), '_', str_replace(groups[2], ' ', '.'),
        ##           '_prot_diffExp_heatmap.pdf',sep=''),
        ##     width=5 + 0.12 * sum(prot.data.combat$Condition. %in% make.names(groups)),
        ##     height=3 + 0.12 * sum(de$adj.P.Val<0.05))
        
        ## scaled.sig <- scale(prot.data.combat[prot.data.combat$Condition. %in% make.names(groups),
        ##                                      rownames(de)[de$adj.P.Val < 0.05]])
        ## range <- min(max(abs(scaled.sig)),3)

        ## annot_df <- data.frame(Diagnosis = factor(diag.map[prot.data.combat[prot.data.combat$Condition.
        ##                                                                     %in% make.names(groups), 'Condition.']], levels=groups),
        ##                        `Infant Sex` = as.factor(prot.data.combat[prot.data.combat$Condition.
        ##                                                               %in% make.names(groups), 'InfSex']),
        ##                        `Gestational Age` = prot.data.combat[prot.data.combat$Condition.
        ##                                                          %in% make.names(groups), 'WksGest'],
        ##                        ## Race = prot.data.combat[prot.data.combat$Condition.
        ##                        ##                         %in% groups, 'Race'],
        ##                        ## PrePregBMI = prot.data.combat[prot.data.combat$Condition.
        ##                        ##                               %in% groups, 'PrePregBMI'],
        ##                        ## Smoking = prot.data.combat[prot.data.combat$Condition.
        ##                        ##                            %in% groups, 'Smoking'],
        ##                        ## DeliveryType = prot.data.combat[prot.data.combat$Condition.
        ##                        ##                            %in% groups, 'FDELTYPE'],
        ##                        row.names=rownames(prot.data.combat[prot.data.combat$Condition.
        ##                                                            %in% make.names(groups),]))

        ## temp.annot.cols <- list()
        ## for (g in groups) {
        ##     temp.annot.cols[['Diagnosis']][g] <- annot.colors[['Diagnosis']][g]
        ## }
        
        ## pheatmap(t(scaled.sig), breaks = seq(-range, range, length.out = 101),
        ##          annotation_col=annot_df, main=paste(groups[1], 'vs.', groups[2]),
        ##          clustering_method='ward.D', clustering_distance_rows = "euclidean",
        ##          clustering_distance_cols = "euclidean", annotation_colors=temp.annot.cols, 
        ##          labels_row=rownames(de)[de$adj.P.Val < 0.05])
        ## dev.off()


        ## png(paste('plots-cond3/', str_replace(groups[1], ' ', '-'), '_', str_replace(groups[2], ' ', '.'),
        ##           '_prot_diffExp_heatmap.png',sep=''),
        ##     width=5 + 0.12 * sum(prot.data.combat$Condition. %in% make.names(groups)),
        ##     height=3 + 0.12 * sum(de$adj.P.Val<0.05), units='in', res=400)
        
        ## pheatmap(t(scaled.sig), breaks = seq(-range, range, length.out = 101),
        ##          annotation_col=annot_df, main=paste(groups[2], 'vs.', groups[1]),
        ##          clustering_method='ward.D', clustering_distance_rows = "euclidean",
        ##          clustering_distance_cols = "euclidean", annotation_colors=temp.annot.cols, 
        ##          labels_row=rownames(de)[de$adj.P.Val < 0.05])
        ## dev.off()
    }

    de.meta <- prot.meta %>% filter(OlinkID %in% prots2keep) %>% filter(Assay %in% rownames(de)) %>% group_by(Assay) %>% mutate_all(function(x) paste(x, collapse=';')) %>% filter(row_number()==1) %>% ungroup() %>% as.data.frame
    rownames(de.meta) <- de.meta$Assay
    ## rownames(de) <- make.names(rownames(de))
    
    de <- cbind(de.meta[rownames(de),], de)
    
    write.csv(de, paste('DE-cond3/diffExpProt_', str_replace(groups[1], ' ', '-'), '_',
                        str_replace(groups[2], ' ', '.'), '.csv', sep=''))
    write.csv(de[de$adj.P.Val<0.05,], paste('DE-cond3/diffExpProt_', str_replace(groups[1], ' ', '-'), '_',
                                            str_replace(groups[2], ' ', '.'), '_FDR05.csv', sep=''))
    write.csv(de[de$adj.P.Val<0.1,], paste('DE-cond3/diffExpProt_', str_replace(groups[1], ' ', '-'), '_',
                                           str_replace(groups[2], ' ', '.'), '_FDR1.csv', sep=''))
}

colnames(prot.data.combat)

prot.data.combat$FLT1overPGF <- prot.data.combat$FLT1 / prot.data.combat$PGF

prot.data.combat$Condition. <- factor(prot.data.combat$Condition.)

FLT1overPGF.tests <- data.frame()

for (groups in comparisons) {
    dataset_group <- prot.data.combat[prot.data.combat$Condition. %in% groups,]
    dataset_group$Condition. <- factor(dataset_group$Condition.)
    FLT1overPGF.tests <- rbind(FLT1overPGF.tests,
                               wilcox_test(dataset_group, FLT1overPGF ~ Condition.))
}

FLT1overPGF.tests <- FLT1overPGF.tests %>% adjust_pvalue(method = 'holm') %>%
    filter(p.adj < 0.05) %>%
    mutate(y.position = c(1.2, 1.3, 1.4, 1.5, 1.6))

library(ggpubr)
ggviolin(prot.data.combat, 'Condition.', 'FLT1overPGF', fill='Condition.', add = "boxplot") +
    stat_pvalue_manual(FLT1overPGF.tests)

ggsave('prot-FLT1overPGF-dists.pdf',
       ggviolin(prot.data.combat, 'Condition.', 'FLT1overPGF',
                fill='Condition.', ylab='FLT1 / PGF', add = "boxplot") +
       stat_pvalue_manual(FLT1overPGF.tests), height=6, width=6)

ggsave('prot-FLT1overPGF-dists.png',
       ggviolin(prot.data.combat, 'Condition.', 'FLT1overPGF',
                fill='Condition.', ylab='FLT1 / PGF', add = "boxplot") +
       stat_pvalue_manual(FLT1overPGF.tests), height=6, width=6)


dataset_group<-prot.data.combat ## [prot.data.combat$Condition. %in% groups,]
print(table(dataset_group$Condition.))
## design a model - control will always be on top as 0
#### VUE all ####
design<-model.matrix(~0 + VUE + WksGest + factor(InfSex) + factor(Race) + PrePregBMI + factor(Smoking) + factor(FDELTYPE) + factor(Labor.initiation), dataset_group)
colnames(design)<-c("VUE0", "VUE1", "GestationalAge", "InfantSex", "RaceOther", "RaceW", "PrePregBMI", "Smoking", "DeliveryType", "LaborInitiation")
## make contrast - what to compare
contrast<- makeContrasts(Diff = "VUE1 - VUE0", levels=design)
print(contrast)
## apply linear model to each protein
## Robust regression provides an alternative to least squares regression that works with less restrictive assumptions. Specifically, it provides much better regression coefficient estimates when outliers are present in the data
fit<-lmFit(t(dataset_group[rownames(design),16:ncol(dataset_group)]), design=design, maxit=1000)
## apply contrast
contrast_fit<-contrasts.fit(fit, contrast)
## apply empirical Bayes smoothing to the SE
ebays_fit<-eBayes(contrast_fit, robust=T)
## summary
print(summary(decideTests(ebays_fit)))
## extract DE results
de<-topTable(ebays_fit, n=ncol(dataset_group), adjust.method="fdr", confint=TRUE)

de <- de[order(de$P.Value),]

de[de$adj.P.Val < 0.05,]

scaled.sig <- scale(prot.data.combat[,rownames(de)[de$adj.P.Val < 0.05]])
range <- min(max(abs(scaled.sig)),3.5)

range

annot_df <- data.frame(VUE = as.factor(prot.data.combat[, 'VUE']),
                       Infant.Sex = as.factor(prot.data.combat[, 'InfSex']),
                       Gestational.Age = prot.data.combat[, 'WksGest'], 
                       row.names=rownames(prot.data.combat[,]))

pdf(paste('VUE_NoVUE_prot_diffExp_heatmap.pdf',sep=''),
    width=5 + 0.12 * nrow(prot.data.combat),
    height=3 + 0.12 * sum(de$adj.P.Val<0.05))

pheatmap(t(scaled.sig), breaks = seq(-range, range, length.out = 101),
         annotation_col=annot_df, main="VUE vs No VUE",
         clustering_method='ward.D', clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean") ## , 

dev.off()

## pheatmap(t(prot.data.combat[,c('HSPB6', 'MIF', 'SCLY', 'WFIKKN2','ICAM-2','PSGL-1','TFF3','WksGest')]), breaks = seq(-3, 3, length.out = 101), scale='row', annotation_col=annot_df, clustering_method='ward.D')

## bidir.nodes <- c()

## for (e in g$edges) {
##     split.edge <- str_split(e, ' ', n=3)[[1]]
##     if (split.edge[2]=='<->') {
##         bidir.nodes <- c(bidir.nodes, split.edge[1], split.edge[3])
##     }
## }

## bidir.nodes <- unique(bidir.nodes)

## bidir.nodes <- bidir.nodes[bidir.nodes %in% colnames(prot.data[,-c(1:9)])]

## length(bidir.nodes)

## pheatmap(t(prot.data.combat[,bidir.nodes]), breaks = seq(-3, 3, length.out = 101), scale='row', annotation_col=annot_df, clustering_method='ward.D')
