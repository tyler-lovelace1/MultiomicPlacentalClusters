library(dplyr)
library(stringr)
library(limma)
library(sva)
library(ggplot2)
library(pheatmap)
library(readxl)
source('utils.R')

## clin <- read.csv('clinical-211205.csv', na.strings='#N/A', row.names=1)

## clin <- read_excel('Primary clinical variables.xlsx', 1) %>% as.data.frame
## colnames(clin) <- make.names(colnames(clin))
## clin <- clin[-1,]
## rownames(clin) <- make.names(clin$Study_ID.ID)
## head(clin)

## clin <- clin %>% mutate_if(function(x) mean(is.na(as.numeric(x))) < 0.5, as.numeric)

## ## rownames(clin) <- str_replace_all(rownames(clin), '-', '.')

## ## head(clin)

## clin$Batch <- ifelse(str_detect(rownames(clin), 'KH'), 'KH', NA)
## clin$Batch <- ifelse(str_detect(rownames(clin), 'MJ'), 'MJ', clin$Batch)
## clin$Batch <- ifelse(str_detect(rownames(clin), 'Mini.DP'), 'Mini-DP', clin$Batch)

## clin$Condition. <- ifelse(clin$Condition. %in% c('PTD','Severe PE','FGR','Control','Control PTD'), clin$Condition., 'FGR+HDP')

## clin$PrePregBMI <- clin$PrePregWt_Kg / clin$Height_Meters.HeightMeters^2

## clin$Smoking <- ifelse(str_detect(clin$Smoke, 'Yes'), 'Yes', NA)
## clin$Smoking <- ifelse(str_detect(clin$Smoke, 'Never'), 'No', clin$Smoking)

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

final.clin <- clin %>% select(Condition., Batch, WksGest, InfSex, Race, PrePregBMI,
                              Smoking, FDELTYPE, Labor.initiation, Labor.no.labor)

final.clin$Race[!final.clin$Race %in% c('W', 'B')] <- 'other'
final.clin$FDELTYPE[final.clin$FDELTYPE==0] <- NA

final.clin %>% mutate_at(c('Race', 'FDELTYPE', 'Labor.initiation', 'Labor.no.labor'), factor)

tail(final.clin)


#### IMPORT AND FILTER METBOLITE DATA ####


metab.meta <- read_excel('UPIT-02-20PHML+ DATA TABLES OB.xlsx', 2) %>% as.data.frame
rownames(metab.meta) <- metab.meta$CHEM_ID
head(metab.meta)

samp.meta <- read_excel('UPIT-02-20PHML+ DATA TABLES OB.xlsx', 3) %>%
    mutate_at('CLIENT_IDENTIFIER', function(x) str_replace_all(x, '-', '.')) %>% as.data.frame
rownames(samp.meta) <- samp.meta$PARENT_SAMPLE_NAME
head(samp.meta[,1:5])

metab <- read_excel('UPIT-02-20PHML+ DATA TABLES OB.xlsx', 8) %>% as.data.frame
## rownames(metab) <- metab$PARENT_SAMPLE_NAME
rownames(metab) <- samp.meta[metab$PARENT_SAMPLE_NAME,'CLIENT_IDENTIFIER']
metab <- metab[,-1]
head(metab[,1:10])

## metab <- metab[,apply(metab, 2, sd)!=0]
metab <- metab %>% select(function(x) length(unique(x)) / length(x) > 0.5)

samp.list <- intersect(rownames(final.clin), rownames(metab))

length(samp.list) == nrow(metab)

cond.totals <- table(factor(final.clin[samp.list,'Condition.']))

cond.totals


pcs <- prcomp(t(scale(metab)))

plot(pcs$sdev[1:50])

summary(pcs)

pc.plot <- pcs$rotation[,1:10]

pc.plot <- cbind(pc.plot, final.clin[rownames(metab),1:9])

ggplot(pc.plot, aes(x=PC1, y=PC2, col=Batch)) +
    geom_point() +
    geom_vline(xintercept=mean(pc.plot$PC1)-4*sd(pc.plot$PC1)) +
    geom_vline(xintercept=mean(pc.plot$PC1)+4*sd(pc.plot$PC1)) +
    geom_hline(yintercept=mean(pc.plot$PC2)-4*sd(pc.plot$PC2)) +
    geom_hline(yintercept=mean(pc.plot$PC2)+4*sd(pc.plot$PC2)) +
    ggtitle('PC1 vs PC2 by Batch') +
    theme_bw()



summary(lm(WksGest ~ 1 + PC1 + PC2 + PC3 + PC4 + PC5 +
               PC6 + PC7 + PC8 + PC9 + PC10, pc.plot))

## Call:
## lm(formula = WksGest ~ 1 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + 
##     PC7 + PC8 + PC9 + PC10, data = pc.plot)

## Residuals:
##      Min       1Q   Median       3Q      Max 
## -12.3903  -1.1703   0.0376   1.3150   7.1121 

## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)  36.8040     0.1138 323.335  < 2e-16 ***
## PC1          -5.6661     2.1019  -2.696  0.00738 ** 
## PC2           5.8857     2.1019   2.800  0.00541 ** 
## PC3         -17.1872     2.1019  -8.177 6.34e-15 ***
## PC4          11.3157     2.1019   5.383 1.39e-07 ***
## PC5         -29.6264     2.1019 -14.095  < 2e-16 ***
## PC6           0.4128     2.1019   0.196  0.84443    
## PC7           2.7376     2.1019   1.302  0.19367    
## PC8          15.8356     2.1019   7.534 4.78e-13 ***
## PC9           4.9126     2.1019   2.337  0.02003 *  
## PC10        -18.5441     2.1019  -8.822  < 2e-16 ***
## ---
## Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## Residual standard error: 2.102 on 330 degrees of freedom
## Multiple R-squared:  0.5777,	Adjusted R-squared:  0.5649 
## F-statistic: 45.14 on 10 and 330 DF,  p-value: < 2.2e-16

summary(glm(Batch=='Mini-DP' ~ 1 + PC1 + PC2 + PC3 + PC4 + PC5 +
                PC6 + PC7 + PC8 + PC9 + PC10, pc.plot, family='binomial'))

## Call:
## glm(formula = Batch == "Mini-DP" ~ 1 + PC1 + PC2 + PC3 + PC4 + 
##     PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family = "binomial", 
##     data = pc.plot)

## Deviance Residuals: 
##     Min       1Q   Median       3Q      Max  
## -2.5022  -0.6703   0.5000   0.7632   1.7739  

## Coefficients:
##             Estimate Std. Error z value Pr(>|z|)    
## (Intercept)   1.2330     0.1479   8.337  < 2e-16 ***
## PC1           6.2868     2.7565   2.281 0.022564 *  
## PC2          -9.3512     2.5441  -3.676 0.000237 ***
## PC3          -3.9076     2.6073  -1.499 0.133947    
## PC4          -2.6501     2.6138  -1.014 0.310620    
## PC5         -12.5444     2.6554  -4.724 2.31e-06 ***
## PC6          -6.8796     2.6763  -2.571 0.010153 *  
## PC7           4.0869     2.7562   1.483 0.138135    
## PC8           4.8851     2.5991   1.880 0.060171 .  
## PC9          -2.3133     2.6700  -0.866 0.386266    
## PC10         -0.8477     2.4940  -0.340 0.733925    
## ---
## Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## (Dispersion parameter for binomial family taken to be 1)

##     Null deviance: 395.64  on 340  degrees of freedom
## Residual deviance: 333.91  on 330  degrees of freedom
## AIC: 355.91

## Number of Fisher Scoring iterations: 4

pdf('uncorrected_metab_PCA_plots.pdf', width=6, height=5)
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
dev.off()

final.clin$Race[!final.clin$Race %in% c('W', 'B')] <- 'other'
final.clin$FDELTYPE[final.clin$FDELTYPE==0] <- NA

mod <- model.matrix(~Condition. + WksGest + InfSex + Race + PrePregBMI + Smoking + factor(FDELTYPE) + factor(Labor.initiation), final.clin[rownames(metab),])

## prot.data.corrected <- removeBatchEffect(t(prot.data[,-c(1:9)]), batch=prot.data$Batch, design=mod)

metab.corrected <- ComBat(t(metab[rownames(mod),]), batch=final.clin[rownames(mod),'Batch'], mod=mod)

pcs <- prcomp(t(scale(t(metab.corrected))))

plot(pcs$sdev[1:50])

summary(pcs)

pc.plot <- pcs$rotation[,1:10]

pc.plot <- cbind(pc.plot, final.clin[colnames(metab.corrected),1:9])



summary(lm(WksGest ~ 1 + PC1 + PC2 + PC3 + PC4 + PC5 +
               PC6 + PC7 + PC8 + PC9 + PC10 + Condition. +
               InfSex + Race + PrePregBMI + Smoking + factor(FDELTYPE),
           pc.plot))

## Call:
## lm(formula = WksGest ~ 1 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + 
##     PC7 + PC8 + PC9 + PC10 + Condition. + InfSex + Race + PrePregBMI + 
##     Smoking + factor(FDELTYPE), data = pc.plot)

## Residuals:
##     Min      1Q  Median      3Q     Max 
## -9.5463 -0.9244  0.1443  1.1361  4.0838 

## Coefficients:
##                         Estimate Std. Error t value Pr(>|t|)    
## (Intercept)            3.802e+01  5.541e-01  68.613  < 2e-16 ***
## PC1                   -4.696e+00  1.738e+00  -2.701 0.007307 ** 
## PC2                   -8.282e+00  1.805e+00  -4.590 6.54e-06 ***
## PC3                    1.055e+01  1.843e+00   5.727 2.49e-08 ***
## PC4                   -6.378e+00  1.769e+00  -3.605 0.000366 ***
## PC5                    1.798e+01  2.049e+00   8.773  < 2e-16 ***
## PC6                    5.809e+00  1.788e+00   3.249 0.001289 ** 
## PC7                   -4.091e+00  1.718e+00  -2.381 0.017886 *  
## PC8                   -1.519e+01  1.861e+00  -8.161 9.36e-15 ***
## PC9                   -1.117e+01  1.727e+00  -6.468 4.06e-10 ***
## PC10                   8.463e+00  1.780e+00   4.755 3.09e-06 ***
## Condition.Control PTD -1.988e+00  5.041e-01  -3.944 0.000100 ***
## Condition.FGR         -1.052e+00  3.500e-01  -3.007 0.002864 ** 
## Condition.FGR+HDP     -3.421e+00  4.648e-01  -7.361 1.78e-12 ***
## Condition.PTD         -3.608e+00  2.859e-01 -12.621  < 2e-16 ***
## Condition.Severe PE   -2.487e+00  3.160e-01  -7.870 6.56e-14 ***
## InfSexM               -4.526e-02  1.967e-01  -0.230 0.818210    
## Raceother              6.594e-01  3.982e-01   1.656 0.098801 .  
## RaceW                  2.654e-01  2.664e-01   0.996 0.319930    
## PrePregBMI             1.168e-02  1.566e-02   0.746 0.456255    
## SmokingYes             1.711e-01  2.094e-01   0.817 0.414452    
## factor(FDELTYPE)2     -2.407e-04  2.170e-01  -0.001 0.999116    
## ---
## Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## Residual standard error: 1.677 on 299 degrees of freedom
## Multiple R-squared:  0.7449,	Adjusted R-squared:  0.727 
## F-statistic: 41.59 on 21 and 299 DF,  p-value: < 2.2e-16


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
## -2.6388  -0.3600   0.4157   0.6464   1.9877  

## Coefficients:
##                        Estimate Std. Error z value Pr(>|z|)   
## (Intercept)           -2.793316   3.605627  -0.775  0.43851   
## PC1                    1.788042   2.964957   0.603  0.54647   
## PC2                    1.175463   3.161441   0.372  0.71003   
## PC3                   -0.159594   3.203111  -0.050  0.96026   
## PC4                    1.120181   2.903114   0.386  0.69960   
## PC5                    0.655921   3.711513   0.177  0.85972   
## PC6                    2.474072   2.977737   0.831  0.40605   
## PC7                   -1.374694   2.855758  -0.481  0.63025   
## PC8                   -0.837607   3.468667  -0.241  0.80918   
## PC9                   -0.247256   2.911226  -0.085  0.93232   
## PC10                  -0.861902   3.025796  -0.285  0.77576   
## Condition.Control PTD -0.998875   0.710201  -1.406  0.15959   
## Condition.FGR         -1.479747   0.540447  -2.738  0.00618 **
## Condition.FGR+HDP     -2.508943   0.807700  -3.106  0.00189 **
## Condition.PTD         -0.027715   0.615131  -0.045  0.96406   
## Condition.Severe PE   -0.115644   0.588502  -0.197  0.84421   
## WksGest                0.112898   0.092364   1.222  0.22159   
## InfSexM               -0.001107   0.329272  -0.003  0.99732   
## Raceother              1.663680   0.721551   2.306  0.02113 * 
## RaceW                  1.039790   0.410690   2.532  0.01135 * 
## PrePregBMI             0.004075   0.025789   0.158  0.87443   
## SmokingYes            -0.430745   0.325598  -1.323  0.18586   
## factor(FDELTYPE)2     -1.104935   0.341342  -3.237  0.00121 **
## ---
## Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## (Dispersion parameter for binomial family taken to be 1)

##     Null deviance: 369.03  on 320  degrees of freedom
## Residual deviance: 273.97  on 298  degrees of freedom
## AIC: 319.97

## Number of Fisher Scoring iterations: 5

## pc.plot <- pc.plot %>% mutate_at('Acetaminophen.y.n', as.factor)

pdf('corrected_metab_PCA_plots.pdf', width=6, height=5)
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

## ggplot(pc.plot, aes(x=PC1, y=PC2, col=Acetaminophen.y.n)) +
##     geom_point() +
##     ggtitle('PC1 vs PC2 by Infant Sex') +
##     theme_bw()

dev.off()

groups <- c('Control PTD', 'Control')

metab.data.combat <- cbind(final.clin[colnames(metab.corrected),],t(metab.corrected))

metab.only.combat <- metab.data.combat[,10:ncol(metab.data.combat)]
colnames(metab.only.combat) <- metab.meta[colnames(metab.only.combat),'CHEMICAL_NAME']

write.csv(metab.only.combat, 'data/combat_metab_data.csv')


## metab.data.combat <- read.csv('data/combat_metab_data.csv', row.names=1)

## metab.data.combat <- cbind(final.clin[rownames(metab.data.combat),], metab.data.combat)



run_limma<-function(groups) {
    ## extract data by group
    dataset_group<-metab.data.combat ## [prot.data.combat$Condition. %in% groups,]
    print(table(dataset_group$Condition.))
    ## design a model - control will always be on top as 0

    #### COND GROUP 1 ####
    ## design<-model.matrix(~0 + as.factor(dataset_group$Condition.) + dataset_group$WksGest + as.factor(dataset_group$InfSex))
    ## colnames(design)<-c("Control", "Control.PTD", "FGR", "FGR.HDP", "PTD", "Severe.PE", "GestationalAge", "InfantSex")

    ## #### COND GROUP 2 ####
    ## design<-model.matrix(~0 + as.factor(dataset_group$Condition.) + dataset_group$WksGest + as.factor(dataset_group$InfSex) + as.factor(dataset_group$Race) + dataset_group$PrePregBMI + as.factor(dataset_group$Smoking) + as.factor(dataset_group$FDELTYPE) + as.factor(dataset_group$Labor.initiation))
    ## colnames(design)<-c("Control", "Control.PTD", "FGR", "FGR.HDP", "PTD", "Severe.PE",
    ##                     "GestationalAge", "InfantSex", "RaceOther", "RaceW", "PrePregBMI",
    ##                     "Smoking", "DeliveryType", "LaborInitiation")

    #### COND GROUP 3 ####
    design<-model.matrix(~0 + as.factor(dataset_group$Condition.) + dataset_group$WksGest + as.factor(dataset_group$InfSex) + as.factor(dataset_group$Race) + dataset_group$PrePregBMI + as.factor(dataset_group$Smoking) + as.factor(dataset_group$FDELTYPE) + as.factor(dataset_group$Labor.initiation) + as.factor(dataset_group$Labor.no.labor))
    colnames(design)<-c("Control", "Control.PTD", "FGR", "FGR.HDP", "PTD", "Severe.PE",
                        "GestationalAge", "InfantSex", "RaceOther", "RaceW", "PrePregBMI",
                        "Smoking", "DeliveryType", "LaborInitiation", "Labor")
    
    
    ## make contrast - what to compare
    contrast<- makeContrasts(Diff = paste(str_replace(groups, '[+ ]', '.'), collapse=" - "), levels=design)
    print(contrast)
    ## apply linear model to each protein
    ## Robust regression provides an alternative to least squares regression that works with less restrictive assumptions. Specifically, it provides much better regression coefficient estimates when outliers are present in the data
    fit<-lmFit(t(dataset_group[11:ncol(dataset_group)]), design=design, maxit=1000)
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

## pdf('eBayes_robust_metab_diffExp_heatmaps.pdf', width=12, height=12)
for (groups in comparisons) {
    
    print(groups)
    de <- run_limma(groups)
    print(paste('DE-cond3/diffExpMetab_', str_replace(groups[1], ' ', '-'), '_',
                str_replace(groups[2], ' ', '.'), '.csv', sep=''))

    if (sum(de$adj.P.Val < 0.05)>1) {

        if (min(nrow(de[!is.na(de$adj.P.Val) & de$adj.P.Val<0.05,]), 500) == 500) {
            col.labels <- metab.meta[rownames(de)[1:500],'CHEMICAL_NAME']
            col.idxs <- rownames(de)[1:500]
        } else {
            col.labels <- metab.meta[rownames(de)[de$adj.P.Val < 0.05],'CHEMICAL_NAME']
            col.idxs <- rownames(de)[de$adj.P.Val < 0.05]
        }

        scaled.sig <- scale(metab.data.combat[metab.data.combat$Condition. %in% make.names(groups), col.idxs])    
        
        range <- max(abs(scaled.sig))

        print(range)

        range <- min(max(abs(scaled.sig)),3)
        
        annot_df <- data.frame(Diagnosis = factor(diag.map[metab.data.combat[metab.data.combat$Condition.
                                                                             %in% make.names(groups), 'Condition.']], levels=groups),
                               `Infant Sex` = as.factor(metab.data.combat[metab.data.combat$Condition.
                                                                        %in% make.names(groups), 'InfSex']),
                               `Gestational Age` = metab.data.combat[metab.data.combat$Condition.
                                                                     %in% make.names(groups), 'WksGest'], 
                               row.names=rownames(metab.data.combat[metab.data.combat$Condition.
                                                                    %in% make.names(groups),]))

        temp.annot.cols <- list()
        for (g in groups) {
            temp.annot.cols[['Diagnosis']][g] <- annot.colors[['Diagnosis']][g]
        }
       
        pdf(paste('plots-cond3/', str_replace(groups[1], ' ', '-'), '_',
                  str_replace(groups[2], ' ', '.'),
                  '_metab_diffExp_heatmap.pdf',sep=''),
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
                  '_metab_diffExp_heatmap.png',sep=''),
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
            col.labels <- metab.meta[rownames(de)[1:25],'CHEMICAL_NAME']
            col.idxs <- rownames(de)[1:25]
        } else {
            col.labels <- metab.meta[rownames(de)[de$adj.P.Val < 0.05],'CHEMICAL_NAME']
            col.idxs <- rownames(de)[de$adj.P.Val < 0.05]
        }

        scaled.sig <- scale(metab.data.combat[metab.data.combat$Condition. %in% make.names(groups), col.idxs])
            
        
        range <- max(abs(scaled.sig))

        print(range)

        range <- min(max(abs(scaled.sig)),3)
        
        pdf(paste('plots-cond3/', str_replace(groups[1], ' ', '-'), '_',
                  str_replace(groups[2], ' ', '.'),
                  '_metab_diffExp_top25_heatmap_unlabeled.pdf',sep=''),
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
                  '_metab_diffExp_top25_heatmap_unlabeled.png',sep=''),
            width=8,
            height=7,
            units='in', res=400)

        pheat <- pheatmap(t(scaled.sig), breaks = seq(-range, range, length.out = 101),
                 annotation_col=annot_df, main=paste(groups[1], 'vs.', groups[2]),
                 ## cluster_rows=sum(de$adj.P.Val < 0.05)>1,
                 clustering_method='ward.D', clustering_distance_rows = "euclidean",
                 clustering_distance_cols = "euclidean", annotation_colors=temp.annot.cols,
                 labels_row=col.labels, show_colnames=F, show_rownames=F)

        write.csv(data.frame(Labels=col.labels[pheat$tree_row[['order']]]),
                  paste('plots-cond3/', str_replace(groups[1], ' ', '-'), '_',
                        str_replace(groups[2], ' ', '.'),
                        '_metab_diffExp_top25_heatmap_row_labels.csv',sep=''),
                  row.names=F)
        
        dev.off()
        
        ## pdf(paste('plots-cond3/', str_replace(groups[1], ' ', '-'), '_', str_replace(groups[2], ' ', '.'),
        ##           '_metab_diffExp_heatmap.pdf',sep=''),
        ##     width=5 + 0.12 * sum(metab.data.combat$Condition. %in% make.names(groups)),
        ##     height=3 + 0.12 * sum(de$adj.P.Val<0.05))
        
        ## scaled.sig <- scale(metab.data.combat[metab.data.combat$Condition. %in% make.names(groups),
        ##                                       rownames(de)[de$adj.P.Val < 0.05]])
        ## range <- min(max(abs(scaled.sig)),3)

        ## annot_df <- data.frame(Diagnosis = factor(diag.map[metab.data.combat[metab.data.combat$Condition.
        ##                                                                      %in% make.names(groups), 'Condition.']], levels=groups),
        ##                        `Infant Sex` = as.factor(metab.data.combat[metab.data.combat$Condition.
        ##                                                                 %in% make.names(groups), 'InfSex']),
        ##                        `Gestational Age` = metab.data.combat[metab.data.combat$Condition.
        ##                                                            %in% make.names(groups), 'WksGest'], 
        ##                        row.names=rownames(metab.data.combat[metab.data.combat$Condition.
        ##                                                             %in% make.names(groups),]))
        ## temp.annot.cols <- list()
        ## for (g in groups) {
        ##     temp.annot.cols[['Diagnosis']][g] <- annot.colors[['Diagnosis']][g]
        ## }
       
        ## pheatmap(t(scaled.sig), breaks = seq(-range, range, length.out = 101),
        ##          annotation_col=annot_df, main=paste(groups[1], 'vs.', groups[2]),
        ##          cluster_rows=sum(de$adj.P.Val < 0.05)>1,
        ##          clustering_method='ward.D', clustering_distance_rows = "euclidean",
        ##          clustering_distance_cols = "euclidean", annotation_colors=temp.annot.cols,
        ##          labels_row=metab.meta[rownames(de)[de$adj.P.Val < 0.05],'CHEMICAL_NAME'])
        ## dev.off()

        ## png(paste('plots-cond3/', str_replace(groups[1], ' ', '-'), '_', str_replace(groups[2], ' ', '.'),
        ##           '_metab_diffExp_heatmap.png',sep=''),
        ##     width=5 + 0.12 * sum(metab.data.combat$Condition. %in% make.names(groups)),
        ##     height=3 + 0.12 * sum(de$adj.P.Val<0.05), units='in', res=400)
        
        ## pheatmap(t(scaled.sig), breaks = seq(-range, range, length.out = 101),
        ##          annotation_col=annot_df, main=paste(groups[1], 'vs.', groups[2]),
        ##          cluster_rows=sum(de$adj.P.Val < 0.05)>1,
        ##          clustering_method='ward.D', clustering_distance_rows = "euclidean",
        ##          clustering_distance_cols = "euclidean", annotation_colors=temp.annot.cols,
        ##          labels_row=metab.meta[rownames(de)[de$adj.P.Val < 0.05],'CHEMICAL_NAME'])
        ## dev.off()
    }

    de <- cbind(metab.meta[rownames(de),], de)
    
    write.csv(de, paste('DE-cond3/diffExpMetab_', str_replace(groups[1], ' ', '-'), '_',
                        str_replace(groups[2], ' ', '.'), '.csv', sep=''))
    write.csv(de[de$adj.P.Val<0.05,], paste('DE-cond3/diffExpMetab_', str_replace(groups[1], ' ', '-'), '_',
                                            str_replace(groups[2], ' ', '.'), '_FDR05.csv', sep=''))
    write.csv(de[de$adj.P.Val<0.1,], paste('DE-cond3/diffExpMetab_', str_replace(groups[1], ' ', '-'), '_',
                        str_replace(groups[2], ' ', '.'), '_FDR1.csv', sep=''))
}
## dev.off()


smoker.metabs <- metab[,as.character(metab.meta[grep('cotinine', metab.meta$CHEMICAL_NAME),'CHEM_ID'])]

colnames(smoker.metabs) <- metab.meta[grep('cotinine', metab.meta$CHEMICAL_NAME),'CHEMICAL_NAME']

clin <- read.csv('doppler_clinical.csv', row.names=1)
colnames(clin) <- make.names(colnames(clin))
rownames(clin) <- make.names(clin$Study_ID.ID)
head(clin)

clin <- clin %>% mutate_if(function(x) mean(is.na(as.numeric(x))) < 0.5, as.numeric)

clin$Condition. <- ifelse(clin$Condition. %in% c('PTD','Severe PE','FGR','Control','Control PTD'), clin$Condition., 'FGR+HDP')

clin$PrePregBMI <- clin$PrePregWt_Kg / clin$Height_Meters.HeightMeters^2

clin$Smoking <- ifelse(str_detect(clin$Smoke, 'Yes'), 'Yes', NA)
clin$Smoking <- ifelse(str_detect(clin$Smoke, 'Never'), 'No', clin$Smoking)

clin$Race[!clin$Race %in% c('W', 'B')] <- 'other'
clin$FDELTYPE[clin$FDELTYPE==0] <- NA

clin <- clin[!is.na(clin$hvar.Clusters),]

samp.ids <- rownames(clin)

smoke.ids <- intersect(samp.ids, rownames(smoker.metabs))

library(reshape2)

plotdf <- melt(cbind(clin[smoke.ids,c('Smoke', 'Smoking','Condition.','hvar.Clusters')], smoker.metabs[smoke.ids,]), id.vars=c('Smoke', 'Smoking','Condition.','hvar.Clusters'), value.name='Expression', variable.name='Metabolite')

plotdf$Clusters <- factor(plotdf$hvar.Clusters)

head(plotdf)

ggplot(plotdf, aes(x=Metabolite, y=Expression, fill=Smoking)) +
    geom_boxplot() +
    stat_compare_means() +
    theme_classic()

ggplot(plotdf, aes(x=Metabolite, y=Expression, fill=Clusters)) +
    geom_boxplot() +
    stat_compare_means() +
    theme_classic()

ggplot(plotdf, aes(x=Metabolite, y=Expression, fill=Condition.)) +
    geom_boxplot() +
    stat_compare_means() +
    theme_classic()


plotdf <- melt(cbind(clin[samp.ids,c('Smoke', 'Smoking','Condition.','hvar.Clusters')], smoker.metabs[samp.ids,]), id.vars=c('Smoke', 'Smoking','Condition.','hvar.Clusters'), value.name='Expression', variable.name='Metabolite')

plotdf$Clusters <- factor(plotdf$hvar.Clusters)

head(plotdf)

library(ggpubr)
library(RColorBrewer)

brewer.colors <- brewer.pal(n=12, "Paired")## )(10)

brewer.colors <- c(brewer.colors[1:2], brewer.colors[11], brewer.colors[5:6], brewer.colors[8], brewer.colors[3:4], brewer.colors[9:10])

pdf('smoking_metabolites_plots.pdf', width=6, height=6)
ggplot(plotdf, aes(x=Metabolite, y=Expression, fill=Smoking)) +
    geom_violin() +
    geom_boxplot(width=0.1, position=position_dodge(width=0.9)) +
    stat_compare_means() +
    theme_classic()

ggplot(plotdf, aes(x=Metabolite, y=Expression, fill=Clusters)) +
    geom_violin() +
    geom_boxplot(width=0.1, position=position_dodge(width=0.9)) +
    stat_compare_means() +
    scale_fill_manual(values=brewer.colors[7:10]) +
    theme_classic()

ggplot(plotdf, aes(x=Metabolite, y=Expression, fill=Condition.)) +
    geom_violin() +
    geom_boxplot(width=0.1, position=position_dodge(width=0.9)) +
    stat_compare_means() +
    scale_fill_manual(values=brewer.colors[1:6]) +
    theme_classic()
dev.off()
