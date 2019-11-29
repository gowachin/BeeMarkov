### R code from vignette source 'Markov_Report_Guerra_Jaunatre.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: packages
###################################################
library(xtable)
library(ggplot2)
library(reshape)
# library(plyr)
# library(MASS)
# library(ResourceSelection)
# library(ggpubr)
# library(tidyr)

if("BeeMarkov" %in% installed.packages()) {
  library(BeeMarkov)
}else{
  # library(devtools)
  # install_github("gowachin/BeeMarkov")
  # library(BeeMarkov)
}

library(seqinr)


###################################################
### code chunk number 2: working (eval = FALSE)
###################################################
## setwd("manuscrit")
## rerun = TRUE


###################################################
### code chunk number 3: transition
###################################################
transition <- function(file, # fichier
                       l_word = 1, # longueur des mots
                       n_seq = 1, # Nombre de sequences a analyser
                       log = TRUE,
                       type = "") {
  seq <-seqinr::read.fasta(file)
  Nseq <- length(seq)
  
  if (Nseq < n_seq) { # check if user want to input too many seq in the matrix learning
    stop(paste("n_seq is larger than the number of sequence in this file ( ", Nseq, " sequences for this file).", sep = ""))
  }
  
  l <- sapply(seq, length)
  if (l_word > min(l) | l_word >= 10) {
    stop("This is really too much for me, abort!!!")
  }

  tmp <-seqinr::count(seq[[1]], l_word) + 1 # add 1 occurence to have at least 1 obs
  
  cat(paste("      ============ Training model M",type,l_word," ============        \n", sep = ""))
  pb <- utils::txtProgressBar(min = 0, max = n_seq, style = 3)
  for (i in 2:n_seq) {
    utils::setTxtProgressBar(pb, i)
    
    tmp <- tmp +seqinr::count(seq[[i]], l_word)
  }
  close(pb)
  
  # # alternativ way, slower but prettier
  # cat(paste("      ============ Training model M",type,l_word," ============        \n", sep = ""))
  # l_count = function(Seq,n = l_word){count(Seq,n)} 
  # tmp <- rowSums(sapply(seq,l_count))
  
  if (l_word>1) {
    i <- 1
    wind = 4 
    for(j in 0:((length(tmp)/wind)-1)){
      tmp[(1+wind*j):(wind+wind*j)] <- tmp[(1+wind*j):(wind+wind*j)] * 4^(i-1) / sum(tmp[(1+wind*j):(wind+wind*j)] )
    }
    cat('      ============ Computing conditionnal probabilities ============        \n')
  } else {
    tmp = tmp / sum(tmp)
  }
  # possibility to compute without log
  if (log) {
    tmp = log(tmp)
  }
  return(tmp)
}

mP <-transition(file = "../raw_data/mus_cpg_app.fa", n_seq = 1160, l_word = 2, log = F)

mP <- matrix(mP, byrow = TRUE, ncol = 4, dimnames = list(unique(substr(names(mP),1,1)),unique(substr(names(mP),1,1))))

matable <- xtable(x = mP , label = "tab:trans")
# Notez les doubles \\ nécessaires dans R, c'est la "double escape rule"
print(matable, file = "fig/tab_trans.tex",size = "tiny", NA.string = "NA",
      table.placement = "!t",
      floating = FALSE,
      caption.placement="top",
      include.rownames = TRUE,
      include.colnames = TRUE,
      latex.environments = "center")


###################################################
### code chunk number 4: threshold (eval = FALSE)
###################################################
## quality <- function(file, # fichier
##                     pos_training = NULL, # file to train with for positive
##                     neg_training = NULL, # file to train with for negative
##                     trans_pos = NULL, #transition matrice pos
##                     trans_neg = NULL, #transition matrice pos
##                     l_word_pos = 1, # word lenght for transition table positif
##                     l_word_neg = 1, # word length for transition table negatif
##                     n_train = 1, # number of sequences to train with
##                     n_seq = 1, # number of sequences to analyse
##                     quiet = FALSE
## ){
##   if(is.null(c(trans_pos,trans_neg))){
##     trans_pos <-transition(file = pos_training, n_seq = n_train, l_word = l_word_pos, type ="+")
##     trans_neg <- transition(file = neg_training, n_seq = n_train, l_word = l_word_neg, type = "-")
##   } 
##   
##   seq <-seqinr::read.fasta(file)
##   Nseq <- length(seq)
##   
##   if(Nseq < n_seq){ # check if user want to input too many seq in the matrix learning
##     stop(paste("n_seq is larger than the number of sequence in this file ( ",Nseq," sequences for this file).", sep = ""))
##   }
##   
##   result <- data.frame(VP = rep(FALSE,n_seq),
##                        FN = rep(TRUE,n_seq),
##                        pos = rep(NA,n_seq),
##                        neg = rep(NA,n_seq)
##   )
##   
##   p_init_pos = log(1/4^l_word_pos)
##   p_init_neg = log(1/4^l_word_neg)
##   
##   if(!quiet) {cat('      ============ Computing sensi and speci for the test sequences ============       \n')
##     pb <- utils::txtProgressBar(min = 0, max = n_seq, style = 3)}
##   for(i in 1:n_seq){
##     if(!quiet) utils::setTxtProgressBar(pb, i)
##     n_word_pos <-seqinr::count(seq[[i]], l_word_pos)
##     n_word_neg <-seqinr::count(seq[[i]], l_word_neg)
##     result$pos[i] <- p_init_pos + sum( trans_pos * n_word_pos )
##     result$neg[i] <- p_init_neg + sum( trans_neg * n_word_neg )
##     if(result$pos[i] > result$neg[i]){
##       result$VP[i] <- TRUE ; result$FN[i] <- FALSE 
##     } else { 
##       result$VP[i] <- FALSE ; result$FN[i] <- TRUE
##     }
##   }
##   if(!quiet) close(pb)
##   
##   tmp <- colSums(result[,1:2])
##   return(tmp)
## }
## 
## threshold <- function(pos_test, # fichier
##                      neg_test,
##                      pos_training, # file to train with for positive
##                      neg_training, # file to train with for negative
##                      pos_seq = c(1:2),
##                      neg_seq = c(1:2),
##                      n_train = 1, # number of sequences to train with
##                      n_seq = 1 # number of sequences to analyse
## ){
##   
##   # choix =utils::menu(c("yes","no"),
##   #              title = "You will launch long computation, do you wish to procede further ?")
##   # if (choix ==1) cat("\n      ============ Go take a good coffee ============       \n\n")
##   # if (choix ==2) stop("You stopped the computations")
##   
##   
##   trans_pos <- list()
##   trans_neg <- list()
##   
##   cat('\n      ============ Training modeles ============       \n\n')
##   for(i in pos_seq){
##     trans_pos[[i]] <- transition(file = pos_training, n_seq = n_train, l_word = i, type = "+")
##   }
##   cat('\n')
##   for(j in neg_seq){
##     trans_neg[[j]] <- transition(file = neg_training, n_seq = n_train, l_word = j, type = "-")
##   }
##   
##   cat('\n      ============ Training modeles ============       \n\n')    
##   sensi <- speci <- matrix(rep(0,length(pos_seq)*length(neg_seq)),ncol = length(pos_seq),nrow = length(neg_seq))
##   for(i in pos_seq){
##     for(j in neg_seq){
##       cat('\n')
##       cat(paste("      ============ Model M+ (",i,"/",j,") ============        \n", sep = ""))
##       pos <- quality(file = pos_test, # fichier
##                      trans_pos = trans_pos[[i]], #transition matrice pos
##                      trans_neg = trans_neg[[j]], #transition matrice neg
##                      l_word_pos = i, # transition table positif
##                      l_word_neg = j, # transition table negatif
##                      n_train = n_train, # number of sequences to train with
##                      n_seq = n_seq, # number of sequences to analyse
##                      quiet = TRUE
##       )
##       cat(paste("      ============ Model M- (",i,"/",j,") ============        \n", sep = ""))
##       neg <- quality(file = neg_test, # fichier
##                      trans_pos = trans_pos[[i]], #transition matrice pos
##                      trans_neg = trans_neg[[j]], #transition matrice neg
##                      l_word_pos = i, # transition table positif
##                      l_word_neg = j, # transition table negatif
##                      n_train = n_train, # number of sequences to train with
##                      n_seq = n_seq, # number of sequences to analyse
##                      quiet = TRUE
##       )
##       
##       sensi[i,j] <- pos[1] / sum(pos)
##       speci[i,j] <- neg[2] / sum(neg)
##       
##     }
##   }
##   
##   cat(paste("      ============ Come back from your coffee ============        \n", sep = ""))
##   colnames(sensi) <- colnames(speci) <- paste("-",neg_seq, sep="")
##   rownames(sensi) <- rownames(speci) <- paste("+",pos_seq, sep="")
##   
##   final <-list(sensi,speci) ; names(final) = c("sensi","speci")
##   
##   return(final)
## }
## 


###################################################
### code chunk number 5: seuil (eval = FALSE)
###################################################
## if("final.RData" %in% list.files("fig/") && !rerun){
##   load("fig/final.RData")
## }else{
##   final <- threshold("../raw_data/mus_cpg_test.fa", # fichier
##   "../raw_data/mus_tem_test.fa",
##   "../raw_data/mus_cpg_app.fa", # file to train with for positive
##   "../raw_data/mus_tem_app.fa", # file to train with for negative
##   pos_seq = c(1:2),
##   neg_seq = c(1:2),
##   n_train = 1160, # number of sequences to train with
##   n_seq = 1163 # number of sequences to analyse
## )
## save(final,file = "fig/final.RData")
## }
## 
## table <- cbind(melt(final$sensi), melt(final$speci)[, 3])
## colnames(table) <- c("M", "m", "Sensi", "Speci")
## table$score <- table$Sensi + table$Speci
## 
## visual <- ggplot(table, aes(M, m)) +
##   geom_raster(aes(fill = score), hjust = 0.5, vjust = 0.5, interpolate = FALSE) +
##   geom_contour(aes(z = score)) + xlab("CpG+") + ylab("CpG-")
## visual
## ggsave('visual.pdf', plot = visual, device = "pdf", path = 'fig/',scale = 3, width = 7, height = 4, units = "cm", limitsize = T)
## rm(visual)
## 
## table[which(table$tot == max(table$tot)),]


###################################################
### code chunk number 6: trans_mod
###################################################
  l_c <- 1 / 1000
  l_nc <- 1 / 125000
  
   trans_mod <- log(matrix(c(
    1 - l_c, l_nc,
    l_c, 1 - l_nc
  ),
  ncol = 2, nrow = 2 , dimnames = list(c("M+","M-"),c("M+","M-"))
  ))

matable <- xtable(x = trans_mod , label = "tab:trans_mod")
# Notez les doubles \\ nécessaires dans R, c'est la "double escape rule"
print(matable, file = "fig/trans_mod.tex",size = "tiny", NA.string = "NA",
      table.placement = "!t",
      floating = FALSE,
      caption.placement="top",
      include.rownames = TRUE,
      include.colnames = TRUE,
      latex.environments = "center")


###################################################
### code chunk number 7: github_install (eval = FALSE)
###################################################
## # NOT RUN
## library(devtools)
## install_github("gowachin/BeeMarkov")
## library(BeeMarkov)


###################################################
### code chunk number 8: data_visual (eval = FALSE)
###################################################
## icu <- read.delim("../raw_data/icu.txt", header = TRUE)
## 
## # transf en facteurs pour mieux analyser
## icu$LOC[which(icu$LOC == 2)] <- 1
## 
## # creating alternative df with factors for binomial variables
## othr <- data.frame(sapply(icu[, -c(1, 3, 11, 12)], as.factor))
## icuf <- cbind(icu[, 1], othr, icu[, c(3, 11, 12)])
## rm(othr)
## 
## s <- summary(icuf[, -c(1, 19:21)])
## s <- as.data.frame(s)
## s$Freq <- as.character(s$Freq)
## s <- separate(s, # tibble, dataframe
##   Freq, # column separated
##   sep = ":", # separator between columns
##   into = c("rm", "count"), # names of new variables
##   remove = TRUE
## )
## s <- s[, -c(1, 3)]
## s$count <- as.numeric(s$count)
## s$count <- replace_na(s$count, 0)
## s$l <- 1:dim(s)[1]
## s <- ddply(s, "Var2", transform, label_count = cumsum(count))
## s <- s[order(s$l), ]
## s <- s[, -3]
## s$label <- as.character(s$count)
## s$label[which(s$label == "0")] <- ""
## 
## visual = ggplot(s, aes(x = Var2,y = count))+
##   geom_bar(stat = "identity", aes(fill = count), colour="black")+ 
##   geom_text(aes(y=label_count, label=label), vjust=1.6, 
##             color="white", size=3.5) +
##   scale_fill_gradient2(low = "#999999", high = "#E69F00",
##                        mid = "#56B4E9", midpoint = 100 , name = "Effectif") +
##   labs(x = "Variables Qualitatives", y = "") +
##   theme_classic() +
##   theme(axis.text.x = element_text(angle = 90)) 
## 
## ggsave('visual.pdf', plot = visual, device = "pdf", path = 'fig/',scale = 3, width = 7, height = 4, units = "cm", limitsize = T)
## # sum(is.na(icu)) # no missing data
## rm(s,visual)


###################################################
### code chunk number 9: quali_cor (eval = FALSE)
###################################################
## sampleData <- function(data, percent.calib)
## {
##   nCalib = nrow(data) * percent.calib/100
##   iCalib = sample(1:nrow(data), nCalib)
##   calib = data[iCalib,] 
##  # test = data[-iCalib,] 
##   return(
##     #list(
##     calib=calib
##     #, test=test)
##     )
## }
## 
## icuvs_1 = sampleData(icu, 70)
## 
## qualits_all = colnames(icu[,-c(1,2,3,5,11,12)])
## 
## chtest <- function(datas, namesvector) {
##   m1 <- matrix(
##     data = NA, nrow = length(namesvector), ncol = length(namesvector), byrow = FALSE,
##     dimnames = list(namesvector, namesvector)
##   )
##   m2 <- matrix(
##     data = NA, nrow = length(namesvector), ncol = length(namesvector), byrow = FALSE,
##     dimnames = list(namesvector, namesvector)
##   )
## 
##   for (i in 1:length(namesvector)) {
##     var1 <- as.matrix(datas[namesvector[i]])
## 
##     for (j in which(namesvector != namesvector[i])) {
##       var2 <- as.matrix(datas[namesvector[j]])
##       tutu <- chisq.test(table(x = var1, y = var2), simulate.p.value = TRUE)
##       # cat(namesvector[i], "VS", namesvector[j], "\n");
##       # print(tutu)
##       # cat(rep("-", getOption("width"))); cat("\n")
##       if (tutu$p.value < 0.05) m1[i, j] <- round(tutu$p.value, 4)
##       m2[i, j] <- tutu$statistic
##     }
##   }
## 
##   m1[lower.tri(m1)] <- "-"
##   m1 <- as.data.frame(m1)
##   m2[lower.tri(m2)] <- "-"
##   m2 <- as.data.frame(m2)
##   return(list("pvalues" = m1, "statistics" = m2))
## }
## 
## chtest_all = chtest(icu, qualits_all)
## 
## #chall_vs_1 = chtest(icuvs_1$calib, qualits_all)
## 
## ##visualisation of all stuff with or without sampling
## #tibble::as.tibble(cbind.data.frame("names" = rownames(chtest_all[["pvalues"]]),chtest_all[["pvalues"]]))
## correl = tibble::as.tibble(cbind.data.frame("names" = rownames(chtest_all[["pvalues"]]),chtest_all[["pvalues"]]))
## 
## correl = as.data.frame(correl[,-1])
## rownames(correl) <- colnames(correl)
## 
## matable <- xtable(x = correl , label = "quali_cor")
## # Notez les doubles \\ nécessaires dans R, c'est la "double escape rule"
## print(matable, file = "fig/quali_cor.tex", size = "tiny", NA.string = "NA",
##       table.placement = "!t",
##       floating = FALSE,
##       caption.placement="top",
##       include.rownames = TRUE,
##       include.colnames = TRUE,
##       latex.environments = "center")
## # On veut des '.' au lieu des des NA


###################################################
### code chunk number 10: quanti_cor (eval = FALSE)
###################################################
## quantit <- c("HRA", "SYS")
## # 1) Check normality -- graphically
## 
## # ggqqplot(icu, x = quantit, combine = TRUE) # compares variable distribution with a theoretical normal distr
## 
## # 2) Check normality -- test SW
## hra_shap <- shapiro.test(icu$HRA)
## sys_shap <- shapiro.test(icu$SYS)
## # confirmé normalité - on peut faire cor.test sous hypothèse distrib paramétrique (Pearson)
## var_corelations <- cor(icu[quantit])
## var_pvalues <- var_corelations ##### pourquoi ces deux lignes?????
## var_pvalues[] <- NA
## for (i in (2:nrow(var_corelations))) {
##   for (j in (1:(i - 1)))
##   {
##     # cat("\ncor.test entre",colnames(var_corelations)[i],"et",colnames(var_corelations)[j])
##     var_pvalues[i, j] <- round(cor.test(icu[, colnames(var_corelations)[i]],
##       icu[, colnames(var_corelations)[j]],
##       method = "kendall",
##       alternative = "two.sided"
##     )$p.value,
##     digits = 4
##     )
##   }
## }
## var_pvalues[2,1] # pas de corrélation significative


###################################################
### code chunk number 11: tri (eval = FALSE)
###################################################
## icu = icu[,c(2:4,6,9,11:14,18)]


###################################################
### code chunk number 12: step (eval = FALSE)
###################################################
## # sample
## glmBase_s <- glm(STA ~ 1, data = icuvs_1, family = "binomial")
## stepwise_s <- stepAIC(glmBase_s, scope = list(
##   upper = ~ AGE * HRA * SYS * INF * SER * GENDER * PRE * TYP * PCO,
##   lower = ~1
## ), direction = "both", trace = F)
## cat("=================[sample data]=================")
## stepwise_s
## 
## # total
## glmBase_t <- glm(STA ~ 1, data = icu, family = "binomial")
## stepwise_t <- stepAIC(glmBase_t, scope = list(
##   upper = ~ AGE * HRA * SYS * INF * SER * GENDER * PRE * TYP * PCO,
##   lower = ~1
## ), direction = "both", trace = F)
## cat("=================[total data]=================")
## stepwise_t


###################################################
### code chunk number 13: model_70% (eval = FALSE)
###################################################
## # avec les variables choisies par le stepwise: TYP, AGE, SYS, AGE:SYS
## mod_s = glm(STA ~ AGE + TYP + SYS, data = icuvs_1, family = "binomial")
## summary(mod_s)
## anova(mod_s, test = "Chisq")
## 
## mod_s1 <- glm(STA ~ AGE, data = icuvs_1, family = "binomial")
## summary(mod_s1)
## anova(mod_s1, test = "Chisq")
## 
## mod_s2 <- glm(STA ~ SYS, data = icuvs_1, family = "binomial")
## summary(mod_s2)
## anova(mod_s2, test = "Chisq")
## 
## mod_s3 <- glm(STA ~ TYP, data = icuvs_1, family = "binomial")
## summary(mod_s3)
## anova(mod_s3, test = "Chisq")
## 
## rs_1 = anova(mod_s, mod_s1, test = "Chisq")
## rs_2 = anova(mod_s, mod_s2, test = "Chisq")
## rs_3 = anova(mod_s, mod_s3, test = "Chisq")
## 
## results1 = data.frame(
##   "AGE" = c(rs_1$`Resid. Dev`[2],rs_1$`Resid. Dev`[1]), 
## "SYS" = c(rs_2$`Resid. Dev`[2],rs_2$`Resid. Dev`[1]), 
## "TYP" = c(rs_3$`Resid. Dev`[2],rs_3$`Resid. Dev`[1]),
## row.names = c("Dév. modèle 1 variables", "Dév. modèle 3 variables"))
## 
## # pour tout le jeu de données
## mod = glm(STA ~ AGE + TYP + SYS, data = icu, family = "binomial")
## summary(mod)
## anova(mod, test = "Chisq")
## 
## mod1 <- glm(STA ~ AGE, data = icu, family = "binomial")
## summary(mod1)
## anova(mod1, test = "Chisq")
## 
## mod2 <- glm(STA ~ SYS, data = icu, family = "binomial")
## summary(mod2)
## anova(mod2, test = "Chisq")
## 
## mod3 <- glm(STA ~ TYP, data = icu, family = "binomial")
## summary(mod3)
## anova(mod3, test = "Chisq")
## 
## r1 = anova(mod, mod1, test = "Chisq")
## r2 = anova(mod, mod2, test = "Chisq")
## r3 = anova(mod, mod3, test = "Chisq")
## 
## results2 = data.frame(
##   "AGE" = c(r1$`Resid. Dev`[2],r1$`Resid. Dev`[1]), 
## "SYS" = c(r2$`Resid. Dev`[2],r2$`Resid. Dev`[1]), 
## "TYP" = c(r3$`Resid. Dev`[2],r3$`Resid. Dev`[1]),
## row.names = c("Dév. modèle 1 variables", "Dév. modèle 3 variables"))
## 
## 
## matable <- xtable(x = results1 , label = "result1")
## # Notez les doubles \\ nécessaires dans R, c'est la "double escape rule"
## print(matable, file = "fig/result1.tex", NA.string = "NA",
##       table.placement = "!t",
##       floating = FALSE,
##       caption.placement="top",
##       include.rownames = TRUE,
##       include.colnames = TRUE,
##       latex.environments = "center")
## 
## 
## matable <- xtable(x = results2 , label = "result2")
## # Notez les doubles \\ nécessaires dans R, c'est la "double escape rule"
## print(matable, file = "fig/result2.tex", NA.string = "NA",
##       table.placement = "!t",
##       floating = FALSE,
##       caption.placement="top",
##       include.rownames = TRUE,
##       include.colnames = TRUE,
##       latex.environments = "center")


###################################################
### code chunk number 14: validation model 70% (eval = FALSE)
###################################################
## # d'abord pour SAMPLE :::::::::::::::::::::::::::::::::::::::::::::::
## mod_s = glm(STA ~ AGE + TYP + SYS, data = icuvs_1, family = "binomial")
## 
## # d'après pdf regression logistique
## rh1_s = hoslem.test(icuvs_1$STA, fitted(mod_s))
## 
## # Analyse des résidus: d'après pdf regression logistique diapo 34 - on confirme distribution normale
## # mais ça ne sert plus à grande chose à part ça
## plot(rstudent(mod_s), type = "p", ylim = c(-3,3))
## abline(h = 2, col = "red")
## abline(h = -2, col = "red")
## 
## 
## # distance de Cook: y a-t-il des points très influents ? 
## plot(cooks.distance(mod_s), xlab = "Individus", ylab = "")
## cooks.distance(mod_s)[cooks.distance(mod_s)>1]


###################################################
### code chunk number 15: performance model 70% (eval = FALSE)
###################################################
## library(pROC)
## # d'abord pour SAMPLE :::::::::::::::::::::::::::::::::::::::::::::::
## pred.test = predict(mod_s)
## roc_test = roc(icuvs_1$STA, pred.test)
## pred.test.bin = ifelse(pred.test>0.5, 1, 0) # transformation en binaire
## 
## counting_s = table(icuvs_1$STA, pred.test.bin)
## 
## plot(roc_test)
## auc_s = roc_test$auc # il y a plus de surface en dessous de la courbe: donc plus d'erreurs... étonnant ¬¬
## 
## 
## # puis pour TOUT JEU DONNÉES  :::::::::::::::::::::::::::::::::::::::::::::::
## pred_tot = predict(mod)
## roc_tot = roc(icu$STA, pred_tot)
## pred_tot_bin = ifelse(pred_tot>0.5, 1, 0) # transformation en binaire
## 
## counting_t = table(icu$STA, pred_tot_bin)
## 
## plot(roc_tot)
## auc_t = roc_tot$auc


###################################################
### code chunk number 16: model pour 100% et graph (eval = FALSE)
###################################################
## # d'abord pour SAMPLE :::::::::::::::::::::::::::::::::::::::::::::::
## mod_GOOD = glm(STA ~ AGE * SYS + TYP, data = icu, family = "binomial")
## 
## # Analyse des résidus: d'après pdf regression logistique diapo 34 - on confirme distribution normale
## # mais ça ne sert plus à grande chose à part ça
## plot(rstudent(mod_GOOD), type = "p", ylim = c(-3,3))
## abline(h = 2, col = "red")
## abline(h = -2, col = "red")
## 
## 
## # distance de Cook: y a-t-il des points très influents ? non plus
## plot(cooks.distance(mod_GOOD))
## cooks.distance(mod_GOOD)[cooks.distance(mod_GOOD)>1]
## 
## 
## # Analyses performance du modèle
## pred_tot_g = predict(mod_GOOD)
## roc_tot_g = roc(icu$STA, pred_tot_g)
## pred_tot_gbin = ifelse(pred_tot_g>0.5, 1, 0) # transformation en binaire
## 
## counting_t_g = table(icu$STA, pred_tot_gbin)
## 
## plot(roc_tot_g)
## auc_t_g = roc_tot_g$auc
## 
## 
## #graphique
## par(mfrow= c(2,1))
## pred.tot2 = round(pred_tot) # complètement à la rache
## plot(icu$STA, col = "green", main = "Sans effets croisés")
## points(pred.tot2, col = "red")
## 
## 
## pred.tot2_g = round(pred_tot_g) # complètement à la rache
## plot(icu$STA, col = "green", main = "Avec effets croisés")
## points(pred.tot2_g, col = "red") 
## # bon... pas terrible
## 


###################################################
### code chunk number 17: table (eval = FALSE)
###################################################
## matable <- xtable(x = summary(icuf[,-c(1,19:21)]) , label = "tset",
## caption = "Lblablabla")
## # Notez les doubles \\ nécessaires dans R, c'est la "double escape rule"
## print(matable, file = "fig/tset.tex", size = "tiny", NA.string = "-")
## # On veut des '.' au lieu des des NA


