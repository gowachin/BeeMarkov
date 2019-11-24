#### Script christelle #####
# library(seqinr)

# Lecture d’un fichier fasta :
# cpg_A=seqinr::read.fasta(file= "raw_data/mus_cpg_app.fa")
# tem_A=seqinr::read.fasta(file= "raw_data/mus_tem_app.fa")
# Donne le nombre de sequences dans le fichier
# length(cpg_A)
# Extraction de la premiere sequence du fichier
# cpg_A1=cpg_A[[1]]
# Fonction compte : seqinr::count the words with a certain number of letter
# seqinr::count(cpg_A1,3)

##### functions (in package) ####

transition <- function(file, # fichier
                       l_word = 1, # longueur des mots
                       n_seq = 1, # Nombre de sequences a analyser
                       # l_need = NULL, # length to downgrade if needed for model comparaison
                       log = TRUE,
                       type = "") {
  seq <- seqinr::read.fasta(file)
  Nseq <- length(seq)

  if (Nseq < n_seq) { # check if user want to input too many seq in the matrix learning
    stop(paste("n_seq is larger than the number of sequence in this file ( ", Nseq, " sequences for this file).", sep = ""))
  }

  l <- sapply(seq, length)
  if (l_word > min(l) | l_word >= 10) {
    stop("This is really too much for me, abort!!!")
  }

  # if (!is.null(l_need)) {
  #   if (l_need > l_word) {
  #     stop("l_need is higher than l_word, you can't downscale higter, by definition...")
  #   }else if(l_need == l_word){
  #     l_need = NULL
  #     warning("l_need is equal than l_word, therefore there is no downscale here")
  #   }else if(l_need < 1){
  #     l_need = NULL
  #     warning("l_need is inferior to 1, therefore there is no downscale here")
  #   }
  # }

  tmp <- seqinr::count(seq[[1]], l_word) + 1 # add 1 occurence to have at least 1 obs

  cat(paste("      ============ Training model M", type, l_word, " ============        \n", sep = ""))
  pb <- utils::txtProgressBar(min = 0, max = n_seq, style = 3)
  for (i in 2:n_seq) {
    utils::settxtProgressBar(pb, i)

    tmp <- tmp + seqinr::count(seq[[i]], l_word)
  }
  close(pb)

  # # alternativ way, slower but prettier
  # cat(paste("      ============ Training model M",type,l_word," ============        \n", sep = ""))
  # l_count = function(Seq,n = l_word){count(Seq,n)}
  # tmp <- rowSums(sapply(seq,l_count))

  # transformation to length needed
  # if (!is.null(l_need)) {
  if (l_word > 1) {
    # gap = l_word - l_need
    # for(i in 1:gap){
    i <- 1
    wind <- 4 # ^i
    for (j in 0:((length(tmp) / wind) - 1)) {
      tmp[(1 + wind * j):(wind + wind * j)] <- tmp[(1 + wind * j):(wind + wind * j)] * 4^(i - 1) / sum(tmp[(1 + wind * j):(wind + wind * j)])
      # print(tmp[(1+wind*j):(wind+wind*j)])
    }
    # cat(paste('compress round = ', i," \n", sep=""))
    # cat(paste('compress wind = ', wind," \n", sep=""))
    cat("      ============ Computing conditionnal probabilities ============        \n")
    # }
  } else {
    tmp <- tmp / sum(tmp)
  }

  # possibility to compute without log
  if (log) {
    tmp <- log(tmp)
  }

  return(tmp)
}
# mP <-transition(file = "raw_data/mus_cpg_app.fa", n_seq = 1160, l_word = 3, log = F) #, l_need = 2)
#  mM <- transition(file = "raw_data/mus_tem_app.fa", n_seq = 1160, l_word = 2, log = F)
#
# mP ; mM
#
# m1 <-transition(file = "raw_data/mus_cpg_app.fa", n_seq = 1160, l_word = 1, log = F)
# m2 <-transition(file = "raw_data/mus_cpg_app.fa", n_seq = 1160, l_word = 2, log = F)
# m3 <-transition(file = "raw_data/mus_cpg_app.fa", n_seq = 1160, l_word = 3, log = F) #, l_need = 2)
#
# m1 ; m2 ; m3

control <- function(file, # fichier
                    pos_training = NULL, # file to train with for positive
                    neg_training = NULL, # file to train with for negative
                    trans_pos = NULL, # transition matrice pos
                    trans_neg = NULL, # transition matrice pos
                    l_word_pos = 1, # word lenght for transition table positif
                    l_word_neg = 1, # word length for transition table negatif
                    n_train = 1, # number of sequences to train with
                    n_seq = 1, # number of sequences to analyse
                    quiet = FALSE) {
  if (is.null(c(trans_pos, trans_neg))) {
    # if(l_word_pos == l_word_neg){
    trans_pos <- transition(file = pos_training, n_seq = n_train, l_word = l_word_pos, type = "+")
    trans_neg <- transition(file = neg_training, n_seq = n_train, l_word = l_word_neg, type = "-")
    # } else {
    # # need to downscale one model
    # l_order <- order(c(l_word_pos,l_word_neg) )
    # if(diff(l_order)>0){
    #   # print("l_word_neg > l_word_pos")
    #   trans_pos <-transition(file = pos_training, n_seq = n_train, l_word = l_word_pos, type ="+")
    #   trans_neg <- transition(file = neg_training, n_seq = n_train, l_word = l_word_neg, type = "-", l_need = l_word_pos)
    # }else{
    #   # print("l_word_neg < l_word_pos")
    #   trans_pos <-transition(file = pos_training, n_seq = n_train, l_word = l_word_pos, type ="+", l_need = l_word_neg)
    #   trans_neg <- transition(file = neg_training, n_seq = n_train, l_word = l_word_neg, type = "-")
    # }

    # in fact, need to downscale every model, but just to l_word_pos -1
    #   trans_pos <-transition(file = pos_training, n_seq = n_train, l_word = l_word_pos, type ="+", l_need = l_word_pos-1)
    #   trans_neg <- transition(file = neg_training, n_seq = n_train, l_word = l_word_neg, type = "-", l_need = l_word_neg-1)
    # }
  }

  seq <- seqinr::read.fasta(file)
  Nseq <- length(seq)

  if (Nseq < n_seq) { # check if user want to input too many seq in the matrix learning
    stop(paste("n_seq is larger than the number of sequence in this file ( ", Nseq, " sequences for this file).", sep = ""))
  }

  # maybe an option with equivalent transition matrice
  # if(is.null(transition) & force){ # check if user want to force the computation with a transition table with equivalent P()
  #   warning("No transition table input. By default, will use a transition table with equivalent values")
  #   transtition = 1
  # } else if(is.null(transition) & !force) {
  #
  # }

  result <- data.frame(
    VP = rep(FALSE, n_seq),
    FN = rep(TRUE, n_seq),
    pos = rep(NA, n_seq),
    neg = rep(NA, n_seq)
  )

  p_init_pos <- log(1 / 4^l_word_pos)
  p_init_neg <- log(1 / 4^l_word_neg)

  if (!quiet) {
    cat("      ============ Computing sensi and speci for the test sequences ============       \n")
    pb <- utils::txtProgressBar(min = 0, max = n_seq, style = 3)
  }
  for (i in 1:n_seq) {
    if (!quiet) utils::settxtProgressBar(pb, i)
    n_word_pos <- seqinr::count(seq[[i]], l_word_pos)
    n_word_neg <- seqinr::count(seq[[i]], l_word_neg)
    result$pos[i] <- p_init_pos + sum(trans_pos * n_word_pos)
    result$neg[i] <- p_init_neg + sum(trans_neg * n_word_neg)
    if (result$pos[i] > result$neg[i]) {
      result$VP[i] <- TRUE
      result$FN[i] <- FALSE
    } else {
      result$VP[i] <- FALSE
      result$FN[i] <- TRUE
    }
  }
  if (!quiet) close(pb)

  tmp <- colSums(result[, 1:2])
  return(tmp)
}

# control(file = "raw_data/mus_cpg_test.fa", # fichier
#         pos_training = "raw_data/mus_cpg_app.fa", # file to train with
#         neg_training = "raw_data/mus_tem_app.fa", # file to train with
#         #trans_pos = 42,
#         l_word_pos = 2, # transition table positif
#         l_word_neg = 3, # transition table negatif
#         n_train = 1160, # number of sequences to train with
#         n_seq = 1163 # number of sequences to analyse
# )

assembly <- function(pos_test, # fichier
                     neg_test,
                     pos_training, # file to train with for positive
                     neg_training, # file to train with for negative
                     pos_seq = c(1:2),
                     neg_seq = c(1:2),
                     n_train = 1, # number of sequences to train with
                     n_seq = 1 # number of sequences to analyse
) {
  choix <- utils::menu(c("yes", "no"),
    title = "You will launch long computation, do you wish to procede further ?"
  )
  if (choix == 1) cat("\n      ============ Go take a good coffee ============       \n\n")
  if (choix == 2) stop("You stopped the computations")


  trans_pos <- list()
  trans_neg <- list()

  cat("      ============ Training modeles ============       \n\n")
  for (i in pos_seq) {
    trans_pos[[i]] <- transition(file = pos_training, n_seq = n_train, l_word = i, type = "+")
  }
  cat("\n")
  for (j in neg_seq) {
    trans_neg[[j]] <- transition(file = neg_training, n_seq = n_train, l_word = j, type = "-")
  }

  cat("
      ============ Training modeles ============       \n\n")
  sensi <- speci <- matrix(rep(0, length(pos_seq) * length(neg_seq)), ncol = length(pos_seq), nrow = length(neg_seq))
  for (i in pos_seq) {
    for (j in neg_seq) {
      cat("\n")
      cat(paste("      ============ Model M+ (", i, "/", j, ") ============        \n", sep = ""))
      pos <- control(
        file = pos_test, # fichier
        trans_pos = trans_pos[[i]], # transition matrice pos
        trans_neg = trans_neg[[j]], # transition matrice neg
        l_word_pos = i, # transition table positif
        l_word_neg = j, # transition table negatif
        n_train = n_train, # number of sequences to train with
        n_seq = n_seq, # number of sequences to analyse
        quiet = TRUE
      )
      cat(paste("      ============ Model M- (", i, "/", j, ") ============        \n", sep = ""))
      neg <- control(
        file = neg_test, # fichier
        trans_pos = trans_pos[[i]], # transition matrice pos
        trans_neg = trans_neg[[j]], # transition matrice neg
        l_word_pos = i, # transition table positif
        l_word_neg = j, # transition table negatif
        n_train = n_train, # number of sequences to train with
        n_seq = n_seq, # number of sequences to analyse
        quiet = TRUE
      )

      sensi[i, j] <- pos[1] / sum(pos)
      speci[i, j] <- neg[2] / sum(neg)

      # cat(paste('sensi=',sensi[i,j],"\n",
      #           'speci=',speci[i,j],"\n",
      #           sep=''))

      # sensi = VP / VP + FN
      # speci = VN / VN + FP
    }
  }

  cat(paste("      ============ Come back from your coffee ============        \n", sep = ""))
  colnames(sensi) <- colnames(speci) <- paste("-", neg_seq, sep = "")
  rownames(sensi) <- rownames(speci) <- paste("+", pos_seq, sep = "")

  final <- list(sensi, speci)
  names(final) <- c("sensi", "speci")

  return(final)
}

#### computations to choose model ####

# takes long, warning
cpg_fin <- assembly("raw_data/mus_cpg_test.fa", # fichier
  "raw_data/mus_tem_test.fa",
  "raw_data/mus_cpg_app.fa", # file to train with for positive
  "raw_data/mus_tem_app.fa", # file to train with for negative
  pos_seq = c(1:6),
  neg_seq = c(1:6),
  n_train = 1160, # number of sequences to train with
  n_seq = 1163 # number of sequences to analyse
)

txt <- "Computation has ended"
GET(paste("https://smsapi.free-mobile.fr/sendmsg?user=17267063&pass=CDZDEAQ49d1q3X&msg=%",
  paste(as.character(charToRaw(txt)), collapse = "%"),
  sep = ""
))

# ça ne marche pas pour les couples où le minimum > 1 avec projet de base
# tout redescendua 1 pour avoir des resultats...a voir si d'un point de vue mathematique ça colle
final <- cpg_fin

library(reshape)

table <- cbind(melt(cpg_fin$sensi), melt(cpg_fin$speci)[, 3])
colnames(table) <- c("M", "m", "Sensi", "Speci")
table$tot <- table$Sensi + table$Speci

library(ggplot2)

ggplot(table, aes(M, m)) +
  geom_raster(aes(fill = tot), hjust = 0.5, vjust = 0.5, interpolate = FALSE) +
  geom_contour(aes(z = tot))

table[which(table$tot == max(table$tot)), ]

#### viterbi functions ####
viterbi <- function(file = "raw_data/mus1.fa",
                    pos_training = "raw_data/mus_cpg_app.fa", # file to train with for positive
                    neg_training = "raw_data/mus_tem_app.fa", # file to train with for negative
                    l_word_pos = 1,
                    l_word_neg = 1,
                    n_train = 1160,
                    n_ana = 1,
                    l_c = 1000, # length coding
                    l_nc = 125000 # length non-coding
) {
  # if(l_word_pos == l_word_neg){
  trans_pos <- transition(file = pos_training, n_seq = n_train, l_word = l_word_pos, type = "+")
  trans_neg <- transition(file = neg_training, n_seq = n_train, l_word = l_word_neg, type = "-")
  # } else {
  #   # need to downscale one model
  #   l_order <- order(c(l_word_pos,l_word_neg) )
  #
  #   trans_pos <-transition(file = pos_training, n_seq = n_train, l_word = l_word_pos, type ="+", l_need = 1)
  #   trans_neg <- transition(file = neg_training, n_seq = n_train, l_word = l_word_neg, type = "-", l_need = 1)
  # }

  seq <- seqinr::read.fasta(file)
  if (n_ana > 1) stop("Analysis for multiple files is not implemented yet")

  raw_seq <- seq[[1]]
  long <- length(raw_seq)

  beg <- max(l_word_pos, l_word_neg)

  # compute p_initial and transition matrix for markov
  pos_init <- neg_init <- log(0.5)
  l_c <- 1 / l_c
  l_nc <- 1 / l_nc

  trans_mod <- log(matrix(c(
    1 - l_c, l_nc,
    l_c, 1 - l_nc
  ),
  ncol = 2, nrow = 2
  ))
  colnames(trans_mod) <- rownames(trans_mod) <- c("c", "nc")

  # initialisation
  ncol <- 7
  proba <- matrix(rep(NA, long * ncol), ncol = ncol)
  colnames(proba) <- c("M+", "M-", "model", "length", "rep_length", "begin", "end")
  # old version is v2 takes too long
  # v1 = function(){
  # trans_pos[which(names(trans_pos)==paste(raw_seq[(beg-l_word_pos+1):beg],collapse = ""))]
  # }
  # v2 = function(){
  # min(count(raw_seq[(beg-l_word_pos+1):beg], l_word_pos) * trans_pos)
  # }
  # system.time(v1())
  # system.time(v2())

  proba[beg, 1] <- trans_pos[which(names(trans_pos) == paste(raw_seq[(beg - l_word_pos + 1):beg], collapse = ""))] + pos_init
  proba[beg, 2] <- trans_neg[which(names(trans_neg) == paste(raw_seq[(beg - l_word_neg + 1):beg], collapse = ""))] + neg_init
  if (proba[beg, 1] > proba[beg, 2]) {
    proba[beg, 3] <- 1
  } else {
    proba[beg, 3] <- 2
  }
  proba[beg, c(4, 5)] <- 1
  tmp <- proba[beg, c(6, 7)] <- beg
  # base before initialisation...what to put??? ####
  # proba[1:beg-1,c(1,2)] <- -42
  proba[1:beg - 1, 5] <- beg - 1
  proba[1:beg - 1, 3] <- 3
  proba[1, 6] <- 1
  proba[beg - 1, c(4, 7)] <- beg - 1

  # head(proba) ; tmp

  # long = 1000

  beg <- beg + 1
  cat("      ============ Viterbi is running ============       \n\n")
  pb <- utils::txtProgressBar(min = beg, max = long, style = 3)
  for (i in beg:long) {
    utils::setTxtProgressBar(pb, i)

    # proba d'avoir la base sous M
    pM <- trans_pos[which(names(trans_pos) == paste(raw_seq[(i - l_word_pos + 1):i], collapse = ""))] + max(
      proba[i - 1, 1] + trans_mod[1, 1],
      proba[i - 1, 2] + trans_mod[2, 1]
    )

    # proba d'avoir la base sous m
    pm <- trans_neg[which(names(trans_neg) == paste(raw_seq[(i - l_word_neg + 1):i], collapse = ""))] + max(
      proba[i - 1, 2] + trans_mod[2, 2],
      proba[i - 1, 1] + trans_mod[1, 2]
    )

    proba[i, 1] <- pM
    proba[i, 2] <- pm
    # print(proba[i,])
    if (proba[i, 1] > proba[i, 2]) {
      proba[i, 3] <- 1
    } else {
      proba[i, 3] <- 2
    }
    # length information
    if (proba[i, 3] == proba[i - 1, 3]) {
      proba[i, 4] <- proba[i - 1, 4] + 1 # increase part length
      proba[i - 1, c(4, 7)] <- NA # erase length in previous ligne
    } else {
      proba[i, 4] <- 1 # initiate new part length
      proba[i - 1, 7] <- i - 1 # put end value of precedent part
      proba[c(tmp:(i - 1)), 5] <- proba[i - 1, 4] # rep value of length for precedent part
      proba[i, 6] <- tmp <- i # put begin value of the actual part
    }
  }
  close(pb)

  # closing table
  proba[i, 7] <- i # put end value of precedent part
  proba[c(tmp:(i)), 5] <- proba[i, 4] # rep value of length for precedent part

  # head(proba)
  # add a column for the line number
  proba <- cbind(c(1:dim(proba)[1]), proba)
  colnames(proba)[1] <- "n"

  return(proba)
}

library(BeeMarkov)
mus1 <- viterbi(
  l_word_pos = 5,
  l_word_neg = 4
)

plot(x = mus1[, 1], y = rep(1, max(mus1[, 1])), col = mus1[, 4])

head(mus1)

mus1[c(735:850), ]
min <- 1
max <- dim(mus1)[1]
plot(x = min:max, y = mus1[c(min:max), 4], col = mus1[c(min:max), 4])

plot(x = min:max, y = log10(mus1[c(min:max), 6]), col = mus1[c(min:max), 4])

# count the length of different parts
seq <- mus1 # as.data.frame(mus1)

smoothing <- function(seq,
                      l_word_pos = 1,
                      l_word_neg = 1,
                      smooth_win = 10) {
  beg <- max(l_word_pos, l_word_neg)

  seq <- cbind(seq, seq[, 4])
  colnames(seq)[9] <- c("smoothed")

  cat("      ============ Smoothing boucle ============       \n")
  pb <- utils::txtProgressBar(min = 1, max = max(seq[, 1]), style = 3)
  for (i in 1:dim(seq)[1]) {
    utils::setTxtProgressBar(pb, i)
    if (seq[i, 6] <= smooth_win) {
      seq[i, 9] <- 3
    }
  }
  close(pb)

  seq <- cbind(seq, seq[, c(5:8)])
  colnames(seq)[10:13] <- paste("S_", colnames(seq[, c(10:13)]), sep = "")

  seq[beg, 10] <- 1
  tmp <- seq[beg, c(12, 13)] <- beg
  seq[-c(1:beg), c(10:13)] <- NA

  head(seq)
  beg <- beg + 1
  cat("      ============ Smoothing length ============      \n")
  pb <- utils::txtProgressBar(min = beg - 1, max = dim(seq)[1], style = 3)
  for (i in beg:dim(seq)[1]) {
    utils::setTxtProgressBar(pb, i)
    if (seq[i, 9] == seq[i - 1, 9]) {
      seq[i, 10] <- seq[i - 1, 10] + 1 # increase part length
      seq[i - 1, c(10, 13)] <- NA # erase length in previous ligne
    } else {
      seq[i, 10] <- 1 # initiate new part length
      seq[i - 1, 13] <- i - 1 # put end value of precedent part
      seq[c(tmp:(i - 1)), 11] <- seq[i - 1, 10] # rep value of length for precedent part
      seq[i, 11] <- tmp <- i # put begin value of the actual part
    }
  }
  close(pb)

  # closing table
  seq[i, 13] <- i # put end value of precedent part
  seq[c(tmp:(i)), 11] <- seq[i, 10] # rep value of length for precedent part
  
  return(seq)
}

s_mus1 <- smoothing(mus1,5,4,20)

s20_mus1 <- s_mus1

s20_mus1[236317,]

mus1[236317,]

summary(as.factor(s6_mus1[,10]))
summary(as.factor(s20_mus1[,10]))

plot(x = s20_mus1[,1], y = log(s20_mus1[,11]), col = s20_mus1[,9])
abline(h = log(20), col = "blue")

seq[1:6,]
seq[c(735:850),]

summary(as.factor(seq[,5]))
summary(as.factor(seq[,10]))

hist(seq[,10],nclass = max(seq[,1]), xlim = c(1,20))

plot(x = seq[,1], y = seq[,])

sum(seq[,9] == seq[,4])/dim(seq)[1]

# time control ####
# try = 6
# l_word = try
# n_seq = 1160 
# seq <- cpg_A
# Nseq <- length(seq)
# 
# v1 = function(){
# test <-count(cpg_A[[1]], l_word)
# l_count = function(seq,n = try){count(seq,n)}
# test = rowSums(sapply(cpg_A,l_count) )
# test = test / sum(test)
# }
# 
# v2 = function(){
#   tmp <-count(seq[[1]], l_word)
#   for(i in 2:n_seq){
#     tmp <- tmp + seqinr::count(seq[[i]], l_word)
#   }
# tmp / sum(tmp)
# }
# 
# system.time(v1())
# system.time(v2())
