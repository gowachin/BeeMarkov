#### Script christelle #####
# library(seqinr)

# Lecture d’un fichier fasta :
# cpg_A=read.fasta(file= "raw_data/mus_cpg_app.fa")
# tem_A=read.fasta(file= "raw_data/mus_tem_app.fa")
# Donne le nombre de séquences dans le fichier
# length(cpg_A)
# Extraction de la première séquence du fichier
# cpg_A1=cpg_A[[1]]
# Fonction compte : count the words with a certain number of letter
# count(cpg_A1,3)

transition <- function(file, # fichier
                       l_word = 1, # longueur des mots
                       n_seq = 1, # Nombre de sequences à analyser
                       l_need = NULL, # length to downgrade if needed for model comparaison
                       log = TRUE,
                       type = "") {
  seq <- read.fasta(file)
  Nseq <- length(seq)

  if (Nseq < n_seq) { # check if user want to input too many seq in the matrix learning
    stop(paste("n_seq is larger than the number of sequence in this file ( ", Nseq, " sequences for this file).", sep = ""))
  }

  l <- sapply(seq, length)
  if (l_word > min(l) | l_word >= 10) {
    stop("This is really too much for me, abort!!!")
  }

  if (!is.null(l_need)) {
    if (l_need > l_word) {
      stop("l_need is higher than l_word, you can't downscale higter, by definition...")
    }else if(l_need == l_word){
      l_need = NULL
      warning("l_need is equal than l_word, therefore there is no downscale here")
    }
  }

  tmp <- count(seq[[1]], l_word) + 1 # add 1 occurence to have at least 1 obs

  cat(paste("      ============ Training model M",type,l_word," ============        \n", sep = ""))
  pb <- txtProgressBar(min = 0, max = n_seq, style = 3)
  for (i in 2:n_seq) {
    setTxtProgressBar(pb, i)

    tmp <- tmp + count(seq[[i]], l_word)
  }
  close(pb)

  # l_count = function(Seq,n = l_word){count(Seq,n)} # alternativ way, slower but pretty
  # tmp <- rowSums(sapply(seq,l_count))
  
  # cat('avant modif')
  # print(tmp)
  
  # transformation to length needed
  if (!is.null(l_need)) {
    gap = l_word - l_need
    for(i in 1:gap){
      wind = 4^i
      for(j in 0:((length(tmp)/wind)-1)){
        tmp[(1+wind*j):(wind+wind*j)] <- tmp[(1+wind*j):(wind+wind*j)] * 4^(i-1) / sum(tmp[(1+wind*j):(wind+wind*j)] )
        # print(tmp[(1+wind*j):(wind+wind*j)])
      }
      cat(paste('compress round = ', i," \n", sep=""))
      cat(paste('compress wind = ', wind," \n", sep=""))
    }
  } else {
    tmp = tmp / sum(tmp)
  }
  
  # cat('apres modif')
  # print(tmp)
  
  # possibility to compute without log
  if (log) {
    tmp = log(tmp)
  }
  
  # cat('final')
  # print(tmp)
  return(tmp)
}
# mP <-transition(file = "raw_data/mus_cpg_app.fa", n_seq = 1160, l_word = 2, log = F)
# mM <- transition(file = "raw_data/mus_tem_app.fa", n_seq = 1160, l_word = 1, log = F)
# 
# mP ; mM
# 
# m1 <-transition(file = "raw_data/mus_cpg_app.fa", n_seq = 1160, l_word = 1, log = F)
# m2 <-transition(file = "raw_data/mus_cpg_app.fa", n_seq = 1160, l_word = 2, log = F)
# m3 <-transition(file = "raw_data/mus_cpg_app.fa", n_seq = 1160, l_word = 3, log = F, l_need = 2)
# 
# m1 ; m2 ; m3

control <- function(file, # fichier
                    pos_training, # file to train with for positive
                    neg_training, # file to train with for negative
                    l_word_pos = 1, # word lenght for transition table positif
                    l_word_neg = 1, # word length for transition table negatif
                    n_train = 1, # number of sequences to train with
                    n_seq = 1 # number of sequences to analyse
){
  if(l_word_pos == l_word_neg){
    trans_pos <-transition(file = pos_training, n_seq = n_train, l_word = l_word_pos, type ="+")
    trans_neg <- transition(file = neg_training, n_seq = n_train, l_word = l_word_neg, type = "-")
  } else {
    # need to downscale one model
    l_order <- order(c(l_word_pos,l_word_neg) )
    # if(diff(l_order)>0){
    #   # print("l_word_neg > l_word_pos")
    #   trans_pos <-transition(file = pos_training, n_seq = n_train, l_word = l_word_pos, type ="+")
    #   trans_neg <- transition(file = neg_training, n_seq = n_train, l_word = l_word_neg, type = "-", l_need = l_word_pos)
    # }else{
    #   # print("l_word_neg < l_word_pos")
    #   trans_pos <-transition(file = pos_training, n_seq = n_train, l_word = l_word_pos, type ="+", l_need = l_word_neg)
    #   trans_neg <- transition(file = neg_training, n_seq = n_train, l_word = l_word_neg, type = "-")
    # }
    
    # in fact, need to downscale every model
    trans_pos <-transition(file = pos_training, n_seq = n_train, l_word = l_word_pos, type ="+", l_need = 1)
    trans_neg <- transition(file = neg_training, n_seq = n_train, l_word = l_word_neg, type = "-", l_need = 1)
  }
  
  seq <- read.fasta(file)
  Nseq <- length(seq)
  
  if(Nseq < n_seq){ # check if user want to input too many seq in the matrix learning
    stop(paste("n_seq is larger than the number of sequence in this file ( ",Nseq," sequences for this file).", sep = ""))
  }
  
  # maybe an option with equivalent transition matrice
  # if(is.null(transition) & force){ # check if user want to force the computation with a transition table with equivalent P()
  #   warning("No transition table input. By default, will use a transition table with equivalent values")
  #   transtition = 1
  # } else if(is.null(transition) & !force) {
  #   
  # }
  
  result <- data.frame(VP = rep(FALSE,n_seq),
                       FN = rep(TRUE,n_seq),
                       pos = rep(NA,n_seq),
                       neg = rep(NA,n_seq)
                       )
  
  p_init_pos = log(1/4^l_word_pos)
  p_init_neg = log(1/4^l_word_neg)
  
  cat('      ============ Computing sensi and speci for the test sequences ============       \n')
  pb <- txtProgressBar(min = 0, max = n_seq, style = 3)
  for(i in 1:n_seq){
    setTxtProgressBar(pb, i)
    n_word_pos <- count(seq[[i]], l_word_pos)
    n_word_neg <- count(seq[[i]], l_word_neg)
    result$pos[i] <- p_init_pos + sum( trans_pos * n_word_pos )
    result$neg[i] <- p_init_neg + sum( trans_neg * n_word_neg )
    if(result$pos[i] > result$neg[i]){
      result$VP[i] <- TRUE ; result$FN[i] <- FALSE 
    } else { 
      result$VP[i] <- FALSE ; result$FN[i] <- TRUE
      }
  }
  close(pb)
  
  tmp <- colSums(result[,1:2])
  return(tmp)
}

# control(file = "raw_data/mus_cpg_test.fa", # fichier
#         pos_training = "raw_data/mus_cpg_app.fa", # file to train with
#         neg_training = "raw_data/mus_tem_app.fa", # file to train with
#         l_word_pos = 2, # transition table positif
#         l_word_neg = 4, # transition table negatif
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
){
  
  choix = menu(c("yes","no"),
               title = "You will launch long computation, do you wish to procede further ?")
  if (choix ==1) cat("\n      ============ Go take a good coffee ============       \n\n")
  if (choix ==2) stop("You stopped the computations")
  
  sensi <- speci <- matrix(rep(0,length(pos_seq)*length(neg_seq)),ncol = length(pos_seq),nrow = length(neg_seq))
  
  for(i in 1:length(pos_seq)){
    for(j in 1:length(neg_seq)){
      
      cat('\n      ============ Computing M+ ============       \n\n')
      pos <- control(file = pos_test, # fichier
                     pos_training = pos_training, # file to train with
                     neg_training = neg_training, # file to train with
                     l_word_pos = pos_seq[i], # transition table positif
                     l_word_neg = neg_seq[j], # transition table negatif
                     n_train = n_train, # number of sequences to train with
                     n_seq = n_seq # number of sequences to analyse
      )
      cat('\n      ============ Computing M- ============       \n')
      neg <- control(file = neg_test, # fichier
                     pos_training = pos_training, # file to train with
                     neg_training = neg_training, # file to train with
                     l_word_pos = pos_seq[i], # transition table positif
                     l_word_neg = neg_seq[j], # transition table negatif
                     n_train = n_train, # number of sequences to train with
                     n_seq = n_seq # number of sequences to analyse
      )

      sensi[i,j] <- pos[1] / sum(pos)
      speci[i,j] <- neg[2] / sum(neg)

      # cat(paste('sensi=',sensi[i,j],"\n",
      #           'speci=',speci[i,j],"\n",
      #           sep=''))
      
      # sensi = VP / VP + FN
      # speci = VN / VN + FP
    }
  }
  
  colnames(sensi) <- colnames(speci) <- paste("-",neg_seq, sep="")
  rownames(sensi) <- rownames(speci) <- paste("+",pos_seq, sep="")
  
  final <-list(sensi,speci) ; names(final) = c("sensi","speci")
  
  return(final)
}


# takes long, warning
cpg_fin = assembly("raw_data/mus_cpg_test.fa", # fichier
         "raw_data/mus_tem_test.fa",
         "raw_data/mus_cpg_app.fa", # file to train with for positive
         "raw_data/mus_tem_app.fa", # file to train with for negative
         pos_seq = c(1:7),
         neg_seq = c(1:7),
         n_train = 1160, # number of sequences to train with
         n_seq = 1163 # number of sequences to analyse
)

txt = "Computation has ended"
GET(paste('https://smsapi.free-mobile.fr/sendmsg?user=17267063&pass=CDZDEAQ49d1q3X&msg=%',
          paste(as.character(charToRaw(txt)), collapse = "%"), sep = ""))

# ça ne marche pas pour les couples où le minimum > 1 avec projet de base
# tout redescendua 1 pour avoir des resultats...a voir si d'un point de vue mathématique ça colle
final

library(reshape)

table = cbind(melt(cpg_fin$sensi), melt(cpg_fin$speci)[,3])
colnames(table) = c("M","m","Sensi","Speci")
table$tot = table$Sensi+table$Speci

library(ggplot2)

ggplot(table, aes(M,m)) +
  geom_raster(aes(fill = tot), hjust=0.5, vjust=0.5, interpolate=FALSE) +
  geom_contour(aes(z = tot))  

table[which(table$tot == max(table$tot)),]

# viterbi ####
viterbi <- function(file = "raw_data/mus1.fa",
                    pos_training = "raw_data/mus_cpg_app.fa", # file to train with for positive
                    neg_training = "raw_data/mus_tem_app.fa", # file to train with for negative
                    l_word_pos = 2,
                    l_word_neg = 1,
                    n_train = 1160,
                    n_ana = 1,
                    l_c = 1000, # length coding
                    l_nc = 125000 # length non-coding
                    ) {
  if(l_word_pos == l_word_neg){
    trans_pos <-transition(file = pos_training, n_seq = n_train, l_word = l_word_pos, type ="+")
    trans_neg <- transition(file = neg_training, n_seq = n_train, l_word = l_word_neg, type = "-")
  } else {
    # need to downscale one model
    l_order <- order(c(l_word_pos,l_word_neg) )
    
    trans_pos <-transition(file = pos_training, n_seq = n_train, l_word = l_word_pos, type ="+", l_need = 1)
    trans_neg <- transition(file = neg_training, n_seq = n_train, l_word = l_word_neg, type = "-", l_need = 1)
  }
  
  seq <- read.fasta(file)
  if(n_ana > 1) stop('Analysis for multiple files is not implemented yet')
  
  raw_seq <- seq[[1]]
  long <- length(raw_seq)
  
  beg <- max(l_word_pos,l_word_neg)
  
  proba <- matrix(rep(NA,long*3),ncol = 3) ; colnames(proba) <- c("M+","M-")
  
  proba[1:l_word_pos-1,1] = 0
  proba[1:l_word_pos-1,2] = 0
  # head(proba)
  
  cat('      ============ Viterbi is running ============       \n')
  pb <- txtProgressBar(min = beg, max = long, style = 3)
  for(i in beg:long){
    setTxtProgressBar(pb, i)
    # # proba d'avoir la base sous M
    # pM <- min( count(raw_seq[(i-l_word_pos+1):i], l_word_pos) * trans_pos)
    # 
    # # proba d'avoir la base sous M
    # pm <- min(count(raw_seq[(i-l_word_neg+1):i], l_word_neg) * trans_neg)
    # 
    # proba[i,1] = pM
    # proba[i,2] = pm
    if(proba[i,1]>proba[i,2]) {proba[i,3]=2} else {proba[i,3]=3}
  }
  close(pb)

  # head(proba)
}

plot(x = 1:long, y = rep(1,long), col = proba[,3])

# # scrap ####
# test = count(cpg_A1,2)
# sum(matrix(test/sum(test),ncol = 4)[,1])
# # il crée bien la matrice comme il faut avec count
# titres = matrix(names(test), ncol = 4)
# titres


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
#     tmp <- tmp + count(seq[[i]], l_word)
#   }
# tmp / sum(tmp)
# }
# 
# system.time(v1())
# system.time(v2())
