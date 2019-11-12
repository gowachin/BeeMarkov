#### Script christelle #####
# library(seqinr)
# Lecture d’un fichier fasta :
cpg_A=read.fasta(file= "raw_data/mus_cpg_app.fa")
tem_A=read.fasta(file= "raw_data/mus_tem_app.fa")
# Donne le nombre de séquences dans le fichier
length(cpg_A)
# Extraction de la première séquence du fichier
cpg_A1=cpg_A[[1]]
# Fonction compte : count the words with a certain number of letter
count(cpg_A1,3)

file <- "raw_data/mus_cpg_app.fa"

transition <- function(file, # fichier
                     l_word = 1, # longueur des mots
                     n_seq = 1 # Nombre de sequences à analyser
                     #log = FALSE
                     ){
  seq <- read.fasta(file)
  Nseq <- length(seq)
  
  if(Nseq < n_seq){ # check if user want to input too many seq in the matrix learning
    stop(paste("n_seq is larger than the number of sequence in this file ( ",Nseq," sequences for this file).", sep = ""))
  }
  
  l <- sapply(seq,length)
  if(l_word > min(l) | l_word >= 10){
    stop("This is really too much for me, abort!!!")
  }
  
  tmp <-count(seq[[1]], l_word) +1 # add 1 occurence to have at least 1 obs
  for(i in 2:n_seq){
    tmp <- tmp + count(seq[[i]], l_word)
  }
  
  # l_count = function(Seq,n = l_word){count(Seq,n)} # alternativ way, slower but pretty
  # tmp <- rowSums(sapply(seq,l_count))

  # if(log){
  #   log(tmp) - log(sum(tmp))
  # } else {
  log(tmp / sum(tmp))
  # }
}
# mP <-transition(file = "raw_data/mus_cpg_app.fa", n_seq = 1160, l_word = 2)
# mM <- transition(file = "raw_data/mus_tem_app.fa", n_seq = 1160, l_word = 1)
# 
# mP ; mM

control <- function(file, # fichier
                    pos_training, # file to train with for positive
                    neg_training, # file to train with for negative
                    l_word_pos = 1, # word lenght for transition table positif
                    l_word_neg = 1, # word length for transition table negatif
                    n_train = 1, # number of sequences to train with
                    n_seq = 1 # number of sequences to analyse
){
  
  trans_pos <-transition(file = pos_training, n_seq = n_train, l_word = l_word_pos)
  trans_neg <- transition(file = neg_training, n_seq = n_train, l_word = l_word_neg)
  
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
  
  for(i in 1:n_seq){
    n_word_pos <- count(seq[[i]], l_word_pos)
    n_word_neg <- count(seq[[i]], l_word_neg)
    result$pos[i] <- p_init_pos + sum( trans_pos * n_word_pos )
    result$neg[i] <- p_init_neg + sum( trans_neg * n_word_neg )
    if(result$pos[i] < result$neg[i]){
      result$VP[i] <- TRUE ; result$FN[i] <- FALSE 
    } else { 
      result$VP[i] <- FALSE ; result$FN[i] <- TRUE
      }
  }
  tmp <- colSums(result[,1:2])
  return(tmp)
}

control(file = "raw_data/mus_tem_test.fa", # fichier
        pos_training = "raw_data/mus_cpg_app.fa", # file to train with
        neg_training = "raw_data/mus_tem_app.fa", # file to train with
        l_word_pos = 2, # transition table positif
        l_word_neg = 2, # transition table negatif
        n_train = 1160, # number of sequences to train with
        n_seq = 1163 # number of sequences to analyse
)

assembly <- function(pos_test, # fichier
                     neg_test,
                     pos_training, # file to train with for positive
                     neg_training, # file to train with for negative
                     pos_seq = c(1:2),
                     neg_seq = c(1:2),
                     n_train = 1, # number of sequences to train with
                     n_seq = 1 # number of sequences to analyse
){
  
  sensi <- speci <- matrix(rep(0,length(pos_seq)*length(neg_seq)),ncol = length(pos_seq),nrow = length(neg_seq))
  
  for(i in 1:length(pos_seq)){
    for(j in 1:length(neg_seq)){
      
      pos <- control(file = pos_test, # fichier
                     pos_training = pos_training, # file to train with
                     neg_training = neg_training, # file to train with
                     l_word_pos = pos_seq[i], # transition table positif
                     l_word_neg = neg_seq[j], # transition table negatif
                     n_train = n_train, # number of sequences to train with
                     n_seq = n_seq # number of sequences to analyse
      )
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
      
      # sensi = VP / VP + FN
      # speci = VN / VN + FP
    }
  }
  return(list(sensi,speci))
}

# doit merdouiller ici, a voir
assembly("raw_data/mus_cpg_test.fa", # fichier
         "raw_data/mus_tem_test.fa",
         "raw_data/mus_cpg_app.fa", # file to train with for positive
         "raw_data/mus_tem_app.fa", # file to train with for negative
         pos_seq = c(1:2),
         neg_seq = c(1:2),
         n_train = 1160, # number of sequences to train with
         n_seq = 1163 # number of sequences to analyse
)

# scrap ####
test = count(cpg_A1,2)
sum(matrix(test/sum(test),ncol = 4)[,1])
# il crée bien la matrice comme il faut avec count
titres = matrix(names(test), ncol = 4)
titres


# time control ####
try = 6
l_word = try
n_seq = 1160 
seq <- cpg_A
Nseq <- length(seq)

v1 = function(){
test <-count(cpg_A[[1]], l_word)
l_count = function(seq,n = try){count(seq,n)}
test = rowSums(sapply(cpg_A,l_count) )
test = test / sum(test)
}

v2 = function(){
  tmp <-count(seq[[1]], l_word)
  for(i in 2:n_seq){
    tmp <- tmp + count(seq[[i]], l_word)
  }
tmp / sum(tmp)
}

system.time(v1())
system.time(v2())