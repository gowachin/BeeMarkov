#' BeeMarkov
#'
#' A package to study markov chains.
#'
#' @docType package
#' @name beemarkov
NULL


#' @import graphics
#' @import seqinr
#' @import utils
NULL

#' transition
#'
#' Compute the transition model for one Markov model on DNA.
#'
#' @param file a file (fasta) to read and use as training for model
#' @param l_word length of words for the model. Equal to the "order of the model + 1" 
#' @param n_seq number of sequences to train the model with
#' @param log boolean information if the transition matrix of the model is in log or not.
#' @param type a indicator ("+" or "-") for printing pretty stuff during computations
#'
#'
#' @author Jaunatre Maxime <maxime.jaunatre@etu.univ-grenoble-alpes.fr>
#'
#' @export
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

#' quality
#'
#' Compute the quality of the model and return a TRUE Positive or FALSE Negative information. 
#' It test the positive and negative model and assign every sequence of the fasta file to one model or the other. 
#' Therefore, the file must countain sequences which are know to be from one model.
#'
#' @param file a file (fasta) to read and use to test model
#' @param pos_training a file (fasta) to read and train the positive model
#' @param neg_training a file (fasta) to read and train the negative model
#' @param trans_pos a transition matrix if it was already computed. Therefore, no need to train models. Warning, it must be in log
#' @param trans_neg a transition matrix if it was already computed. Therefore, no need to train models. Warning, it must be in log
#' @param l_word_pos length of words for the model. Equal to the "order of the model + 1" 
#' @param l_word_neg length of words for the model. Equal to the "order of the model + 1" 
#' @param n_train  number of sequences to train with
#' @param n_seq number of sequences to analyse
#' @param quiet if some informations are print or not (boolean)
#'
#'
#' @author Jaunatre Maxime <maxime.jaunatre@etu.univ-grenoble-alpes.fr>
#'
#' @export
quality <- function(file, # fichier
                    pos_training = NULL, # file to train with for positive
                    neg_training = NULL, # file to train with for negative
                    trans_pos = NULL, #transition matrice pos
                    trans_neg = NULL, #transition matrice pos
                    l_word_pos = 1, # word lenght for transition table positif
                    l_word_neg = 1, # word length for transition table negatif
                    n_train = 1, # number of sequences to train with
                    n_seq = 1, # number of sequences to analyse
                    quiet = FALSE
){
  if(is.null(c(trans_pos,trans_neg))){
    trans_pos <-transition(file = pos_training, n_seq = n_train, l_word = l_word_pos, type ="+")
    trans_neg <- transition(file = neg_training, n_seq = n_train, l_word = l_word_neg, type = "-")
  } 
  
  seq <-seqinr::read.fasta(file)
  Nseq <- length(seq)
  
  if(Nseq < n_seq){ # check if user want to input too many seq in the matrix learning
    stop(paste("n_seq is larger than the number of sequence in this file ( ",Nseq," sequences for this file).", sep = ""))
  }
  
  result <- data.frame(VP = rep(FALSE,n_seq),
                       FN = rep(TRUE,n_seq),
                       pos = rep(NA,n_seq),
                       neg = rep(NA,n_seq)
  )
  
  p_init_pos = log(1/4^l_word_pos)
  p_init_neg = log(1/4^l_word_neg)
  
  if(!quiet) {cat('      ============ Computing sensi and speci for the test sequences ============       \n')
    pb <- utils::txtProgressBar(min = 0, max = n_seq, style = 3)}
  for(i in 1:n_seq){
    if(!quiet) utils::setTxtProgressBar(pb, i)
    n_word_pos <-seqinr::count(seq[[i]], l_word_pos)
    n_word_neg <-seqinr::count(seq[[i]], l_word_neg)
    result$pos[i] <- p_init_pos + sum( trans_pos * n_word_pos )
    result$neg[i] <- p_init_neg + sum( trans_neg * n_word_neg )
    if(result$pos[i] > result$neg[i]){
      result$VP[i] <- TRUE ; result$FN[i] <- FALSE 
    } else { 
      result$VP[i] <- FALSE ; result$FN[i] <- TRUE
    }
  }
  if(!quiet) close(pb)
  
  tmp <- colSums(result[,1:2])
  return(tmp)
}

#' threshold
#'
#' Compute the quality for multiple set of values in vectors. 
#' Compute TRUE Positive or FALSE Negative information and FALSE Positive and TRUE Negative. 
#' All of this allow to compute sensitivoty and specificity of every model.
#'  
#' Compute the quality of the model and return a TRUE Positive or FALSE Negative information. 
#' It test the positive and negative model and assign every sequence of the fasta file to one model or the other. 
#' Therefore, the file must countain sequences which are know to be from one model.
#'
#' @param pos_test a file (fasta) to read and use to test positive model
#' @param neg_test a file (fasta) to read and use to test negative model
#' @param pos_training a file (fasta) to read and train the positive model
#' @param neg_training a file (fasta) to read and train the negative model
#' @param pos_seq a vector for lengths of words for the model. Equal to the "order of the model + 1" 
#' @param neg_seq a vector fof lengths of words for the model. Equal to the "order of the model + 1" 
#' @param n_train  number of sequences to train with
#' @param n_seq number of sequences to analyse
#'
#'
#' @author Jaunatre Maxime <maxime.jaunatre@etu.univ-grenoble-alpes.fr>
#'
#' @export
threshold <- function(pos_test, # fichier
                     neg_test,
                     pos_training, # file to train with for positive
                     neg_training, # file to train with for negative
                     pos_seq = c(1:2),
                     neg_seq = c(1:2),
                     n_train = 1, # number of sequences to train with
                     n_seq = 1 # number of sequences to analyse
){
  
  choix =utils::menu(c("yes","no"),
               title = "You will launch long computation, do you wish to procede further ?")
  if (choix ==1) cat("\n      ============ Go take a good coffee ============       \n\n")
  if (choix ==2) stop("You stopped the computations")
  
  
  trans_pos <- list()
  trans_neg <- list()
  
  cat('\n      ============ Training modeles ============       \n\n')
  for(i in pos_seq){
    trans_pos[[i]] <- transition(file = pos_training, n_seq = n_train, l_word = i, type = "+")
  }
  cat('\n')
  for(j in neg_seq){
    trans_neg[[j]] <- transition(file = neg_training, n_seq = n_train, l_word = j, type = "-")
  }
  
  cat('\n      ============ Training modeles ============       \n\n')    
  sensi <- speci <- matrix(rep(0,length(pos_seq)*length(neg_seq)),ncol = length(pos_seq),nrow = length(neg_seq))
  for(i in pos_seq){
    for(j in neg_seq){
      cat('\n')
      cat(paste("      ============ Model M+ (",i,"/",j,") ============        \n", sep = ""))
      pos <- quality(file = pos_test, # fichier
                     trans_pos = trans_pos[[i]], #transition matrice pos
                     trans_neg = trans_neg[[j]], #transition matrice neg
                     l_word_pos = i, # transition table positif
                     l_word_neg = j, # transition table negatif
                     n_train = n_train, # number of sequences to train with
                     n_seq = n_seq, # number of sequences to analyse
                     quiet = TRUE
      )
      cat(paste("      ============ Model M- (",i,"/",j,") ============        \n", sep = ""))
      neg <- quality(file = neg_test, # fichier
                     trans_pos = trans_pos[[i]], #transition matrice pos
                     trans_neg = trans_neg[[j]], #transition matrice neg
                     l_word_pos = i, # transition table positif
                     l_word_neg = j, # transition table negatif
                     n_train = n_train, # number of sequences to train with
                     n_seq = n_seq, # number of sequences to analyse
                     quiet = TRUE
      )
      
      sensi[i,j] <- pos[1] / sum(pos)
      speci[i,j] <- neg[2] / sum(neg)
      
    }
  }
  
  cat(paste("      ============ Come back from your coffee ============        \n", sep = ""))
  colnames(sensi) <- colnames(speci) <- paste("-",neg_seq, sep="")
  rownames(sensi) <- rownames(speci) <- paste("+",pos_seq, sep="")
  
  final <-list(sensi,speci) ; names(final) = c("sensi","speci")
  
  return(final)
}
