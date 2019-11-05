#### Script christelle #####
library(seqinr)
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
  if(length(seq) < n_seq){ # check if user want to input too many seq in the matrix learning
    warning(paste("n_seq is larger than the number of sequence in this file ( ",length(seq)," sequences for this file).", sep = ""))
  }
  if(l_word > 10){
    warning("This is really too much for me, abort!!!")
  }
  tmp <-count(seq[[1]], l_word)
  for(i in 2:n_seq){
    tmp <- tmp + count(seq[[i]], l_word)
  }
  # if(log){
  #   log(tmp) - log(sum(tmp))
  # } else {
    tmp / sum(tmp)
  # }
}

test = count(cpg_A1,2)
sum(matrix(test/sum(test),ncol = 4)[,1])
# il crée bien la matrice comme il faut avec count
titres = matrix(names(test), ncol = 4)
titres