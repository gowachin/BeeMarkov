#### compiling ####
setwd("manuscrit")
Sweave("Markov_Report_Guerra_Jaunatre.Rnw")
Stangle("Markov_Report_Guerra_Jaunatre.Rnw")
tools::texi2pdf("Markov_Report_Guerra_Jaunatre.tex") 
file.remove(list.files()[ which(! list.files() %in% c("child_test.Rnw", 
                                                      "fig", "tex",
                                                      "Markov_Report_Guerra_Jaunatre.Rnw",
                                                      "Markov_Report_Guerra_Jaunatre.R",
                                                      "Markov_Report_Guerra_Jaunatre.pdf",
                                                      "Markov_Report_Guerra_Jaunatre_oldversion.Rnw") )]
)
setwd("../")
system('ls')
system('xdg-open manuscrit/Markov_Report_Guerra_Jaunatre.pdf ') # works on ubuntu


############################
# Exo de matrices #
############################
P = matrix(c(0.3,0.2,0.7,0.8), ncol = 2) ; P
Pi = matrix(c(0.9,0.1), ncol = 2) ; Pi
P1 = Pi%*%P ; P1

# Exo 2
P = matrix(c(0.75,0,0.2,
             0.25,0.25,0,
             0,0.75,0.8), ncol = 3) ; P
Pi = matrix(c(0.9,0.1,0), ncol = 3) ; Pi
P1 = Pi%*%P ; P1
P2 = P1%*%P ; P2

# Exo 3
PA = matrix(c(0, 0.3, 0, 0.7,
             0.6, 0, 0.4, 0,
             0, 0.8, 0, 0.2,
             0.3, 0, 0.7, 0), ncol = 4) ; PA
Pi = matrix(c(1/8,0,7/8,0), ncol = 4) ; Pi
P1 = Pi%*%PA ; P1

PB = matrix(c(0.3, 0, 0.7, 0,
              0, 0.8, 0, 0.2,
              0.1, 0, 0.9, 0,
              0, 0.4, 0, 0.6), ncol = 4) ; PB
P1 = Pi%*%PB ; P1
