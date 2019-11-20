##						-*- R -*-
## empty Rprofile.site for R on Debian
##
## Copyright (C) 2008 - 2018  Dirk Eddelbuettel and GPL'ed

# ## Example of .Rprofile
#options(width=65, digits=5)
.First <- function() cat(paste(c("\n   Welcome to R!\n\n",
                                 " This is the Markov Project\n\n these packages are loaded by default :",
                                 pack,
                                 "\n\n")
                               , collapse = " "))

pack = c("beepr","seqinr","httr") #ajouter des packages par default

options(defaultPackages=c(getOption("defaultPackages"),
                          pack)) 

#cat(paste(c("these packages are loaded by default :",pack), collapse = " "))

.Last <- function()  cat("\n   Goodbye!\n\n")

## We set the cloud mirror, which is 'network-close' to everybody, as default
#local({
#  r <- getOption("repos")
#  r["CRAN"] <- "https://cloud.r-project.org"
#  options(repos = r)
#})
