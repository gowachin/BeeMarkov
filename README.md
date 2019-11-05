# Markov

chaine de markov pour travailler sur l'adn, predire ilot cpg


## model M1

dependance la base precedente

besoin d'une matrice de transition et d'une distribution initiale (souvent 1/4)
somme des lihgnes de la matrice de transition == 1, car de ligne vers colonne

modeles très sequentiels. mais les modèles sont constant le long de la séquence

pour apprentissage d'un modele de type Mm, il faut compter les frequences des mots de taille m+1

## ilots CPG

travail sur un seul brin
probleme si region GC+ et GC- sur la meme sequence, faire sur une fenetre glissante (1ere proposition historique)

construire 2 modeles qui 
- caracterise les ilots cpg 
- caracterise les modles hors cpg

2 fichiers app pour construire les modeles "...-app" pour faire l'apprentissage
2 fichiers test pour tester les modeles "...-test" 

besoin library `seqinr` pour lire les fasta

modeles a faire de taille 1 a 6, a nous de choisir

indication stat : calculer les sensibilité et specificité


## To do : 


