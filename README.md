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

## td 2

m+ et m- 
utilisation de fenetre glissante pour passer de l'une a l'autre, mais question sur la taille de la fenetre, le taux de recouvrement entre deux fenetre...

pourquoi ne pas regrouper les deux modeles en un seul.
processus observable qui identifie les lettres
processus caché qui identifie la sucession de régions (modele de markov d'ordre 1)

### segmentation de la sequence

differentes ecoles : 
approche supervisée, apprentissage a partir de sequences déjà connues, algorithme de Viterbi *a faire*
approche non supervisée, premier modele aleatoire, puis apprentissage, algorithme EM

risque d'avoir des zones de changement entre deux états m+ ou m-, a lisser sans doute, ou alors juste catégoriser comme mal définis
on evite les alternance pour une base

loi géometrique qui dicte selon qu'on reste plus ou moins longtemps dans un etat
caractérisées par l'esperance de X, qui vaut 1/P

## To do : 

function qui cree le model
function qui calcule specificité et sensibilité
function qui test la vraisemblance
function qui segmente les trois chromosomes (algo viterbi)
