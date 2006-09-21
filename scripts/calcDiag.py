#!/usr/bin/python2.4

# Ce script cherche les diagonales parmi les gens dans la matrice, puis les diagonales de diagonales ...
#    de maniere recursive et s'arrete quand aucun nouvel arrangement n'est cree
# Apres, il essaie de rassembler les regions assez proches, puis reitere l'operation, jusqu'a
#    arriver a un etat "stable"
# Entree : fichier d'orthologues a 12 champs

# INITIALISATION #

# Librairies
import sys
import os

sys.path.append(os.environ['HOME'] + "/work/scripts/utils")
import myTools
import myOrthos
import myMaths


# Arguments
(noms_fichiers, options) = myTools.checkArgs( \
	["genesList.conf", "genesAncestraux.list"],\
	[("trouMaxDiag", int, 0), ("buildByCase", bool, True), ("espece1", str, ""), ("espece2", str, ""), ("ancestralGenome", bool, False)], \
	"Calcule la longueur moyenne des diagonales" \
)


geneBank = myOrthos.GeneBank(noms_fichiers[0], [options["espece1"], options["espece2"]])

if len(geneBank.dicEspeces) == 0 or (len(geneBank.dicEspeces) == 1 and not options["ancestralGenome"]):
	print >> sys.stderr, "Can't retrieve -%(espece1)s- and -%(espece2)s- in the species list" % options
	sys.exit(1)

genesAnc = myOrthos.AncestralGenome(noms_fichiers[1], options["ancestralGenome"])
if options["espece1"] != "":
	genome1 = geneBank.dicEspeces[options["espece1"]]
else:
	genome1 = genesAnc

if options["espece2"] != "":
	genome2 = geneBank.dicEspeces[options["espece2"]]
else:
	genome2 = genesAnc

print >> sys.stderr, "Extraction des orthologues ",
nbDiag = 0.
nbDiag2 = 0.
lenDiag = 0.
lenDiag2 = 0.

orthos1 = dict( [(K1,[]) for K1 in genome1.lstChr] )
orthos2 = dict( [(K2,[]) for K2 in genome2.lstChr] )

ii = 0
for c in genesAnc.lstChr:
	for g in genesAnc.lstGenes[c]:
		g1 = [x for x in g.names if x in genome1.dicGenes]
		g2 = [x for x in g.names if x in genome2.dicGenes]
		for x in g1:
			(K1,i1) = genome1.dicGenes[x]
			for y in g2:
				(K2,i2) = genome2.dicGenes[y]
				orthos1[K1].append( (i1,K2,i2,ii) )
				orthos2[K2].append( (i2,K1,i1,ii) )
				ii += 1
	sys.stderr.write(".")

for K1 in genome1.lstChr:
	orthos1[K1].sort()
sys.stderr.write(".")
for K2 in genome2.lstChr:
	orthos2[K2].sort()
print >> sys.stderr, "OK"
	

print >> sys.stderr, "Extraction des diagonales ",
for K1 in genome1.lstChr:
	for K2 in genome2.lstChr:
		
		t1 = [x[-1] for x in orthos1[K1] if x[1] == K2]
		t2 = [x[-1] for x in orthos2[K2] if x[1] == K1]
		x = [len(d) for d in myMaths.extractDiags(t1, t2, options["buildByCase"], options["trouMaxDiag"])]
		for u in x:
			print u
		y = [l for l in x if l >= 2]
		nbDiag += len(x)
		nbDiag2 += len(y)
		lenDiag += sum(x)
		lenDiag2 += sum(y)
	sys.stderr.write(".")

print >> sys.stderr, " OK"
print >> sys.stderr, nbDiag, nbDiag2, lenDiag, lenDiag2
print >> sys.stderr, lenDiag/nbDiag
print >> sys.stderr, lenDiag2/nbDiag2

