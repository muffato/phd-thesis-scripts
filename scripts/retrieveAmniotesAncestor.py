#! /usr/bin/python2.4


##################
# INITIALISATION #
##################

# Librairies
import sys
import math
import random
import os
import utils.myGenomes
import utils.myTools
import utils.myMaths


########
# MAIN #
########


# Arguments
(noms_fichiers, options) = utils.myTools.checkArgs( \
	["genesAncestraux"], \
	[], \
	"Reconstruit le genome de l'ancetre de 1 et 2 a partir de l'outgroup et des genes de cet ancetre" \
)

genesAnc = utils.myGenomes.loadGenome(noms_fichiers["genesAncestraux"])

# On lit toutes les donnees en filtrant celles qui nous permettent d'arriver a l'ancetre des mammiferes
# Il faut integrer les diagonales et les coller bout a bout

assoc = utils.myTools.myCombinator([])

nb = 0
for ligne in sys.stdin:
	c = ligne.split("\t")
	try:
		# On a un couple de genes
		i = int(c[0])
		j = int(c[1])
		OK = eval(c[2])
		NO = eval(c[3])
		if len(OK) >= 3: # and 'C' in OK and 'T' in OK:
			assoc.addLink([i,j])
	except Exception:
		# On a une diagonale
		if c[0] not in ["AO"]:
			d = eval(c[1])
			if len(d) >= 2:
				assoc.addLink(d)
	nb += 1
	if (nb % 500000) == 0:
		print >> sys.stderr, nb, len(assoc.getGrp())


res = assoc.getGrp()

c = 0
n = 0
for a in res:
	#if len(a) < 100:
	#	#print >> sys.stderr, len(a),
	#	continue
	n += len(a)
	c += 1
	print len(a)
	#print a
	continue
	for i in a:
		#print i
		print chr(96+c), " ".join(genesAnc.lstGenes[utils.myGenomes.AncestralGenome.defaultChr][i].names)

print >> sys.stderr
print >> sys.stderr, n, "genes repartis sur", c, "chromosomes"

