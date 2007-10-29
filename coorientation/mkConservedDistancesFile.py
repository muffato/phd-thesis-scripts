#! /users/ldog/muffato/python -OO

__doc__ = """
Lit un fichier de diagonales et extrait les genes synteniques d'une espece (pour un ancetre specifique)
Renvoie les distances inter-genes correspondantes
"""

import sys
import utils.myGenomes
import utils.myTools


(noms_fichiers,options) = utils.myTools.checkArgs(["diagsFile", "genesFile"], [("ancestr",str,""), ("species",str,"")], __doc__)

print >> sys.stderr, "Chargement des diagonales ...",
f = utils.myTools.myOpenFile(noms_fichiers["diagsFile"], "r")
allDiags = []
for ligne in f:
	if not ligne.startswith(options["ancestr"]):
		continue
	t = ligne.split('\t')
	if t[2] == options["species"]:
		allDiags.append(t[4])
	elif t[5] == options["species"]:
		allDiags.append(t[7])
f.close()
print >> sys.stderr, len(allDiags), "OK"

genome = utils.myGenomes.Genome(noms_fichiers["genesFile"])
allPairs = set()
print >> sys.stderr, "Extraction des paires ...",
for diag in allDiags:
	genes = diag.split()
	for (g1,g2) in utils.myTools.myIterator.slidingTuple(genes):
		(c1,i1) = genome.dicGenes[g1]
		(c2,i2) = genome.dicGenes[g2]
		# Meme chromosme, j'espere !
		if c1 != c2:
			print >> sys.stderr, g1, g2
			continue
		# Les deux genes
		g1 = genome.lstGenes[c1][i1]
		g2 = genome.lstGenes[c2][i2]
		# Ordonnes selon 5' -> 3'
		if i1 > i2:
			(g1,g2) = (g2,g1)
		# Eviter les doublons
		if (g1,g2) in allPairs:
			continue
		allPairs.add( (g1,g2) )
		# Genes imbriques ou chevaichants
		if g1.end >= g2.beginning:
			continue
		if g1.strand == g2.strand:
			print "U", g2.beginning - g1.end
		elif g1.strand == 1:
			print "C", g2.beginning - g1.end
		else:
			print "D", g2.beginning - g1.end
print >> sys.stderr, len(allPairs), "OK"


