#! /users/ldog/muffato/python -OO

__doc__ = """
Affiche les chromosomes orthologues a partir d'un seuil en pourcentage de nombre de genes
"""

##################
# INITIALISATION #
##################

# Librairies
import string
import sys
import os
import utils.myGenomes
import utils.myTools
import utils.myMaths
import utils.myPsOutput


def getOrthologuesChromosomes(genome1, genome2, genesAnc):

	res = {}
	for c in genome1.lstChr:

		lst = []
		for gene in genome1.lstGenes[c]:
			
			tg = gene.names
			if genesAnc != None:
				tg = utils.myMaths.flatten([genesAnc.lstGenes[cc][ii].names for (cc,ii) in [genesAnc.dicGenes[g] for g in tg if g in genesAnc.dicGenes]])
			
			for g in tg:
				if g in genome2.dicGenes:
					col = genome2.dicGenes[g][0]
					lst.append(col)
					break

		count = [(lst.count(x),x) for x in set(lst)]
		count.sort()
		count.reverse()
		nb = (len(lst)*options["minHomology"])/100
		
		res[c] = []
		for (n,c2) in count:
			res[c].append( (c2,n) )
			nb -= n
			if nb <= 0:
				break
	return res


########
# MAIN #
########

# Arguments
(noms_fichiers, options) = utils.myTools.checkArgs( \
	["GenomeADessiner", "GenomeReference"], \
	[("minHomology",int,90), ("orthologuesList",str,"")], \
	__doc__ \
)

genome1 = utils.myGenomes.loadGenome(noms_fichiers["GenomeADessiner"])
genome2 = utils.myGenomes.loadGenome(noms_fichiers["GenomeReference"])
if options["orthologuesList"] != "":
	genesAnc = utils.myGenomes.AncestralGenome(options["orthologuesList"], False, False)
else:
	genesAnc = None


res = getOrthologuesChromosomes(genome1, genome2, genesAnc)

for c in res:
	print "%s\t%s" % (c,"   ".join(["%s (%d)" % x for x in res[c]]))

