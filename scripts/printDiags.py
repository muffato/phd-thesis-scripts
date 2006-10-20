#! /usr/bin/python2.4

__doc__ = """
Extrait toutes les diagonales de genes entre deux especes
"""

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
	["genome1", "genome2", "orthologuesList"], \
	[("keepGaps",bool,True), ("fusionThreshold",int,-1)], \
	__doc__ \
)


# 1. On lit tous les fichiers
genesAnc = utils.myGenomes.AncestralGenome(noms_fichiers[2], False)
genome1 = utils.myMaths.translateGenome(utils.myGenomes.loadGenome(noms_fichiers[0]), "genome1", genesAnc)
genome2 = utils.myMaths.translateGenome(utils.myGenomes.loadGenome(noms_fichiers[1]), "genome2", genesAnc)

nbDiag = 0
nbDiag2 = 0
lenDiag = 0
lenDiag2 = 0

def mkStats(c1, c2, nbTot, d):
	global nbDiag, nbDiag2, lenDiag, lenDiag2
	print " ".join([str(x) for x in d])
	nbDiag += 1
	lenDiag += len(d)
	if len(d) >= 2:
		nbDiag2 += 1
		lenDiag2 += len(d)
	#for g in d:
	#	print 1, " ".join(genesAnc.lstGenes[utils.myGenomes.AncestralGenome.defaultChr][g].names)

print >> sys.stderr, "Extraction des diagonales ...",
utils.myMaths.iterateDiags(genome1, genome2, options["keepGaps"], options["fusionThreshold"], mkStats)
print >> sys.stderr, "OK"
print >> sys.stderr, nbDiag, nbDiag2, lenDiag, lenDiag2
print >> sys.stderr, float(lenDiag)/nbDiag
print >> sys.stderr, float(lenDiag2)/nbDiag2

