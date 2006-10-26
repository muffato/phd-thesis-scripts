#! /usr/bin/python2.4

__doc__ = """
Extrait toutes les diagonales de genes entre deux especes
"""

##################
# INITIALISATION #
##################

# Librairies
import os
import sys
import utils.myGenomes
import utils.myTools
import utils.myDiags
import cPickle

########
# MAIN #
########

# Arguments
(noms_fichiers, options) = utils.myTools.checkArgs( \
	["ancGenesFile"], \
	[("ancetre",str,""), ("output",str,""), ("fusionThreshold",int,-1), ("minimalLength",int,1), \
	("ancGenesFile",str,"/users/ldog/muffato/work/data/ancGenes/ancGenes.%s.list.bz2")], \
	__doc__ \
)

genesAnc = utils.myGenomes.AncestralGenome(noms_fichiers["ancGenesFile"], False)

nb = 1
for l in sys.stdin:
	c = l.split()
	for i in c:
		print nb, " ".join(genesAnc.lstGenes[utils.myGenomes.AncestralGenome.defaultChr][int(i)].names)
	nb += 1


