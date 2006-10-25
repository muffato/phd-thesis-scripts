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


########
# MAIN #
########

# Arguments
(noms_fichiers, options) = utils.myTools.checkArgs( \
	["genesList.conf"], \
	[("ancetre",str,""), ("output",str,""), ("fusionThreshold",int,-1), ("minimalLength",int,1), \
	("ancGenesFile",str,"/users/ldog/muffato/work/data/ancGenes/ancGenes.%s.list.bz2")], \
	__doc__ \
)

#genesAnc = utils.myGenomes.AncestralGenome(noms_fichiers["ancGenesFile"], False)
geneBank = utils.myGenomes.GeneBank(noms_fichiers["genesList.conf"])

for l in sys.stdin:
	c = l.split()
	for i in range(len(c)):
		for j in range(i):
			(e1,_,_) = geneBank.dicGenes[c[i]]
			(e2,_,_) = geneBank.dicGenes[c[j]]
	os.system("bzgrep %s /users/ldog/muffato/work/data/orthologs/orthos.%s.%s.list.bz2" % (c[i],e1,e2))


