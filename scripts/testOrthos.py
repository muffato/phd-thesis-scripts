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
	["genesAnc"], \
	[("ancetre",str,""), ("output",str,""), ("fusionThreshold",int,-1), ("minimalLength",int,1), \
	("ancGenesFile",str,"/users/ldog/muffato/work/data/ancGenes/ancGenes.%s.list.bz2")], \
	__doc__ \
)


f = open(noms_fichiers["genesAnc"], 'r')

lst = f.readlines()[7:-1]
f.close()

for i in range(len(lst)):
	c = lst[i].split()
	for j in range(len(c)):
		if float(c[j]) > 0.5 and float(c[j]) < 1.5:
			print i,j+i+1
		


sys.exit(0)


genesAnc = utils.myGenomes.AncestralGenome(noms_fichiers["genesAnc"], False)

nb = 1
for l in sys.stdin:
	c = l.split()
	for i in c:
		print nb, " ".join(genesAnc.lstGenes[utils.myGenomes.AncestralGenome.defaultChr][int(i)].names)
	nb += 1

sys.exit(0)
phylTree = utils.myBioObjects.PhylogeneticTree(noms_fichiers["phylTree.conf"])
nbEsp = len(phylTree.getSpecies(phylTree.root))

for anc in phylTree.items:
	groupes = [phylTree.getSpecies(e) for (e,_) in phylTree.items[anc]]
	fils = utils.myMaths.flatten(groupes)

	nbO = nbEsp-len(fils)
	nbA = len(groupes[0])
	nbB = len(groupes[1])
	print anc, nbA*nbB + nbA*nbO + nbB*nbO
