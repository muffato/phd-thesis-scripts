#! /users/ldog/muffato/python -OO

__doc__ = """
"""

##################
# INITIALISATION #
##################

# Librairies
import sys
import utils.myGenomes
import utils.myTools
import utils.myDiags


#############
# FONCTIONS #
#############



########
# MAIN #
########

# Arguments
(noms_fichiers, options) = utils.myTools.checkArgs( ["studiedGenome", "referenceGenome", "orthologuesList"], [], __doc__)


# Chargement des fichiers
genome1 = utils.myGenomes.loadGenome(noms_fichiers["studiedGenome"])
genome2 = utils.myGenomes.loadGenome(noms_fichiers["referenceGenome"])
genesAnc = utils.myGenomes.loadGenome(noms_fichiers["orthologuesList"])

def printDiag( ((e1,c1,d1), (e2,c2,d2), da) ):
	global memeSens, seg1, seg2
	if (d1[-1]-d1[0])*(d2[-1]-d2[0]) > 0:
		memeSens += len(da)
		seg1 += 1
	else:
		seg2 += 1


chr = genome1.lstChr[:]

for c in chr:
	memeSens = 0
	seg1 = 0
	seg2 = 0
	genome1.lstChr = [c]
	genome2.lstChr = [c]
	utils.myDiags.calcDiags(genome1.nom, genome2.nom, genome1, genome2, genesAnc, printDiag, 1, -1, False, False)
	length = len(genome1.lstGenes[c])
	
	print "%d\t%d\t%.2f\t%d\t%d/%d" % (c,length,100.*max(memeSens,length-memeSens)/length,seg1+seg2-1,seg1,seg2)
