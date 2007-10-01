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
(noms_fichiers, options) = utils.myTools.checkArgs( ["referenceGenome","studiedGenome", "orthologuesList"], [], __doc__)


# Chargement des fichiers
genome1 = utils.myGenomes.Genome(noms_fichiers["referenceGenome"])
genome2 = utils.myGenomes.Genome(noms_fichiers["studiedGenome"])
genesAnc = utils.myGenomes.Genome(noms_fichiers["orthologuesList"])

chr = genome2.lstChr[:]
totMemeSens = 0.
totLength = 0.
totNbPaires = 0.
totPairesOK = 0.
totNbSegments = 0.

for c in chr:
	memeSens = 0
	seg1 = 0
	seg2 = 0
	#genome1.lstChr = [c]
	genome2.lstChr = [c]
	for ((c1,d1), (c2,d2), da) in utils.myDiags.calcDiags(genome1, genome2, genesAnc, 1, -1, False, False):
		if (d1[-1]-d1[0])*(d2[-1]-d2[0]) > 0:
			memeSens += len(d1)
			seg1 += 1
		else:
			seg2 += 1

	length = len(genome2.lstGenes[c])

	nbPerfectPairs = 0.
	for i in xrange(length-1):
		try:
			(c1,i1) = genome1.dicGenes[genome2.lstGenes[c][i].names[0]]
			(c2,i2) = genome1.dicGenes[genome2.lstGenes[c][i+1].names[0]]
			if abs(i1-i2) == 1:
				nbPerfectPairs += 1.
		except KeyError:
			continue

	memeSens = max(memeSens,length-memeSens)
	#print "%d\t%d\t%.4f\t%d\t%d/%d\t%.4f\t%.4f" % (c,length,100.*memeSens/length,seg1+seg2-1,seg1,seg2,float(length)/(seg1+seg2), 100.*nbPerfectPairs/(length-1))
	totMemeSens += memeSens
	totLength += length
	totNbPaires += length-1
	totPairesOK += nbPerfectPairs
	totNbSegments += seg1+seg2

print totMemeSens
print totLength
print totNbPaires
print totPairesOK
print totNbSegments
print "%.4f\t%.4f\t%.4f" % (100.*totMemeSens/totLength,100.*totPairesOK/totNbPaires,totLength/totNbSegments)

