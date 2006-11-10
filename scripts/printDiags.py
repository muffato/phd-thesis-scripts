#! /users/ldog/muffato/python

__doc__ = """
Extrait toutes les diagonales de genes entre deux especes
"""

##################
# INITIALISATION #
##################

# Librairies
import sys
import utils.myGenomes
import utils.myTools
import utils.myDiags


########
# MAIN #
########


# Arguments
(noms_fichiers, options) = utils.myTools.checkArgs( \
	["genome1", "genome2", "orthologuesList"], \
	[("onAncGenes",bool,False), ("fusionThreshold",int,-1), ("sameStrand",bool,True)], \
	__doc__ \
)


# 1. On lit tous les fichiers
genome1 = utils.myGenomes.loadGenome(noms_fichiers["genome1"])
genome2 = utils.myGenomes.loadGenome(noms_fichiers["genome2"])
genesAnc = utils.myGenomes.loadGenome(noms_fichiers["orthologuesList"])
lstGenesAnc = genesAnc.lstGenes[utils.myGenomes.AncestralGenome.defaultChr]

print >> sys.stderr, "Preparation ...",
transGenome1 = utils.myDiags.translateGenome(genome1, genesAnc)
locations = [[] for x in lstGenesAnc]
for ianc in xrange(len(lstGenesAnc)):
	for g in lstGenesAnc[ianc].names:
		if g not in genome2.dicGenes:
				continue
		(c,i) = genome2.dicGenes[g]
		locations[ianc].append( (c,i,genome2.lstGenes[c][i].strand) )

nbDiag = 0
nbDiag2 = 0
lenDiag = 0
lenDiag2 = 0

def mkStats(c1, c2, d1, d2):
	global nbDiag, nbDiag2, lenDiag, lenDiag2
	d = d1
	
	if options["onAncGenes"]:
		d = [transGenome1[c1][i] for i in d1]
	else:
		d = [genome1.lstGenes[c1][i].names[0] for i in d1]
	print " ".join([str(x) for x in d])
	
	nbDiag += 1
	lenDiag += len(d)
	if len(d) >= 2:
		nbDiag2 += 1
		lenDiag2 += len(d)

print >> sys.stderr, "Extraction des diagonales ...",
utils.myDiags.iterateDiags(transGenome1, locations, options["fusionThreshold"], options["sameStrand"], mkStats)
print >> sys.stderr, "OK"
print >> sys.stderr, nbDiag, nbDiag2, lenDiag, lenDiag2
print >> sys.stderr, float(lenDiag)/nbDiag
print >> sys.stderr, float(lenDiag2)/nbDiag2

