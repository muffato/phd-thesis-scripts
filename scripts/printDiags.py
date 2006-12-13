#! /users/ldog/muffato/python -OO

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



def translateGenome(genome):
	newGenome = {}
	for c in genome.lstChr:
		newGenome[c] = [(genesAnc.dicGenes.get(g.names[0], (0,-1))[1],g.strand) for g in genome.lstGenes[c]]
		if options["keepOnlyOrthos"]:
			newGenome[c] = [x for x in newGenome[c] if x[0] != -1]
	return newGenome


########
# MAIN #
########


# Arguments
(noms_fichiers, options) = utils.myTools.checkArgs( \
	["genome1", "genome2", "orthologuesList"], \
	[("fusionThreshold",int,-1), ("sameStrand",bool,True), ("keepOnlyOrthos",bool,False)], \
	__doc__ \
)


# 1. On lit tous les fichiers
genome1 = utils.myGenomes.loadGenome(noms_fichiers["genome1"])
genome2 = utils.myGenomes.loadGenome(noms_fichiers["genome2"])
genesAnc = utils.myGenomes.loadGenome(noms_fichiers["orthologuesList"])
lstGenesAnc = genesAnc.lstGenes[utils.myGenomes.Genome.defaultChr]

print >> sys.stderr, "Preparation ...",
transGenome1 = translateGenome(genome1)
transGenome2 = translateGenome(genome2)
del genesAnc
locations = [[] for x in lstGenesAnc]
for c in genome2.lstChr:
	for i in xrange(len(transGenome2[c])):
		(ianc,s) = transGenome2[c][i]
		if ianc != -1:
			locations[ianc].append( (c,i,s) )

nbDiag = 0
nbDiag2 = 0
lenDiag = 0
lenDiag2 = 0

def mkStats(c1, c2, d1, d2):
	global nbDiag, nbDiag2, lenDiag, lenDiag2
	
	da = [str(transGenome1[c1][i][0]) for i in d1]
	da1 = [genome1.lstGenes[c1][i].names[0] for i in d1]
	da2 = [genome2.lstGenes[c2][i].names[0] for i in d2]
	print "%d\t%s\t%s\t%s\t%s\t%s" % (len(da), c1, c2, " ".join(da), " ".join(da1), " ".join(da2))
	
	nbDiag += 1
	lenDiag += len(da)
	if len(da) >= 2:
		nbDiag2 += 1
		lenDiag2 += len(da)

print >> sys.stderr, "Extraction des diagonales ...",
utils.myDiags.iterateDiags(transGenome1, locations, options["fusionThreshold"], options["sameStrand"], mkStats)
print >> sys.stderr, "OK"

print >> sys.stderr, nbDiag, nbDiag2, lenDiag, lenDiag2
print >> sys.stderr, float(lenDiag)/nbDiag
print >> sys.stderr, float(lenDiag2)/nbDiag2

