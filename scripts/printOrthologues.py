#! /users/ldog/muffato/python -OO

__doc__ = """
	Compare deux genomes et renvoie la liste des couples de genes orthologues avec les details sur leurs positions
"""

##################
# INITIALISATION #
##################

# Librairies
import sys
import utils.myGenomes
import utils.myTools



########
# MAIN #
########

# Arguments
(noms_fichiers, options) = utils.myTools.checkArgs( \
	["studiedGenome", "referenceGenome"], \
	[("orthologuesList",str,""), ("includeGaps",bool,False), ("includeScaffolds",bool,False), ("includeRandoms",bool,False), ("reverse",bool,False)], \
	__doc__
)


# Chargement des fichiers
genome1 = utils.myGenomes.Genome(noms_fichiers["studiedGenome"])
genome2 = utils.myGenomes.Genome(noms_fichiers["referenceGenome"])
if options["reverse"]:
	x = genome1
	genome1 = genome2
	genome2 = x
if options["orthologuesList"] != "":
	genesAnc = utils.myGenomes.Genome(options["orthologuesList"])
else:
	genesAnc = genome2

# Les chromosomes a etudier
chr1 = genome1.lstChr
chr2 = genome2.lstChr
if options["includeScaffolds"]:
	chr1.extend(genome1.lstScaff)
	chr2.extend(genome2.lstScaff)
if options["includeRandoms"]:
	chr1.extend(genome1.lstRand)
	chr2.extend(genome2.lstRand)

table12 = genome1.buildOrthosTable(chr1, genome2, chr2, options["includeGaps"], genesAnc)

# Fichier avec les noms des paires de genes orthologues et leurs coordonnees
for c1 in chr1:
	for (i1,t) in table12[c1]:
		g1 = genome1.lstGenes[c1][i1]
		for (c2,i2) in t:
			g2 = genome2.lstGenes[c2][i2]
			print utils.myTools.printLine( [c1,g1.beginning,g1.end,g1.strand,"/".join(g1.names), c2,g2.beginning,g2.end,g2.strand,"/".join(g2.names)] )
			#r = (c1, g1.beginning, g1.end, g1.strand, "/".join(g1.names), \
			#	c2, g2.beginning, g2.end, g2.strand, "/".join(g2.names))
			#print "\t".join([str(x) for x in r])

