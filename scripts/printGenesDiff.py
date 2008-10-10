#! /users/ldog/muffato/python

__doc__ = """
	Compare deux genomes et renvoie les changements de familles de genes
"""

# Librairies
import sys
import utils.myTools
import utils.myGenomes


def doAnalysis(genome1, genome2, txtDiff, txtDupDiff):
	for g2 in genome2:
		lg1 = genome1.getPosition(g2.names)
		if len(lg1) == 0:
			print txtDiff, " ".join(g2.names)
			all.add(g2.names)
		elif len(lg1) > 1:
			print txtDupDiff, "/".join([" ".join(genome1.lstGenes[c1][i1].names) for (c1,i1) in lg1])
			all.update([genome1.lstGenes[c1][i1].names for (c1,i1) in lg1])
			all.add(g2.names)

# Arguments
arguments = utils.myTools.checkArgs( [("studiedGenome",file), ("referenceGenome",file)], [], __doc__)

# Chargement des fichiers
genome1 = utils.myGenomes.Genome(arguments["studiedGenome"])
genome2 = utils.myGenomes.Genome(arguments["referenceGenome"])

all = set()

doAnalysis(genome1, genome2, "+", "--")
doAnalysis(genome2, genome1, "-", "++")

for g2 in genome2:
	if g2.names not in all:
		print "=", " ".join(g2.names)

