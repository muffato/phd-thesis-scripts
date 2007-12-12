#! /users/ldog/muffato/python -OO

__doc__ = """
	Compare deux genomes et renvoie les changements de familles de genes
"""

# Librairies
import sys
import utils.myGenomes
import utils.myTools


def doAnalysis(genome1, genome2, txtDiff, txtDupDiff):
	for g2 in genome2:
		lg1 = genome1.getPosition(g2.names)
		if len(lg1) == 0:
			print txtDiff, " ".join(g2.names)
			all.add( tuple(sorted(g2.names)) )
		elif len(lg1) > 1:
			print txtDupDiff, "/".join([" ".join(genome1.lstGenes[c1][i1].names) for (c1,i1)  in lg1])
			all.update([ tuple(sorted(genome1.lstGenes[c1][i1].names)) for (c1,i1) in lg1])
			all.add( tuple(sorted(g2.names)) )

# Arguments
(noms_fichiers, options) = utils.myTools.checkArgs( ["studiedGenome", "referenceGenome"], [], __doc__)

# Chargement des fichiers
genome1 = utils.myGenomes.Genome(noms_fichiers["studiedGenome"])
genome2 = utils.myGenomes.Genome(noms_fichiers["referenceGenome"])

all = set()

doAnalysis(genome1, genome2, "+", "--")
doAnalysis(genome2, genome1, "-", "++")

for g2 in genome2:
	if tuple(sorted(g2.names)) not in all:
		print "=", " ".join(g2.names)

