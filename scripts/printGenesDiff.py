#!/usr/bin/env python2

__doc__ = """
	Compare deux genomes et renvoie les changements de familles de genes
"""

# Librairies
import sys
import utils.myTools
import utils.myGenomes


def getGeneTxt(g):
	return "/".join(g.names) + ":%s:%d-%d:%d" % (g.chromosome, g.beginning, g.end, g.strand)


def doAnalysis(genome1, genome2, txtDiff, txtDupDiff):
	for g2 in genome2:
		lg1 = genome1.getPosition(g2.names)
		if len(lg1) == 0:
			#print txtDiff, " ".join(g2.names)
			print txtDiff, getGeneTxt(g2)
			all.add(g2.names)
		elif len(lg1) > 1:
			print txtDupDiff, (" %s" % txtDupDiff).join(getGeneTxt(genome1.lstGenes[c1][i1]) for (c1,i1) in lg1)
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
		#print "=", " ".join(g2.names)
		print "=", getGeneTxt(g2)

