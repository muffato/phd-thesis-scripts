#!/usr/bin/env python2
# Copyright (c) 2010 IBENS/Dyogen : Matthieu MUFFATO, Alexandra LOUIS, Hugues ROEST CROLLIUS
# Mail : hrc@ens.fr or alouis@biologie.ens.fr
# Licences GLP v3 and CeCILL v2

__doc__ = """
	Renvoie specificite/sensibilite en chromosomes/genes
"""

import sys
import utils.myGenomes
import utils.myTools
import utils.myDiags
import utils.myPsOutput


# Fabrique la liste des orthologues entre les deux genomes
def buildOrthosTable(genome1, chr1, genome2, chr2):

	# Tous les orthologues entre pour les chromosomes OK du genome 1
	res = {}
	for c1 in chr1:
		res[c1] = {}
		for g1 in genome1.lstGenes[c1]:
			tmp = set()
			for s in g1.names:
				if s in genome2.dicGenes:
					tmp.add(genome2.dicGenes[s])
			# On ne garde que les chromsomes OK du genome 2
			tmp = [(c,i) for (c,i) in tmp if c in chr2]
			res[c1][genome1.dicGenes[s][1]] = [(c,i) for (c,i) in tmp if c in chr2]
			
	nbPairesOK = 0.
	nbPairesTotal = 0.
	nbGenesOK = 0.
	nbGenesTotal = 0.

	for c1 in chr1:
		score = dict.fromkeys(chr2, 0)
		for i1 in res[c1]:
			for (c2,i2) in res[c1][i1]:
				score[c2] += 1

		for c2 in chr2:
			nbPairesOK += (score[c2]*(score[c2]-1))/2.
		nbPairesTotal += (len(res[c1])*(len(res[c1])-1))/2.
		
		nbGenesOK += max(score.values())
		nbGenesTotal += len(res[c1])

	if nbPairesTotal != 0:
		print "%.4f\t%.4f" % (100.*nbPairesOK/nbPairesTotal, 100.*nbGenesOK/nbGenesTotal),
	else:
		print "%.4f\t%.4f" % (0., 0.),


# Arguments
arguments = utils.myTools.checkArgs( [("referenceGenome",file)], [], __doc__)

genome2 = utils.myGenomes.Genome(arguments["referenceGenome"])

for l in sys.stdin:
	l = l.replace('\n', '')
	# Chargement des fichiers
	genome1 = utils.myGenomes.Genome(l)
	print l,
	buildOrthosTable(genome1, genome1.lstChr, genome2, genome2.lstChr)
	buildOrthosTable(genome2, genome2.lstChr, genome1, genome1.lstChr)
	print

