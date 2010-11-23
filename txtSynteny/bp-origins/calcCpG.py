#!/usr/bin/env python2
# Copyright (c) 2010 IBENS/Dyogen : Matthieu MUFFATO, Alexandra LOUIS, Hugues ROEST CROLLIUS
# Mail : hrc@ens.fr or alouis@biologie.ens.fr
# Licences GLP v3 and CeCILL v2

__doc__ = """
	Lit un genome, charge la sequence et calcule la densite en CpG
"""

import os
import sys
import bisect
import operator
import itertools
import collections

import utils.myFile
import utils.myMaths
import utils.myTools
import utils.myGenomes

arguments = utils.myTools.checkArgs([("genesList",file)], [("chromSequence",str,"/ressources/seq/HUMAN/NCBI36_HG18/dna/UNMASKED/Homo_sapiens.NCBI36.46.dna.chromosome.%s.fa")], __doc__)

genome = utils.myGenomes.Genome(arguments["genesList"])

for c in genome.lstGenes:
	try:
		print >> sys.stderr, "Loading DNA sequence for chromosome %s ..." % c,
		(_,seq) = utils.myGenomes.myFASTA.iterFile(arguments["chromSequence"] % c).next()
		print >> sys.stderr, len(seq), "OK"
		for g in genome.lstGenes[c]:
	
			if g.strand > 0:
				x1 = g.beginning - 1500
				x2 = g.beginning + 1500
			else:
				x1 = g.end - 1500
				x2 = g.end + 1500
			(n,gc) = utils.myGenomes.getMonoNuc(seq, x1, x2, "GC")
			(n1,cpg) = utils.myGenomes.getDiNuc(seq, x1, x2, ["CG"])
			print g.names[0], n, gc, n1, cpg,
			if (n1 != 0) and (n != 0):
				print (float(cpg)/n1) / (((float(gc)/n)/2)**2)
			else:
				print None
		
	except IOError:
		print >> sys.stderr, "NO"

