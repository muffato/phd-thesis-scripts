#!/usr/bin/env python2
# Copyright (c) 2010 IBENS/Dyogen : Matthieu MUFFATO, Alexandra LOUIS, Hugues ROEST CROLLIUS
# Mail : hrc@ens.fr or alouis@biologie.ens.fr
# Licences GLP v3 and CeCILL v2

__doc__ = """
	Lit les CNE et transforme les positions genomiques en positions indexees
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

arguments = utils.myTools.checkArgs([("list",file), ("humanChr",file)], [], __doc__)
#str,"/ressources/seq/HUMAN/NCBI36_HG18/dna/UNMASKED/Homo_sapiens.NCBI36.46.dna.chromosome.%s.fa")], __doc__)

print >> sys.stderr, "Loading DNA sequence ...",
(chrom,seq) = utils.myGenomes.myFASTA.iterFile(arguments["humanChr"]).next()
chrom = chrom.split()[0]
print >> sys.stderr, chrom, len(seq), "OK"

f = utils.myFile.openFile(arguments["list"], "r")
for l in f:
	l = l[:-1]
	t = l.split("\t")
	if t[2] == chrom:
		x1 = int(t[3])
		x2 = int(t[4])
		(n,gc) = utils.myGenomes.getMonoNuc(seq, x1, x2, "GC")
		if n != 0:
			print l + "\t%d\t%d\t%.2f" % (n, gc, (100.*gc)/n)
		else:
			print l + "\t0\t0\tNone"
	else:
		print l
f.close()

