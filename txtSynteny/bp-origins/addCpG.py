#!/usr/bin/env python2

__doc__ = """
	Renvoie le ration CpG observe/attendu pour tous les intervalles du chromosome charge
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
		(n1,cpg) = utils.myGenomes.getDiNuc(seq, x1, x2, ["CG"])
		if (n1 != 0) and (n != 0):
			print l + "\t%.2f" % ( (float(cpg)/n1) / (((float(gc)/n)/2)**2) )
		else:
			print l + "\tNone"
	else:
		print l #+ "\t" + "0\t0\tNone"
f.close()

