#! /users/ldog/muffato/python

import sys
import itertools
import collections

import utils.myFile
import utils.myMaths
import utils.myTools
import utils.myGenomes
import utils.myPhylTree

arguments = utils.myTools.checkArgs( [("fastaSeq",file)], [("startPos",int,0), ("endPos",int,0)], "Compte les codons")

count = collections.defaultdict(int)

for (_,l) in utils.myGenomes.iterFile(arguments["fastaSeq"])
	if len(l)%3 != 0:
		continue
	assert (len(l)%3) == 0
	for i in xrange(arguments["startPos"], len(l)/3-arguments["endPos"]):
		count[l[3*i:3*i+3]] += 1

nb = float(sum(count.values()))/100

for codon in itertools.product("ACGT", "ACGT", "ACGT"):
	codon = ''.join(codon)
	print utils.myFile.myTSV.printLine([codon, utils.myGenomes.codon2aa[codon], count[codon], "%.4f" % (count[codon]/nb)])
