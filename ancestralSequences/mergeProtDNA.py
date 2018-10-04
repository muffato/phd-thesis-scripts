#!/usr/bin/env python2

import sys
import itertools
import collections

import utils.myFile
import utils.myMaths
import utils.myTools
import utils.myGenomes
import utils.myPhylTree


arguments = utils.myTools.checkArgs( [("reconsSeqProt",file), ("reconsSeqDNA",file), ("simulSeqDNA",file)], [], "Renvoie la meilleure sequence ancestrale ADN compte tenu des acides amines")


def loadSeq(filename):
	name = None
	dicReconsSeq = {}
	f = utils.myFile.openFile(filename, "r")
	for l in f:
		l = l[:-1]
		if l.startswith("NAME"):
			name = l[5:]
			dicReconsSeq[name] = {}
		elif l.startswith("PROBA"):
			dicReconsSeq[name][l[6]] = [float(x) for x in l[8:].split()]
	f.close()
	return dicReconsSeq

seqDNA = loadSeq(arguments["reconsSeqDNA"])
seqProt = loadSeq(arguments["reconsSeqProt"])

dicRealSeq = {}
f = utils.myFile.openFile(arguments["simulSeqDNA"], "r")
for l in f:
	if l.startswith("node"):
		#dicRealSeq[int(l[4:10])] = l[10:-1].replace(" ", "")
		dicRealSeq[int(l[4:10])] = ''.join(utils.myGenomes.codon2aa[x] for x in l[10:-1].split())
f.close()
print >> sys.stderr, dicRealSeq.keys()


lcodons = [''.join(codon) for codon in itertools.product("ACGT", "ACGT", "ACGT")]

def prod(it):
	p = 1
	for x in it:
		p *= x
	return p

for anc in seqProt:
	assert anc in seqDNA, anc
	refAnc = int(anc[5:]) + 22
	assert refAnc in dicRealSeq
	for i in xrange(len(seqProt[anc].values()[0])):
		bestProt = max((seqProt[anc][p][i],p) for p in seqProt[anc])[1]
		pc = collections.defaultdict(float)
		for codon in lcodons:
			pc[utils.myGenomes.codon2aa[codon]] += prod(seqDNA[anc][codon[j]][3*i+j] for j in xrange(3))
		bestCodon = max((pc[p],p) for p in seqProt[anc])[1]
		tmpDNA = ''.join(max((seqDNA[anc][n][3*i+j],n) for n in "ACGT")[1] for j in xrange(3))
		bestDNA = utils.myGenomes.codon2aa[tmpDNA]
		print refAnc, dicRealSeq[refAnc][i], bestDNA, bestProt, bestCodon, prod(seqDNA[anc][tmpDNA[j]][3*i+j] for j in xrange(3)), "%.5g" % seqProt[anc][bestProt][i], "%.5g" % pc[bestCodon]

