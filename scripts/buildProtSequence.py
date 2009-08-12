#! /users/ldog/muffato/python

import sys
import itertools
import collections

import utils.myFile
import utils.myMaths
import utils.myTools
import utils.myGenomes
import utils.myPhylTree


arguments = utils.myTools.checkArgs( [("reconsSeqDNA",file)], [], "Renvoie la meilleure sequence ancestrale ADN compte tenu des acides amines")


lcodons = [''.join(codon) for codon in itertools.product("ACGT", "ACGT", "ACGT")]
#lcodons = [c for c in lcodons if utils.myGenomes.codon2aa[c] != '*']
laa = sorted(utils.myGenomes.aa2codon)
#laa.remove('*')

def prod(it):
	p = 1
	for x in it:
		p *= x
	return p


dicReconsSeq = {}

def proceed():
	res = dict((aa,"") for aa in laa)
	res[None] = ""
	length = len(dicReconsSeq['A'])
	assert length%3 == 0
	length /= 3
	nbGap = 0
	nbUnknown = 0
	for i in xrange(length):
		pc = collections.defaultdict(float)
		for codon in lcodons:
			pc[utils.myGenomes.codon2aa[codon]] += prod(dicReconsSeq[codon[j]][3*i+j] for j in xrange(3))
		for aa in laa:
			res[aa] += " %g" % pc[aa]
			#res[aa].append( "%g" % pc[aa] )
		lpos = sorted([(pc[aa],aa) for aa in laa], reverse = True)
		if lpos[0][0] <= 0:
			res[None] += '-'
			#res[None].append('-')
			nbGap += 1
		elif lpos[0][0] == lpos[1][0]:
			res[None] += 'X'
			#res[None].append('X')
			nbUnknown += 1
		else:
			res[None] += lpos[0][1]
			#res[None].append(lpos[0][1])
	print "LENGTH\t%d\t%d\t%d" % (length, length-nbGap, length-nbGap-nbUnknown)
	#print "SEQ\t%s" % (''.join(res[None]))
	print "SEQ\t%s" % res[None]
	for aa in laa:
		#print "PROBA-%s\t%s" % (aa, ' '.join(res[aa]))
		print "PROBA-%s\t%s" % (aa, res[aa][1:])
	print

f = utils.myFile.openFile(arguments["reconsSeqDNA"], "r")
for l in f:
	l = l[:-1]
	if l.startswith("NAME") or l.startswith("SPECIES"):
		if len(dicReconsSeq) != 0:
			proceed()
			dicReconsSeq = {}
		print l
	elif l.startswith("PROBA"):
		dicReconsSeq[l[6]] = [float(x) for x in l[8:].split()]
f.close()

proceed()

