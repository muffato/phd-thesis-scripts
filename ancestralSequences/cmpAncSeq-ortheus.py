#! /users/ldog/muffato/python

import sys
import itertools

import utils.myFile
import utils.myMaths
import utils.myTools
import utils.myGenomes
import utils.myPhylTree

arguments = utils.myTools.checkArgs( [("ExpectedAligment",file), ("ComputedAligment",file), ("tree",file)], [], "Compare la sequence reconstruite a la sequence attendue")

tree = utils.myPhylTree.PhylogeneticTree(arguments["tree"])

shiftID = len(tree.listSpecies) + 1
for i in xrange(len(tree.listAncestr)):
	print >> sys.stderr, i, i + shiftID, tree.species["NAME_%d" % i]

dicRealSeq = {}
f = utils.myFile.openFile(arguments["ExpectedAligment"], "r")
for l in f:
	if l.startswith("node"):
		#dicRealSeq[int(l[4:10])] = l[10:-1].replace(" ", "")
		dicRealSeq[int(l[4:10])] = ''.join(utils.myGenomes.codon2aa[x] for x in l[10:-1].split())
f.close()
print >> sys.stderr, dicRealSeq.keys()

tmp = utils.myGenomes.myFASTA.loadFile(arguments["ComputedAligment"])
print >> sys.stderr, tmp.keys()
tmp2 = dict([(tree.lastCommonAncestor(x.split("_")),s.replace("-","")) for (x,s) in tmp.iteritems()])
dicReconsSeq = dict([(x,''.join(utils.myGenomes.codon2aa[s[3*i:3*i+3]] for i in xrange(len(s)/3))) for (x,s) in tmp2.iteritems()])
print >> sys.stderr, dicReconsSeq.keys()

for anc in xrange(len(tree.listAncestr)):
	anc1 = "NAME_%d" % anc
	anc2 = anc + shiftID
	assert anc1 in dicReconsSeq, (anc,anc1)
	assert anc2 in dicRealSeq, (anc,anc2)
	assert "-" not in dicReconsSeq[anc1], dicReconsSeq[anc1]
	assert len(dicRealSeq[anc2]) == len(dicReconsSeq[anc1]), anc
	for (i,(expected,recons)) in enumerate(itertools.izip(dicRealSeq[anc2], dicReconsSeq[anc1])):
		#if (i % 3) == 2:
		#	continue

		print anc2, expected,recons

