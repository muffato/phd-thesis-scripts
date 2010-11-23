#!/usr/bin/env python2
# Copyright (c) 2010 IBENS/Dyogen : Matthieu MUFFATO, Alexandra LOUIS, Hugues ROEST CROLLIUS
# Mail : hrc@ens.fr or alouis@biologie.ens.fr
# Licences GLP v3 and CeCILL v2

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
		dicRealSeq[int(l[4:10])] = l[10:-1].replace(" ", "")
		#dicRealSeq[int(l[4:10])] = ''.join(utils.myGenomes.codon2aa[x] for x in l[10:-1].split())
f.close()
print >> sys.stderr, dicRealSeq.keys()

name = None
dicReconsSeq = {}
f = utils.myFile.openFile(arguments["ComputedAligment"], "r")
for l in f:
	l = l[:-1]
	if l.startswith("NAME"):
		name = l[5:]
		dicReconsSeq[name] = {}
	elif l.startswith("PROBA"):
		#print >> sys.stderr, "!%s!" % l[6], name
		if l[6] != '*':
			#print >> sys.stderr, "OK"
			dicReconsSeq[name][l[6]] = [float(x) for x in l[8:].split()]
f.close()

print >> sys.stderr, dicReconsSeq.keys()

lN = "ACGT"
#lN = "ACDEFGHIKLMNPQRSTVWY"
idN = dict((x,i) for (i,x) in enumerate(lN))

for anc in xrange(len(tree.listAncestr)):
	anc1 = "NAME_%d" % anc
	anc2 = anc + shiftID
	assert anc1 in dicReconsSeq, (anc,anc1)
	assert anc2 in dicRealSeq, (anc,anc2)
	#print >> sys.stderr, anc, anc1, anc2, dicReconsSeq[anc1].keys()
	assert set(dicReconsSeq[anc1]) == set(lN), (anc1,set(dicReconsSeq[anc1]),set(lN))
	assert set(len(dicReconsSeq[anc1][x]) for x in lN) == set([len(dicRealSeq[anc2])]), ([len(dicReconsSeq[anc1][x]) for x in lN], len(dicRealSeq[anc2]))
	for (i,(expected,proba)) in enumerate(zip(dicRealSeq[anc2], zip(*(dicReconsSeq[anc1][x] for x in lN)))):
		#if (i % 3) != 2:
		#	continue
		ipredicted = max((p,i) for (i,p) in enumerate(proba))[1]
		print anc2, expected, proba[idN[expected]], lN[ipredicted], proba[ipredicted]

