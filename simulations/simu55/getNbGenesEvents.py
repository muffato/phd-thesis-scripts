#!/usr/bin/env python2

__doc__ = """
	Lit les arbres et renvoie le nombre d'evenements de chaque type par million d'annees
"""

import sys
import math
import random
import itertools
import collections

import utils.myMaths
import utils.myTools
import utils.myPhylTree


# Arguments
arguments = utils.myTools.checkArgs([("phylTree.conf",file), ("geneStats",str)], [], "A partir du detail de l'evolution des genes, renvoie les stats globales")

# L'arbre phylogenetique
phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])
sumBranchLength = sum(l for (_,l) in phylTree.parent.itervalues())
print "sumBranchLength", sumBranchLength

nbDup = 0
nbGain = 0
nbLoss = 0
dupNumber = 0

f = utils.myFile.openFile(arguments["geneStats"], "r")
f.readline()
for l in f:
	(lst,new) = eval(l)
	nbGain += len(new)
	nbLoss += sum(1 for x in lst.itervalues() if len(x) == 0)
	nbDup += sum(1 for x in lst.itervalues() if len(x) >= 2)
	dupNumber += sum(len(x) for x in lst.itervalues() if len(x) >= 2)
f.close()

print "nbGain", nbGain, float(nbGain)/sumBranchLength
print "nbLoss", nbLoss, float(nbLoss)/sumBranchLength
print "nbDup", nbDup, float(nbDup)/sumBranchLength
print "dupNumber", dupNumber, float(dupNumber)/nbDup

