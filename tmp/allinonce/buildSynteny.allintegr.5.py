#!/usr/bin/env python2

__doc__ = """
	A partir de diagonales pair-wise, construit des versions integrees simultanement pour chaque ancetre
"""


import sys
import itertools
import collections

import utils.myPhylTree
import utils.myGenomes
import utils.myFile
import utils.myTools
import utils.myMaths

import utils.myDiags


# Arguments
arguments = utils.myTools.checkArgs( \
	[("phylTree.conf",file), ("diagGroups",file) ], \
	[("OUT.ancDiags",str,"anc/diags.%s.list.bz2"), \
	("ancGenesFile",str,"~/work/data/ancGenes/ancGenes.%s.list.bz2")], \
	__doc__ \
)

class ancContigs:

	def __init__(self):
		self.succ = {}
		self.pred = {}
		self.allsucc = collections.defaultdict(set)
		self.allpred = collections.defaultdict(set)
		self.age = {}

	def isCompatible(self, (g1,s1), (g2,s2), age):
		if (g1,s1) in self.succ:
			if (g2,s2) in self.pred:
				return False
			if (g2,s2) in self.succ:
				return False
			if self.age[(g1,s1)] <= age:
				return False
		if (g2,s2) in self.allpred:
			if (g1,s1) in self.pred:
				return False
			if (g1,s1) in self.succ:
				return False
			if self.age[(g2,-s2)] <= age:
				return False
		return True

	def removeFollowing(self, (g1,s1)):
		(g2,s2) = self.succ[(g1,s1)]
		allsucc = self.allsucc[(g1,s1)].copy()
		allpred = self.allpred[(g2,s2)].copy()
		for x in allpred:
			self.allsucc[x].difference_update(allsucc)
		for x in allsucc:
			self.allpred[x].difference_update(allpred)
		assert len(self.allsucc[(g1,s1)]) == 0
		assert len(self.allpred[(g2,s2)]) == 0
		del self.succ[(g1,s1)]
		del self.pred[(g2,s2)]
		del self.succ[(g2,-s2)]
		del self.pred[(g1,-s1)]
		del self.allsucc[(g1,s1)]
		del self.allpred[(g2,s2)]
		del self.allsucc[(g2,-s2)]
		del self.allpred[(g1,-s1)]
		del self.age[(g1,s1)]
		del self.age[(g2,-s2)]

	def addLink(self, (g1,s1), (g2,s2), age):
		assert (g1,s1) not in self.succ
		assert (g2,s2) not in self.pred

		allsucc = set([(g2,s2)])
		if (g2,s2) in self.allsucc:
			allsucc.update(self.allsucc[(g2,s2)])
		self.succ[(g1,s1)] = (g2,s2)

		allpred = set([(g1,s1)])
		if (g1,s1) in self.allpred:
			allpred.update(self.allpred[(g1,s1)])
		self.pred[(g2,s2)] = (g1,s1)

		for x in allpred:
			self.allsucc[x].update(allsucc)
		for y in allsucc:
			self.allpred[y].update(allpred)

		self.age[(g1,s1)] = age

	def extractContigs(self):
		extr = set(self.succ).difference(self.pred)
		seen = set()

		for x in extr:
			if x in seen:
				continue
			l = [x]
			while x in self.succ:
				x = self.succ[x]
				l.append(x)
			seen.add((x[0],-x[1]))
			yield l

phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])

all = []
f = utils.myFile.openFile(arguments["diagGroups"], "r")
for l in f:
	lst = [x for x in eval(l) if x[0] in phylTree.listAncestr]
	dic2 = collections.defaultdict(int)
	for (a,g1,g2,s1,s2) in lst:
		dic2[a] += 1
	if (len(dic2) > 0) and (max(dic2.values()) == 1):
		all.append( (len(lst), lst) )
f.close()

print >> sys.stderr, len(all)


all.sort(reverse=True)

cont = collections.defaultdict(ancContigs)

for (s,l) in all:
	age = max(phylTree.ages[x[0]] for x in l)
	for (anc,g1,g2,s1,s2) in l:
		if cont[anc].isCompatible((g1,s1), (g2,s2), age) != cont[anc].isCompatible((g2,-s2), (g1,-s1), age):
			print (anc,g1,g2,s1,s2)
			print cont[anc].succ.get((g1,s1)), cont[anc].age.get((g1,s1))
			print cont[anc].pred.get((g2,s2)), cont[anc].age.get((g2,-s2))
			print cont[anc].succ.get((g2,-s2)), cont[anc].pred.get((g1,-s1))
			print cont[anc].allsucc.get((g1,s1)), cont[anc].allpred.get((g2,s2))
			print cont[anc].allsucc.get((g2,-s2)), cont[anc].allpred.get((g1,-s1))
			print
			break
		if not cont[anc].isCompatible((g1,s1), (g2,s2), age):
			print "notwanted", s, anc, (g1,s1), (g2,s2)
			break
	else:
		for (anc,g1,g2,s1,s2) in l:
			if (g1,s1) in cont[anc].succ:
				(g3,s3) = cont[anc].succ[(g1,s1)]
				cont[anc].removeFollowing((g1,s1))
				cont[anc].addLink((g2,s2), (g3,s3), age)
				cont[anc].addLink((g3,-s3), (g2,-s2), age)
			elif (g2,s2) in cont[anc].pred:
				(g3,s3) = cont[anc].pred[(g2,s2)]
				cont[anc].removeFollowing((g3,s3))
				cont[anc].addLink((g3,s3), (g1,s1), age)
				cont[anc].addLink((g1,-s1), (g3,-s3), age)
			cont[anc].addLink((g1,s1), (g2,s2), age)
			cont[anc].addLink((g2,-s2), (g1,-s1), age)

print >> sys.stderr, sum(len(x.pred) for x in cont.itervalues()), utils.myMaths.myStats.txtSummary([len(x.pred) for x in cont.itervalues()])

for anc in cont:
	#notseen = set(xrange(len(utils.myGenomes.Genome(arguments["ancGenesFile"] % anc).lstGenes[None])))
	notseen = set()
	f = utils.myFile.openFile(arguments["OUT.ancDiags"] % anc, "w")
	for x in cont[anc].extractContigs():
		notseen.difference_update(g for (g,s) in x)
		print >> f, utils.myFile.myTSV.printLine([anc, len(x), utils.myFile.myTSV.printLine([g for (g,s) in x], delim=" "), utils.myFile.myTSV.printLine([s for (g,s) in x], delim=" ")])
	for g in notseen:
		print >> f, utils.myFile.myTSV.printLine([anc, 1, g, 0])
	print >> sys.stderr, anc, len(cont[anc].succ)/2, utils.myMaths.myStats.txtSummary([len(x) for x in cont[anc].extractContigs()])
	f.close()


