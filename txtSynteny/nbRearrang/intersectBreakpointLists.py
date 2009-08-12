#! /users/ldog/muffato/python

__doc__ = """
	Parcourt les deux genomes et extrait les intervalles conserves, les points de cassure
	  en tenant compte des evenements de genes (gain/perte/duplication)
"""

import sys
import collections
import itertools

import utils.myFile
import utils.myDiags
import utils.myMaths
import utils.myTools
import utils.myGenomes
import utils.myPhylTree


# Argument:
modes = ["diags","interv","sdnds","intcons"]
arguments = utils.myTools.checkArgs( [("humanGenome",file), ("list1",file), ("list2",file)], [("type1",str,modes),("type2",str,modes)], __doc__)

genome = utils.myGenomes.Genome(arguments["humanGenome"])

def loadRearrangFromDiags(s):
	f = utils.myFile.openFile(s, "r")
	lstBP = []
	for l in f:
		if not l.startswith("GAIN"):
			continue
		t = l.split()
		tmp = []
		for (g1,g2) in itertools.product([t[9],t[11]], [t[17],t[19]]):
			(c1,i1) = genome.dicGenes[g1]
			(c2,i2) = genome.dicGenes[g2]
			assert c1 == c2
			if i1 < i2:
				tmp.append( (abs(i1-i2),g1,g2) )
			else:
				tmp.append( (abs(i1-i2),g2,g1) )
		lstBP.append( min(tmp)[1:] )
		#print "INI-diags", lstBP[-1]
	f.close()
	return lstBP

def loadRearrangFromIntervals(s):
	f = utils.myFile.openFile(s, "r")
	lstBP = []
	for l in f:
		if l.startswith("-"):
			continue
		t = l.split()
		if (t[6] == "SAME_ORIENT") or (t[6] == "DUPLICATES") or (t[6] == "TWO_ENDS"):
			continue
		lstBP.append( (t[1],t[2]) )
		#print "INI-interv", lstBP[-1]
	f.close()
	return lstBP

def loadRearrangFromBP(s, i1, i2):
	f = utils.myFile.openFile(s, "r")
	lstBP = []
	for l in f:
		if l.startswith("LOST"):
			continue
		t = l.split()
		#lstBP.append( (t[i1],t[i2]) )
		(c1,j1) = genome.dicGenes[t[i1]]
		(c2,j2) = genome.dicGenes[t[i2]]
		assert c1 == c2
		if j1 < j2:
			lstBP.append( (t[i1],t[i2]) )
		else:
			lstBP.append( (t[i2],t[i1]) )

		#print "INI-bp", lstBP[-1]
	f.close()
	return lstBP

if arguments["type1"] == modes[0]:
	bp1 = loadRearrangFromDiags(arguments["list1"])
elif arguments["type1"] == modes[1]:
	bp1 = loadRearrangFromIntervals(arguments["list1"])
elif arguments["type1"] == modes[2]:
	bp1 = loadRearrangFromBP(arguments["list1"], 0, 1)
else:
	bp1 = loadRearrangFromBP(arguments["list1"], 2, 3)
print >> sys.stderr, len(bp1)
if arguments["type2"] == modes[0]:
	bp2 = loadRearrangFromDiags(arguments["list2"])
elif arguments["type2"] == modes[1]:
	bp2 = loadRearrangFromIntervals(arguments["list2"])
elif arguments["type2"] == modes[2]:
	bp2 = loadRearrangFromBP(arguments["list2"], 0, 1)
else:
	bp2 = loadRearrangFromBP(arguments["list2"], 2, 3)
print >> sys.stderr, len(bp2)

def expand(lst):
	newlst = collections.defaultdict(list)
	for (g1,g2) in lst:
		(c1,i1) = genome.dicGenes[g1]
		(c2,i2) = genome.dicGenes[g2]
		assert c1 == c2
		assert i1 < i2
		newlst[c1].append( (i1,i2,[(i,i+1) for i in xrange(i1,i2)]) )
		#print "TRANSFORM", g1, g2, i1, i2
	#for c in newlst:
	#	newlst[c].sort()
	return newlst

inside1 = expand(bp1)
print >> sys.stderr, len(inside1)
inside2 = expand(bp2)
print >> sys.stderr, len(inside2)

def txt(c, l):
	return "+".join("%d>%d(%s>%s)" % (g1, g2, genome.lstGenes[c][g1].names[0], genome.lstGenes[c][g2].names[0]) for (g1,g2) in l)

for c in set(inside1).intersection(inside2):

	grp = utils.myTools.myCombinator()
	for (_,_,x1) in inside1[c]:
		grp.addLink(x1)
	for (_,_,x2) in inside2[c]:
		grp.addLink(x2)
	
	d1 = collections.defaultdict(list)
	for (i1a,i1b,x1) in inside1[c]:
		d1[grp.dic[x1[0]]].append( (i1a,i1b) )
	d2 = collections.defaultdict(list)
	for (i2a,i2b,x2) in inside2[c]:
		d2[grp.dic[x2[0]]].append( (i2a,i2b) )
	
	for i in set(d1).intersection(d2):
		print "BOTH", c, txt(c, d1[i]), txt(c, d2[i])
	
	for i in set(d1).difference(d2):
		print "ONLY1", c, txt(c, d1[i])
	
	for i in set(d2).difference(d1):
		print "ONLY2", c, txt(c, d2[i])

for c in set(inside1).difference(inside2):
	for x in inside1[c]:
		print "ONLY1", c, txt(c, [x[:2]])

for c in set(inside2).difference(inside1):
	for x in inside2[c]:
		print "ONLY2", c, txt(c, [x[:2]])

