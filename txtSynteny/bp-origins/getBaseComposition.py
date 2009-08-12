#! /users/ldog/muffato/python

__doc__ = """
	Lit une liste de coordonnees et renvoie la composition pour chacune des positions
"""

import sys
import collections

import utils.myFile
import utils.myTools
import utils.myGenomes

arguments = utils.myTools.checkArgs([("list",file)], [("chromFile",str,""), ("addReverseComplement",bool,True), ("halfWindowSize",int,5000)], __doc__)

lst = collections.defaultdict(list)
f = utils.myFile.openFile(arguments["list"], "r")
for l in f:
	t = l.split()
	lst[t[0]].append( (int(t[1]),int(t[2])) )
f.close()
print >> sys.stderr, "Coords loaded"

def gen():
	return [0]*(2*arguments["halfWindowSize"]+1)
res = {"A": gen(), "C": gen(), "G": gen(), "T": gen()}


fw = {"A": "A", "T": "T", "C": "C", "G": "G"}
rev = {"A": "T", "T": "A", "C": "G", "G": "C"}

def addSegment(chrom, x1, x2, dir):
	x1 = max(x1, 0)
	x2 = min(x2, len(chrom)-1)
	(dic,it) = (fw,xrange(x1, x2+1)) if dir > 0 else (rev,xrange(x2, x1-1, -1))
	for (i,n) in enumerate(it):
		if chrom[n] in "ACGT":
			res[dic[chrom[n]]][i] += 1

for c in lst:
	try:
		chrom = utils.myGenomes.myFASTA.loadFile(arguments["chromFile"] % c).values()[0]
		print >> sys.stderr, "Chrom %s loaded" % c
	except IOError:
		print >> sys.stderr, "Chrom %s not found" % c
		continue
	for (pos,dir) in lst[c]:
		addSegment(chrom, pos-arguments["halfWindowSize"], pos+arguments["halfWindowSize"], dir)
		if arguments["addReverseComplement"]:
			addSegment(chrom, pos-arguments["halfWindowSize"], pos+arguments["halfWindowSize"], -dir)
	print >> sys.stderr, "Chrom %s processed" % c

for i in xrange(2*arguments["halfWindowSize"]+1):
	print i-arguments["halfWindowSize"], res["A"][i], res["C"][i], res["G"][i], res["T"][i]

