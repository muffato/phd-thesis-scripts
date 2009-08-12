#! /users/ldog/muffato/python

__doc__ = """
	Lit les CNE et transforme les positions genomiques en positions indexees
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

arguments = utils.myTools.checkArgs( \
	[("list",file), ("cneFile",file)], \
	[], \
	__doc__ \
)


print >> sys.stderr, "Loading breakpoints ...",
sets = {}
f = utils.myFile.myTSV.reader(arguments["list"])
for t in f.csvobject:
	if t[0] not in sets:
		sets[t[0]] = {}
	if t[1] not in sets[t[0]]:
		sets[t[0]][t[1]] = set()
	sets[t[0]][t[1]].add( (t[-5],int(t[-4]),int(t[-3])) )
f.file.close()
print >> sys.stderr, "OK"

print >> sys.stderr, "Creating arrays ...",
res = {}
for anc in sets:
	res[anc] = {}
	for typ in sets[anc]:
		res[anc][typ] = (collections.defaultdict(list),{})
		n = 0
		for (c,x1,x2) in sets[anc][typ]:
			if x1 < x2:
				res[anc][typ][0][c].append( (x1,x2) )
				n += 1
		for x in res[anc][typ][0]:
			res[anc][typ][0][x].sort()
			res[anc][typ][1][x] = [0] * len(res[anc][typ][0][x])
		print >> sys.stderr, anc, typ, n
print >> sys.stderr, "OK"


def loadFile(name, fields, (chrom,x1,x2), title):
	print >> sys.stderr, "Loading %s ..." % title,
	for line in utils.myFile.myTSV.readTabular(name, fields):
		assert line[x1] <= line[x2]
		for anc in res:
			for typ in res[anc]:
				r = res[anc][typ]
				i0 = bisect.bisect(r[0][line[chrom]], (line[x1],line[x2]))
				i = i0 - 1
				while i >= 0:
					ref = r[0][line[chrom]][i]
					assert ref[0] <= line[x1]
					if ref[1] < line[x1]:
						break
					y2 = min(line[x2], ref[1])
					r[1][line[chrom]][i] += (y2-line[x1]+1)
					i -= 1
				i = i0
				while i < len(r[0][line[chrom]]):
					ref = r[0][line[chrom]][i]
					assert line[x1] <= ref[0]
					if ref[0] > line[x2]:
						break
					y2 = min(line[x2], ref[1])
					r[1][line[chrom]][i] += (y2-ref[0]+1)
					i += 1
	print >> sys.stderr, "OK"

	for anc in res:
		for typ in res[anc]:
			#print >> sys.stderr, anc, title, typ
			#print >> sys.stderr, "%s:" % typ, (100. * sum(res[anc][typ][1])) / sum(res[anc][typ][0]), sum(res[anc][typ][1]), sum(res[anc][typ][0])
			for c in res[anc][typ][0]:
				if c not in res[anc][typ][1]:
					assert len(res[anc][typ][0][c]) == 0
					continue
				for ((x1,x2),ncne) in zip(res[anc][typ][0][c], res[anc][typ][1][c]):
					print anc, title, typ, x2-x1+1, ncne
		#print >> sys.stderr



# Les CNEs
loadFile(arguments["cneFile"], [str,str,int,int,str,str,int,str], (1,2,3), "CNE")


