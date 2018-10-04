#!/usr/bin/env python2

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

arguments = utils.myTools.checkArgs( [("list",file), ("trackFile",file), ("chrom",int), ("x1",int),("x2",int)], [], __doc__)


print >> sys.stderr, "Loading track ...",
track = {}
f = utils.myFile.myTSV.reader(arguments["trackFile"])
for t in f.csvobject:
	if t[arguments["chrom"]] not in track:
		track[t[arguments["chrom"]]] = []
	track[t[arguments["chrom"]]].append( (int(t[arguments["x1"]]),int(t[arguments["x2"]])) )
f.file.close()

for x in track.itervalues():
	x.sort()
print >> sys.stderr, "OK"

used = {}
for c in track:
	used[c] = [0] * len(track[c])
	#print >> sys.stderr, c, len(track[c])

def getDensity(chrom, x1, x2):
	if (x1 > x2) or (chrom not in track):
		return "0\t0.00"
	nTrack = 0

	#print >> sys.stderr, "new:", chrom, x1, x2
	i0 = bisect.bisect(track[chrom], (x1,x2))
	i = i0 - 1
	while i >= 0:
		ref = track[chrom][i]
		assert ref[0] <= x1
		if ref[1] < x1:
			break
		y2 = min(x2, ref[1])
		#print >> sys.stderr, "%s/%d: adding %d in %d %s" % (chrom,i, y2-x1+1, ref[1]-ref[0]+1, ref)
		nTrack += (y2-x1+1)
		#used[chrom][i] += (y2-x1+1)
		#if used[chrom][i] > (track[chrom][i][1] - track[chrom][i][0] + 1):
			#print >> sys.stderr, "Depassement en %s/%d : %d sur %d" % (chrom,i, used[chrom][i], (track[chrom][i][1] - track[chrom][i][0] + 1))
		i -= 1
	i = i0
	while i < len(track[chrom]):
		ref = track[chrom][i]
		assert x1 <= ref[0]
		if ref[0] > x2:
			break
		y2 = min(x2, ref[1])
		#print >> sys.stderr, "%s/%d: adding %d in %d %s" % (chrom,i, y2-ref[0]+1, ref[1]-ref[0]+1, ref)
		nTrack += (y2-ref[0]+1)
		#used[chrom][i] += (y2-ref[0]+1)
		#if used[chrom][i] > (track[chrom][i][1] - track[chrom][i][0] + 1):
			#print >> sys.stderr, "Depassement en %s/%d : %d sur %d" % (chrom,i, used[chrom][i], (track[chrom][i][1] - track[chrom][i][0] + 1))
		i += 1

	return "%d\t%.2f" % (nTrack, (100.*nTrack)/(x2-x1+1))


f = utils.myFile.openFile(arguments["list"], "r")
for l in f:
	l = l[:-1]
	t = l.split("\t")
	x1 = int(t[3])
	x2 = int(t[4])
	print l + "\t" + getDensity(t[2], x1, x2)
f.close()

#for chrom in track:
#	for (l,(x1,x2)) in zip(used[chrom],track[chrom]):
#		print >> sys.stderr, l, x2-x1+1


