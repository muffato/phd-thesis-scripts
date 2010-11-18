#!/usr/bin/env python2

import sys
import itertools
import collections

import utils.myFile
import utils.myMaths
import utils.myTools
import utils.myGenomes
import utils.myPhylTree

arguments = utils.myTools.checkArgs( [("allRootSeqs",file)], [], "Calcule les frequences de chaque codon")

count = collections.defaultdict(int)

f = utils.myFile.openFile(arguments["allRootSeqs"], "r")
for l in f:
	l = l[:-1]
	assert (len(l)%3) == 0
	for i in xrange(len(l)/3):
		count[l[3*i:3*i+3]] += 1
f.close()

count["TGA"] = 0
count["TAA"] = 0
count["TAG"] = 0

print count

nb = float(sum(count.values()))

print nb

for l in itertools.product("TCAG", "TCAG"):
	for b in "TCAG":
		codon = ''.join(l) + b
		print "%.6f" % (count[codon]/nb),
		#print codon, count[codon], "%.6f" % (count[codon]/nb)
	print

