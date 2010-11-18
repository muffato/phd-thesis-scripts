#!/usr/bin/env python2

__doc__ = """
	Realise des tirages aleatoires et mesure les moyennes des ensembles.
	En deduit une p-value pour une moyenne observee
"""

import sys
import math
import random
import collections

import utils.myFile
import utils.myTools
import utils.myMaths

# Arguments
arguments = utils.myTools.checkArgs( [("dataFile",file), ("sampleSize",int), ("testedValue",float)], [("nbSamples",int,1000000)], __doc__)

data = []
f = utils.myFile.openFile(arguments["dataFile"], "r")
for l in f:
	data.extend([float(x) for x in l.split()])
f.close()

size = arguments["sampleSize"]
n1 = 0
res = []
for i in xrange(arguments["nbSamples"]):
	if (i % 10000) == 0:
		print >> sys.stderr, "%d / %d done" % (i, arguments["nbSamples"])
	s = sum(random.sample(data, size)) / float(size)
	res.append(s)
	if s < arguments["testedValue"]:
		n1 += 1

print utils.myMaths.myStats.txtSummary(res)
print n1, arguments["nbSamples"]-n1, (100.*n1)/arguments["nbSamples"], 100.-(100.*n1)/arguments["nbSamples"]

