#!/usr/bin/env python2

__doc__ = """
	Selectionne les points de cassure dont le score est superieure a la moyenne de tous les scores rencontres
"""

import sys
import collections

import utils.myFile
import utils.myTools


arguments = utils.myTools.checkArgs( [("bpFile",file)], [], __doc__)

f = utils.myFile.openFile(arguments["bpFile"], "r")
dic = collections.defaultdict(list)
for l in f:
	dic[l.split()[0]].append(l)
f.close()

for key in ["LOST", "GAIN"]:
	data = []
	for l in dic[key]:
		for x in l.split():
			if (x[0] == "[") and (x[-1] == "]"):
				data.append(int(x[1:-1]))
	mean = sum(data)/float(len(data))
	n = 0
	for l in dic[key]:
		t = l.split()
		if int(t[4][6:]) >= mean:
			n += 1
			print l,
	print >> sys.stderr, key, mean, "%d/%d" % (n, len(dic[key]))

