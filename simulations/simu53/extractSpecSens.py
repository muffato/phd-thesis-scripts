#!/usr/bin/env python2

import sys

for l in sys.stdin:
	t = l.split()
	if t[0] == "ALLPAIRS":
		anc = t[1]
		scores = [int(x) for x in t[2:]]
	elif t[0] == "REFPAIRS":
		tot = int(t[1])
		s = 0
		txt = anc
		for x in scores[1:]:
			s += x
			txt = txt + "\t%.2f %.2f" % ((100.*s)/scores[0], (100.*s)/tot)
			#print "SPEC", anc, (100.*s)/scores[0]
			#print "SENS", anc, (100.*s)/tot
		print txt

