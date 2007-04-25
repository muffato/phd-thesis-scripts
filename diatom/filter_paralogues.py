#! /users/ldog/muffato/python -OO


# Librairies
import sys
import math

ev = float(sys.argv[1])
id = float(sys.argv[2])

for l in sys.stdin:
	c = l.split()
	evalue = float(c[10])
	if evalue <= 0 or math.log10(evalue) <= ev:
		a = int(c[0])
		b = int(c[1][3:])
		if a != b and float(c[2]) >= id:
			#print l,
			#print "Pt%d Pt%d" % (a,b)
			print "Tp%d Tp%d" % (a,b)


