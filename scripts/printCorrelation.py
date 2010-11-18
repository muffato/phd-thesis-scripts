#!/usr/bin/env python2

__doc__ = """
	Renvoie la correlation entre les deux jeux de donnees
"""

import sys

import utils.myFile
import utils.myMaths


lx1 = []
lx2 = []
spearman = False

if "-rank" in sys.argv:
	sys.argv.remove("-rank")
	spearman = True
f = utils.myFile.openFile(sys.argv[1], "r") if len(sys.argv) > 1 else sys.stdin

for l in f:
	l = l.strip()
	if l.startswith("#"):
		continue
	t = l.replace('\n', '').split()
	if len(t) == 2:
		try:
			d1 = float(t[0])
			d2 = float(t[1])
			lx1.append(d1)
			lx2.append(d2)
		except ValueError:
			pass
f.close()

if spearman:
	def transform(l):
		size = len(l)
		ls = sorted((x,i) for (i,x) in enumerate(l))
		newl = []
		i = 0
		while i < size:
			currval = ls[i][0]
			s = i
			j = i + 1
			while (j < size) and (currval == ls[j][0]):
				s += j
				j += 1
			for k in xrange(i, j):
				newl.append( (ls[i][1], s/float(j-i)) )
			i = j
		return [newi for (oldi,newi) in sorted(newl)]

	lx1 = transform(lx1)
	lx2 = transform(lx2)

print utils.myMaths.myStats.correlation(lx1, lx2)

