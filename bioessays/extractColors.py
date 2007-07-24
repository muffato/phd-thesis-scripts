#! /users/ldog/muffato/python -OO

import sys
import operator
import itertools

tab = [range(1028) for i in xrange(1284)]

for l in sys.stdin:
	c = l.split()
	x = int(c[0])
	y = int(c[1])
	color = (int(c[2][0:2], 16), int(c[2][2:4], 16), int(c[2][4:6], 16))
	tab[x][y] = color


#chrom = {"a": (0,182,220), "b":(111,190,78), "c":(208,144,143), "d":(253,202,43), "e":(136,144,150), "f":(241,210,227), "g":(80,199,238), "h":(210,199,51), "i":(0,175,81), "j":(237,86,258), "k":(253,239,132), "l":(38,123,191), "n":(234,50,50), "-":(255,255,255)}
chrom = {"a": (0,176,215), "b":(51,181,78), "c":(199,41,146), "d":(248,201,34), "e":(138,143,147), "f":(241,187,185), "g":(1,190,236), "h":(200,191,55), "i":(1,169,80), "j":(239,97,167), "k":(252,234,43), "l":(52,56,160), "n":(234,48,47), "-":(255,255,255)}

#chrom = {"a": (0,176,215), "b":(51,181,78), "c":(199,41,146), "d":(248,201,34), "e":(138,143,147), "f":(241,187,185), "g":(1,190,236), "h":(200,191,55), "i":(1,169,80), "j":(239,97,167), "k":(252,234,43), "l":(52,56,160), "n":(234,48,47), "white":(255,255,255)}
#for c in chrom:
#	s = "#%02X:%02X:%02X" % chrom[c]
#	print "| sed 's/%s/%s/g'" % (s,c),

#sys.exit(0)

# MED3 = Homme
#for ch in xrange(24):
#	currChrom = []
#	x = int(255 + 34.5*ch)
#	for y in xrange(80+9*ch,960):
#		col = [0.,0.,0.]
#		for xx in xrange(x-7,x+7):


# MED4 = Tetraodon
#for ch in xrange(29):
#	currChrom = []
#	x = int(63 + 42*ch)
#	for y in xrange(170,877):
#		col = [0.,0.,0.]
#		for xx in xrange(x-10,x+10):

# MED5 = Medaka
#for ch in xrange(29):
#	currChrom = []
#	x = int(65 + 42*ch)
#	for y in xrange(135,910):
#		col = [0.,0.,0.]
#		for xx in xrange(x-10,x+10):

# MED6 = Zebrafish
#for ch in xrange(29):
#	currChrom = []
#	x = int(62 + 42*ch)
#	for y in xrange(177,872):
#		col = [0.,0.,0.]
#		for xx in xrange(x-10,x+10):


d = {}
for ch in xrange(29):
	currChrom = []
	x = int(62 + 42*ch)
	for y in xrange(177,872):
		col = [0.,0.,0.]
		for xx in xrange(x-10,x+10):
			for i in xrange(3):
				col[i] += tab[xx][y][i]
		for i in xrange(3):
			col[i] /= 20.
		for c in chrom:
			d[c] = sum( [ (chrom[c][i]-col[i])**2 for i in xrange(3) ] )
		it = d.items()
		it.sort(key = operator.itemgetter(1))
		print >> sys.stderr, ch+1, it[0][0], tab[x][y], it
		currChrom.append(it[0][0])
	i0 = None
	for (i,val) in enumerate(currChrom):
		if val != "-":
			if i0 == None:
				i0 = i
			i1 = i
	for (val,items) in itertools.groupby(currChrom[i0:i1+1]):
		print "#%02X:%02X:%02X" % chrom[val],
		print " %f ;" % (len(list(items))*.1),
	print

