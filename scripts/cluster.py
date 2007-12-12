#! /users/ldog/muffato/python -OO

import os
import sys
import utils.myTools

(_,options) = utils.myTools.checkArgs([], [("compr",str,"gzip --best")], "Clustering by compression")

lst = []
for l in sys.stdin:
	lst.append(l[:-1].replace("'", "\\'"))
nb = len(lst)
print >> sys.stderr, nb

def getIniSize(s):
	return os.stat(s)[6]

def getCompSize(s):
	(stdin,stdout) = os.popen2("cat %s | %s | wc -c" % (s,options["compr"]))
	stdin.close()
	n = long(stdout.readline())
	stdout.close()
	return n

def getCompSize2(s1, s2):
	(stdin,stdout) = os.popen2("cat %s %s | %s | wc -c" % (s1,s2,options["compr"]))
	stdin.close()
	n = long(stdout.readline())
	stdout.close()
	return n

#def getDist(i1, i2):
#	return 1. - float(sizesComp[i1] + sizesComp[i2] - sizesComp2[i1][i2]) / max(sizesComp[i1],sizesComp[i2])

#sizesIni = [getIniSize(x) for x in lst]

sizesComp = [getCompSize(x) for x in lst]

#sizesComp2 = [[getCompSize2(x,y) for x in lst] for y in lst]

print nb

#dist = [[getDist(x,y) for x in xrange(len(lst))] for y in xrange(len(lst))]

for (x,nx) in enumerate(lst):
	print nx
	cx = sizesComp[x]
	for (y,ny) in enumerate(lst):
		if x == y:
			print 0
		else:
			cy = sizesComp[y]
			d = 1. - float(cx + cy - getCompSize2(nx,ny)) / max(cx,cy)
			print d
	print

