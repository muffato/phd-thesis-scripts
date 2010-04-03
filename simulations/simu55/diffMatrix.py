#! /users/ldog/muffato/python

__doc__ = """
	Calcule la matrice des differences relatives entre les deux matrices
"""

import sys
import itertools
import collections

import utils.myFile
import utils.myMaths
import utils.myTools

#v1/rearrang-RandomTrees/4/log.pairwise:Extraction des diagonales entre Canis familiaris et Homo sapiens [Boreoeutheria] ... 2 [5/9/17] [11/19/31] 87 [12.47/10.91-1866] OK


# Arguments
arguments = utils.myTools.checkArgs( [("matrix1",file), ("matrix2",file)], [], __doc__, showArgs=False )

def loadMatrix(s):
	f = utils.myFile.openFile(s, "r")
	data = [l[:-1].split("\t") for l in f]
	f.close()
	cols = data[0][1:]
	rows = [x[0] for x in data][1:]
	val = [[float(y) if y != "None" else None for y in x[1:]] for x in data[1:]]
	return (cols,rows,val)

(cols1,rows1,mat1) = loadMatrix(arguments["matrix1"])
(cols2,rows2,mat2) = loadMatrix(arguments["matrix2"])

assert cols1 == cols2, (cols1,cols2)
assert rows1 == rows2, (rows1,rows2)

diff = []
print utils.myFile.myTSV.printLine([""] + cols1)
for (y,r) in enumerate(rows1):
	res = [r]
	for (x,c) in enumerate(cols1):
		x1 = mat1[y][x]
		x2 = mat2[y][x]
		if (x1 is None) or (x2 is None):
			res.append(None)
		else:
			res.append(round(100*(x2-x1)/x1,3))
	diff.extend(x for x in res[1:] if x is not None)
	print utils.myFile.myTSV.printLine(res)

print >> sys.stderr, utils.myMaths.myStats.txtSummary(diff, withN50=False)

