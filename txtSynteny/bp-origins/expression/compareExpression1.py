#!/usr/bin/env python2

__doc__ = """
	Charge les liens GO et renvoie les niveaux d'intersection
"""

import sys
import math
import itertools
import collections

import utils.myFile
import utils.myTools
import utils.myMaths
import utils.myGenomes

# Arguments
arguments = utils.myTools.checkArgs( [("exp_data",file), ("gene_pairs",file)], [], __doc__)

exp_data = {}
f = utils.myFile.openFile(arguments["exp_data"], "r")
for l in f:
	t = l.split()
	exp_data[t[1]] = [float(x) for (i,x) in enumerate(t[3:]) if (i+1) not in [1,17,70,71,72,74,75]]
f.close()

f = utils.myFile.openFile(arguments["gene_pairs"], "r")
for l in f:
	(g1,g2) = l.split()
	if (g1 not in exp_data) or (g2 not in exp_data):
		continue
	l = [abs(math.log(e1/e2, 2)) for (e1,e2) in zip(exp_data[g1], exp_data[g2]) if e1<4000 and e2<4000]
	if len(l) == 0:
		continue
	score = round(sum(l)/len(l), 4)
	#score = sum([abs(math.log(e1/e2, 2)) for (e1,e2) in zip(exp_data[g1], exp_data[g2])])
	print g1, g2, score
f.close()


