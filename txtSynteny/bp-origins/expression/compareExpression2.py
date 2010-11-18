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
exp_data_prop = {}
f = utils.myFile.openFile(arguments["exp_data"], "r")
for l in f:
	t = l.split()
	exp_data[t[1]] = [float(x) for (i,x) in enumerate(t[3:]) if (i+1) not in [1,17,70,71,72,74,75]]
	s = sum(exp_data[t[1]])
	exp_data_prop[t[1]] = [x/s for x in exp_data[t[1]]]
f.close()

f = utils.myFile.openFile(arguments["gene_pairs"], "r")
for l in f:
	(g1,g2) = l.split()
	if (g1 not in exp_data) or (g2 not in exp_data):
		continue
	
	# abs(log2(meanExp1/meanExp2))
	scoreA = abs(math.log(sum(exp_data[g1])/sum(exp_data[g2]), 2))
	
	# mean( abs(log2(exp1.i/exp2.i)) on tissue i )
	scoreB = sum([abs(math.log(e1/e2, 2)) for (e1,e2) in zip(exp_data[g1], exp_data[g2])]) / len(exp_data[g1])
	
	# mean( abs(log2((exp1.i/meanExp1)/(exp2.i/meanExp2))) on tissue i )
	scoreC = sum([abs(math.log(e1/e2, 2)) for (e1,e2) in zip(exp_data_prop[g1], exp_data_prop[g2])]) / len(exp_data[g1])
	
	print g1, g2, round(scoreA,4), round(scoreB,4), round(scoreC,4)

f.close()

