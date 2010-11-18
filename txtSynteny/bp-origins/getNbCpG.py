#!/usr/bin/env python2

__doc__ = """
	Dessine les densites de coorientation
"""

import sys
import math
import collections

import utils.myFile
import utils.myTools
import utils.myMaths
import utils.myGenomes
import utils.myPsOutput

# Arguments
arguments = utils.myTools.checkArgs( [("dataFile",file), ("cpglist",file)], [], __doc__)

allcpg = set()
f = utils.myFile.openFile(arguments["cpglist"], "r")
for l in f:
	t = l[:-1].split()
	if float(t[5]) > 0.375:
		allcpg.add(t[0])
f.close()

f = utils.myFile.openFile(arguments["dataFile"], "r")
for l in f:
	l = l[:-1]
	t = l.split('\t')
	
	(g1,g2) = t[1].split('/')
	(s1,s2) = t[5].split('/')
	n = 0
	if (s1 == "-1") and (g1 in allcpg):
		n += 1
	if (s2 == "1") and (g2 in allcpg):
		n += 1
	
	print l + ("\t%d" % n)
	
	#if (s1 != s2) and (n == 2):
	#	print l,
		

f.close()

