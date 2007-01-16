#! /users/ldog/muffato/python -OO

__doc__ = """
A partir de toutes les diagonales extraites entre les especes,
  reconstruit les chromosomes (ou scaffold) de chaque ancetre.
"""


##################
# INITIALISATION #
##################

# Librairies
import math
import sys
import utils.myGenomes
import utils.myTools
import utils.myMaths


nb1_o2o = 0
nb1_many = 0
for l in utils.myTools.myOpenFile(sys.argv[1], 'r'):
	c = l.split()
	if "one2one" in c[6]:
		nb1_o2o += 1
	elif "many" in c[6]:
		nb1_many += 1
	else:
		print "????"


s = set([])
for l in utils.myTools.myOpenFile(sys.argv[2], 'r'):
	c = l.split()
	if "one2one" in c[4]:
		s.add(c[0])

nb2_o2o = 0
nb2_many = 0
for l in utils.myTools.myOpenFile(sys.argv[2], 'r'):
	c = l.split()
	if c[0] in s:
		nb2_o2o += 1
	else:
		nb2_many += 1

print "*one2one.V1=%d"%nb1_o2o, "*many.V1=%d"%nb1_many, "WITH*one2one.V2=%d" % len(s), "LINES*one2one.V2=%d" % nb2_o2o, "ONLY*many.V2=%d"%nb2_many

