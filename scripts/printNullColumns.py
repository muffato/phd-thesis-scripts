#! /users/ldog/muffato/python -OO

__doc__ = """
	Affiche les profils de colonne NULL
"""

import sys
import utils.myTools

count = utils.myTools.defaultdict(int)
for l in sys.stdin:
	t = l.replace('\n','').split('\t')
	signatures = tuple([x!="\N" for x in t])
	count[signatures] += 1

for (x,c) in count.iteritems():
	print x, c

