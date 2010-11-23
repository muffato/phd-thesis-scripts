#!/usr/bin/env python2
# Copyright (c) 2010 IBENS/Dyogen : Matthieu MUFFATO, Alexandra LOUIS, Hugues ROEST CROLLIUS
# Mail : hrc@ens.fr or alouis@biologie.ens.fr
# Licences GLP v3 and CeCILL v2

__doc__ = """
	Affiche les profils de colonne NULL
"""

import sys
import utils.myTools

count = collections.defaultdict(int)
for l in sys.stdin:
	t = l.replace('\n','').split('\t')
	signatures = tuple([x!="\N" for x in t])
	count[signatures] += 1

for (x,c) in count.iteritems():
	print x, c

