#!/usr/bin/env python2
# Copyright (c) 2010 IBENS/Dyogen : Matthieu MUFFATO, Alexandra LOUIS, Hugues ROEST CROLLIUS
# Mail : hrc@ens.fr or alouis@biologie.ens.fr
# Licences GLP v3 and CeCILL v2

__doc__ = """
	Compare les genes successifs dans le genome moderne et renvoie les scores de simultaneite
	Renvoie les memes scores pour les memes genes pris dans des ordres aleatoires
"""

import sys
import math
import random
import itertools
import collections

import utils.myMaths
import utils.myTools
import utils.myGenomes
import utils.myPhylTree

newd = lambda : collections.defaultdict(int)
run = collections.defaultdict(newd)

tot = 0
for l in sys.stdin:
	t = l.split()
	i = int(t[0])
	s = int(10*float(t[1]))
	run[i][s] += 1
	tot += 1

nb = dict( (i,sum(run[i].itervalues())) for i in run)
nbX = sum(nb[i] for i in run if i > 0) / (len(run)-1.)
factor = nb[0]/nbX
for s in [-10] + range(11):
	l = [run[i][s]*factor for i in run if i > 0]
	m = utils.myMaths.myStats.mean(l)
	print float(s)/10, m, 1.96*utils.myMaths.myStats.stddev(l, m=m)/math.sqrt(len(l))

