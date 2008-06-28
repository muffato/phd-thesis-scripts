#! /users/ldog/muffato/python -OO

__doc__ = """
	Lit des donnees dans un fichier et affiche les valeurs de l'histogramme
"""

import sys
import math
import utils.myTools
import utils.myMaths

# Arguments
arguments = utils.myTools.checkArgs([("fichier",file)], [("type",str,"float"), ("nbBins",int,1000)], __doc__)

# Lit un flot de nombres et affiche les stats
lst = []
t = eval(arguments["type"])
f = utils.myTools.myOpenFile(arguments["fichier"], 'r')
for l in f:
	c = l.split()
	lst.extend(t(x) for x in c)
f.close()

nb = arguments["nbBins"]
mn = min(lst)
mx = max(lst)
res = [0] * nb
for x in lst:
	i = int(nb * (float(x-mn)/(mx-mn)))
	if i == nb:
		i = nb-1
	res[i] += 1
for (i,x) in enumerate(res):
	print float(i*(mx-mn))/nb+mn, x

