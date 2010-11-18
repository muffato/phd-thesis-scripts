#!/usr/bin/env python2

__doc__ = """
	Lit des donnees dans un fichier et affiche les valeurs de l'histogramme
"""

import sys

import utils.myFile
import utils.myTools

arguments = utils.myTools.checkArgs( \
	[("files",utils.myTools.FileList(0))], \
	[("nbBins",int,100), ("freq",bool,False), ("min",float,None), ("max",float,None)], \
	__doc__, \
showArgs=False)

# Les fichiers a regrouper (par defaut, on lit l'entree standard)
files = arguments["files"]
if len(files) == 0:
	files.append("-")

lst = []
for f in files:
	# Ouverture du fichier
	if f == "-":
		f = sys.stdin
	else:
		f = utils.myFile.openFile(f, 'r')
	# Lecture & regroupement
	for l in f:
		c = l.split()
		for x in c:
			try:
				x = int(x)
			except ValueError:
				x = float(x)
			lst.append(x)
	# Fermeture
	f.close()

c = float(len(lst)) if arguments["freq"] else 1.
nbBins = arguments["nbBins"]
if nbBins < 0:
	import collections
	d = collections.defaultdict(int)
	for x in lst:
		d[x] += 1
	for x in sorted(d):
		print x, d[x]/c
else:
	mn = min(lst) if arguments["min"] is None else arguments["min"]
	mx = max(lst) if arguments["max"] is None else arguments["max"]
	res = [0] * nbBins
	for x in lst:
		i = int(nbBins * (float(x-mn)/(mx-mn)))
		if i == nbBins:
			i = nbBins-1
		res[i] += 1
	for (i,x) in enumerate(res):
		print float(i*(mx-mn))/nbBins+mn, x/c

