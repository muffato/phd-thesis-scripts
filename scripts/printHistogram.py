#! /users/ldog/muffato/python

__doc__ = """
	Lit des donnees dans un fichier et affiche les valeurs de l'histogramme
"""

import sys

import utils.myFile

# Les fichiers a regrouper (par defaut, on lit l'entree standard)
files = set(sys.argv[1:])

nbBins = 100
withFreq = False
for f in list(files):
	if f.startswith("-nbBins="):
		nbBins = int(f[8:])
		files.discard(f)
	elif f.startswith("-freq"):
		withFreq = True
		files.discard(f)

if len(files) == 0:
	files.add("-")

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

c = float(len(lst)) if withFreq else 1.

if nbBins < 0:
	import collections
	d = collections.defaultdict(int)
	for x in lst:
		d[x] += 1
	for x in sorted(d):
		print x, d[x]/c
else:
	mn = min(lst)
	mx = max(lst)
	res = [0] * nbBins
	for x in lst:
		i = int(nbBins * (float(x-mn)/(mx-mn)))
		if i == nbBins:
			i = nbBins-1
		res[i] += 1
	for (i,x) in enumerate(res):
		print float(i*(mx-mn))/nbBins+mn, x/c

