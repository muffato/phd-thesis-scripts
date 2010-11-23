#!/usr/bin/env python2
# Copyright (c) 2010 IBENS/Dyogen : Matthieu MUFFATO, Alexandra LOUIS, Hugues ROEST CROLLIUS
# Mail : hrc@ens.fr or alouis@biologie.ens.fr
# Licences GLP v3 and CeCILL v2

__doc__ = """
	Lit des donnees dans un fichier et les regroupe
"""

import sys
import operator
import collections

import utils.myFile
import utils.myTools


arguments = utils.myTools.checkArgs([("files",utils.myTools.FileList(1))], [("neededProp",float,1.)], __doc__, showArgs=False)

# Impression des diagonales ancestrales de Homininae ... 2 [3/11/30] [28/59/103] 314 [24.18/35.19-785] + 3496 singletons OK

data = collections.defaultdict(dict)
for s in arguments["files"]:
	f = utils.myFile.openFile(s, 'r')
	for l in f:
		if "ancestrales" in l:
			t = l.split()
			data[t[5]][s] = float(t[11][1:].split('/')[0])
	f.close()

reffile = arguments["files"][0]
files = arguments["files"][1:]

for anc in data:
	ref = data[anc][reffile]
	files = [f for f in files if data[anc][f] >= ref*arguments["neededProp"]]

print >> sys.stderr, files
files.insert(0, reffile)

for anc in data:
	val = [(data[anc][x],x) for x in files]
	valm = max(val)[0]
	print anc, [x for x in val if x[0] == valm][0][1]
	#print anc, max([(data[anc][x],x) for x in files])[1]

