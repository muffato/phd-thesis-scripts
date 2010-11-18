#!/usr/bin/env python2

__doc__ = """
	Affiche une matrice recapitulant toutes les comparaisons pair-wise
"""

import sys
import collections
import itertools

import utils.myFile
import utils.myMaths
import utils.myTools
import utils.myPhylTree

#Extraction des diagonales entre Canis familiaris et Homo sapiens [Boreoeutheria] ... 2 [5/9/17] [11/19/31] 87 [12.47/10.91-1866] OK


# Arguments
arguments = utils.myTools.checkArgs( [("genomes",utils.myTools.FileList(1))], [("index",int,8), ("remove2XSpecies",str,"")], __doc__ )

alldata = []
allesp = []
for fs in arguments["genomes"]:
	data = {}
	f = utils.myFile.openFile(fs, "r")
	for l in f:
		if l.startswith("Extraction"):
			t = l.replace("[", " ").replace("/", " ").replace("]", " ").replace("-", " ").split()
			if int(t[-2]) == 0:
				continue
			data[(t[4]+" "+t[5],t[7]+" "+t[8])] = [float(t[x]) for x in range(-12,-3)+[-2]]
	f.close()
	alldata.append(data)
	allesp.append(frozenset(e for (e,_) in data) | frozenset(e for (_,e) in data))

# On verifie qu'il y a les memes comparaisons dans tous les fichiers
allespS = set(allesp)
assert len(allespS) == 1

if len(arguments["remove2XSpecies"]) != 0:
	phylTree = utils.myPhylTree.PhylogeneticTree(arguments["remove2XSpecies"])
	lesp = sorted(phylTree.lstEsp6X | phylTree.lstEspFull)
else:
	lesp = sorted(allesp[0])

print "\t".join([""] + lesp)
for e1 in lesp:
	res = [e1]
	for e2 in lesp:
		l = [d.get((e1,e2), d.get((e2,e1))) for d in alldata]
		res.append(utils.myMaths.myStats.mean([x[arguments["index"]] for x in l if x is not None]))
	print utils.myFile.myTSV.printLine(res)


