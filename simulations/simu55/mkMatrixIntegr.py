#! /users/ldog/muffato/python

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

#Impression des diagonales ancestrales de Boreoeutheria ... 2 [5/173/566] [531/704/1194] 1874 [344.68/409.68-69] + 119 singletons OK


# Arguments
arguments = utils.myTools.checkArgs( [("genomes",utils.myTools.FileList(1))], [("remove2XSpecies",str,"")], __doc__ )

alldata = []
allesp = []
for fs in arguments["genomes"]:
	data = {}
	f = utils.myFile.openFile(fs, "r")
	for l in f:
		if "ancestrales" in l:
			t = l.replace("[", " ").replace("/", " ").replace("]", " ").replace("-", " ").split()
			if len(t) > 15:
				data[" ".join(t[5:-16])] = [float(t[x]) for x in range(-15,-6)+[-5,-3]]
	f.close()
	alldata.append(data)
	allesp.append(frozenset(data))

# On verifie qu'il y a les memes comparaisons dans tous les fichiers
allespS = set(allesp)
assert len(allespS) == 1

if len(arguments["remove2XSpecies"]) != 0:
	phylTree = utils.myPhylTree.PhylogeneticTree(arguments["remove2XSpecies"])
	#lesp = sorted(anc for anc in phylTree.listAncestr if min([len(set(phylTree.species[e]).difference(phylTree.lstEsp2X)) for (e,_) in phylTree.items[anc]]) > 0)
	lesp = sorted(anc for anc in phylTree.listAncestr if len([e for (e,_) in phylTree.items[anc] if  len(set(phylTree.species[e]).difference(phylTree.lstEsp2X)) > 0]) >= 2)
else:
	lesp = sorted(allesp[0])

print utils.myFile.myTSV.printLine(["", "min", "25%", "50%", "75%", "N75", "N50", "N25", "max", "mean", "nb", "singletons"])
for e in lesp:
	res = [e]
	for i in xrange(11):
		res.append(utils.myMaths.myStats.mean([d[e][i] for d in alldata]))
	print utils.myFile.myTSV.printLine(res)


