#!/usr/bin/env python2

__doc__ = """
	Verifie que toutes les contraintes d'integrite sont respectees
"""

import sys
import utils.myTools
import utils.myPhylTree


arguments = utils.myTools.checkArgs( [("databaseFiles",str), ("phylTree.conf",file)], [], __doc__)
phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])
inputFiles = {}
for x in ["especes","genes","arbres","syntenies"]:
	inputFiles[x] = utils.myTools.myOpenFile(arguments["databaseFiles"] % x, 'r')


# ESPECES
##########
dataEspece = list(utils.myFile.myTSV.readTabular(arguments["databaseFiles"] % "especes", [int,str,int,str]))
dicEspID = {}
for x in dataEspece:
	dicEspID[x[0]] = x[1]
	# Le nom existe
	assert x[1] in phylTree.officialName
	# L'age est correct
	assert x[2] == phylTree.ages[x[1]]
# unicite de l'ID
assert len(dataEspece) == len(dicEspID)


# ARBRES
#########
dataArbres = list(utils.myFile.myTSV.readTabular(inputFiles["arbres"], [int,int,int,str,str,int,int]))
allgenesID = set(x[0] for x in dataArbres)
# unicite de l'ID
assert len(dataArbres) == len(allgenesID)
# unicite de gauche
assert len(dataArbres) == len(set(x[1] for x in dataArbres))
# unicite de droite
assert len(dataArbres) == len(set(x[2] for x in dataArbres))
# double unicite
assert (2*len(dataArbres)) == len(set(x[1] for x in dataArbres) | set(x[2] for x in dataArbres))
# Parent refere a un ID
assert set(int(x[3]) for x in dataArbres if x[3] != "\N").issubset(allgenesID)
# Distance est NULL uniquement si Parent est NULL
assert len([x for x in dataArbres if (x[3] == "\N") ^ (x[4] == "\N")]) == 0
# Duplication vaut 0/1/2
assert set(x[5] for x in dataArbres).issubset(range(3))
# Espece refere a un ID
assert set(x[6] for x in dataArbres).issubset(dicEspID)
# Test de la structure d'arbre par rapport aux gauche/droite et aux especes
dicArbre = {}
for x in dataArbres:
	if x[3] == "\N":
		par = None
	else:
		par = int(x[3])
	dicArbre[x[0]] = (x[1],x[2],par,x[5],dicEspID[x[6]])
for (x,(g,d,par,_,esp)) in dicArbre.iteritems():
	if par in dicArbre:
		assert (dicArbre[par][0] < g) and (dicArbre[par][1] > d)
		esp2 = dicArbre[par][4]
		assert phylTree.isChildOf(esp, esp2)
		#print esp == esp2, dicArbre[x][3], dicArbre[par][3]
		#if esp == esp2:
		#	assert dicArbre[par][3] != 0
		#else:
		#assert (dicArbre[par][3] == 0) ^ (esp == esp2), (x,dicArbre[x],dicArbre[par])


# SYNTENIES
############
dataSynt = list(utils.myFile.myTSV.readTabular(inputFiles["syntenies"], [int]*4))
# Orientation
assert set(x[2] for x in dataSynt).issubset([1,-1])
# Gene refere a un ID
assert set(x[1] for x in dataSynt).issubset(allgenesID)
# Unicite (bloc,pos)
assert len(set((x[0],x[3]) for x in dataSynt)) == len(dataSynt)

dicSynt = collections.defaultdict(set)
for x in dataSynt:
	dicSynt[x[0]].add(x[3])
for x in dicSynt.itervalues():
	assert x == set(range(len(x)))


# GENES
########
dataGenes = list(utils.myFile.myTSV.readTabular(inputFiles["genes"], [int,int,str,str,int,int,str,str,str]))
# Unicite (esp,chrom,pos)
assert len(set((x[1],x[3],x[4]) for x in dataGenes)) == len(dataGenes)
# Unicite ID
assert len(set(x[0] for x in dataGenes)) == len(dataGenes)
# gene refere a un ID
assert set(x[0] for x in dataGenes).issubset(allgenesID)
# esp refere a un ID
assert set(x[1] for x in dataGenes).issubset(dicEspID)
for x in dataGenes:
	if (x[6] != "\N") and (x[7] != "\N"):
		assert int(x[6]) < int(x[7])

for f in inputFiles.itervalues():
	f.close()


