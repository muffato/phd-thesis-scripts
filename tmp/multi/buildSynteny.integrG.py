#!/usr/bin/env python2
# Copyright (c) 2010 IBENS/Dyogen : Matthieu MUFFATO, Alexandra LOUIS, Hugues ROEST CROLLIUS
# Mail : hrc@ens.fr or alouis@biologie.ens.fr
# Licences GLP v3 and CeCILL v2

__doc__ = """
	A partir de diagonales pair-wise, construit des versions integrees qui correspondent a des segments de chromosomes ancestraux
"""


import sys
import collections

import utils.myPhylTree
import utils.myGenomes
import utils.myFile
import utils.myTools
import utils.myMaths

import myDiags


# Arguments
arguments = utils.myTools.checkArgs( \
	[("phylTree.conf",file), ("target",str), ("diagGroups",file)], \
	[("minimalWeight",int,1), ("power",float,10.), \
	("IN.ancDiags",str,""), \
	("OUT.ancDiags",str,"anc/diags.%s.list.bz2"), \
	("ancGenesFile",str,"~/work/data/ancGenes/ancGenes.%s.list.bz2")], \
	__doc__ \
)


# L'arbre phylogenetique
phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])

targets = set(myDiags.getTargets(phylTree, arguments["target"])[1])

graphs = {}
ancGenomes = {}
for anc in targets:
	ancGenomes[anc] = utils.myGenomes.Genome(arguments["IN.ancDiags"] % phylTree.fileName[anc])
	graphs[anc] = myDiags.WeightedDiagGraph()

# Est-ce que la paire est vue dans l'ancetre
def isAncOK(anc, (g1,s1), (g2,s2)):
	(c1,i1) = ancGenomes[anc].dicGenes[g1]
	(c2,i2) = ancGenomes[anc].dicGenes[g2]
	if c1 != c2:
		return False
	if i2 != (i1 + s1*ancGenomes[anc].lstGenes[c1][i1].strand):
		return False
	if s1*ancGenomes[anc].lstGenes[c1][i1].strand != s2*ancGenomes[anc].lstGenes[c2][i2].strand:
		return False
	return True


# Remplissage des graphes et calcul du poids
f = utils.myFile.openFile(arguments["diagGroups"], "r")
for l in f:
	lst = eval(l)
	nbTot = 0
	nbOK = 0
	for x in lst:
		if x[0] in targets:
			nbTot += 1
			if isAncOK(x[0], (str(x[1]),x[3]), (str(x[2]),x[4])):
				nbOK += 1
	if nbTot > 0:
		weight = arguments["power"] ** (float(nbOK)/nbTot)
		print "weight", weight
		for x in lst:
			if x[0] in targets:
				graphs[x[0]].addDiag([(x[1],x[3]),(x[2],x[4])], weight=weight)
f.close()

def do(anc):

	genesAnc = utils.myGenomes.Genome(arguments["ancGenesFile"] % phylTree.fileName[anc])
	nbGenesAnc = len(genesAnc.lstGenes[None])

	print >> sys.stderr, "Blocs integres de %s ..." % anc,
	print "new anc", anc
	print "nb genes", nbGenesAnc
	
	# Coupe du graphe
	graphs[anc].printIniGraph()
	graphs[anc].cleanGraphTopDown(arguments["minimalWeight"])

	f = utils.myFile.myTSV.writer(arguments["OUT.ancDiags"] % phylTree.fileName[anc])
	s = []
	singles = set(xrange(nbGenesAnc))

	# Extraction des blocs integres
	for (d,dw) in graphs[anc].getBestDiags():
		ds = [x[1] for x in d]
		da = [x[0] for x in d]

		if len(da) > 1:
			singles.difference_update(da)
			s.append( len(da) )

		res = [anc, len(da), utils.myFile.myTSV.printLine(da," "), utils.myFile.myTSV.printLine(ds, " "), utils.myFile.myTSV.printLine(dw," ")]
		f.csvobject.writerow(res)

	# Singletons restants
	for x in singles:
		res = [anc, 1, x, 1, ""]
		f.csvobject.writerow(res)

	f.file.close()

	print >> sys.stderr, utils.myMaths.myStats.txtSummary(s), "+ %d singletons OK" % len(singles)

# Traitement final
for anc in targets:
	do(anc)


