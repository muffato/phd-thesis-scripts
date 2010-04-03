#! /users/ldog/muffato/python

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
	[("minimalWeight",int,1), \
	("OUT.ancDiags",str,"anc/diags.%s.list.bz2"), \
	("ancGenesFile",str,"~/work/data/ancGenes/ancGenes.%s.list.bz2")], \
	__doc__ \
)


# L'arbre phylogenetique
phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])

targets = set(myDiags.getTargets(phylTree, arguments["target"])[1])

graphs = {}
for anc in targets:
	graphs[anc] = myDiags.WeightedDiagGraph()


# Remplissage des graphes et calcul du poids
f = utils.myFile.openFile(arguments["diagGroups"], "r")
for l in f:
	lst = eval(l)
	for x in lst:
		if x[0] in targets:
			graphs[x[0]].addDiag([(x[1],x[3]),(x[2],x[4])], weight=len(lst))
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


