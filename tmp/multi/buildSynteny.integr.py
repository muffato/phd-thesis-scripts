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
	[("phylTree.conf",file), ("target",str), ("step1PairwiseDiags",file), ("diagGroups",file)], \
	[("minimalWeight",int,1), ("power",float,10.), \
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

# Liste des paires de genes (ancestrales + modernes) bloquees
blocked = set()
f = utils.myFile.openFile(arguments["step1PairwiseDiags"], "r")
for l in f:
	t = l[:-1].split("\t")
	assert len(t) == 9, (l,t)
	assert t[0] == "2"
	strands = tuple(int(x) for x in t[7].split())
	group = [(t[1],) + tuple(t[3].split()) + strands, (t[4],) + tuple(t[6].split()) + strands]
	for x in t[8].split("|"):
		(anc,_,la) = x.partition("=")
		group.append( (anc,) + tuple(int(x) for x in la.split()) + strands )
	blocked.update(group)
	blocked.update((x[0],x[2],x[1],-x[4],-x[3]) for x in group)
f.close()
print >> sys.stderr, len(blocked)

# Remplissage des graphes et calcul du poids
f = utils.myFile.openFile(arguments["diagGroups"], "r")
for l in f:
	lst = eval(l)
	nblocked = sum(1 for x in lst if x in blocked)
	weight = arguments["power"] ** (float(nblocked) / len(lst))
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


