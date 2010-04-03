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
	[("phylTree.conf",file), ("target",str), ("pairwiseDiags",file)], \
	[("minimalWeight",int,1), ("minimalLength",int,2), \
	("OUT.ancDiags",str,"anc/diags.%s.list.bz2"), \
	("ancGenesFile",str,"~/work/data/ancGenes/ancGenes.%s.list.bz2")], \
	__doc__ \
)


# L'arbre phylogenetique
phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])

targets = set(myDiags.getTargets(phylTree, arguments["target"])[1])

print >> sys.stderr, "Chargement des diagonales pairwise ...",
pairwiseDiags = collections.defaultdict(list)
n = 0
f = utils.myFile.openFile(arguments["pairwiseDiags"], "r")
for l in f:
	t = l[:-1].split("\t")
	assert len(t) == 9
	diagS = [int(x) for x in t[7].split()]
	if len(diagS) < arguments["minimalLength"]:
		continue
	for x in t[8].split("|"):
		n += 1
		(anc,x) = x.split("=")
		if anc in targets:
			pairwiseDiags[anc].append((x,diagS))
f.close()
print >> sys.stderr, n, "OK"

def do(anc):

	genesAnc = utils.myGenomes.Genome(arguments["ancGenesFile"] % phylTree.fileName[anc])

	graph = myDiags.WeightedDiagGraph()

	print >> sys.stderr, "Blocs pairwise de %s ..." % anc,
	s = []
	for (da,ds) in pairwiseDiags[anc]:
		d = zip([int(x) for x in da.split()], ds)
		s.append(len(d))
		graph.addDiag(d)
	print >> sys.stderr, utils.myMaths.myStats.txtSummary(s), "OK"

	print >> sys.stderr, "Blocs integres de %s ..." % anc,
	print "NEWANC", anc
	
	# Coupe du graphe
	graph.printIniGraph()
	graph.cleanGraphTopDown(arguments["minimalWeight"])

	f = utils.myFile.myTSV.writer(arguments["OUT.ancDiags"] % phylTree.fileName[anc])
	s = []
	singles = set(xrange(len(genesAnc.lstGenes[None])))

	# Extraction des blocs integres
	for (d,dw) in graph.getBestDiags():
		ds = [x[1] for x in d]
		da = [x[0] for x in d]

		if len(d) > 1:
			s.append( len(d) )

		singles.difference_update(da)
		res = [anc, len(da), utils.myFile.myTSV.printLine(da," "), utils.myFile.myTSV.printLine(ds, " "), utils.myFile.myTSV.printLine(dw," ")]
		f.csvobject.writerow(res)

	for x in singles:
		res = [anc, 1, x, 1, ""]
		f.csvobject.writerow(res)
	f.file.close()

	print >> sys.stderr, utils.myMaths.myStats.txtSummary(s), "+ %d singletons OK" % len(singles)

# Traitement final
for anc in targets:
	do(anc)


