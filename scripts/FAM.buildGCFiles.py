#! /users/ldog/muffato/python -OO

__doc__ = """
	Lit une serie d'arbres resultats de NHML.
	Extrait les GC ancestraux et les place pour chaque famille
"""

import sys
import utils.myTools
import utils.myPhylTree

# Arguments
(noms_fichiers, options) = utils.myTools.checkArgs( ["phylTree.conf", "res"], [("ancGenesFile",str,""), ("outputGCFile",str,"")], __doc__ )

# L'arbre phylogenetique
phylTree = utils.myPhylTree.PhylogeneticTree(noms_fichiers["phylTree.conf"])

utils.myTools.mkDir(options["outputGCFile"])

# 1. On charge toutes les familles et on garde en leurs coordonnees
dicFamilles = {}
dicGC = {}
for esp in phylTree.allNames:
	nb = 0
	for gene in utils.myGenomes.Genome(options["ancGenesFile"] % phylTree.fileName[esp]):
		names = frozenset(gene.names)
		dicFamilles[names] = (esp, gene.beginning)
		nb += 1
	dicGC[esp] = [None] * nb

# 2. On charge les arbres, reconstruit les familles et range les GC
f = utils.myTools.myOpenFile(noms_fichiers["res"], "r")
for l in f:
	
	t = l.split('\t')
	
	names = frozenset(eval(t[1]))
	GC = float(t[3])
	(esp,i) = dicFamilles.get( frozenset(names), (None,0) )
	if esp != None:
		dicGC[esp][i] = GC
	else:
		print >> sys.stderr, "ERROR", names

f.close()

for (esp,tab) in dicGC.iteritems():
	f = utils.myTools.myOpenFile(options["outputGCFile"] % phylTree.fileName[esp], "w")
	print >> f, utils.myTools.printLine(tab, "\n")
	f.close()

