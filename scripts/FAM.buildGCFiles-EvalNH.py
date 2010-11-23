#!/usr/bin/env python2
# Copyright (c) 2010 IBENS/Dyogen : Matthieu MUFFATO, Alexandra LOUIS, Hugues ROEST CROLLIUS
# Mail : hrc@ens.fr or alouis@biologie.ens.fr
# Licences GLP v3 and CeCILL v2

__doc__ = """
	Lit une serie d'arbres resultats de NHML.
	Extrait les GC ancestraux et les place pour chaque famille
"""

import sys
import utils.myTools
import utils.myPhylTree

# Arguments
arguments = utils.myTools.checkArgs( [("phylTree.conf",file), ("trees",file)], [("ancGenesFile",str,""), ("outputGCFile",str,"")], __doc__ )

# L'arbre phylogenetique
phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])

utils.myTools.mkDir(arguments["outputGCFile"])

# 1. On charge toutes les familles et on garde en leurs coordonnees
dicFamilles = {}
dicGC = {}
for esp in phylTree.allNames:
	nb = 0
	for gene in utils.myGenomes.Genome(arguments["ancGenesFile"] % phylTree.fileName[esp]):
		names = frozenset(gene.names)
		dicFamilles[names] = (esp, gene.beginning)
		nb += 1
	dicGC[esp] = [None] * nb

# 2. On charge les arbres, reconstruit les familles et range les GC
f = utils.myTools.myOpenFile(arguments["trees"], "r")
for arbre in f:

	# Lit les nb prochains caracteres de l'arbre
	def readStr(nb):
		global pos
		ret = arbre[pos:pos+nb]
		pos += nb
		return ret

	def readTree():
		#print >> sys.stderr, "new readTree (%d) *%s*" % (pos,arbre[pos:])
		if arbre[pos] == '(':
			#assert (readStr(1) == '(')
			readStr(1) # '('
			e1 = readTree()
			#print >> sys.stderr, "e1", e1
			#assert (readStr(2) == ', ')
			readStr(2) # ', '
			e2 = readTree()
			#print >> sys.stderr, "e2", e2
			#assert (readStr(1) == ')')
			readStr(1) # ')'
			elt = (e1,e2)
		else:
			x = arbre.index(' ', pos)
			elt = readStr(x-pos)
			#print >> sys.stderr, "elt", elt
			#assert (readStr(1) == ' ')
			readStr(1) # ' '
		
		x = arbre.find(':', pos)
		if x < 0:
			return (elt,None)

		gc = float(readStr(x-pos))
		#assert (readStr(1) == ':')
		(readStr(1) == ':')
		while arbre[pos] in "0123456789.":
			readStr(1)
		#print >> sys.stderr, "gc", gc
		return (elt,gc)
	pos = 0
	tree = readTree()
	print tree

	def writeGC( (subTree,GC) ):
		if type(subTree) == tuple:
			(e1,e2) = subTree
			names = writeGC(e1) + writeGC(e2)
		else:
			names = [subTree]
		#print frozenset(names), GC
		(esp,i) = dicFamilles.get( frozenset(names), (None,0) )
		if esp != None:
			print frozenset(names), esp, i, GC
			dicGC[esp][i] = GC
		return names

	writeGC(tree)

f.close()

for (esp,tab) in dicGC.iteritems():
	f = utils.myTools.myOpenFile(arguments["outputGCFile"] % phylTree.fileName[esp], "w")
	print >> f, utils.myFile.myTSV.printLine(tab, "\n")
	f.close()

