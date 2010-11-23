#!/usr/bin/env python2
# Copyright (c) 2010 IBENS/Dyogen : Matthieu MUFFATO, Alexandra LOUIS, Hugues ROEST CROLLIUS
# Mail : hrc@ens.fr or alouis@biologie.ens.fr
# Licences GLP v3 and CeCILL v2

__doc__ = """
	Renvoie pour chaque gene ancestral la liste des branches sur lesquelles il a subi un evenement
"""

import re
import sys

import utils.myMaths
import utils.myTools
import utils.myGenomes
import utils.myPhylTree

arguments = utils.myTools.checkArgs( [("phylTree.conf",file), ("genesFile",str), ("ancGenesFile",str)], [("countDup",bool,True), ("countLoss",bool,True)], __doc__)

phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])

# Chargement des tous les fichiers
genes = {}
todo = {}
for e in phylTree.listSpecies:
	genes[e] = utils.myGenomes.Genome(arguments["genesFile"] % phylTree.fileName[e])
for a in phylTree.listAncestr:
	genes[a] = utils.myGenomes.Genome(arguments["ancGenesFile"] % phylTree.fileName[a])
	todo[a] = set(g.names[0] for g in genes[a])

allnames = set()
for a in phylTree.listAncestr:
	allnames.update(todo[a])
for s in sorted(phylTree.listAncestr, key=lambda a: phylTree.ages[a], reverse=True):
	continue
	new = set(todo[s]).difference(allnames)
	maybedup = allnames.difference(todo[s])
	l = [x for x in new if not any(x.startswith(y) for y in maybedup)]
	print >> sys.stderr, "On %s: initially %d trees = %d with different names, %d to add = %d with different names, %d new" % (s, len(allnames), len(maybedup), len(todo[s]), len(new), len(l))
	allnames.update(l)
evol = dict( (g,[None]*len(phylTree.indNames)) for g in allnames )
print >> sys.stderr, len(evol)

def countEvents(node, c, i):
	r = {}
	for (e,_) in phylTree.items.get(node, []):
		l = genes[e].getPosition(genes[node].lstGenes[c][i].names)
		#r[e] = int(len(l) > 1 if e in phylTree.lstEsp2X else len(l) != 1)
		r[e] = int( ((len(l) > 1) and arguments["countDup"]) or (((len(l) == 0) and (e not in phylTree.lstEsp2X)) and arguments["countLoss"]) )
		for (cc,ii) in l:
			r.update(countEvents(e, cc, ii))
	for (b,x) in r.iteritems():
		evol[gene][phylTree.indNames[b]] = x
	return r

for gene in allnames:
	def lookupBeginning(node):
		if node in phylTree.listAncestr:
			if gene in todo[node]:
				(c,i) = genes[node].dicGenes[gene]
				countEvents(node, c, i)
				evol[gene][phylTree.indNames[node]] = 1
			else:
				for (e,_) in phylTree.items[node]:
					lookupBeginning(e)
	lookupBeginning(phylTree.root)
	print gene, evol[gene]

