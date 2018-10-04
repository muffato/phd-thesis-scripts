#!/usr/bin/env python2

__doc__ = """
	Lit les arbres de proteines et cree les fichiers de genes ancestraux
"""

import sys
import collections

import utils.myFile
import utils.myTools
import utils.myPhylTree
import utils.myProteinTree

# Arguments
arguments = utils.myTools.checkArgs( [("phylTree.conf",file), ("proteinTree",file)], [("OUT.ancGenesFile",str,"")], __doc__ )

phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])

dupCount = collections.defaultdict(int)
def futureName(name, dup):
	if dup >= 2:
		dupCount[name] += 1
		return name + utils.myProteinTree.getDupSuffix(dupCount[name], False)
	else:
		return name

###########################################
# Sauvegarde toutes les familles de genes #
###########################################
def extractGeneFamilies(node, baseName, previousAnc, lastWrittenAnc):

	newAnc = info[node]['taxon_name']
	(toWrite,newLastWritten,isroot) = utils.myProteinTree.getIntermediateAnc(phylTree, previousAnc, lastWrittenAnc, newAnc, info[node]['Duplication'] >= 2)

	info[node]['family_name'] = baseName

	# Les genes des descendants
	if node in data:
		allGenes = []
		for (g,_) in data[node]:
			allGenes.extend( extractGeneFamilies(g, futureName(baseName, info[node]['Duplication']), newAnc, newLastWritten) )
	else:
		allGenes = [ info[node]["gene_name"] ]

	for a in toWrite:
		geneFamilies[a].append( [baseName] + allGenes )

	return allGenes


print >> sys.stderr, "Mise en forme des arbres ...",
nb = 0
geneFamilies = collections.defaultdict(list)
for (r,data,info) in utils.myProteinTree.loadTree(arguments["proteinTree"]):
	nb += 1
	extractGeneFamilies(r, info[r]["tree_name"], None, None)
	utils.myProteinTree.printTree(sys.stdout, data, info, r)
print >> sys.stderr, "%d arbres OK" % nb


for (anc,lst) in geneFamilies.iteritems():
	print >> sys.stderr, "Ecriture des familles de %s ..." % anc,
	f = utils.myFile.openFile(arguments["OUT.ancGenesFile"] % phylTree.fileName[anc], "w")
	for gg in lst:
		print >> f, " ".join(gg)
	f.close()
	print >> sys.stderr, len(lst), "OK"

