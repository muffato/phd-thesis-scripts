#! /users/ldog/muffato/python -OO

__doc__ = """
	Lit les arbres de proteines et cree les fichiers de genes ancestraux
"""


# Librairies
import os
import sys
import utils.myTools
import utils.myPhylTree
import utils.myProteinTree

# Arguments
arguments = utils.myTools.checkArgs( [("phylTree.conf",file), ("proteinTree",file)], [("OUT.ancGenesFile",str,"")], __doc__ )

phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])


####################################################################
# On considere que les duplications 'dubious' ne sont pas valables #
####################################################################
def isDuplicatedNode(inf):
	return (inf['Duplication'] != 0) and ('dubious_duplication' not in inf)


###########################################
# Sauvegarde toutes les familles de genes #
###########################################
def extractGeneFamilies(node, baseName, previousAnc, lastWrittenAnc):

	newAnc = info[node]['taxon_name']
	(toWrite,newLastWritten,isroot) = utils.myProteinTree.getIntermediateAnc(phylTree, previousAnc, lastWrittenAnc, newAnc, isDuplicatedNode(info[node]))

	info[node]['family_name'] = baseName

	# Les genes des descendants
	if node in data:
		allGenes = []
		for (i,(g,d)) in enumerate(data[node]):
			if isDuplicatedNode(info[node]):
				newName = baseName + (".%d" % (i+1))
			else:
				newName = baseName
			allGenes.extend( extractGeneFamilies(g, newName, newAnc, newLastWritten) )
	else:
		allGenes = [ info[node]["gene_name"] ]

	for a in toWrite:
		geneFamilies[a].append( [baseName] + allGenes )

	return allGenes


print >> sys.stderr, "Mise en forme des arbres ...",
nb = 0
geneFamilies = utils.myTools.defaultdict(list)
for (r,data,info) in utils.myProteinTree.loadTree(arguments["proteinTree"]):
	nb += 1
	extractGeneFamilies(r, "FAM%d" % nb, None, None)
	utils.myProteinTree.printTree(sys.stdout, data, info, r)
print >> sys.stderr, "%d arbres OK" % nb


utils.myTools.mkDir(arguments["OUT.ancGenesFile"])
for (anc,lst) in geneFamilies.iteritems():
	print >> sys.stderr, "Ecriture des familles de %s ..." % anc,
	f = utils.myTools.myOpenFile(arguments["OUT.ancGenesFile"] % phylTree.fileName[anc], "w")
	for gg in lst:
		print >> f, " ".join(gg)
	f.close()
	print >> sys.stderr, len(lst), "OK"

