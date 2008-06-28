#! /users/ldog/muffato/python -OO

__doc__ = """
	Lit les arbres de proteines et cree les arbres au format Newick
"""


# Librairies
import os
import sys
import utils.myTools
import utils.myPhylTree
import utils.myProteinTree

# Arguments
arguments = utils.myTools.checkArgs( [("phylTree.conf",file), ("proteinTree",file)], [("treeAncestor",str,"")], __doc__ )

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

	if isroot:
		trueRoots.append(node)
		global nbA
		nbA += 1
		#baseName = "FAM%d" % nbA
		return
	else:
		for (g,_) in data.get(node,[]):
			extractGeneFamilies(g, baseName, newAnc, newLastWritten)
	
	return
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
		allGenes = [ links[node][0] ]

	for a in toWrite:
		geneFamilies[a].append( [baseName] + allGenes )

	return allGenes


if arguments["treeAncestor"] == "":
	treeAncestor = phylTree.root
else:
	treeAncestor = arguments["treeAncestor"]

print >> sys.stderr, "Mise en forme des arbres ...",
nb = 0
nbA = 0
for (r,data,info) in utils.myProteinTree.loadTree(arguments["proteinTree"]):
	trueRoots = []
	nb += 1
	extractGeneFamilies(r, "NONAME", None, None)
	for r in trueRoots:
		utils.myProteinTree.printTree(sys.stdout, data, info, r)
	continue
	# Restreint les arbres a la racine que l'on demande
	todo = [r]
	while len(todo) > 0:
		r = todo.pop()
		if phylTree.isChildOf(info[r]['taxon_name'], treeAncestor):
			# Extrait les genes ancestraux et renvoie les vraies racines des arbres
			trueRoots = []
			nb += 1
			extractGeneFamilies(r, "NONAME", None, None)
			for r in trueRoots:
				utils.myProteinTree.printTree(sys.stdout, data, info, r)
				#printNewickTree(ft5, r)
		else:
			todo.extend( [x for (x,_)  in data[r]] )
print >> sys.stderr, "%d (%d) arbres OK" % (nb,nbA)
#ft5.close()
sys.exit(0)

utils.myTools.mkDir(arguments["OUT.ancGenesFile"])
for (anc,lst) in geneFamilies.iteritems():
	print >> sys.stderr, "Ecriture des familles de %s ..." % anc,
	f = utils.myTools.myOpenFile(arguments["OUT.ancGenesFile"] % phylTree.fileName[anc], "w")
	for gg in lst:
		print >> f, " ".join(gg)
	f.close()
	print >> sys.stderr, len(lst), "OK"

