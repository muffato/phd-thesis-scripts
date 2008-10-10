#! /users/ldog/muffato/python

import utils.myTools
import utils.myPhylTree

arguments = utils.myTools.checkArgs([("phylTree.conf",file)],[],"")

# Renvoie l'arbre au format avec des parentheses
def convertToFlatFile(anc):

	a = phylTree.fileName[anc] # anc.replace(' ', '.')
	if anc in phylTree.listSpecies:
		return a
	else:
		return "(" + ",".join([convertToFlatFile(e) + ":" + str(l) for (e,l) in phylTree.items[anc]]) + ")%s|%d" % (a,phylTree.ages[anc])


phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])
print convertToFlatFile(phylTree.root), ";"

