#! /users/ldog/muffato/python -OO

import utils.myTools
import utils.myPhylTree

(noms,_) = utils.myTools.checkArgs(["phylTree.conf"],[],"")

# Renvoie l'arbre au format avec des parentheses
def convertToFlatFile(anc):

	a = anc.replace(' ', '.')
	if anc in phylTree.listSpecies:
		return a
	else:
		return "(" + ",".join([convertToFlatFile(e) + ":" + str(l) for (e,l) in phylTree.items[anc]]) + ")%s|%d" % (a,phylTree.ages[anc])


phylTree = utils.myPhylTree.PhylogeneticTree(noms["phylTree.conf"])
print convertToFlatFile(phylTree.root), ";"

