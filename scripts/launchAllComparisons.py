#! /users/ldog/muffato/python -OO

# Librairies
import os
import sys
import utils.myTools
import utils.myPhylTree

# Arguments
(noms_fichiers, options) = utils.myTools.checkArgs( ["phylTree.conf"], [("command",str,""), ("ancGenomesFile",str,""), ("genesFile",str,""), ("ancGenesFile",str,"")], "Lance la commande pour toutes les branches de l'arbre")

phylTree = utils.myPhylTree.PhylogeneticTree(noms_fichiers["phylTree.conf"])

if options["ancGenomesFile"] == "":
	options["ancGenomesFile"] = "%s"
if options["ancGenesFile"] == "":
	options["ancGenesFile"] = "%s"
if options["genesFile"] == "":
	options["genesFile"] = "%s"

def do(node):

	for (e,_) in phylTree.items.get(node, []):

		nodeF = options["ancGenomesFile"] % phylTree.fileName[node]
		nodeA = options["ancGenesFile"] % phylTree.fileName[node]
		if e in phylTree.listSpecies:
			eF = options["genesFile"] % phylTree.fileName[e]
		else:
			eF = options["ancGenomesFile"] % phylTree.fileName[e]
		
		command = options["command"].replace("=A=", nodeF).replace("=G=", nodeA).replace("=E=", eF)
		print >> sys.stderr, node, e
		print "%s\t%s\t" % (node,e),
		sys.stdout.flush()
		os.system(command)
		sys.stdout.flush()
		
		do(e)

do(phylTree.root)

