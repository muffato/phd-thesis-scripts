#! /users/ldog/muffato/python -OO

import os
import sys
import utils.myTools
import utils.myPhylTree

(noms_fichiers, options) = utils.myTools.checkArgs( ["phylTree.conf"], [("mode",str,["species","branches"]), ("ancGenomesFile",str,""), ("genesFile",str,""), ("ancGenesFile",str,"")], "Compare les jeux de genes de chaque espece")

phylTree = utils.myPhylTree.PhylogeneticTree(noms_fichiers["phylTree.conf"])

def selectOutput(command, init):
	print >> sys.stderr, "Computing %s/%s ..." % (init[0], init[1]),
	(stdin,stdout,stderr) = os.popen3(command)
	stdin.close()
	res = [None] * 8
	for l in stdout:
		pass
	stdout.close()
	for l in stderr:
		t = l.split()
		if l.startswith("Evolution"):
			res[0] = t[6]
			res[1] = t[13]
		elif "perdus" in l:
			res[2] = t[3]
		elif "nouveaux" in l:
			res[3] = t[3]
		elif "duplications" in l:
			if res[4] == None:
				res[4] = t[2]
			else:
				res[5] = t[3]
		elif "cassures" in l:
			res[6] = t[4]
		elif "fusions" in l:
			res[7] = t[4]
			break
	stderr.close()
	print utils.myTools.printLine(init + res)
	print >> sys.stderr, "OK"


def do(node):

	for (e,l) in phylTree.items.get(node, []):

		nodeA = options["ancGenomesFile"] % phylTree.fileName[node]
		nodeG = options["ancGenesFile"] % phylTree.fileName[node]
		if e in phylTree.listSpecies:
			eE = options["genesFile"] % phylTree.fileName[e]
		else:
			eE = options["ancGenomesFile"] % phylTree.fileName[e]
		command = "/users/ldog/muffato/work/scripts/printOrthologousChr.py %s %s -orthologuesList=%s +includeScaffolds" % (nodeA,eE,nodeG)
		selectOutput(command, [node,e,l])
		do(e)

if options["mode"] == "species":

	for (e1,e2) in utils.myTools.myIterator.tupleOnStrictUpperList(phylTree.listSpecies):
		par = phylTree.dicParents.get(e1).get(e2)
		fE1 = options["genesFile"] % phylTree.fileName[e1]
		fE2 = options["genesFile"] % phylTree.fileName[e2]
		fA = options["ancGenesFile"] % phylTree.fileName[par]
		command = "/users/ldog/muffato/work/scripts/printOrthologousChr.py %s %s -orthologuesList=%s +includeScaffolds" % (fE1,fE2,fA)
		selectOutput(command, [e1,e2,2*phylTree.ages.get(par)])

else:
	do(phylTree.root)
