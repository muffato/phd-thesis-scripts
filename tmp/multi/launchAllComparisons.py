#!/usr/bin/env python2
# Copyright (c) 2010 IBENS/Dyogen : Matthieu MUFFATO, Alexandra LOUIS, Hugues ROEST CROLLIUS
# Mail : hrc@ens.fr or alouis@biologie.ens.fr
# Licences GLP v3 and CeCILL v2

import os
import sys

import utils.myFile
import utils.myTools
import utils.myPhylTree

arguments = utils.myTools.checkArgs( [("phylTree.conf",file)], [("ancGenomesFile",str,""), ("genesFile",str,""), ("ancGenesFile",str,"")], "Compare les jeux de genes de chaque espece")

phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])


def do(node):

	for (e,l) in phylTree.items.get(node, []):

		nodeA = arguments["ancGenomesFile"] % phylTree.fileName[node]
		nodeG = arguments["ancGenesFile"] % phylTree.fileName[node]
		if e in phylTree.listSpecies:
			eE = arguments["genesFile"] % phylTree.fileName[e]
		else:
			eE = arguments["ancGenomesFile"] % phylTree.fileName[e]
		command = "/users/ldog/muffato/workspace/scripts/printOrthologousChr.py %s %s -orthologuesList=%s +includeScaffolds" % (nodeA,eE,nodeG)
		command = "/users/ldog/muffato/workspace/txtSynteny/nbRearrang/cmpIntervals.py %s %s" % (nodeA,eE)
		#command = "/bin/ls %s %s " % (nodeA,eE)
		os.system(command)
		do(e)

do(phylTree.root)
