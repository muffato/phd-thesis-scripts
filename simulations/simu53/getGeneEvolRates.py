#!/usr/bin/env python2
# Copyright (c) 2010 IBENS/Dyogen : Matthieu MUFFATO, Alexandra LOUIS, Hugues ROEST CROLLIUS
# Mail : hrc@ens.fr or alouis@biologie.ens.fr
# Licences GLP v3 and CeCILL v2

import sys
import subprocess

import utils.myFile
import utils.myTools
import utils.myPhylTree

arguments = utils.myTools.checkArgs( [("phylTree.conf",file)], [("genesFile",str,""), ("ancGenesFile",str,"")], "Compare les jeux de genes de chaque espece")

phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])


def do(node):

	nodeA = arguments["ancGenesFile"] % phylTree.fileName[node]
	for (e,age) in phylTree.items.get(node, []):
		if e in phylTree.listSpecies:
			eE = arguments["genesFile"] % phylTree.fileName[e]
		else:
			eE = arguments["ancGenesFile"] % phylTree.fileName[e]
		
		print >> sys.stderr, "Computing %s/%s ..." % (node,e),
		p = subprocess.Popen(["/users/ldog/muffato/workspace/scripts/printGenesDiff.py",nodeA,eE], shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE, close_fds=True)
		gains = 0
		pertes = 0
		dup1 = 0
		dup2 = 0
		for l in p.stdout:
			t = l.split()
			if t[0] == "+":
				gains += 1
			if t[0] == "-":
				pertes += 1
			if t[0] == "++":
				dup1 += 1
				dup2 += len(t)-1
		p.stdout.close()
		p.stderr.close()
		print utils.myFile.myTSV.printLine([node, e, age, gains, pertes, dup1, dup2])
		print >> sys.stderr, "OK"

		do(e)

do(phylTree.root)

