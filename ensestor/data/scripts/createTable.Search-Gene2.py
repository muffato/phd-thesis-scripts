#!/usr/bin/env python2
# Copyright (c) 2010 IBENS/Dyogen : Matthieu MUFFATO, Alexandra LOUIS, Hugues ROEST CROLLIUS
# Mail : hrc@ens.fr or alouis@biologie.ens.fr
# Licences GLP v3 and CeCILL v2

__doc__ = """
	Cree la base Genomicus (table Search & Gene modifiee pour les noms & descriptions)
	Necessite Gene pour en creer la nouvelle version
	Lit les fichiers les listes des noms & descriptions
"""

import sys

import utils.myFile
import utils.myTools
import utils.myPhylTree


arguments = utils.myTools.checkArgs(
	[("phylTree.conf",file), ("geneTable",file)],
	[("namesFile",str,""), ("descrFile",str,""), ("output",str,"")],
	__doc__
)

phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])

dicNames = {}
dicDescr = {}
for esp in phylTree.listSpecies:
	print >> sys.stderr, "Loading names & descr of", esp, "...",
	for t in utils.myFile.myTSV.readTabular(arguments["namesFile"] % phylTree.fileName[esp], (str,str)):
		if len(t[1]) > 0:
			assert t[0] not in dicNames
			dicNames[t[0]] = t[1]
	for t in utils.myFile.myTSV.readTabular(arguments["descrFile"] % phylTree.fileName[esp], (str,str)):
		if len(t[1]) > 0:
			assert t[0] not in dicDescr
			dicDescr[t[0]] = t[1]
	print >> sys.stderr, "OK"

print >> sys.stderr, "Transforming Gene ...",
fG = utils.myFile.openFile(arguments["output"] % "Gene", "w")
fS = utils.myFile.openFile(arguments["output"] % "Search", "w")
for t in utils.myFile.myTSV.readTabular(arguments["geneTable"], [str]*10):
	line = list(t)
	print >> fS, utils.myFile.myTSV.MySQLFileWriter([line[2], line[0]])
	if line[2] in dicDescr:
		line[9] = dicDescr[line[2]]
	if line[2] in dicNames:
		print >> fS, utils.myFile.myTSV.MySQLFileWriter([dicNames[line[2]], line[0]])
		line[2] = "%s (%s)" % (dicNames[line[2]],line[2])
	print >> fG, utils.myFile.myTSV.MySQLFileWriter(line)
fG.close()
fS.close()
print >> sys.stderr, "OK"

