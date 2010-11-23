#!/usr/bin/env python2
# Copyright (c) 2010 IBENS/Dyogen : Matthieu MUFFATO, Alexandra LOUIS, Hugues ROEST CROLLIUS
# Mail : hrc@ens.fr or alouis@biologie.ens.fr
# Licences GLP v3 and CeCILL v2

import sys

import utils.myFile
import utils.myTools

arguments = utils.myTools.checkArgs( [("data",file), ("output",str)], [], "Cree les fichiers separes de sequences ancestrales" )

f = utils.myFile.openFile(arguments["data"], "r")
output = None
for l in f:
	if l.startswith("ID"):
		if output is not None:
			print >> output
			output.close()
		output = utils.myFile.openFile(arguments["output"] % l.split()[1], "a")
		continue
	#if l.startswith("NAME"):
	#	isanc = l.split()[1].startswith("NAME_")
	#if isanc:
	print >> output, l,
f.close()
if output is not None:
	print >> output
	output.close()

