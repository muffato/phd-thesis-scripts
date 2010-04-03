#! /users/ldog/muffato/python

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

