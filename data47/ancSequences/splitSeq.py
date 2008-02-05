#! /users/ldog/muffato/python -OO

import sys
import utils.myMaths
import utils.myTools
import utils.myPhylTree

(noms_fichiers, options) = utils.myTools.checkArgs( ["InputSeq"], [("outputSeq",str,"")], "!")

currID = None
currF = None
backup = None

f = utils.myTools.myOpenFile(noms_fichiers["InputSeq"], "r")
for l in f:
	if l.startswith("ID"):
		newID = int(l[:-1].split()[1])
		currF = backup
		if newID != currID:
			if currF != None:
				currF.close()
			print currID
			currID = newID
			backup = currF = utils.myTools.myOpenFile(options["outputSeq"] % currID, "w")
		#l = ""
	#elif l.startswith("SPECIES"):
	#	x = l[:-1].split('\t')[1].split()
	#	if len(x) == 1:
	#		currF = utils.myTools.null
	#	l = ("ID\t%d\n" % currID) + l
	#elif l.startswith("PROBA"):
	#	lst = [" %g" % float(x) for x in l[:-1].split()[1:]]
	#	l = l.split()[0] + ''.join(lst) + '\n'
	#if '\t' not in l:
	#	l = l.replace(' ', '\t', 1)
	currF.write(l)
	#print >> currF, l,
f.close()

