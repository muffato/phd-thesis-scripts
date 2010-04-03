#! /users/ldog/muffato/python

__doc__ = """
	Enleve les intervalles dont une des bornes s'est dupliquee
"""

import utils.myGenomes
import utils.myFile
import utils.myTools


arguments = utils.myTools.checkArgs( [("ancGenes",file), ("modernGenome",str), ("intervals",file)], [], __doc__)

ancGenes = utils.myGenomes.Genome(arguments["ancGenes"])
modernGenome = {}

def nbCopies(s):
	(c,i) = ancGenes.dicGenes[s]
	return len([x for x in ancGenes.lstGenes[c][i].names if x in modernGenome[species].dicGenes])

f = utils.myFile.openFile(arguments["intervals"], "r")
for l in f:
	if ("Set1:none" in l) or ("Set1:region" in l):
		species = l.split()[0]
		if species not in modernGenome:
			modernGenome[species] = utils.myGenomes.Genome(arguments["modernGenome"] % species)
		(f1,f2) = l.split("\t")[1].split()
		if (nbCopies(f1) == 1) and (nbCopies(f2) == 1):
			print l,
		else:
			print l.replace("Set2:pair", "Set2:duppair"),
	else:
		print l,
f.close()

