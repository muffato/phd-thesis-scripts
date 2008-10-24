#! /users/ldog/muffato/python

__doc__ = """
	Separe les valeurs d'un fichier de scores selon si les diagonales etaient sur le meme chromosome ancestral
"""

import sys
import utils.myTools
import utils.myPhylTree


# Argument
arguments = utils.myTools.checkArgs( [("sim",int)], [("IN.genomeFile",str,""), ("IN.scoresFile",str,""), ("IN.diagsFile",str,"")], __doc__ )

def doSim(sim):

	# On charge le genome pour savoir le chromosome de chaque diagonale
	chromAssoc = {}
	for (c,diag,_) in utils.myFile.myTSV.readTabular(arguments["IN.genomeFile"] % sim, [int,str,str]):
		for i in diag.split():
			chromAssoc[int(i)] = c

	# On charge les diagonales
	numDiagsOK = set()
	for (i,(_,_,diag,_,_,_)) in enumerate(utils.myFile.myTSV.readTabular(arguments["IN.diagsFile"] % sim, [str,int,str,str,str,str])):
		if len(set([chromAssoc[int(g)] for g in diag.split()])) == 1:
			numDiagsOK.add(i)
			maxdiag = i
	
	# On charge le fichier de scores et on selectionne les bons scores des mauvais
	for (i1,i2,score) in utils.myFile.myTSV.readTabular(arguments["IN.scoresFile"] % sim, [int,int,str], ' '):
		if (i1 not in numDiagsOK) or (i2 not in numDiagsOK):
			continue
		if chromAssoc[i1] == chromAssoc[i2]:
			print score
		else:
			print >> sys.stderr, score

doSim(arguments["sim"])

