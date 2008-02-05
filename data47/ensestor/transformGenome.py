#! /users/ldog/muffato/python -OO


# Librairies
import sys
import utils.myGenomes
import utils.myTools
import utils.myMaths
import utils.myPhylTree

# Arguments
############
(noms_fichiers, options) = utils.myTools.checkArgs( ["blocsList", "genome"], [], "transforme un genome avec les indices des blocs")

lstBlocs = {}
f = utils.myTools.myOpenFile(noms_fichiers["blocsList"], "r")
for l in f:
	(i,_,bloc) = l[:-1].partition('\t')
	lstBlocs[bloc] = i
f.close()

f = utils.myTools.myOpenFile(noms_fichiers["genome"], "r")
for l in f:
	t = l[:-1].split('\t')
	bloc = t[2] + "\t" + t[3]
	if len(t[2].split()) == 1:
		t[1] = "0"
	print t[0] + "\t" + t[1] + "\t" + lstBlocs[bloc]
f.close()

