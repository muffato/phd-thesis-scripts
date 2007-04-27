#! /users/ldog/muffato/python -OO


# INITIALISATION #

# Librairies
import sys
import os

sys.path.append(os.environ['HOME'] + "/work/scripts/utils")
import myOrthos
import myTools
import myMaths

(noms_fichiers, options) = myTools.checkArgs(["GENOME_ANCESTRAL"], [("filtreOK",str,""), ("filtreNO",str,"")], "")

genesAnc = myOrthos.AncestralGenome(noms_fichiers[0], True, False)

lst = set([])
oui = set([x for x in options["filtreOK"]])
non = set([x for x in options["filtreNO"]])

for s in sys.stdin:
	t = set([])
	for g in s.split():
		if g in genesAnc.dicGenes:
			t.add(genesAnc.dicGenes[g])
	
	ss = [c for (c,i) in t]
	
	if len(oui) != 0 and len(oui.intersection(ss)) == 0:
		continue
	
	if len(non.intersection(ss)) != 0:
		continue
	
	for (c,i) in t:
		print " ".join(genesAnc.lstGenes[c][i][-1]),
	print s,

