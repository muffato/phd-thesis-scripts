#! /usr/bin/python2.4


# INITIALISATION #

# Librairies
import sys
import os

sys.path.append(os.environ['HOME'] + "/work/scripts/utils")
import myOrthos
import myTools
import myMaths

(noms_fichiers, options) = myTools.checkArgs(["GENOME_ANCESTRAL"], [], "")

genesAnc = myOrthos.AncestralGenome(noms_fichiers[0], True)

chroms = { "ALPHA":['n','y','z'], "BETA":['u','l','i'], "GAMMA":['j','m','k'], \
"DELTA":['g','c'], "EPSILON":['f','d'], "PHI":['e','p','q','s','w','x'] }

lstChroms = chroms.keys()
lstChroms.sort()

for s in sys.stdin:
	t = set([])
	ss = s.split()[1:]
	for g in ss:
		if g in genesAnc.dicGenes:
			t.add(genesAnc.dicGenes[g])

	if len(t) < 2:
		continue
		
	s = set([c for (c,i) in t])
	for x in lstChroms:
		if len(s.intersection(chroms[x])) != 0:
			res = x
			break
	else:
		continue
		
	print res,
	for (c,i) in t:
		print " ".join(genesAnc.lstGenes[c][i]),
	print
