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

genesAnc = myOrthos.AncestralGenome(noms_fichiers[0], True, False)

chromsMain = { "ALPHA":['n','y','z'], "BETA":['u','l','i'], "GAMMA":['j','m','k'], \
"DELTA":['g','c'], "EPSILON":['f','d'], "PHI":['w','x'] }

chromsAltern = { "ALPHA":['s','t','k'], "BETA":['o','v','r'], "GAMMA":['q'], \
"DELTA":['q','b','a','p'], "EPSILON":['e','p','v'], "PHI":['e','s','q'] }

lstChroms = chromsMain.keys()
lstChroms.sort()

for s in sys.stdin:
	t = set([])
	ss = s.split()
	for g in ss:
		if g in genesAnc.dicGenes:
			t.add(genesAnc.dicGenes[g])

	if len(t) < 2:
		continue

	score = dict([(c,0) for c in lstChroms])
	for (c,_) in t:
		for x in lstChroms:
			if c in chromsMain[x]:
				score[x] += 2
			elif c in chromsAltern[x]:
				score[x] += 1
	
	res = myMaths.sortDict(score)[0]
		
	print res,
	for (c,i) in t:
		print " ".join(genesAnc.lstGenes[c][i].names),
	print s,

