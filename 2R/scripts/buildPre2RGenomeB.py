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

newGenome = dict( [(c,[]) for c in genesAnc.lstChr] )

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
		
	for (c,i) in t:
		newGenome[c].append( (i,res) )


for c in genesAnc.lstChr:
	lastP = 0
	lastV = 0
	newGenome[c].sort()
	for (i,res) in newGenome[c]:
		if lastP != 0:
			for j in range(lastP, (lastP+i)/2):
				print lastV, " ".join(genesAnc.lstGenes[c][j].names)
		else:
			lastP = -i
		for j in range((lastP+i)/2, i):
			print res, " ".join(genesAnc.lstGenes[c][j].names)
		lastP = i
		lastV = res
	for j in range(lastP, len(genesAnc.lstGenes[c])):
		print lastV, " ".join(genesAnc.lstGenes[c][j].names)



