#! /users/ldog/muffato/python -OO


# INITIALISATION #

# Librairies
import sys
import operator
import utils.myTools
import utils.myGenomes

(noms_fichiers, options) = utils.myTools.checkArgs(["GENOME_ANCESTRAL"], [("addOtherGenes",bool,False)], "")

genesAnc = utils.myGenomes.loadGenome(noms_fichiers["GENOME_ANCESTRAL"])


chromsMain = { "ALPHA":['n','y','z'], "BETA":['u','l','i'], "GAMMA":['j','m','k'], \
"DELTA":['g','c'], "EPSILON":['f','d'], "PHI":['w','x'] }

chromsAltern = { "ALPHA":['s','t','k'], "BETA":['o','v','r'], "GAMMA":['q'], \
"DELTA":['q','b','a','p'], "EPSILON":['e','p','v'], "PHI":['e','s','q'] }



chromsMain = { "ALPHA":[22,32,8,9], "BETA":[7,29,17], "GAMMA":[21,6,20], \
"DELTA":[4,19,13,10], "EPSILON":[28,11], "PHI":[25,5] }

chromsAltern = { "ALPHA":[], "BETA":[1], "GAMMA":[1], \
"DELTA":[], "EPSILON":[3,23], "PHI":[3,23] }


newGenome = dict( [(c,[]) for c in genesAnc.lstChr] )

for s in sys.stdin:
	t = set([])
	ss = s.split()
	for g in ss:
		if g in genesAnc.dicGenes:
			t.add(genesAnc.dicGenes[g])

	if len(t) < 2:
		continue

	score = dict.fromkeys(chromsMain, 0)
	for (c,_) in t:
		c = int(c)
		for x in chromsMain:
			if c in chromsMain[x]:
				score[x] += 2
			elif c in chromsAltern[x]:
				score[x] += 1
	
	res = sorted(score.items(), key=operator.itemgetter(1), reverse=True)[0][0]
		
	if options["addOtherGenes"]:
		print res,s,
	for (c,i) in t:
		newGenome[c].append( (i,res) )
		if options["addOtherGenes"]:
			print " ".join(genesAnc.lstGenes[c][i].names),
if not options["addOtherGenes"]:
	sys.exit(0)

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




