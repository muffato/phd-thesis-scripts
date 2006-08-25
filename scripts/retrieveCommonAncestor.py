#! /usr/bin/python2.4


##################
# INITIALISATION #
##################

# Librairies
import string
import sys
import math
import random
import os

sys.path.append(os.environ['HOME'] + "/work/scripts/utils")
import myOrthos
import myTools
import myMaths


########
# MAIN #
########


# Arguments
(noms_fichiers, options) = myTools.checkArgs( \
	["genesList.conf", "listeOrthologues"], \
	[], \
	"" \
)


geneBank = myOrthos.MyGeneBank(noms_fichiers[0])
orthologues = myOrthos.AncestralGenome(noms_fichiers[1], False)

genomeH = geneBank.dicEspeces['H']
genomeD = geneBank.dicEspeces['D']
genomeM = geneBank.dicEspeces['M']
genomeO = geneBank.dicEspeces['O']
genomeC = geneBank.dicEspeces['C']


def compareGenomes(genome1, genome2):

	lst = []
	for c1 in genome1.lstChr:
		
		l = []
		lastC2 = 0
		
		for (_,_,_,g1) in genome1.lstGenes[c1]:
		
			if g1 not in orthologues.dicGenes:
				continue
			(c,i) = orthologues.dicGenes[g1]
			t = orthologues.lstGenes[c][i]

			newC2 = [genome2.dicGenes[g][0] for g in t if g in genome2.dicGenes]
			
			if len(newC2) == 0:
				continue
			
			if newC2[0] != lastC2:
				if len(l) >= 1:
					lst.append(l)
				l = [i]
				lastC2 = newC2[0]
			else:
				l.append(i)

		lst.append(l)
	return lst

superLst = []
for x in geneBank.dicEspecesNonDup:
	for y in geneBank.dicEspecesNonDup:
		if x != y:
			l = compareGenomes(geneBank.dicEspeces[x], geneBank.dicEspeces[y])
			superLst.append(l)
			print x, y, len(l), sum([len(p) for p in l])

chrom = superLst[0]

for lst in superLst[1:]:
	
	dic = {}
	for i in range(len(chrom)):
		for g in chrom[i]:
			dic[g] = i

	for l in lst:
		l0 = []
		for g in l:
			if g in dic:
				i = dic[g]
				l0.extend(chrom[i])
				chrom[i] = []
			else:
				l0.append(g)
		chrom.append(l0)

newChrom = [x for x in chrom if len(x) > 0]
print len(dic), len(newChrom), len([x for x in chrom if len(x) >= 10]), len([x for x in chrom if len(x) >= 50]), len([x for x in chrom if len(x) >= 100])
l = [len(x) for x in newChrom]
print sum(l), l

sys.exit(0)











for cH in genomeH.lstChr:
	
	#cH = 1
	print cH,
	
	l = 0
	lastCD = 0
	lastCM = 0
	lastCO = 0
	
	for (_,_,_,gH) in genomeH.lstGenes[cH]:
	
		if gH not in orthologues.dicGenes:
			continue
		(c,i) = orthologues.dicGenes[gH]
		t = orthologues.lstGenes[c][i]

		newCD = [genomeD.dicGenes[g][0] for g in t if g in genomeD.dicGenes]
		if len(newCD) != 0:
			newCD = newCD[0]
		else:
			newCD = lastCD

		newCO = [genomeO.dicGenes[g][0] for g in t if g in genomeO.dicGenes]
		if len(newCO) != 0:
			newCO = newCO[0]
		else:
			newCO = lastCO

		newCM = [genomeM.dicGenes[g][0] for g in t if g in genomeM.dicGenes]
		if len(newCM) != 0:
			newCM = newCM[0]
		else:
			newCM = lastCM
		
		nbChgt = [lastCD != newCD, lastCO != newCO, lastCM != newCM].count(True)

		if nbChgt >= 1:
			if l >= 1:
				print "*%d (%s-%s-%s)" % (l, lastCD, lastCM, lastCO),
				#print "%s %s %s" % (lastCD, lastCM, lastCO)
			l = 1
			lastCD = newCD
			lastCO = newCO
			lastCM = newCM
		else:
			l += 1


	print "*%d (%s-%s-%s)" % (l, lastCD, lastCM, lastCO)
	#print "%s %s %s" % (lastCD, lastCM, lastCO)
	#sys.exit(1)



sys.exit(0)
if len(geneBank.dicEspeces) != 2:
	print >> sys.stderr, "Can't retrieve %s and %s in %s" % (options["espece1"], options["espece2"], noms_fichiers[0])
	sys.exit(1)
genomeAnc = myOrthos.AncestralGenome(noms_fichiers[1], True)

#options["seuilIdentiteMin"] = float(options["seuilIdentiteMin"]) / 100.

lstTout = []
# 2. On cree les blocs ancestraux tries et on extrait les diagonales
for chrAnc in genomeAnc.lstChr:

	n = len(genomeAnc.lstGenes[chrAnc])
	
	print >> sys.stderr, "Chromosome %s (%d genes) ... " % (chrAnc, n),
	
	if options["seuilLongueurMin"] >= 1:
		tailleMin = options["seuilLongueurMin"]
	else:
		tailleMin = options["seuilLongueurMin"] * n
	
	(blocsAnc, _) = genomeAnc.splitChr(geneBank, chrAnc)
	
	blocs1 = dict([(y,[x[2] for x in blocsAnc[options["espece1"]][y]]) for y in blocsAnc[options["espece1"]]])
	blocs2 = dict([(y,[x[2] for x in blocsAnc[options["espece2"]][y]]) for y in blocsAnc[options["espece2"]]])
	
	bl = {}
	bl2 = {}
	# On extrait les blocs avec s% d'identite
	for x in blocs1:
		xx = set(blocs1[x])
		if len(xx) < tailleMin:
			continue
		for y in blocs2:
			yy = set(blocs2[y])
			if len(yy) < tailleMin:
				continue
			zz = xx & yy
			print float(len(zz))*100./float(len(xx)), float(len(zz))*100./float(len(yy))
			if len(zz) > options["seuilIdentiteMin"]*len(xx) and len(zz) > options["seuilIdentiteMin"]*len(yy):
			#if len(zz) > options["seuilIdentiteMin"]:
				if x in bl:
					bl[x].add(y)
				else:
					bl[x] = set([y])
				if y in bl2:
					bl2[y].add(x)
				else:
					bl2[y] = set([x])
	
	
	# On rassemble les liens	
	lst = []
	while len(bl) > 0:
		(c,v) = bl.items()[0]
		l1 = set([c])
		l2 = v
		for j in range(100):
			for x in l1:
				l2.update(bl[x])
			for y in l2:
				l1.update(bl2[y])
		for x in l1:
			bl.pop(x)
		lst.append( (chrAnc, l1,l2) )
		print >> sys.stderr, l1, l2

	for j in range(00):
		# On place les blocs vides
		for x in blocs1:
			s = 0
			best = -1
			xx = set(blocs1[x])
			for (_,l1,l2) in lst:
				if x in l1:
					break
				r = set([])
				for y in l2:
					r.update(xx & set(blocs2[y]))
				if len(r) > s:
					s = len(r)
					best = l1
			else:
				if best != -1:
					print >> sys.stderr, "Ajout de", x, "H a", best
					best.add(x)
				elif j == 99:
					print >> sys.stderr, "Rien trouve pour", x, "H", blocs1[x]

		# On place les blocs vides
		for y in blocs2:
			s = 0
			best = -1
			yy = set(blocs2[y])
			for (_,l1,l2) in lst:
				if y in l2:
					break
				r = set([])
				for x in l1:
					r.update(yy & set(blocs1[x]))
				if len(r) > s:
					s = len(r)
					best = l2
			else:
				if best != -1:
					print >> sys.stderr, "Ajout de", y, "P a", best
					best.add(y)
				elif j == 99:
					print >> sys.stderr, "Rien trouve pour", y, "P", blocs2[y]
				
	for (c,l1,l2) in lst:
		xH = set([])
		for cH in l1:
			xH.update(blocs1[cH])
		xP = set([])
		for cP in l2:
			xP.update(blocs2[cP])

		lstTout.append( (c,l1,l2,xH,xP) )
	#lstTout.extend(lst)
	#print lst
	print >> sys.stderr, len(lst), "blocs chez l'ancetre %(espece1)s/%(espece2)s" % options

#print lstTout
print >> sys.stderr, len(lstTout)
print >> sys.stderr, len(lstTout), "blocs chez l'ancetre %(espece1)s/%(espece2)s" % options

lst = []
for (c,l1,l2,g1,g2) in lstTout:
	for (cc,x,y,gg1,gg2) in lst:
		if (x.issubset(l1) or x.issuperset(l1)) and (y.issubset(l2) or y.issuperset(l2)):
			x.update(l1)
			y.update(l2)
			cc.append(c)
			gg1.append(g1)
			gg2.append(g2)
			break
	else:
		lst.append( ([c],l1,l2,[g1],[g2]) )
print >> sys.stderr, len(lst), "chromosomes chez l'ancetre %(espece1)s/%(espece2)s" % options

cc = 0
for (l,l1,l2,g1,g2) in lst:
	print >> sys.stderr, l,l1,l2
	continue
	for i in range(len(l)):
		c = l[i]
		g = g1[i] | g2[i]
		for x in g:
			print chr(97+cc),
			for y in genomeAnc.lstGenes[c][x]:
				print y,
			print
	cc += 1

