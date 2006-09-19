#! /usr/bin/python2.4


##################
# INITIALISATION #
##################

# Librairies
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
	["genesList.conf", "genomeOutgroup"], \
	[("seuilLongueurMin", float, 0.1), ("seuilIdentiteMin", float, 33), ("espece1", str, 'H'), ("espece2", str, 'C')], \
	"" \
)


# 1. On lit tous les fichiers
geneBank = myOrthos.GeneBank(noms_fichiers[0], [options["espece1"], options["espece2"]])
#if len(geneBank.dicEspeces) != 2:
#	print >> sys.stderr, "Can't retrieve %s and %s in %s" % (options["espece1"], options["espece2"], noms_fichiers[0])
#	sys.exit(1)
genomeAnc = myOrthos.loadGenome(noms_fichiers[1])

#options["seuilIdentiteMin"] = float(options["seuilIdentiteMin"]) / 100.

lstTout = []
# 2. On cree les blocs ancestraux tries et on extrait les diagonales
for chrAnc in genomeAnc.lstChr:

	n = len(genomeAnc.lstGenes[chrAnc])
	
	print >> sys.stderr, "Chromosome %s (%d genes) ... " % (chrAnc, n),
	
	#if options["seuilLongueurMin"] >= 1:
	#	tailleMin = options["seuilLongueurMin"]
	#else:
	#	tailleMin = options["seuilLongueurMin"] * n
	
	(blocsAnc, _) = genomeAnc.splitChr(geneBank, chrAnc)

	print dict([ (x,len(blocsAnc['C'][x])) for x in blocsAnc['C']])
	continue
	sys.exit(0)
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
			#print float(len(zz))*100./float(len(xx)), float(len(zz))*100./float(len(yy))
			#if len(zz) > options["seuilIdentiteMin"]*len(xx) and len(zz) > options["seuilIdentiteMin"]*len(yy):
			if len(zz) > options["seuilIdentiteMin"]:
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

	for j in range(100):
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
	print >> sys.stderr, len(lst), "blocs chez l'ancetre %(espece1)s/%(espece2)s" % options

sys.exit(0)
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
	nb = 0
	for i in range(len(l)):
		c = l[i]
		g = g1[i] | g2[i]
		for x in g:
			nb += 1
			print chr(97+cc),
			for y in genomeAnc.lstGenes[c][x]:
				print y,
			print
	print >> sys.stderr, l, nb, l1,l2
	cc += 1

