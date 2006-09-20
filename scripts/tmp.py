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
	["genome1", "genome2", "genomeOutgroup", "orthologuesList"], \
	[("seuilLongueurMin", float, 0.1), ("seuilIdentiteMin", float, 33), ("espece1", str, 'H'), ("espece2", str, 'C')], \
	"Reconstruit le genome de l'ancetre de 1 et 2 a partir de l'outgroup et des genes de cet ancetre" \
)


# 1. On lit tous les fichiers
genome1 = myOrthos.loadGenome(noms_fichiers[0])
genome2 = myOrthos.loadGenome(noms_fichiers[1])
genomeOutgroup = myOrthos.loadGenome(noms_fichiers[2])
lstGenes = myOrthos.AncestralGenome(noms_fichiers[3], False)


# La table d'association geneAncestral -> genetOugroup
assocGeneOutgroup = {}
for (i,_,_,ts)in lstGenes.lstGenes[myOrthos.AncestralGenome.defaultChr]:
	a = set([])
	for s in ts:
		if s in genomeOutgroup.dicGenes:
			a.add(genomeOutgroup.dicGenes[s])

	if len(a) >= 1:
		assocGeneOutgroup[i] = a.pop()


# On construit les blocs ancestraux
# Pour chaque chromosome de l'outgroup, c'est la liste des genes du genome qui se trouvent dessus
# Ces genes sont regroupes selon leur chromosome dans l'espece
# Les genes sont des indices dans lstGenes

def makeAncRegions(genome, outgroup, lstGenes):
	blocsAncs = dict([(ca, dict([(c,set([])) for c in genome.lstChr])) for ca in genomeOutgroup.lstChr])
	for c in genome.lstChr:
		for (_,_,_,(s,)) in genome.lstGenes[c]:
			if s not in lstGenes.dicGenes or s not in outgroup.dicGenes:
				continue
			(_,i) = lstGenes.dicGenes[s]
			(ca,_) = outgroup.dicGenes[s]
			blocsAncs[ca][c].add(i)
	return blocsAncs

blocsAnc1 = makeAncRegions(genome1, genomeOutgroup, lstGenes)
blocsAnc2 = makeAncRegions(genome2, genomeOutgroup, lstGenes)












#sys.exit(0)
lstTout = []
# 2. On cree les blocs ancestraux tries et on extrait les diagonales

#print "graph {"
for chrAnc in genomeOutgroup.lstChr:

	n = len(genomeOutgroup.lstGenes[chrAnc])
	
	print >> sys.stderr, "Chromosome %s (%d genes) ... " % (chrAnc, n),
	
	if options["seuilLongueurMin"] >= 1:
		tailleMin = options["seuilLongueurMin"]
	else:
		tailleMin = options["seuilLongueurMin"] * n
	
	blocs1 = blocsAnc1[chrAnc]
	blocs2 = blocsAnc2[chrAnc]
	
	bl = {}
	bl2 = {}
	#print "{"
	# On extrait les blocs avec s% d'identite
	for x in blocs1:
		xx = blocs1[x]
		if len(xx) < tailleMin:
			continue
		for y in blocs2:
			yy = blocs2[y]
			if len(yy) < tailleMin:
				continue
			zz = xx & yy
			if len(zz) == 0:
				continue
			#print len(zz)
			#continue
			#print float(len(zz))*100./float(len(xx)), float(len(zz))*100./float(len(yy))
			#if len(zz) >= options["seuilIdentiteMin"]*len(xx) and len(zz) >= options["seuilIdentiteMin"]*len(yy):
			if len(zz) >= options["seuilIdentiteMin"]:
				#print "\"%s.H.%s (%d)\" -- \"%s.P.%s (%d)\" [label=\"%d\"]" % (chrAnc,x,len(xx), chrAnc,y,len(yy), len(zz))
				if x in bl:
					bl[x].add(y)
				else:
					bl[x] = set([y])
				if y in bl2:
					bl2[y].add(x)
				else:
					bl2[y] = set([x])
	
	#print "}"
	#break
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
				#print >> sys.stderr, "Rien trouve pour", x, "H", blocs1[x]
				#continue
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
				#print >> sys.stderr, "Rien trouve pour", y, "P", blocs2[y]
				#continue
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
	#print >> sys.stderr, len(lst), "blocs chez l'ancetre %(espece1)s/%(espece2)s" % options

print >> sys.stderr, len(lstTout), "blocs chez l'ancetre %(espece1)s/%(espece2)s" % options
#print "}"
#sys.exit(0)

lst = []
for (c,l1,l2,g1,g2) in lstTout:
	for (cc,x,y,gg1,gg2) in lst:
		#if (x.issubset(l1) or x.issuperset(l1)) and (y.issubset(l2) or y.issuperset(l2)):
		if False:
			x.update(l1)
			y.update(l2)
			cc.append(c)
			gg1.append(g1)
			gg2.append(g2)
			break
	else:
		lst.append( ([c],l1,l2,[g1],[g2]) )
print >> sys.stderr, len(lst), "chromosomes chez l'ancetre %(espece1)s/%(espece2)s" % options
#sys.exit(1)
cc = 0
for (l,l1,l2,g1,g2) in lst:
	nb = 0
	for i in range(len(l)):
		c = l[i]
		g = g1[i] | g2[i]
		for x in g:
			nb += 1
			#print chr(97+cc),
			#for y in genomeOutgroup.lstGenes[c][x]:
			#	print y,
			print x
	print >> sys.stderr, l, nb, l1,l2
	cc += 1

