#! /users/ldog/muffato/python


##################
# INITIALISATION #
##################

# Librairies
import sys
import math
import random
import os
import utils.myOrthos
import utils.myTools
import utils.myMaths


########
# MAIN #
########


# Arguments
(noms_fichiers, options) = utils.myTools.checkArgs( \
	["genome1", "genome2", "genomeOutgroup", "orthologuesList"], \
	[("seuilLongueurMin", float, 0.1), ("seuilIdentiteMin", float, 33), ("seuilIdentiteMin2", float, 33), ("espece1", str, 'H'), ("espece2", str, 'C')], \
	"Reconstruit le genome de l'ancetre de 1 et 2 a partir de l'outgroup et des genes de cet ancetre" \
)


# 1. On lit tous les fichiers
genome1 = utils.myOrthos.loadGenome(noms_fichiers[0])
genome2 = utils.myOrthos.loadGenome(noms_fichiers[1])
genomeOutgroup = utils.myOrthos.loadGenome(noms_fichiers[2])
genesAnc = utils.myOrthos.AncestralGenome(noms_fichiers[3], False)


# La table d'association geneAncestral -> genetOugroup
assocGeneOutgroup = {}
for g in genesAnc.lstGenes[utils.myOrthos.AncestralGenome.defaultChr]:
	a = set([])
	for s in g.names:
		if s in genomeOutgroup.dicGenes:
			a.add(genomeOutgroup.dicGenes[s])

	if len(a) >= 1:
		assocGeneOutgroup[g.beginning] = a.pop()


# On construit les blocs ancestraux
# Pour chaque chromosome de l'outgroup, c'est la liste des genes du genome qui se trouvent dessus
# Ces genes sont regroupes selon leur chromosome dans l'espece
# Les genes sont des indices dans genesAnc

def makeAncRegions(genome, outgroup, genesAnc):
	blocsAncs = dict([(ca, dict([(c,set([])) for c in genome.lstChr])) for ca in genomeOutgroup.lstChr])
	for c in genome.lstChr:
		for g in genome.lstGenes[c]:
			s = g.names[0]
			if s not in genesAnc.dicGenes or s not in outgroup.dicGenes:
				continue
			(_,i) = genesAnc.dicGenes[s]
			(ca,_) = outgroup.dicGenes[s]
			blocsAncs[ca][c].add(i)
	return blocsAncs

blocsAnc1 = makeAncRegions(genome1, genomeOutgroup, genesAnc)
blocsAnc2 = makeAncRegions(genome2, genomeOutgroup, genesAnc)












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
			#print (100*len(zz))/len(xx)
			#print (100*len(zz))/len(yy)
			#continue
			#print float(len(zz))*100./float(len(xx)), float(len(zz))*100./float(len(yy))
			if options["seuilIdentiteMin"] >= 1:
				test = (len(zz) >= options["seuilIdentiteMin"])
			else:
				test = (len(zz) >= options["seuilIdentiteMin"]*len(xx)) and (len(zz) >= options["seuilIdentiteMin"]*len(yy))
			if test:
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
		aaa = xH | xP
		#print " ".join(myMaths.flatten([genesAnc.lstGenes[myOrthos.AncestralGenome.defaultChr][u].names for u in aaa]))
	#print >> sys.stderr, len(lst), "blocs chez l'ancetre %(espece1)s/%(espece2)s" % options

print >> sys.stderr, len(lstTout), "blocs chez l'ancetre %(espece1)s/%(espece2)s" % options
#print "}"
#sys.exit(0)

def rassembleSyntenies(g1, g2, genesAnc):

	# On construit les couleurs
	res = {}
	for c in g1.lstChr:

		res[c] = []
		for i in range(len(g1.lstGenes[c])):
		
			tg = myMaths.flatten([genesAnc.lstGenes[cc][ii].names for (cc,ii) in [genesAnc.dicGenes[g] for g in g1.lstGenes[c][i].names if g in genesAnc.dicGenes]])
			
			for g in tg:
				if g in g2.dicGenes:
					res[c].append( (g2.dicGenes[g][0],genesAnc.dicGenes[g][1]) )
					break

	res2 = []
	for c in g1.lstChr:
		
		last = ""
		curr = set([])
		
		for (col,tg) in res[c]:
			if col != last:
				if len(curr) != 0:
					res2.append(curr)
				last = col
				curr = set([])
			curr.add(tg)
		res2.append(curr)
	return res2
	
ens1 = rassembleSyntenies(genome1, genome2, genesAnc)
print >> sys.stderr, len(ens1), sum([len(a) for a in ens1])
ens2 = rassembleSyntenies(genome2, genome1, genesAnc)
print >> sys.stderr, len(ens2), sum([len(a) for a in ens2])

newEns = []
tout = set([])
i = 0
dic = {}
for x in ens1:
	continue
	for y in ens2:
		a = x & y
		if len(a) < 2:
			continue
		newEns.append(a)
		tout.update(a)
		#print >> sys.stderr,  " ", len(a)
		for u in a:
			dic[u] = i
		i += 1
		#print " ".join(myMaths.flatten([genesAnc.lstGenes[myOrthos.AncestralGenome.defaultChr][u].names for u in a]))

newEns = ens1 + ens2
print >> sys.stderr, len(newEns), sum([len(a) for a in newEns])
print >> sys.stderr, len(tout)

#print lstTout

lst = []
for (c,l1,l2,g1,g2) in lstTout:
	g = g1 | g2
	for (cc,x,y,gg1,gg2) in lst:
		gg = gg1 | gg2
		for u in newEns:
			z = g & u
			zz = gg & u
			test = False
			#print (100*len(z))/len(g)
			#print (100*len(zz))/len(gg)
			#print len(z)
			#print len(zz)
			if options["seuilIdentiteMin2"] >= 1:
				test = (len(z) >= options["seuilIdentiteMin2"])
				test = test and (len(zz) >= options["seuilIdentiteMin2"])
			else:
				test = (len(z) >= options["seuilIdentiteMin2"]*len(g))
				test = test and (len(zz) >= options["seuilIdentiteMin2"]*len(gg))
			if test:
				break
		#if (x.issubset(l1) or x.issuperset(l1)) and (y.issubset(l2) or y.issuperset(l2)):
		if test:
			x.update(l1)
			y.update(l2)
			cc.append(c)
			gg1.update(g1)
			gg2.update(g2)
			break
	else:
		lst.append( ([c],l1,l2,g1,g2) )
print >> sys.stderr, len(lst), "presque-chromosomes chez l'ancetre %(espece1)s/%(espece2)s" % options
#sys.exit(1)
res = []
for (l,l1,l2,g1,g2) in lst:
	a = g1 | g2
	if len(a) < 150:
		continue
	res.append( (len(a), l, l1, l2, a) )

print >> sys.stderr, len(res), "presque-chromosomes chez l'ancetre %(espece1)s/%(espece2)s" % options
res.sort()
res.reverse()
c = 0
nbTMP = 0
for (nb,l,l1,l2,a) in res:
	print >> sys.stderr, l, nb, l1,l2
	for i in a:
		#continue
		print chr(97+c), " ".join(genesAnc.lstGenes[myOrthos.AncestralGenome.defaultChr][i].names)
	nbTMP += nb
	c += 1
print >> sys.stderr, "Total nb genes:", nbTMP
