#!/usr/bin/env python2

__doc__ = """
	Lit les CNE et transforme les positions genomiques en positions indexees
"""

import sys
import bisect
import operator
import itertools
import collections

import utils.myTools
import utils.myGenomes

arguments = utils.myTools.checkArgs( \
	[("cneFile",file), ("intervalFile",file), ("strandFile",file)], \
	[("windowSize",int,10)], \
	__doc__ \
)



# Chargement du genome
print >> sys.stderr, "Loading intervals ...",
genome = collections.defaultdict(list)
f = utils.myFile.openFile(arguments["intervalFile"], "r")
for l in f:
	l = eval(l)
	genome[l[0]].append(l[1:])
f.close()

dicGene2Pos = collections.defaultdict(list)
for c in genome:
	assert genome[c] == sorted(genome[c])
	for (x,y) in utils.myTools.myIterator.slidingTuple(genome[c]):
		assert y[0] == x[1]+1
	for (i,x) in enumerate(genome[c]):
		for g in x[3]:
			dicGene2Pos[g].append( (c,i) )
print >> sys.stderr, "OK"

# Orientations des genes
print >> sys.stderr, "Loading strands ...",
strand = {}
f = utils.myFile.openFile(arguments["strandFile"], "r")
for l in f:
	t = l.split()
	assert len(t) == 2
	strand[t[1]] = int(t[0])
f.close()
print >> sys.stderr, "OK"

# Chargement des diagonales
print >> sys.stderr, "Loading pseudo-synteny blocks ...",
diags = collections.defaultdict(list)
fin = {}
for x in dicGene2Pos:

	if isinstance(x, frozenset):
		continue

	(cmin,imin) = min(dicGene2Pos[x])
	(cmax,imax) = max(dicGene2Pos[x])
	assert cmin == cmax
	c = cmin
	xmin = genome[c][imin][0]
	xmax = genome[c][imax][1]
	diags[c].append( xmin )
	fin[(c,xmin)] = (xmax,strand[x])

for x in diags.itervalues():
	x.sort()
print >> sys.stderr, "OK"


step = arguments["windowSize"] * 1000
def addInterv(count, c, xmin, xmax):

	i = bisect.bisect_right(diags[c], xmin )
	
	if (i >= 1) and (xmin < fin[(c,diags[c][i-1])][0]):
		# dans diagonale precedente
		choix1 = diags[c][i-1]
		choix2 = fin[(c,diags[c][i-1])]
		if choix2[1] > 0:
			choix2 = sys.maxint
		else:
			choix1 = -sys.maxint
			choix2 = choix2[0]
		into = True
	else:
		# entre diagonales
		choix1 = fin[(c,diags[c][i-1])] if i >= 1 else -sys.maxint
		if choix1[1] > 0:
			choix1 = -sys.maxint
		else:
			choix1 = choix1[0]
		choix2 = diags[c][i] if i < len(diags[c]) else sys.maxint
		test = fin[(c,diags[c][i])][1] if i < len(diags[c]) else 0
		if test < 0:
			choix2 = sys.maxint
		if (choix1 == -sys.maxint) and (choix2 == sys.maxint):
			return
		into = False

	assert choix1 <= xmin <= xmax <= choix2, (choix1,xmin,xmax,choix2)

	global to
	to = 0
	nb = 0
	def add(size):
		if into:
			n = nb
		else:
			n = -1-nb
		count[n] += size
		global to
		to = to + size
	l = xmax-xmin+1
	while xmin <= xmax:
		
		choix2 -= step
		choix1 += step
		
		#print >> sys.stderr, l, xmin, xmax, xmax-xmin+1, to, xmax-xmin+1+to, choix1, choix2
		if choix1 >= choix2:
			# L'intervalle entier est au meme niveau
			add(xmax-xmin+1)
			break
			
		if xmax < choix1:
			# Objet entierement a gauche
			add(xmax-xmin+1)
			break
		elif xmin < choix1:
			# Objet a cheval a gauche
			add(choix1-xmin)
			xmin = choix1
		else:
			# Objet plus profond a gauche
			pass
		
		if xmin > choix2:
			# Objet entierement a droite
			add(xmax-xmin+1)
			break
		elif xmax > choix2:
			# Objet a cheval a droite
			add(xmax-choix2)
			xmax = choix2
		else:
			# Objet plus profond a droite
			pass
		nb += 1
	assert to == l, (to,l)
	
print >> sys.stderr, "Computing genomic sizes ...",
# Pour la normalisation
tot1 = collections.defaultdict(int)
tot2 = collections.defaultdict(int)
tot3 = collections.defaultdict(int)
for c in genome:
	# L'interieur des diagonales
	for x in genome[c]:
		if x[2] == "intron":
			addInterv(tot1, c, x[0], x[1])
		elif x[2] == "exon":
			addInterv(tot2, c, x[0], x[1])
		else:
			addInterv(tot3, c, x[0], x[1])
print >> sys.stderr, "intron=", sum(tot1.values()), "exon=", sum(tot2.values()), "intergenic=", sum(tot3.values()), "OK"

# Positionnement des CNEs
print >> sys.stderr, "Loading CNEs ...",
res = collections.defaultdict(int)
res1 = collections.defaultdict(int)
res2 = collections.defaultdict(int)
res3 = collections.defaultdict(int)
for line in utils.myFile.myTSV.readTabular(arguments["cneFile"], [str,str,int,int,str,str,int,str]):
	c = utils.myGenomes.commonChrName(line[1])
	i = bisect.bisect(genome[c], (line[2],line[3])) - 1
	if i >= 0:
		ref = genome[c][i]
		if ref[0] <= line[2] <= line[3] <= ref[1]:
			if ref[2] == "intron":
				addInterv(res1, c, line[2], line[3])
			elif ref[2] == "exon":
				addInterv(res2, c, line[2], line[3])
			else:
				addInterv(res3, c, line[2], line[3])
		else:
			print >> sys.stderr, "-",
	else:
		print >> sys.stderr, "_",
print >> sys.stderr, "OK"


print >> sys.stderr, "Printing results ...",
for x in sorted(set(res1.keys()+res2.keys()+res3.keys())):
	lst = [x]
	def pr(x, y):
		if y > 0:
			lst.append( round(100.*float(x)/float(y),2) )
		else:
			lst.append(0)
		lst.append(y)
	pr(res1[x], tot1[x])
	pr(res2[x], tot2[x])
	pr(res3[x], tot3[x])
	pr(res1[x]+res3[x], tot1[x]+tot3[x])
	pr(res1[x]+res2[x]+res3[x], tot1[x]+tot2[x]+tot3[x])
	print utils.myFile.myTSV.printLine(lst)

print >> sys.stderr, "OK"

