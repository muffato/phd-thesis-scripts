#!/usr/bin/env python2
# Copyright (c) 2010 IBENS/Dyogen : Matthieu MUFFATO, Alexandra LOUIS, Hugues ROEST CROLLIUS
# Mail : hrc@ens.fr or alouis@biologie.ens.fr
# Licences GLP v3 and CeCILL v2

__doc__ = """
	Lit les CNE et transforme les positions genomiques en positions indexees
"""

import sys
import math
import bisect
import operator
import itertools
import collections

import utils.myTools
import utils.myGenomes

arguments = utils.myTools.checkArgs( \
	[("cneFile",file), ("intervalFile",file), ("refSpecies_latin",str), ("diagsFile",file)], \
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
	#assert genome[c] == sorted(genome[c])
	#for (x,y) in utils.myTools.myIterator.slidingTuple(genome[c]):
	#	assert y[0] == x[1]+1
	for (i,x) in enumerate(genome[c]):
		for g in x[3]:
			dicGene2Pos[g].append( (c,i) )
print >> sys.stderr, "OK"


# Chargement des diagonales
print >> sys.stderr, "Loading pseudo-synteny blocks ...",
startGenes = collections.defaultdict(list)
endGenes = {}
for x in dicGene2Pos:

	if isinstance(x, frozenset):
		continue

	(cmin,imin) = min(dicGene2Pos[x])
	(cmax,imax) = max(dicGene2Pos[x])
	assert cmin == cmax
	c = cmin
	xmin = genome[c][imin][0]
	xmax = genome[c][imax][1]
	startGenes[c].append( xmin )
	endGenes[(c,xmin)] = xmax

for x in startGenes.itervalues():
	x.sort()
print >> sys.stderr, "OK"


# Chargement des diagonales
print >> sys.stderr, "Loading synteny blocks ...",
startDiags = collections.defaultdict(list)
endDiags = {}
ft = utils.myFile.myTSV.reader(arguments["diagsFile"])
for t in ft.csvobject:

	# Les genes de la diagonale
	if arguments["refSpecies_latin"] in t[2]:
		c = utils.myGenomes.commonChrName(t[3])
		do = t[4]
	else:
		assert arguments["refSpecies_latin"] in t[5]
		c = utils.myGenomes.commonChrName(t[6])
		do = t[7]

	# Les positions de ces genes
	genes = do.split()
	lg = len(genes)
	ind = dict([(g,i) for (i,g) in enumerate(genes)])
	(cmin,imin) = min([min(dicGene2Pos[g]) for g in genes])
	(cmax,imax) = max([max(dicGene2Pos[g]) for g in genes])
	assert c == cmin
	assert c == cmax

	xmin = genome[c][imin][0]
	xmax = genome[c][imax][1]
	startDiags[c].append( xmin )
	endDiags[(c,xmin)] = xmax
ft.file.close()

for x in startDiags.itervalues():
	x.sort()
print >> sys.stderr, "OK"


step = arguments["windowSize"] * 1000
def addInterv(diags, fin, count, c, xmin, xmax):

	i = bisect.bisect_right(diags[c], xmin )
	
	if (i >= 1) and (xmin < fin[(c,diags[c][i-1])]):
		# dans diagonale precedente
		choix1 = diags[c][i-1]
		choix2 = fin[(c,diags[c][i-1])]
		into = True
	else:
		# entre diagonales
		choix1 = fin[(c,diags[c][i-1])] if i >= 1 else -sys.maxint
		choix2 = diags[c][i] if i < len(diags[c]) else sys.maxint
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
tot1b = collections.defaultdict(int)
tot2b = collections.defaultdict(int)
tot3b = collections.defaultdict(int)
for c in genome:
	# L'interieur des diagonales
	for x in genome[c]:
		if x[2] == "intron":
			addInterv(startDiags, endDiags, tot1, c, x[0], x[1])
			addInterv(startGenes, endGenes, tot1b, c, x[0], x[1])
		elif x[2] == "exon":
			addInterv(startDiags, endDiags, tot2, c, x[0], x[1])
			addInterv(startGenes, endGenes, tot2b, c, x[0], x[1])
		else:
			addInterv(startDiags, endDiags, tot3, c, x[0], x[1])
			addInterv(startGenes, endGenes, tot3b, c, x[0], x[1])
print >> sys.stderr, "intron=", sum(tot1.values()), "exon=", sum(tot2.values()), "intergenic=", sum(tot3.values()), "OK"

# Positionnement des CNEs
print >> sys.stderr, "Loading CNEs ...",
res = collections.defaultdict(int)
res1 = collections.defaultdict(int)
res2 = collections.defaultdict(int)
res3 = collections.defaultdict(int)
res1b = collections.defaultdict(int)
res2b = collections.defaultdict(int)
res3b = collections.defaultdict(int)
for line in utils.myFile.myTSV.readTabular(arguments["cneFile"], [str,str,int,int,str,str,int,str]):
	c = utils.myGenomes.commonChrName(line[1])
	i = bisect.bisect(genome[c], (line[2],line[3])) - 1
	if i >= 0:
		ref = genome[c][i]
		if ref[0] <= line[2] <= line[3] <= ref[1]:
			if ref[2] == "intron":
				addInterv(startDiags, endDiags, res1, c, line[2], line[3])
				addInterv(startGenes, endGenes, res1b, c, line[2], line[3])
			elif ref[2] == "exon":
				addInterv(startDiags, endDiags, res2, c, line[2], line[3])
				addInterv(startGenes, endGenes, res2b, c, line[2], line[3])
			else:
				addInterv(startDiags, endDiags, res3, c, line[2], line[3])
				addInterv(startGenes, endGenes, res3b, c, line[2], line[3])
		else:
			print >> sys.stderr, "-",
	else:
		print >> sys.stderr, "_",
print >> sys.stderr, "OK"


print >> sys.stderr, "Printing results ...",
for x in sorted(set(res1.keys()+res2.keys()+res3.keys())):
	lst = [x]
	def pr(x, y, xx, yy):
		if y > 0:
			lst.append( round(100.*float(x)/float(y),2) )
		else:
			lst.append(0)
		lst.append(y)
		if (y > 0) and (yy > 0) and (x > 0) and (xx > 0):
			r = 100.*float(x)/float(y)
			rr = 100.*float(xx)/float(yy)
			lst.append( round(rr, 2) )
			lst.append( round(math.log(r/rr, 2), 2) )
		else:
			lst.append(0)
			lst.append(0)
	pr(res1[x], tot1[x], res1b[x], tot1b[x])
	pr(res2[x], tot2[x], res2b[x], tot2b[x])
	pr(res3[x], tot3[x], res3b[x], tot3b[x])
	pr(res1[x]+res3[x], tot1[x]+tot3[x], res1b[x]+res3b[x], tot1b[x]+tot3b[x])
	pr(res1[x]+res2[x]+res3[x], tot1[x]+tot2[x]+tot3[x], res1b[x]+res2b[x]+res3b[x], tot1b[x]+tot2b[x]+tot3b[x])
	print utils.myFile.myTSV.printLine(lst)

print >> sys.stderr, "OK"

