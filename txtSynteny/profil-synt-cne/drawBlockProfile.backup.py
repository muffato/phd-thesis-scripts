#!/usr/bin/env python2
# Copyright (c) 2010 IBENS/Dyogen : Matthieu MUFFATO, Alexandra LOUIS, Hugues ROEST CROLLIUS
# Mail : hrc@ens.fr or alouis@biologie.ens.fr
# Licences GLP v3 and CeCILL v2

__doc__ = """
	Lit les CNE et transforme les positions genomiques en positions indexees
"""

import os
import sys
import bisect
import operator
import itertools
import collections

import utils.myMaths
import utils.myTools
import utils.myGenomes

arguments = utils.myTools.checkArgs( \
	[("cneFile",file), ("intervalFile",file), ("refSpecies_latin",str), ("diagsFile",file), ("strandFile",file)], \
	[("nbBins",int,1000), ("outputName",str,"")], \
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
print >> sys.stderr, "Loading synteny blocks ...",
diags = collections.defaultdict(list)
fin = {}
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
	#if lg >= 8:
	#	continue
	ind = dict([(g,i) for (i,g) in enumerate(genes)])
	(cmin,imin) = min([min(dicGene2Pos[g]) for g in genes])
	(cmax,imax) = max([max(dicGene2Pos[g]) for g in genes])
	assert c == cmin
	assert c == cmax

	xmin = genome[c][imin][0]
	xmax = genome[c][imax][1]
	diags[c].append( xmin )
	fin[(c,xmin)] = xmax
ft.file.close()

for x in diags.itervalues():
	x.sort()
print >> sys.stderr, "OK"


invStep = arguments["nbBins"]-1

def addInterv(countD, countI, cat, c, xmin, xmax, val1, val2):

	i = bisect.bisect_right(diags[c], xmin )
	
	if (i >= 1) and (xmin < fin[(c,diags[c][i-1])]):
		# dans diagonale precedente
		choix1 = diags[c][i-1]
		choix2 = fin[(c,diags[c][i-1])]
		count = countD[cat]
	else:
		# entre diagonales
		choix1 = fin[(c,diags[c][i-1])] if i >= 1 else None
		choix2 = diags[c][i] if i < len(diags[c]) else None
		if (choix1 is None) or (choix2 is None):
			return
		count = countI[cat]
	
	if not (choix1 <= xmin <= xmax <= choix2):
		print >> sys.stderr, "o",
		return

	# Les nouvelles positions
	choix2 += 1
	l = xmax-xmin+1
	step = (choix2-choix1)/float(arguments["nbBins"])
	nmin = int((xmin-choix1)/step)
	nmax = int((xmax-choix1)/step)
	xmax += 1

	global to
	to = 0
	def add(n, size):
		count[n][0].append( (val1, 1./(nmax-nmin+1)) )
		count[invStep-n][0].append( (val2, 1./(nmax-nmin+1)) )
		count[n][1] += size
		count[invStep-n][1] += size
		global to
		to = to + size

	# La quantite d'objet (pour les densites
	if nmin != nmax:
		for i in xrange(nmin+1, nmax):
			add(i, step)
		add(nmin, step*(nmin+1) - (xmin-choix1))
		add(nmax, (xmax-choix1) - step*nmax)
	else:
		add(nmin, xmax-xmin)

	#assert to == l, (to,l)
	assert abs(to-l)<1e-7, (to,l)



def mkDic():
	d = {}
	for x in ["exon", "intron", "intergenic", "gene"]:
		d[x] = range(arguments["nbBins"])
		for i in xrange(arguments["nbBins"]):
			d[x][i] = [[], 0]
	return d

print >> sys.stderr, "Inserting intervals ...",
totD = mkDic()
totI = mkDic()
for c in genome:
	# L'interieur des diagonales
	for x in genome[c]:
		if x[2] == "intergenic":
			val = []
			for (g1,g2) in itertools.product(x[3][0], x[3][1]):
				val.append((strand[g1]*strand[g2]+1)/2)
			val = sum(val)/float(len(val))
			val1 = val2 = val
		else:
			val = []
			for g in x[3]:
				val.append((strand[g]+1)/2)
			val = sum(val)/float(len(val))
			val1 = val
			val2 = 1-val
		addInterv(totD, totI, x[2], c, x[0], x[1], val1, val2)

# Les genes
for g in dicGene2Pos:
	if isinstance(g, frozenset):
		continue
	(cmin,imin) = min(dicGene2Pos[g])
	(cmax,imax) = max(dicGene2Pos[g])
	assert cmax == cmin
	c = cmin

	xmin = genome[c][imin][0]
	xmax = genome[c][imax][1]

	addInterv(totD, totI, "gene", c, xmin, xmax, 0, 0)

print >> sys.stderr, "OK"

# Les quantites d'espace dans chaque categorie, la somme totale vaut 2*genome
print >> sys.stderr, " ".join("%s=%d" % (x, sum(y[1] for y in totD[x])) for x in totD), "OK"
print >> sys.stderr, " ".join("%s=%d" % (x, sum(y[1] for y in totI[x])) for x in totI), "OK"


# Positionnement des CNEs
print >> sys.stderr, "Loading CNEs ...",
cneD = mkDic()
cneI = mkDic()
for line in utils.myFile.myTSV.readTabular(arguments["cneFile"], [str,str,int,int,str,str,int,str]):
	c = utils.myGenomes.commonChrName(line[1])
	i = bisect.bisect(genome[c], (line[2],line[3])) - 1
	if i >= 0:
		ref = genome[c][i]
		if ref[0] <= line[2] <= line[3] <= ref[1]:
			addInterv(cneD, cneI, ref[2], c, line[2], line[3], 0, 0)
		else:
			print >> sys.stderr, "-",
	else:
		print >> sys.stderr, "_",
print >> sys.stderr, "OK"


print >> sys.stderr, "Printing results ...",
txt = []
for x in xrange(arguments["nbBins"]):
	lst = [(100.*x)/arguments["nbBins"]]
	
	# Affiche le rapport entre les deux valeurs
	def pr(a, b):
		if b > 0:
			lst.append( round(100.*float(a)/float(b),4) )
		else:
			lst.append(0)
	
	# Affiche la moyenne et la mediane de la liste de valeurs
	def prL(l):
		if len(l) > 0:
			s = 0.
			n = 0.
			s2 = 0.
			n2 = 0.
			for (a,b) in l:
				s += a*b
				n += b
				s2 += a
				n2 += 1.
			l.sort()
			n3 = 0
			for (a,b) in l:
				n3 += b
				if n3 > n/2:
					break
			lst.append( round(s/n,4) )
			lst.append( a )
			lst.append( round(s2/n2,4) )
			lst.append( l[len(l)/2][0] )
		else:
			lst.append(0)
			lst.append(0)
			lst.append(0)
			lst.append(0)

	# Affiche longueur moyenne/mediane de l'objet puis la densite par rapport a la reference
	def do(count, ref, cats, refcats):
		l = []
		for c in cats:
			l.extend(count[c][x][0])
		prL(l)
		if ref is not None:
			pr(sum(count[c][x][1] for c in cats), sum(ref[c][x][1] for c in refcats))

	# Stats des objets geniques
	for i in ["gene", "intron", "exon", "intergenic"]:
		do(totD, totD, [i], ["exon", "intron", "intergenic"])
		do(totI, totI, [i], ["exon", "intron", "intergenic"])
	
	# Stats des CNEs
	for i in [["intron"], ["intergenic"], ["intron", "intergenic"]]:
		do(cneD, totD, i, i)
		do(cneI, totI, i, i)
	
	txt.append(lst)
	print utils.myFile.myTSV.printLine(lst)

def do(f, graph):
	for (ng,cols) in enumerate(graph):
		for (ns,i) in enumerate(cols):
			print >> f, "@target G%d.S%d" % (ng,ns)
			print >> f, "@type xy"
			for l in txt:
				print >> f, l[0], l[i]
			print >> f, "&"

import shutil
names = ["meanI", "medianI", "mean1", "median1", "gen"]
for i in xrange(4):
	shutil.copy("header1.xmgr", arguments["outputName"] % names[i])
	f = utils.myFile.openFile(arguments["outputName"] % names[i], "a")
	do(f, [[1+i], [6+i], [11+i], [16+i], [21+i], [26+i], [31+i], [36+i], [41+i,51+i,61+i], [46+i,56+i,66+i]])
	f.close()

shutil.copy("header2.xmgr", arguments["outputName"] % names[-1])
f = utils.myFile.openFile(arguments["outputName"] % names[-1], "a")
do(f, [[5,15,25,35], [10,20,30,40], [45,55,65], [50,60,70]])
f.close()

print >> sys.stderr, "OK"

