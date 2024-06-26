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
	[("repeatFile",file), ("intervalFile",file), ("refSpecies_latin",str), ("diagsFile",file)], \
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

def addInterv(countD, countI, cat, c, xmin, xmax):

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
	xmax += 1
	l = xmax-xmin
	step = (choix2-choix1)/float(arguments["nbBins"])
	nmin = int((xmin-choix1)/step)
	nmax = int((xmax-choix1)/step)

	global to
	to = 0
	nb = 0
	def add(n, size):
		count[n][0].append( (l, 1./(nmax-nmin+1)) )
		count[invStep-n][0].append( (l, 1./(nmax-nmin+1)) )
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

	assert to == l, (to,l)



def mkDic():
	#cats = ["LINE", "SINE", "DNA", "LTR", "exon", "intron", "intergenic", "gene"]
	#cats = ['Alu', 'MIR', 'MER1_type', 'MER2_type', "exon", "intron", "intergenic", "gene"]
	#cats = ['CR1', 'L1', 'L2', 'RTE', "exon", "intron", "intergenic", "gene"]
	cats = ['ERV1', 'ERVK', 'ERVL', 'MaLR', "exon", "intron", "intergenic", "gene"]
	d = {}
	for x in cats:
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
		addInterv(totD, totI, x[2], c, x[0], x[1])
print >> sys.stderr, "OK"

# Les quantites d'espace dans chaque categorie, la somme totale vaut 2*genome
print >> sys.stderr, " ".join("%s=%d" % (x, sum(y[1] for y in totD[x])) for x in totD), "OK"
print >> sys.stderr, " ".join("%s=%d" % (x, sum(y[1] for y in totI[x])) for x in totI), "OK"


# Positionnement des CNEs
print >> sys.stderr, "Loading CNEs ...",
repeatD = mkDic()
repeatI = mkDic()
for line in utils.myFile.myTSV.readTabular(arguments["repeatFile"], [str,int,int,str]):
	c = utils.myGenomes.commonChrName(line[0])
	i = bisect.bisect(genome[c], (line[1],line[2])) - 1
	if i >= 0:
		ref = genome[c][i]
		if ref[0] <= line[1] <= line[2] <= ref[1]:
			addInterv(repeatD, repeatI, line[3], c, line[1], line[2])
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

	#for i in ["LINE", "SINE", "DNA", "LTR"]:
	#for i in ['Alu', 'MIR', 'MER1_type', 'MER2_type']:
	#for i in ['CR1', 'L1', 'L2', 'RTE']:
	for i in ['ERV1', 'ERVK', 'ERVL', 'MaLR']:
		do(repeatD, totD, [i], ["exon", "intron", "intergenic"])
		do(repeatI, totI, [i], ["exon", "intron", "intergenic"])

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
	do(f, [[1+i], [6+i], [11+i], [16+i], [21+i], [26+i], [31+i], [36+i]])
	f.close()

shutil.copy("header2.xmgr", arguments["outputName"] % names[-1])
f = utils.myFile.openFile(arguments["outputName"] % names[-1], "a")
do(f, [[5,15,25,35], [10,20,30,40]])
f.close()

print >> sys.stderr, "OK"

