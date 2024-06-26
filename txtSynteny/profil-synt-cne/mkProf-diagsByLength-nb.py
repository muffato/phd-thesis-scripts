#!/usr/bin/env python2
# Copyright (c) 2010 IBENS/Dyogen : Matthieu MUFFATO, Alexandra LOUIS, Hugues ROEST CROLLIUS
# Mail : hrc@ens.fr or alouis@biologie.ens.fr
# Licences GLP v3 and CeCILL v2

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
	[("cneFile",file), ("intervalFile",file), ("refSpecies_latin",str), ("diagsFile",file)], \
	[("flankLength",int,10)], \
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
annot = {}
lstCNEs = {}
diagLength = {}
for c in genome:
	#assert genome[c] == sorted(genome[c])
	#for (x,y) in utils.myTools.myIterator.slidingTuple(genome[c]):
	#	assert y[0] == x[1]+1
	annot[c] = [None] * len(genome[c])
	lstCNEs[c] = [None] * len(genome[c])
	diagLength[c] = [None] * len(genome[c])
	for (i,x) in enumerate(genome[c]):
		lstCNEs[c][i] = []
		diagLength[c][i] = []
		for g in x[3]:
			dicGene2Pos[g].append( (c,i) )
print >> sys.stderr, "OK"


# Chargement des diagonales
print >> sys.stderr, "Loading synteny blocks ...",
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
	assert cmin == cmax
	assert c == cmax
	
	last = None
	for i in xrange(imin, imax+1):
		diagLength[c][i] = [lg]
		if genome[c][i][2] == "intergenic":
			if middle:
				annot[c][i] = (last-1, last, "intergenic")
			else:
				annot[c][i] = (last,last+1, "intergenic")
		else:
			for g in genome[c][i][3]:
				if g in ind:
					j = min(ind[g], lg-1-ind[g])
					middle = (j == lg-1-ind[g])
					annot[c][i] = last = j + 1
					break
	
	def travel(c, start, step):
		i = start + step
		nb = 1
		while (nb >= -arguments["flankLength"]) and (i >= 0) and (i < len(genome[c])):
			if genome[c][i][2] == "intergenic":
				nb -= 1
				res = (nb, nb+1, "intergenic")
			else:
				res = nb
			if (annot[c][i] is None) or (res >= annot[c][i]):
				annot[c][i] = res
				diagLength[c][i].append(lg)
			i += step
	
	travel(c, imax, 1)
	travel(c, imin, -1)

ft.file.close()
print >> sys.stderr, "OK"


# Positionnement des CNEs
print >> sys.stderr, "Loading CNEs ...",
for line in utils.myFile.myTSV.readTabular(arguments["cneFile"], [str,str,int,int,str,str,int,str]):
	c = utils.myGenomes.commonChrName(line[1])
	i = bisect.bisect(genome[c], (line[2],line[3])) - 1
	if i >= 0:
		ref = genome[c][i]
		if ref[0] <= line[2] <= line[3] <= ref[1]:
			if ref[2] == "exon":
				print >> sys.stderr, "?",
			else:
				lstCNEs[c][i].append( line[3]-line[2]+1)
		else:
			print >> sys.stderr, "-",
	else:
		print >> sys.stderr, "_",
print >> sys.stderr, "OK"


print >> sys.stderr, "Printing results ...",
count = {}
bygene = {}
for c in annot:
	for (coords,t,l,ld) in itertools.izip(genome[c], annot[c], lstCNEs[c], diagLength[c]):
		for x in ld:
			print utils.myFile.myTSV.printLine(["{%d}" % (x,), "[%s]" % (t,), coords, l])
			if x not in count:
				count[x] = collections.defaultdict(list)
				bygene[x] = collections.defaultdict(list)
			d = coords[1]-coords[0]+1
			if coords[2] in ["exon", "intron"]:
				for g in coords[3]:
					bygene[x][g].append( (t,d,l,coords[2]) )
				if coords[2] == "intron":
					count[x][(t,None,"intron")].append( (d,l) )
			else:
				count[x][t].append( (d,l) )

# Groupement des introns/exons de chaque gene
for y in bygene:
	for (g,gg) in bygene[y].iteritems():
		nb = []
		le = []
		li = []
		de = 0
		di = 0
		ne = 0
		for x in gg:
			nb.append(x[0])
			if x[-1] == "exon":
				le.extend(x[2])
				de += x[1]
				ne += 1
			else:
				li.extend(x[2])
				di += x[1]
		t = max([(nb.count(x),x) for x in set(nb)])[1]
		def do(name, d, l):
			if d > 0:
				print utils.myFile.myTSV.printLine(["{%d}" % (y,), "[%s]" % ((t,None,name),), (d,None,name,frozenset([g])), l])
				count[y][(t,None,name)].append( (d,l) )
		do("gene-introns", di, li)
		do("gene-%d-introns" % ne, di, li)
		do("gene-exons", de, le)
		do("gene-%d-exons" % ne, de, le)
		do("gene", di+de, li+le)
		do("gene-%d" % ne, di+de, li+le)


for y in sorted(count):
	for x in sorted(count[y]):
		lreg = []
		lcnes = []
		for (d,l) in count[y][x]:
			lreg.append(d)
			lcnes.extend(l)

		lreg.sort()
		lcnes.sort()
		l1 = len(lreg)
		l2 = len(lcnes)
		s1 = sum(lreg)
		s2 = sum(lcnes)

		print utils.myFile.myTSV.printLine([y, x, l1,s1,(s1/l1)/1000.,lreg[l1/2], l2,s2,round(s2/float(l2),2) if l2>0 else 0,lcnes[l2/2] if l2 > 0 else 0, round(l2/float(l1),3), round(1000000.*l2/s1,2), round(100.*s2/s1,3)])
print >> sys.stderr, "OK"

