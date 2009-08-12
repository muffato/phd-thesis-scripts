#! /users/ldog/muffato/python

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
	[("cneFile",file), ("intervalFile",file)], \
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
for c in genome:
	#assert genome[c] == sorted(genome[c])
	#for (x,y) in utils.myTools.myIterator.slidingTuple(genome[c]):
	#	assert y[0] == x[1]+1
	annot[c] = [None] * len(genome[c])
	lstCNEs[c] = [None] * len(genome[c])
	for (i,x) in enumerate(genome[c]):
		lstCNEs[c][i] = []
		for g in x[3]:
			dicGene2Pos[g].append( (c,i) )
print >> sys.stderr, "OK"


# Chargement des diagonales
print >> sys.stderr, "Loading pseudo-synteny blocks ...",
for x in dicGene2Pos:

	if isinstance(x, frozenset):
		continue

	# Les positions de ces genes
	genes = [x]
	lg = len(genes)
	ind = dict([(g,i) for (i,g) in enumerate(genes)])
	(cmin,imin) = min([min(dicGene2Pos[g]) for g in genes])
	(cmax,imax) = max([max(dicGene2Pos[g]) for g in genes])
	assert cmin == cmax
	c = cmin
	
	last = None
	for i in xrange(imin, imax+1):
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
			i += step
	
	travel(c, imax, 1)
	travel(c, imin, -1)

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
count = collections.defaultdict(list)
bygene = collections.defaultdict(list)
for c in annot:
	for (coords,t,l) in itertools.izip(genome[c], annot[c], lstCNEs[c]):
			print utils.myFile.myTSV.printLine(["[%s]" % (t,), coords, l])
			d = coords[1]-coords[0]+1
			if coords[2] in ["exon", "intron"]:
				for g in coords[3]:
					bygene[g].append( (t,d,l,coords[2]) )
				if coords[2] == "intron":
					count[(t,None,"intron")].append( (d,l) )
			else:
				count[t].append( (d,l) )

# Groupement des introns/exons de chaque gene
for (g,gg) in bygene.iteritems():
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
				print utils.myFile.myTSV.printLine(["[%s]" % ((t,None,name),), (d,None,name,frozenset([g])), l])
				count[(t,None,name)].append( (d,l) )
		do("gene-introns", di, li)
		do("gene-%d-introns" % ne, di, li)
		do("gene-exons", de, le)
		do("gene-%d-exons" % ne, de, le)
		do("gene", di+de, li+le)
		do("gene-%d" % ne, di+de, li+le)

for x in sorted(count):
		lreg = []
		lcnes = []
		for (d,l) in count[x]:
			lreg.append(d)
			lcnes.extend(l)

		lreg.sort()
		lcnes.sort()
		l1 = len(lreg)
		l2 = len(lcnes)
		s1 = sum(lreg)
		s2 = sum(lcnes)

		print utils.myFile.myTSV.printLine([x, l1,s1,(s1/l1)/1000.,lreg[l1/2], l2,s2,round(s2/float(l2),2) if l2>0 else 0,lcnes[l2/2] if l2 > 0 else 0, round(l2/float(l1),3), round(1000000.*l2/s1,2), round(100.*s2/s1,3)])
print >> sys.stderr, "OK"

