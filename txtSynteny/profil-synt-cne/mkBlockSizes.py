#! /users/ldog/muffato/python

__doc__ = """
	Lit les CNE et transforme les positions genomiques en positions indexees
"""

import sys
import bisect
import operator
import itertools
import collections

import utils.myMaths
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


# Initialisation des comptages
nbCNE = {}
longCNE = {}
for (c,x) in genome.lstGenes.iteritems():
	d = {}
	lx = len(x)
	for i in xrange(lx):
		d[i] = 0
		d[(i-1,i)] = 0
	d[(lx-1,lx)] = 0

	for j in xrange(2+1):
		d[-1-j] = 0
		d[(-2-j,-1-j)] = 0
		d[lx+j] = 0
		d[(lx+j,lx+1+j)] = 0
	nbCNE[c] = d.copy()
	longCNE[c] = d.copy()


# Positionnement des CNEs
for line in utils.myFile.myTSV.readTabular(arguments["KX_HDMC_especes.list"], [str,str,int,int,str,str,int,str]):
	
	if line[0] == arguments["refSpecies_common"]:
		
		c = utils.myGenomes.commonChrName(line[1])

		i = bisect.bisect_left(genome.lstGenes[c], (line[2],line[3]))
		i = bisect.bisect_left(genome.lstGenes[c], utils.myGenomes.Gene(c, line[2], line[3], int(line[4]+"1"), line[5]))

		if (i >= 1) and (genome.lstGenes[c][i-1].end > line[2]):
			# intronique
			pos = i-1
		else:
			# inter genique
			pos = (i-1,i)
		
		nbCNE[c][pos] += 1
		longCNE[c][pos] += line[3]-line[2]+1


# Chargement des diagonales
sizes = collections.defaultdict(list)
sizesCNE = collections.defaultdict(list)
sizesCNEl = collections.defaultdict(list)
ft = utils.myFile.myTSV.reader(arguments["diagsFile"])
for t in ft.csvobject:

	# Les genes de la diagonale
	if arguments["refSpecies_latin"] in t[2]:
		do = t[4]
	else:
		assert arguments["refSpecies_latin"] in t[5]
		do = t[7]

	# Les positions de ces genes
	pos = [genome.dicGenes[g] for g in do.split()]
	l = len(pos)
	assert len(set([c for (c,i) in pos])) == 1
	assert abs(pos[-1][1]-pos[0][1]) == (l-1)
	(c,imin) = min(pos)
	assert max(pos) == (c,imin+l-1)

	d = []
	dcne = []
	lcne = []
	for i in xrange(imin, imin+l):
		d.append(genome.lstGenes[c][i].beginning - (genome.lstGenes[c][i-1].end if i >= 0 else 0))
		dcne.append( nbCNE[c][(i-1,i)] )
		lcne.append( longCNE[c][(i-1,i)] )
		d.append(genome.lstGenes[c][i].end - genome.lstGenes[c][i].beginning)
		dcne.append( nbCNE[c][i] )
		lcne.append( longCNE[c][i] )
	if imin+l < len(genome.lstGenes[c]):
		d.append(genome.lstGenes[c][imin+l].beginning - genome.lstGenes[c][imin+l-1].end)
	else:
		d.append(0)
	dcne.append( nbCNE[c][(imin+l-1,imin+l)] )
	lcne.append( longCNE[c][(imin+l-1,imin+l)] )

	d = [x/1000. for x in d]
	lcne = [None if x == 0 else lcne[i]/x for (i,x) in enumerate(dcne)]
	sizes[l].append( d )
	sizesCNE[l].append( dcne )
	sizesCNEl[l].append( lcne )

ft.file.close()

nbDiags = {}
def addRev(dic):
	for (l,ll) in dic.iteritems():
		nbDiags[l] = len(ll)
		ll.extend( [list(reversed(x)) for x in ll] )
		dic[l] = list(itertools.izip(*ll))

addRev(sizes)
addRev(sizesCNE)
addRev(sizesCNEl)

for (l,ll) in sizes.iteritems():

	#d = [list(itertools.ifilter(lambda y: y>0, x)) for x in itertools.izip(*ll)]
	#dcne = list(itertools.izip(*sizesCNE[l]))
	#dcnep = [list(itertools.ifilter(lambda y: y>0, x)) for x in itertools.izip(*sizesCNE[l])]
	#lcne = [list(itertools.ifilter(lambda y: y>0, x)) for x in itertools.izip(*sizesCNEl[l])]
	
	d = [list(itertools.ifilter(lambda y: y>0, x)) for x in sizes[l]]
	dcne = sizesCNE[l]
	dcnep = [list(itertools.ifilter(lambda y: y>0, x)) for x in sizesCNE[l]]
	lcne = [list(itertools.ifilter(lambda y: y>0, x)) for x in sizesCNEl[l]]

	def myapply(f, d):
		lg = utils.myMaths.flatten([d[2*i+1] for i in xrange(l)])
		li = utils.myMaths.flatten([d[2*i+2] for i in xrange(l-1)])
		return [f(lg), f(li)] + [f(x) for x in d]
	
	def pr(d):
		print utils.myFile.myTSV.printLine([l, nbDiags[l]] + d[:len(d)/2+2], func=mystr)
	
	def mystr(x):
		if isinstance(x, float):
			return "%.2f" % x
		elif isinstance(x, int):
			return "%d" % x
		else:
			return str(x)
	
	# Longueur des regions
	md = myapply(lambda x: utils.myMaths.myStats.mean(x), d)
	mdd = myapply(lambda x: utils.myMaths.myStats.getValue(sorted(x), 50), d)
	pr(md)
	pr(mdd)
	
	# Nombre de CNEs
	mdcne = myapply(lambda x: utils.myMaths.myStats.mean(x), dcne)
	mdcnep = myapply(lambda x: utils.myMaths.myStats.mean(x), dcnep)
	mddcnep = myapply(lambda x: utils.myMaths.myStats.getValue(sorted(x), 50), dcnep)
	pr(mdcne)
	pr(mdcnep)
	pr(mddcnep)
	
	# Longueur des CNEs
	mlcne = myapply(lambda x: utils.myMaths.myStats.mean(x), lcne)
	pr(mlcne)
	# inutile car pareil que la moyenne
	#mdlcne = myapply(lambda x: utils.myMaths.myStats.getValue(sorted(x), 50), lcne)
	#print utils.myFile.myTSV.printLine([l, nbDiags[l]] + mdlcne)
	#pr(mdlcne)
	
	# Densite en CNEs
	#pcne = myapply(lambda x: utils.myMaths.myStats.mean(x), [xn/xl for (xl,xn)  in zip(sizes[l], sizesCNE[l]) if xl > 0])
	#pr(pcne)
	pcne = [(None if (xn is None) or (xl is None) else 1000.*xn/xl) for (xl,xn) in zip(md,mdcne)]
	pr(pcne)

	#pcne = [[1000.*xn/xl for (xl,xn) in zip(lxl,lxn) if (xn is not None) and (xl is not None) and (xl != 0) and (xn > 0)] for (lxl,lxn) in zip(sizes[l], sizesCNE[l])]
	#mpcne =  myapply(lambda x: utils.myMaths.myStats.mean(x), pcne)
	#pr(mpcne)
	#mdpcne =  myapply(lambda x: utils.myMaths.myStats.getValue(sorted(x), 50), pcne)
	#pr(mdpcne)

	pcne = [(None if (xn is None) or (xl is None) or (xc is None) else xn*xc/(10.*xl)) for (xl,xn,xc) in zip(md,mdcne,mlcne)]
	pr(pcne)
	
	print	


