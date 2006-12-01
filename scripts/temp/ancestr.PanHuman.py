#! /users/ldog/muffato/python

##################
# INITIALISATION #
##################

# Librairies
import os
import sys
import utils.myGenomes
import utils.myTools
import utils.myMaths
import utils.myDiags
import utils.myPsOutput
import cPickle



def isOK(filsOK, filsNO, outgroupOK, outgroupNO):
	return max(filsNO) == 0

lst = utils.myDiags.DiagRepository()
for l in sys.stdin:
	d = l.split('\t')
	c = [int(x) for x in d[1].split()]
	if len(c) < 3:
		continue
	a = [tuple(x.split('/')) for x in d[2].split()]
	lst.addDiag(c, a)
	
print >> sys.stderr, "lecture OK", lst.nbRealDiag

one2one2one = []
sameChrFils = {}
nearly = []
poubelle = []
for i in xrange(len(lst.lstDiags)):
	a = lst.lstApp[i]
	cH = [c for (e,c) in a if e == 'Human']
	cC = [c for (e,c) in a if e == 'Chimp']
	cM = [c for (e,c) in a if e == 'Macaque']
	if len(cH) == 1 and len(cC) == 1:
		cH = cH[0]
		cC = cC[0]
		if (cH,cC) not in sameChrFils:
			sameChrFils[(cH,cC)] = []
		sameChrFils[(cH,cC)].append( (i,cM) )
	else:
		poubelle.append(i)

outgroup = {}
for x in sameChrFils:
	s = set(utils.myMaths.flatten([y[1] for y in sameChrFils[x]]))
	outgroup[x] = s
for x in outgroup:
	s = [y[0] for y in sameChrFils[x]]
	n = sum([len(lst.lstDiags[y]) for y in s])
	if n < 20:
		print >> sys.stderr, "suppression de", x, n, s
		poubelle.extend(s)
		del sameChrFils[x]
print >> sys.stderr, sameChrFils.keys()
print >> sys.stderr, len(poubelle)

poubelle2 = []
for i in poubelle:
	a = lst.lstApp[i]
	cH = [c for (e,c) in a if e == 'Human']
	cC = [c for (e,c) in a if e == 'Chimp']
	cM = [c for (e,c) in a if e == 'Macaque']
	poss = []
	if len(cH) == 1:
		cH = cH[0]
		poss.extend([x for x in sameChrFils if x[0] == cH])
	if len(cC) == 1:
		cC = cC[0]
		poss.extend([x for x in sameChrFils if x[1] == cC])
	if len(set(poss)) == 1:
		x = poss[0]
		sameChrFils[x].append( (i,cM) )
	else:
		poubelle2.append(i)

print >> sys.stderr, len(poubelle2)
print >> sys.stderr, [(x,len(sameChrFils[x])) for x in sameChrFils]

chromosomes = {}
for x in sameChrFils:
	chromosomes[x] = set(utils.myMaths.flatten([lst.lstDiags[l]  for (l,_) in sameChrFils[x]]))

# TODO: Filtre sur les genes en double
for c in chromosomes:
	for d in chromosomes:
		if c == d:
			continue
		i = chromosomes[c].intersection(chromosomes[d])
		chromosomes[c].difference_update(i)
		chromosomes[d].difference_update(i)


# TODO: Filtre sur les genes sans assignation
tout = set(utils.myMaths.flatten(chromosomes.values()))
genesAnc = utils.myGenomes.loadGenome(sys.argv[1])
lstGenesAnc = genesAnc.lstGenes[utils.myGenomes.AncestralGenome.defaultChr]
gH = utils.myGenomes.loadGenome("~/work/data/genes/genesHuman.list.bz2")
gC = utils.myGenomes.loadGenome("~/work/data/genes/genesChimp.list.bz2")

poubelle3 = []
for i in xrange(len(lstGenesAnc)):
	g = lstGenesAnc[i]
	if i in tout:
		continue
	cH = set([str(gH.dicGenes[x][0]) for x in g.names if x in gH.dicGenes])
	cC = set([str(gC.dicGenes[x][0]) for x in g.names if x in gC.dicGenes])

	poss = []
	if len(cH) == 1:
		cH = cH.pop()
		poss.extend([x for x in sameChrFils if x[0] == cH])
	if len(cC) == 1:
		cC = cC.pop()
		if cC == "2a" or cC == "2b":
			cC = "2"
		poss.extend([x for x in sameChrFils if x[1] == cC])
	if len(set(poss)) == 1:
		x = poss[0]
		chromosomes[x].add(i)
	else:
		print >> sys.stderr, cH, cC, poss
		poubelle3.append(i)

print >> sys.stderr, len(poubelle3)

#genesAnc = utils.myGenomes.loadGenome(sys.argv[1])
#lstGenesAnc = genesAnc.lstGenes[utils.myGenomes.AncestralGenome.defaultChr]
#chromosomes = sameChrFils
for c in chromosomes:
	for i in chromosomes[c]:
		print "%s.%s" % c, " ".join(lstGenesAnc[i].names)
	#print "PanHuman\t%s\t%s" % (" ".join([str(x) for x in chromosomes[c]]), c)


sys.exit(0)
s = ['Human','Chimp','Macaque']
#for (d,_,l) in lst:
for i in poubelle2:
	d = lst.lstDiags[i]
	l = lst.lstApp[i]
	print "\t%s\t%s" % (" ".join([str(x) for x in d]), " ".join(["%s/%s" % (e,c) for (e,c) in l if e in s]))
sys.exit(0)

lst.combinDiags([['Human'],['Chimp']], isOK)
print >> sys.stderr, "combin OK", lst.nbRealDiag

