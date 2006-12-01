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
	cM = [c for (e,c) in a if e == 'Mouse']
	cR = [c for (e,c) in a if e == 'Rat']
	cH = [c for (e,c) in a if e == 'Human']
	if len(cM) == 1 and len(cR) == 1:
		cM = cM[0]
		cR = cR[0]
		if (cM,cR) not in sameChrFils:
			sameChrFils[(cM,cR)] = []
		sameChrFils[(cM,cR)].append( (i,cH) )
	else:
		poubelle.append(i)


print >> sys.stderr, [(x,len(sameChrFils[x])) for x in sameChrFils]
outgroup = {}
for x in sameChrFils:
	s = utils.myMaths.flatten([y[1]*len(lst.lstDiags[y[0]]) for y in sameChrFils[x]])
	ss = set(s)
	outgroup[x] = [(u,s.count(u)) for u in ss]
	print "outgroup", x, outgroup[x]
for c1 in outgroup:
	for c2 in outgroup:
		if c1 == c2:
			continue
		if ((c1[0] == c2[0]) or (c1[1] == c2[1])) and len(outgroup[c1].intersection(outgroup[c2])) > 0:
			print "fusion", c1, c2

for x in outgroup:
	s = [y[0] for y in sameChrFils[x]]
	n = sum([len(lst.lstDiags[y]) for y in s])
	if n < 20:
		print >> sys.stderr, "suppression de", x, n, s
		poubelle.extend(s)
		del sameChrFils[x]
print >> sys.stderr, sameChrFils.keys()
print >> sys.stderr, len(poubelle)
sys.exit(0)
poubelle2 = []
for i in poubelle:
	a = lst.lstApp[i]
	cM = [c for (e,c) in a if e == 'Mouse']
	cR = [c for (e,c) in a if e == 'Rat']
	cH = [c for (e,c) in a if e == 'Human']
	poss = []
	if len(cM) == 1:
		cM = cM[0]
		poss.extend([x for x in sameChrFils if x[0] == cM])
	if len(cR) == 1:
		cR = cR[0]
		poss.extend([x for x in sameChrFils if x[1] == cR])
	if len(set(poss)) == 1:
		x = poss[0]
		sameChrFils[x].append( (i,cH) )
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
gM = utils.myGenomes.loadGenome("~/work/data/genes/genesMouse.list.bz2")
gR = utils.myGenomes.loadGenome("~/work/data/genes/genesRat.list.bz2")

poubelle3 = []
for i in xrange(len(lstGenesAnc)):
	g = lstGenesAnc[i]
	if i in tout:
		continue
	cM = set([str(gM.dicGenes[x][0]) for x in g.names if x in gM.dicGenes])
	cR = set([str(gR.dicGenes[x][0]) for x in g.names if x in gR.dicGenes])

	poss = []
	if len(cM) == 1:
		cM = cM.pop()
		poss.extend([x for x in sameChrFils if x[0] == cM])
	if len(cR) == 1:
		cR = cR.pop()
		poss.extend([x for x in sameChrFils if x[1] == cR])
	if len(set(poss)) == 1:
		x = poss[0]
		chromosomes[x].add(i)
	else:
		print >> sys.stderr, cM, cR, poss
		poubelle3.append(i)

print >> sys.stderr, len(poubelle3)

#genesAnc = utils.myGenomes.loadGenome(sys.argv[1])
#lstGenesAnc = genesAnc.lstGenes[utils.myGenomes.AncestralGenome.defaultChr]
#chromosomes = sameChrFils
for c in chromosomes:
	for i in chromosomes[c]:
		print "%s.%s" % c, " ".join(lstGenesAnc[i].names)
	#print "PanMouse\t%s\t%s" % (" ".join([str(x) for x in chromosomes[c]]), c)


sys.exit(0)
s = ['Mouse','Rat','Human']
#for (d,_,l) in lst:
for i in poubelle2:
	d = lst.lstDiags[i]
	l = lst.lstApp[i]
	print "\t%s\t%s" % (" ".join([str(x) for x in d]), " ".join(["%s/%s" % (e,c) for (e,c) in l if e in s]))
sys.exit(0)

lst.combinDiags([['Mouse'],['Rat']], isOK)
print >> sys.stderr, "combin OK", lst.nbRealDiag


def isOK(filsOK, filsNO, outgroupOK, outgroupNO):
	#return max(filsNO) == 0
	return (max(filsNO) == 0) or (outgroupNO == 0 and outgroupOK >= 5)

lst = utils.myDiags.DiagRepository()
for l in sys.stdin:
	d = l.split('\t')
	c = [int(x) for x in d[1].split()]
	a = [tuple(x.split('/')) for x in d[2].split()]
	lst.addDiag(c, a)
	
print >> sys.stderr, "lecture OK", lst.nbRealDiag

#lst.combinDiags([['Mouse'],['Rat']], isOK)
print >> sys.stderr, "combin OK", lst.nbRealDiag

#s = ['Mouse', 'Rat', 'Human', 'Mouse', 'Dog']
#for (d,_,l) in lst:
#genesAnc = utils.myGenomes.loadGenome(sys.argv[1])
#for c in chromosomes:
	#for i in chromosomes[c]:
	#	print c, " ".join(genesAnc.lstGenes[utils.myGenomes.AncestralGenome.defaultChr][i].names)
#	print "Rodents\t%s\t%s" % (" ".join([str(x) for x in d]), " ".join(["%s/%s" % (e,c) for (e,c) in l if e in s]))

#sys.exit(0)

# 22 + X + Y + MT + 2b
mapping = {}
chromosomes = {}
for i in xrange(1, 22+1):
	s = str(i)
	mapping[('Mouse',s)] = s
	mapping[('Rat',s)] = s
	chromosomes[s] = set([])
mapping[('Mouse','X')] = 'X'
mapping[('Rat','X')] = 'X'
chromosomes['X'] = set([])
mapping[('Mouse','Y')] = 'Y'
mapping[('Rat','Y')] = 'Y'
chromosomes['Y'] = set([])
mapping[('Mouse','MT')] = 'MT'
mapping[('Rat','MT')] = 'MT'
chromosomes['MT'] = set([])
mapping[('Rat','2a')] = '2'
mapping[('Rat','2b')] = '2'


for (d,_,l) in lst:
	hum = [c for (e,c) in l if e =='Mouse']
	chi = [c for (e,c) in l if e =='Rat']
	poss = set([])
	if len(hum) == 1:
		c = ('Mouse',hum[0])
		if c in mapping:
			poss.add(mapping[c])
	if len(chi) == 1:
		c = ('Rat',chi[0])
		if c in mapping:
			poss.add(mapping[c])
	if len(poss) == 1:
		sol = poss.pop()
		#print "on met %s dans %s" % ( " ".join(["%s/%s" % (e,c) for (e,c) in l if e in s]), sol)
		chromosomes[sol].update(d)
	#else:
		#print "PanMouse\t%s\t%s" % (" ".join([str(x) for x in d]), " ".join(["%s/%s" % (e,c) for (e,c) in l if e in s]))
	#print len(d), " ".join(["%s/%s" % (e,c) for (e,c) in l if e in s])


# TODO: Filtre sur les genes en double
for c in chromosomes:
	for d in chromosomes:
		if c == d:
			continue
		i = chromosomes[c].intersection(chromosomes[d])
		chromosomes[c].difference_update(i)
		chromosomes[d].difference_update(i)

tout = set(utils.myMaths.flatten(chromosomes.values()))

# TODO: Filtre sur les genes sans assignation




# Impression
#s = ['Mouse', 'Rat', 'Human', 'Mouse', 'Dog']
genesAnc = utils.myGenomes.loadGenome(sys.argv[1])
for c in chromosomes:
	for i in chromosomes[c]:
		print c, " ".join(genesAnc.lstGenes[utils.myGenomes.AncestralGenome.defaultChr][i].names)
	#print "PanMouse\t%s\t%s" % (" ".join([str(x) for x in chromosomes[c]]), c)

