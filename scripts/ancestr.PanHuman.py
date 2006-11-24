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
	a = [tuple(x.split('/')) for x in d[2].split()]
	lst.addDiag(c, a)
	
print >> sys.stderr, "lecture OK", lst.nbRealDiags()

#lst.combinDiags([['Human'],['Chimp']], isOK)
print >> sys.stderr, "combin OK", lst.nbRealDiags()

# 22 + X + Y + MT + 2b
mapping = {}
chromosomes = {}
for i in xrange(1, 22+1):
	s = str(i)
	mapping[('Human',s)] = s
	mapping[('Chimp',s)] = s
	chromosomes[s] = set([])
mapping[('Human','X')] = 'X'
mapping[('Chimp','X')] = 'X'
chromosomes['X'] = set([])
mapping[('Human','Y')] = 'Y'
mapping[('Chimp','Y')] = 'Y'
chromosomes['Y'] = set([])
mapping[('Human','MT')] = 'MT'
mapping[('Chimp','MT')] = 'MT'
chromosomes['MT'] = set([])
mapping[('Chimp','2a')] = '2'
mapping[('Chimp','2b')] = '2'


for (d,_,l) in lst:
	hum = [c for (e,c) in l if e =='Human']
	chi = [c for (e,c) in l if e =='Chimp']
	poss = set([])
	if len(hum) == 1:
		c = ('Human',hum[0])
		if c in mapping:
			poss.add(mapping[c])
	if len(chi) == 1:
		c = ('Chimp',chi[0])
		if c in mapping:
			poss.add(mapping[c])
	if len(poss) == 1:
		sol = poss.pop()
		#print "on met %s dans %s" % ( " ".join(["%s/%s" % (e,c) for (e,c) in l if e in s]), sol)
		chromosomes[sol].update(d)
	#else:
		#print "PanHuman\t%s\t%s" % (" ".join([str(x) for x in d]), " ".join(["%s/%s" % (e,c) for (e,c) in l if e in s]))
	#print len(d), " ".join(["%s/%s" % (e,c) for (e,c) in l if e in s])


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


for i in xrange(len(lstGenesAnc)):
	g = lstGenesAnc[i]
	if i in tout:
		continue
	hum = [gH.dicGenes[x][0] for x in g.names if x in gH.dicGenes]
	chi = [gC.dicGenes[x][0] for x in g.names if x in gC.dicGenes]
	poss = set([])
	if len(hum) == 1:
		c = ('Human',hum[0])
		if c in mapping:
			poss.add(mapping[c])
	if len(chi) == 1:
		c = ('Chimp',chi[0])
		if c in mapping:
			poss.add(mapping[c])
	if len(poss) == 1:
		sol = poss.pop()
		chromosomes[sol].add(i)
	#else:
	#	print g.names


#sys.exit(0)

# Impression
#s = ['Human', 'Chimp', 'Macaque', 'Mouse', 'Dog']
for c in chromosomes:
	for i in chromosomes[c]:
		print c, " ".join(lstGenesAnc[i].names)
	#print "PanHuman\t%s\t%s" % (" ".join([str(x) for x in chromosomes[c]]), c)

