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
	#return max(filsNO) == 0
	return (max(filsNO) == 0) or (outgroupNO == 0 and outgroupOK >= 5)

lst = utils.myDiags.DiagRepository()
for l in sys.stdin:
	d = l.split('\t')
	c = [int(x) for x in d[1].split()]
	a = [tuple(x.split('/')) for x in d[2].split()]
	lst.addDiag(c, a)
	
print >> sys.stderr, "lecture OK", lst.nbRealDiags()

#lst.combinDiags([['Mouse'],['Rat']], isOK)
print >> sys.stderr, "combin OK", lst.nbRealDiags()

#s = ['Human', 'Chimp', 'Macaque', 'Mouse', 'Dog']
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
#s = ['Mouse', 'Rat', 'Macaque', 'Mouse', 'Dog']
genesAnc = utils.myGenomes.loadGenome(sys.argv[1])
for c in chromosomes:
	for i in chromosomes[c]:
		print c, " ".join(genesAnc.lstGenes[utils.myGenomes.AncestralGenome.defaultChr][i].names)
	#print "PanMouse\t%s\t%s" % (" ".join([str(x) for x in chromosomes[c]]), c)

