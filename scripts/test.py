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



def f1():
	tabGenesAnc = []

	for nom in sys.argv[2:]:
		tmp = utils.myGenomes.loadGenome(nom)
		del tmp.lstGenes
		tabGenesAnc.append(tmp)

	genomeAnc = utils.myGenomes.loadGenome(sys.argv[1])
	del genomeAnc.dicGenes
	print "\t\t%s" % "\t".join([str(x) for x in genomeAnc.lstChr])
	for g in genomeAnc:
		score = dict([(c,0) for c in genomeAnc.lstChr])
		for gen in tabGenesAnc:
			s = g.names[0]
			if s in gen.dicGenes:
				score[gen.dicGenes[s][0]] += 1
		print "%s\t%d\t%s" % (g.chromosome, g.beginning, "\t".join([str(score[x]) for x in score]))

			
def checkNewGenes():
	genesAnc1 = utils.myGenomes.loadGenome(sys.argv[1])
	genesAnc2 = utils.myGenomes.loadGenome(sys.argv[2])
	genesAnc3 = utils.myGenomes.loadGenome(sys.argv[3])
	for c1 in genesAnc1.lstChr:
		ens = set([])
		for g1 in genesAnc1.lstGenes[c1]:
			for s in g1.names:
				if s in genesAnc3.dicGenes:
					ens.add(genesAnc3.dicGenes[s])
		ens2 = []
		for (c,i) in ens:
			for s in genesAnc3.lstGenes[c][i].names:
				if s in genesAnc2.dicGenes:
					ens2.append(genesAnc2.dicGenes[s][0])
		ens2s = set(ens2)
		ens2s.discard("Un_random")
		if len(ens2s) > 1:
			print c1, len(genesAnc1.lstGenes[c1]), ens2
			
	return
	for l in sys.stdin:
		c = l.split()
		if c[0] not in genesAnc1.dicGenes:
			print l,
	





def translateDiagToChrom():

	genesAnc = utils.myGenomes.loadGenome(sys.argv[1])

	nb = 1
	for l in sys.stdin:
		#c = [int(x) for x in l.split('\t')[1].split()]
		c = [int(x) for x in l.split()]
		for i in c:
			print nb, " ".join(genesAnc.lstGenes[utils.myGenomes.AncestralGenome.defaultChr][i].names)
		nb += 1


def buildGraph():
	
	lst = utils.myDiags.DiagRepository()
	for l in sys.stdin:
		c = [int(x) for x in l.split('\t')[1].split()]
		#c = [int(x) for x in l.split()]
		lst.addDiag(c, [])

	combin = utils.myTools.myCombinator([])
	seuil = int(sys.argv[1])
	lst.buildOverlapScores()
	
	print "graph {"
	for i in xrange(len(lst.lstDiags)):
		#print '%d [label="%d.%d"]' % (i,i,len(lst.lstDiags[i]))
		for j in lst.overlapScores[i]:
			nb = lst.overlapScores[i][j]
			if nb >= seuil and i < j:
				combin.addLink([i,j])
				print '%d -- %d [label="%d"]' % (i,j,nb)
				print '%d [label="%d.%d"]' % (i,i,len(lst.lstDiags[i]))
				print '%d [label="%d.%d"]' % (j,j,len(lst.lstDiags[j]))
	print "}"

	for g in combin:
		s = set(utils.myMaths.flatten([lst.lstDiags[i] for i in g]))
		#t = [len(set(vois[x])) for x in s]
		#if max(t) >= 3:
		#print sys.argv[2], " ".join([str(x) for x in s])
		

def tryNewOverlap():
	lst = {}
	phylTree = utils.myBioObjects.PhylogeneticTree(sys.argv[1])
	for anc in phylTree.items:
		lst[anc] = []
	for l in sys.stdin:
		ct = l.split('\t')
		anc = ct[0]
		l = int(ct[1])
		d = [int(x) for x in ct[2].split()]
		esp = [tuple(x.split('.')) for x in ct[3].split()]
		lst[anc].append( (l, d, esp) )
		
	print >> sys.stderr, "lecture OK", " ".join(["%s:%d" % (anc,len(lst[anc])) for anc in diagEntry])
	
	def combin(diags, fils1, fils2, outgroup):
		print >> sys.stderr, "Combinaison de %s.%s^%s ..." % (fils1,fils2,outgroup)
		combin = utils.myTools.myCombinator([])
		for (i,j) in utils.myTools.myMatrixIterator(len(diags), len(diags), utils.myTools.myMatrixIterator.StrictUpperMatrix):
			commun = set([e for (e,_) in diagFinal[i][2].intersection(diagFinal[j][2])])
			fils1OK = len(commun.intersection(fils1)) > 0
			fils2OK = len(commun.intersection(fils2)) > 0
			outgroupOK = len(commun.intersection(outgroup)) > 0
			if (fils1OK and fils2OK) or ((fils1OK or fils2OK) and outgroupOK):
				combin.addLink([i,j])
		for g in combin:
			r = set(utils.myMaths.flatten([diagFinal[i][1] for i in g]))
			print " ".join([str(x) for x in r])

		print >> sys.stderr, "OK"
	
	def recCombin(anc):
		(c1,c2) = phylTree.getBranchesSpecies(anc)
		outgroup = lstEspeces.difference(c1+c2)
		combin(lst[anc], c1, c2, outgroup)
	
	lstEspeces = set(phylTree.getSpecies(phylTree.root))
	recCombin(phylTree.root)

#buildGraph()
translateDiagToChrom()
#buildExtendedDiags()
#buildCliques()
#tryNewOverlap()
