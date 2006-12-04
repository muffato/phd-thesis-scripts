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
		

def buildCliques():
	lst = utils.myDiags.DiagRepository()
	for l in sys.stdin:
		#c = [int(x) for x in l.split('\t')[1].split()]
		c = [int(x) for x in l.split()]
		lst.addDiag(c, [])
	print >> sys.stderr, "lecture OK", lst.nbRealDiag
	seuil = int(sys.argv[1])
	combin = utils.myTools.myCombinator([])
	
	nb = 0
	while nb != lst.nbRealDiag:
		nb = lst.nbRealDiag
		lst.buildCliques()
		if len(lst.cliquesList) <= seuil:
			continue
		cl = lst.cliquesList[-1]
		print >> sys.stderr, "L", [len(x) for x in lst.cliquesList],
		if len(cl) == 0:
			continue
		s = set([])
		for i in cl.pop():
			s.update(lst.lstDiagsSet[i])
		lst.addDiag(list(s), [])
		print >> sys.stderr, lst.nbRealDiag

	for (s,_,_) in lst:
		print " ".join([str(x) for x in s])
	return
	#lst.buildOverlap(75)
	#print >> sys.stderr, lst.nbRealDiag
	#for (s,_,_) in lst:
	#	print " ".join([str(x) for x in s])
	#return
	#nb = lst.nbRealDiag
	#while nb == lst.nbRealDiag:
	#	nb = lst.nbRealDiag
	for i in xrange(len(lst.lstDiags)):
		combin.addLink([i])
	lst.buildCliques()
	print >> sys.stderr, "L", [len(x) for x in lst.cliquesList],
	#return
	#print >> sys.stderr, "L", len(cl3)
	#cl3 = lst.cliquesList[-1]
	cl3b = []
	for cl3 in lst.cliquesList[4:]:
		for c in cl3:
			combin.addLink(c)
	print >> sys.stderr, "combin",
	for c in combin:
		s = set([])
		for i in c:
			s.update(lst.lstDiagsSet[i])
		cl3b.append(s)
	print >> sys.stderr, "mkSuperCliques=",len(cl3b),
	for s in cl3b:
		lst.addDiag(s, [])
		#	cl3b.append(s)
		#print >> sys.stderr, ".",
		#lst2 = utils.myDiags.DiagRepository()
		#lst2.addRepository(lst)
		#print >> sys.stderr, ".",
		#for s in cl3b:
		#	lst2.addDiag(list(s), [])
	print >> sys.stderr, lst.nbRealDiag
	for (s,_,_) in lst:
		print " ".join([str(x) for x in s])
		#lst = lst2
			
def tryNewOverlap():
	lst = []
	for l in sys.stdin:
		#c = [int(x) for x in l.split('\t')[1].split()]
		#a = l.split('/')
		c = [x for x in l.split()]
		lst.append( c )
		#lst.append( (c[0],c[1:]) )
		#lst.append(a[2].split(".") + a[5].split("."))
		#lst.addDiag(c, [])
	print >> sys.stderr, "lecture OK", len(lst)
	dic = {}
	for i in xrange(len(lst)):
		#(e,d) = lst[i]
		d = lst[i]
		#if e not in dic:
		#	dic[e] = {}
		for s in d:
			if s not in dic:
				dic[s] = []
			dic[s].append(i)
	print >> sys.stderr, "dic OK"
	combin = utils.myTools.myCombinator([[x] for x in xrange(len(lst))])
	for s in dic:
		combin.addLink(dic[s])
	print >> sys.stderr, "combin OK"
	genesAnc = utils.myGenomes.loadGenome(sys.argv[1])
	for g in combin:
	#for g in xrange(len(lst)):
		ens = set([])
		for i in g:
			ens.update([str(genesAnc.dicGenes.get(s,("",""))[1]) for s in lst[i]])
			#ens.update(lst[i])
		print " ".join(ens)
	print >> sys.stderr, "print OK"
	



#buildGraph()
translateDiagToChrom()
#buildExtendedDiags()
#buildCliques()
#tryNewOverlap()
