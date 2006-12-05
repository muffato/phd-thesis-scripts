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
		

def tryNewOverlap():
	lst = []
	genesAnc = utils.myGenomes.loadGenome(sys.argv[1])
	for l in sys.stdin:
		ct = l.split('\t')
		a1 = ct[4].split()
		a2 = ct[7].split()
		c = a1+a2
		a1 = [genesAnc.dicGenes.get(s,("",""))[1] for s in a1]
		a2 = [genesAnc.dicGenes.get(s,("",""))[1] for s in a2]
		if "" in a1:
			anc = a2
		else:
			anc = a1
		lst.append( (c, anc, (ct[2],ct[3],a1), (ct[5],ct[6],a2)) )
		
	print >> sys.stderr, "lecture OK", len(lst)
	
	dic = {}
	for i in xrange(len(lst)):
		d = lst[i][0]
		for s in d:
			if s not in dic:
				dic[s] = []
			dic[s].append(i)
	print >> sys.stderr, "dic OK"
	
	combin = utils.myTools.myCombinator([[x] for x in xrange(len(lst))])
	for s in dic:
		combin.addLink(dic[s])
	print >> sys.stderr, "combin OK"

	diagFinal = []
	for g in combin:
		for res in utils.myDiags.getLongestPath([lst[i][1] for i in g]):
			# TODO Projeter sur les chromosomes modernes
			# on teste chaque diagonale (inclusion d'au moins deux elements dans le resultat)
			ok = set([])
			for i in g:
				d = lst[i][1]
				flag = False
				for j in xrange(len(d)-1):
					if (d[j] not in res[0]) or (d[j+1] not in res[0]):
						continue
					if abs(res[0].index(d[j])-res[0].index(d[j+1])) == 1:
						flag = True
						break
				if flag:
					ok.add( (lst[i][2][0],lst[i][2][1]) )
					ok.add( (lst[i][3][0],lst[i][3][1]) )
			print "%d\t%s\t%s" % (len(res[0]), " ".join([str(x) for x in res[0]]), " ".join(["%s.%s" % (e,c) for (e,c) in ok]))
			diagFinal.append( (len(res[0]), res[0], ok) )
			
	print >> sys.stderr, "print OK"
	return
	combin = utils.myTools.myCombinator([])
	fils = [['Human'], ['Chimp']]
	for (i,j) in utils.myTools.myMatrixIterator(len(diagFinal), len(diagFinal), utils.myTools.myMatrixIterator.StrictUpperMatrix):
		commun = set([e for (e,_) in diagFinal[i][2].intersection(diagFinal[j][2])])
		diff = set([e for (e,_) in diagFinal[i][2].symmetric_difference(diagFinal[j][2])])
		filsOK = [len(commun.intersection(x)) for x in fils]
		filsNO = [len(diff.intersection(x)) for x in fils]
		outgroupOK = len(commun) - sum(filsOK)
		outgroupNO = len(diff) - sum(filsNO)
		#print diagFinal[i][2], diagFinal[j][2]
		#print commun, diff, filsOK, filsNO, outgroupOK, outgroupNO
		if max(filsNO) == 0:
		#if outgroupOK >= 1 and min(filsOK) >= 1 and max(filsOK) >= 1:
			combin.addLink([i,j])
	for g in combin:
		r = utils.myMaths.flatten([diagFinal[i][1] for i in g])
		print " ".join([str(x) for x in r])

	print >> sys.stderr, "combin OK"

#buildGraph()
#translateDiagToChrom()
#buildExtendedDiags()
#buildCliques()
tryNewOverlap()
