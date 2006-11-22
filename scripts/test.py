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



def f2():
	f = open(sys.argv[1], 'r')

	lst = f.readlines()[7:-1]
	f.close()

	for i in range(len(lst)):
		continue
		c = lst[i].split()
		for j in range(len(c)):
			if float(c[j]) > 0.5 and float(c[j]) < 1.5:
				print i,j+i+1
			


def translateDiagToChrom():

	genesAnc = utils.myGenomes.loadGenome(sys.argv[1])

	nb = 1
	for l in sys.stdin:
		c = l.split()
		for i in c:
			print nb, " ".join(genesAnc.lstGenes[utils.myGenomes.AncestralGenome.defaultChr][int(i)].names)
		nb += 1

def comptePhyloNodes():
	phylTree = utils.myBioObjects.PhylogeneticTree(sys.argv[1])
	nbEsp = len(phylTree.getSpecies(phylTree.root))

	for anc in phylTree.items:
		groupes = [phylTree.getSpecies(e) for (e,_) in phylTree.items[anc]]
		fils = utils.myMaths.flatten(groupes)

		nbO = nbEsp-len(fils)
		nbA = len(groupes[0])
		nbB = len(groupes[1])
		print anc, nbA*nbB + nbA*nbO + nbB*nbO
	
def compteNbChangements():

	genome1 = utils.myOrthos.loadGenome(sys.argv[1])
	genome2 = utils.myOrthos.loadGenome(sys.argv[2])
	#genesAnc = utils.myOrthos.AncestralGenome(sys.argv[3], False)

	# On construit les couleurs
	res = {}
	somme = 0
	total = 0
	for c in genome1.lstChr:

		res[c] = []
		for gene in genome1.lstGenes[c]:
			
			tg = gene.names
			#tg = utils.myMaths.flatten([genesAnc.lstGenes[cc][ii].names for (cc,ii) in [genesAnc.dicGenes[g] for g in tg if g in genesAnc.dicGenes]])
			
			for g in tg:
				if g in genome2.dicGenes:
					col = genome2.dicGenes[g][0]
					break
			else:
				col = -1
			res[c].append(col)
			total += 1
			if col == c:
				somme += 1

	print >> sys.stderr, somme, total


def buildGraph():
	lst = []
	lst2 = []
	vois = {}
	for l in sys.stdin:
		#c = [int(x) for x in l.split('\t')[1].split()]
		c = [int(x) for x in l.split()]
		lst.append(c)
		lst2.append(set(c))
		if len(c) == 0:
			continue
		if len(c) == 1 and c[0] not in vois:
			vois[c[0]] = []
		else:
			for i in xrange(len(c)-1):
				x = c[i]
				y = c[i+1]
				if x not in vois:
					vois[x] = []
				if y not in vois:
					vois[y] = []
				vois[x].append(y)
				vois[y].append(x)
		
	print "graph {"
	combin = utils.myTools.myCombinator([])
	seuil = int(sys.argv[1])
	for i in xrange(len(lst)):
		for j in xrange(i):
			s = lst2[i].intersection(lst2[j])
			if len(s) == 0:
				continue
			nb = 0
			for x in s:
				nb += min(lst[i].count(x), lst[j].count(x))
			if nb > seuil:
				combin.addLink([i,j])
				print '%d -- %d [label="%d"]' % (i,j,nb)
				print '%d [label="%d.%d"]' % (i,i,len(lst[i]))
				print '%d [label="%d.%d"]' % (j,j,len(lst[j]))
	print "}"
	for g in combin:
		s = set(utils.myMaths.flatten([lst[i] for i in g]))
		#t = [len(set(vois[x])) for x in s]
		#if max(t) >= 3:
		#print sys.argv[2], " ".join([str(x) for x in s])
		

def buildCliques():
	lst = utils.myDiags.DiagRepository()
	combin = utils.myTools.myCombinator([])
	for l in sys.stdin:
		#c = [int(x) for x in l.split('\t')[1].split()]
		c = [int(x) for x in l.split()]
		lst.addDiag(c, [])
	print >> sys.stderr, "lecture OK"
	#lst.buildOverlap(75)
	#print >> sys.stderr, lst.nbRealDiags()
	#for (s,_,_) in lst:
	#	print " ".join([str(x) for x in s])
	#return
	#nb = lst.nbRealDiags()
	#while nb == lst.nbRealDiags():
	#	nb = lst.nbRealDiags()
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
	print >> sys.stderr, lst.nbRealDiags()
	for (s,_,_) in lst:
		print " ".join([str(x) for x in s])
		#lst = lst2
			

#buildGraph()
translateDiagToChrom()
#buildExtendedDiags()
#buildCliques()
