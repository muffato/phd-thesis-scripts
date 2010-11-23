#!/usr/bin/env python2
# Copyright (c) 2010 IBENS/Dyogen : Matthieu MUFFATO, Alexandra LOUIS, Hugues ROEST CROLLIUS
# Mail : hrc@ens.fr or alouis@biologie.ens.fr
# Licences GLP v3 and CeCILL v2

import sys
import numpy
import scipy.linalg.matfuncs
import itertools

import utils.myTools
import utils.myGenomes
import utils.myPhylTree

arguments = utils.myTools.checkArgs( [("phylTree",str), ("alignment-FASTA",str)], [("batch",str,""), ("alphabet",str,"ACGTN")], "Reconstruit la sequence ancestrale")

allBases = arguments["alphabet"][:-1]
allBasesS = set(allBases)
nbBases = len(allBases)
indBases = range(nbBases)

mutations = numpy.identity(nbBases)
dicRates = {}

unknownBase = arguments["alphabet"][-1]
unknownBaseCost = [1. / nbBases] * nbBases


def calcWeightedValue(tree, values, notdefined):
	
	d = len(tree.allNames) * nbBases
	matriceA = numpy.identity(d)
	matriceB = numpy.empty((d,))
	# Par defaut, les resultats seront "notdefined"
	matriceB.fill(notdefined)

	def addFixedValue(e, v):
		for b in indBases:
			matriceB[e+b] = v[b]

	def addLink(e1, e2, l):
		mat = dicRates[l]
		for (b1,b2) in itertools.product(indBases, indBases):
			matriceA[e1+b1][e2+b2] = mat[b1][b2]

	# Construit recursivement la matrice
	def recBuild(ianc, father, length):
		anc = tree.allNames[ianc]

		# Appels recursifs: on selectionne les fils OK (en supposant qu'on le soit nous meme)
		items = [(e*nbBases,p) for (e,p) in tree.numItems[ianc] if recBuild(e, ianc, p)]
		# Le pere si disponible
		if father != None:
			items.append( (father*nbBases,length) )
		ianc *= nbBases

		if anc in values:
			# Si on a une valeur, on met "x = val"
			addFixedValue(ianc, values[anc])
			return True

		elif len(items) >= 2:
			# S'il a suffisament de voisins, on ecrit l'equation
			s = 0.
			for (e,p) in items:
				addLink(ianc, e, p)
				s += 1./max(p,0.00001)
			for b in indBases:
				matriceA[ianc+b][ianc+b] = -s
				matriceB[ianc+b] = 0
			return True
		else:
			return False
	
	# Construction de la matrice
	if len(values) == 0:
		return matriceB
		
	rootNode = tree.indNames[tree.lastCommonAncestor(values.keys())]
	recBuild(rootNode, None, 0)
	# Resolution de l'equation
	res = numpy.linalg.solve(matriceA, matriceB)
	return res



def doWork(tree, seq):

	# Hack amphioxus
	for (n,s) in seq.iteritems():
		if n.startswith("JGI.Brafl1") or n.startswith("ENSTSYT0") or n.startswith("ENSCHOT0") or n.startswith("ENSVPAT0"):
			seq[n] = "-" * len(s)

	dicRates.clear()
	def do(node):
		for (e,l) in tree.items.get(node, []):
			dicRates[l] = scipy.linalg.matfuncs.expm(l*mutations)/max(l,0.00001)
			dicRates[l] = numpy.identity(nbBases)/max(l,0.00001)
			do(e)
	do(tree.root)

	res = []
	n = len(seq.values()[0])
	for i in xrange(n):

		values = {}
		for (e,s) in seq.iteritems():
			c = s[i]
			if c == unknownBase:
				values[e] = unknownBaseCost
			elif c in allBasesS:
				values[e] = [float(c == b) for b in allBases]
		
		res.append(calcWeightedValue(tree, values, -1.))

	for e in tree.listAncestr:
		ie = tree.indNames[e] * nbBases
		seq = ""
		nbM = 0
		nbN = 0
		proba = [""] * nbBases
		for r in res:
			rf = []
			for (ib,b) in enumerate(allBases):
				x = r[ie+ib]
				if x == int(x):
					x = int(x)
				rf.append( (x,b) )
				if x >= 0:
					proba[ib] += " %g" % x
			l = sorted(rf, reverse = True)
			if l[0][0] <= 0:
				# Pas de base alignee ici
				seq += "-"
				nbM += 1
			elif l[0][0] == l[1][0]:
				# Plusieurs nucleotides a probabilites egales
				seq += unknownBase
				nbN += 1
			else:
				# Une seule possibilite
				seq += l[0][1]

		print "NAME\t%s" % e
		print "SPECIES\t%s" % " ".join(tree.species[e])
		print "LENGTH\t%d\t%d\t%d" % (n,n-nbM,n-nbM-nbN)
		print "SEQ\t%s" % seq
		for (ib,b) in enumerate(allBases):
			print "PROBA-%s\t%s" % (b,proba[ib][1:])
		print

if len(arguments["batch"]) == 0:
	tree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree"])
	seq = utils.myGenomes.myFASTA.loadFile(arguments["alignment-FASTA"])
	doWork(tree, seq)
	print >> sys.stderr, "OK"
else:
	for treeID in utils.myTools.getRange(arguments["batch"]):
		print >> sys.stderr, treeID, "...",
		tree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree"] % treeID)
		seq = utils.myGenomes.myFASTA.loadFile(arguments["alignment-FASTA"] % treeID)
		print "ID\t%d" % treeID
		doWork(tree, seq)
	print >> sys.stderr, "OK"

