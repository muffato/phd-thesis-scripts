#! /users/ldog/muffato/python -OO

import sys
import utils.myMaths
import utils.myTools
import utils.myPhylTree

(noms_fichiers, options) = utils.myTools.checkArgs( [], [("start",int,0), ("end",int,0), ("phylTree",str,""), ("gc3.mase",str,"")], "Retrouve le GC ancestral" )

allBases = "ACGT"

for treeID in xrange(options["start"], options["end"]+1):
	
	print >> sys.stderr, treeID

	tree = utils.myPhylTree.PhylogeneticTree(options["phylTree"] % treeID, buildLinks=False)

	f = utils.myTools.myOpenFile(options["gc3.mase"] % treeID, "r")
	nom = None
	seq = {}
	for l in f:
		if l[0] == ";":
			continue
		if nom == None:
			nom = l[:-1]
		else:
			seq[nom] = l[:-1]
			n = len(seq[nom])
			nom = None
	f.close()

	res = []
	for i in xrange(n):
		proba = []
		for base in allBases:
			values = {}
			for s in seq:
				if i < len(seq[s]):
					c = seq[s][i].upper()
				else:
					print >> sys.stderr, "!1", s, i, n
					c = "-"
				if c == base:
					values[s] = 1
				elif c in allBases:
					values[s] = 0
			proba.append(tree.calcWeightedValue(values, -1, None, None))
		res.append(proba)
	for (ie,e) in enumerate(tree.allNames):
		# Pas d'interet
		if len(tree.species[e]) == 1:
			continue
		seq = ""
		nbM = 0
		nbN = 0
		proba = [""] * len(allBases)
		for r in res:
			rf = []
			for (ib,b) in enumerate(allBases):
				x = r[ib][ie]
				if x == int(x):
					x = int(x)
				rf.append( (x,b) )
				proba[ib] += "\t" + str(x)
			l = sorted(rf, reverse = True)
			if l[0][0] <= 0:
				# Pas de base alignee ici
				seq += "-"
				nbM += 1
			elif l[0][0] == l[1][0]:
				# Plusieurs nucleotides a probabilites egales
				seq += "N"
				nbN += 1
			else:
				# Une seule possibilite
				seq += l[0][1]
		
		print utils.myTools.printLine(["ID", treeID])
		print utils.myTools.printLine(["SPECIES", tree.species[e]])
		print utils.myTools.printLine(["LENGTH", n, n-nbM, n-nbM-nbN])
		print utils.myTools.printLine(["SEQ", seq])
		for (ib,b) in enumerate(allBases):
			print ("PROBA-%s" % b) + proba[ib]
		print

