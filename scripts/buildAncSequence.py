#! /users/ldog/muffato/python -OO

import sys
import utils.myMaths
import utils.myTools
import utils.myPhylTree

(noms_fichiers, options) = utils.myTools.checkArgs( [], [("start",int,0), ("end",int,0), ("phylTree",str,""), ("alignment-FASTA",str,"")], "Retrouve le GC ancestral" )

allBases = "ACGT"

for treeID in xrange(options["start"], options["end"]+1):
	
	print >> sys.stderr, treeID, "...",

	tree = utils.myPhylTree.PhylogeneticTree(options["phylTree"] % treeID, buildLinks=False)

	seq = utils.myGenomes.loadFastaFile(options["alignment-FASTA"] % treeID)

	res = []
	n = len(seq.values()[0])
	for i in xrange(n):
		proba = []
		for base in allBases:
			values = {}
			for (e,s) in seq.iteritems():
				if i < len(s):
					c = s[i].upper()
				else:
					print >> sys.stderr, "!1", e, s, i, n
					c = "-"
				if c == base:
					values[e] = 1
				elif c in allBases:
					values[e] = 0
			proba.append(tree.calcWeightedValue(values, -1, None, None))
		res.append(proba)

	for (ie,e) in enumerate(tree.allNames):
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
				if x >= 0:
					proba[ib] += " %g" % x
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
	
		print "ID\t%d" % treeID
		print "SPECIES\t%s" % " ".join(tree.species[e])
		print "LENGTH\t%d\t%d\t%d" % (n,n-nbM,n-nbM-nbN)
		print "SEQ\t%s" % seq
		for (ib,b) in enumerate(allBases):
			print "PROBA-%s\t%s" % (b,proba[ib][1:])
		print

	print >> sys.stderr, "OK"
