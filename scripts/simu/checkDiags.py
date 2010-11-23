#!/usr/bin/env python2
# Copyright (c) 2010 IBENS/Dyogen : Matthieu MUFFATO, Alexandra LOUIS, Hugues ROEST CROLLIUS
# Mail : hrc@ens.fr or alouis@biologie.ens.fr
# Licences GLP v3 and CeCILL v2

__doc__ = """
Prend un genome ancestral et la liste des diagonales inferees a partir des especes modernes
Renvoie le pourcentage de qualite des diagonales.
"""

import sys
import enum
import itertools

import utils.myFile
import utils.myTools
import utils.myGenomes
import utils.myPhylTree


# Arguments
arguments = utils.myTools.checkArgs( \
	[("phylTree.conf",file)], \
	[("sameStrand",bool,True), ("diagsFile",str,""), ("range",str,""), \
	("genesFile",str,"~/work/data/genes/genes.%s.list.bz2"), \
	("ancGenesFile",str,"~/work/data/ancGenes/ancGenes.%s.list.bz2")], \
	__doc__ \
)

phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])
cat = enum.Enum('Nb', 'Parfaites', 'Voisines', 'MemeChr', 'DiffChr')


def do(genome, ancGenes, lstDiags):
	
	genome = utils.myGenomes.Genome(genome)
	ancGenes = utils.myGenomes.Genome(ancGenes)
	print >> sys.stderr, "Opening diags from", lstDiags
	lstDiags = utils.myFile.openFile(lstDiags, "r")
	
	scoresPaires = dict.fromkeys(cat, 0)
	scoresDiags = dict.fromkeys(cat, 0)
	refPairs = sum(len(x)-1 for x in genome.lstGenes.itervalues())
	
	for l in lstDiags:
		t = l.split("\t")
		diag = [int(x) for x in t[2].split()]
		if arguments["sameStrand"]:
			strand = [int(x) for x in t[3].split()]
			score = [int(x) for x in t[4].split()]
		else:
			strand = [0] * len(diag)
			score = [int(x) for x in t[3].split()]

		if len(diag) == 1:
			continue

		scoreDiag = dict.fromkeys(cat, 0)
		for (((g1,s1),(g2,s2)),s) in itertools.izip(utils.myTools.myIterator.slidingTuple(zip(diag,strand)),score):
			
			l1 = genome.getPosition(ancGenes.lstGenes[None][g1].names)
			assert len(l1) == 1
			(c1,i1) = l1.pop()
			t1 = genome.lstGenes[c1][i1].strand
			
			l2 = genome.getPosition(ancGenes.lstGenes[None][g2].names)
			assert len(l2) == 1
			(c2,i2) = l2.pop()
			t2 = genome.lstGenes[c2][i2].strand

			scoreDiag[cat.Nb] += 1
			if c1 != c2:
				res = cat.DiffChr
			elif abs(i1-i2) == 1:
				if (t1*t2 == s1*s2) and (i2 == i1+s1*t1):
					res = cat.Parfaites
				else:
					res = cat.Voisines
			else:
				print "GAP", anc, abs(i1-i2), s
				res = cat.MemeChr
			scoreDiag[res] += 1
			print "PAIR", anc, res.index, s

		for x in xrange(4, 0, -1):
			if scoreDiag[cat[x]] != 0:
				res = cat[x]
				break
		for x in cat:
			scoresPaires[x] += scoreDiag[x]
		
		scoresDiags[cat.Nb] += 1
		scoresDiags[res] += 1

		print "DIAG", anc, " ".join([str(scoreDiag[cat[x]]) for x in xrange(5)])

	lstDiags.close()
	return (scoresDiags, scoresPaires, refPairs)

def printStatus(prefix, scoresDiags, scoresPaires, refPairs):
	print prefix + "DIAGS", anc, " ".join([str(scoresDiags[cat[x]]) for x in xrange(5)])
	print prefix + "PAIRS", anc, " ".join([str(scoresPaires[cat[x]]) for x in xrange(5)])
	print prefix + "REFPAIRS", anc, refPairs

def addScore(s1, s2):
	for x in cat:
		s1[x] += s2[x]


if len(arguments["range"]) == 0:
	for anc in phylTree.listAncestr:
		res = do(arguments["genesFile"] % phylTree.fileName[anc], \
			arguments["ancGenesFile"] % phylTree.fileName[anc], \
			arguments["diagsFile"] % phylTree.fileName[anc])
		printStatus("", *res)

else:
	f = utils.myFile.openFile(arguments["range"], "r")
	todo = []
	for l in f:
		todo.extend(l.split())
	f.close()
	for anc in phylTree.listAncestr:
		scoresDiags = dict.fromkeys(cat, 0)
		scoresPaires = dict.fromkeys(cat, 0)
		refPairs = 0

		for i in todo:

			res = do(arguments["genesFile"] % (i,phylTree.fileName[anc]), \
				arguments["ancGenesFile"] % (i,phylTree.fileName[anc]), \
				arguments["diagsFile"] % (i,phylTree.fileName[anc]))

			printStatus("", *res)
			
			addScore(scoresDiags, res[0])
			addScore(scoresPaires, res[1])
			refPairs += res[2]

		printStatus("ALL", scoresDiags, scoresPaires, refPairs)


