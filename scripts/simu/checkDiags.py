#! /users/ldog/muffato/python

__doc__ = """
Prend un genome ancestral et la liste des diagonales inferees a partir des especes modernes
Renvoie le pourcentage de qualite des diagonales.
"""

import sys
import enum
import itertools

import utils.myFile
import utils.myDiags
import utils.myTools
import utils.myGenomes
import utils.myPhylTree


# Arguments
arguments = utils.myTools.checkArgs( \
	[("phylTree.conf",file)], \
	[("diagsFile",str,""), ("range",str,""), \
	("genesFile",str,"~/work/data/genes/genes.%s.list.bz2"), \
	("ancGenesFile",str,"~/work/data/ancGenes/ancGenes.%s.list.bz2")], \
	__doc__ \
)

phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])
cat = enum.Enum('Nb', 'Parfaites', 'Voisines', 'MemeChr', 'DiffChr')


def do(genome, ancGenes, lstDiags):
	for l in lstDiags:
		t = l.split("\t")
		diag = [int(x) for x in t[2].split()]
		strand = [int(x) for x in t[3].split()]
		score = [int(x) for x in t[4].split()]

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




if len(arguments["range"]) == 0:
	for anc in phylTree.listAncestr:
		scoresPaires = dict.fromkeys(cat, 0)
		scoresDiags = dict.fromkeys(cat, 0)
		refPairs = 0

		genome = utils.myGenomes.Genome(arguments["genesFile"] % phylTree.fileName[anc])
		ancGenes = utils.myGenomes.Genome(arguments["ancGenesFile"] % phylTree.fileName[anc])
		lstDiags = utils.myFile.openFile(arguments["diagsFile"] % phylTree.fileName[anc], "r")

		refPairs += sum(len(x)-1 for x in genome.lstGenes.itervalues())
		do(genome, ancGenes, lstDiags)

		lstDiags.close()

		print "ALLDIAGS", anc, " ".join([str(scoresDiags[cat[x]]) for x in xrange(5)])
		print "ALLPAIRS", anc, " ".join([str(scoresPaires[cat[x]]) for x in xrange(5)])
		print "REFPAIRS", anc, refPairs

else:
	todo = utils.myTools.getRange(arguments["range"])
	for anc in phylTree.listAncestr:
		scoresPaires = dict.fromkeys(cat, 0)
		scoresDiags = dict.fromkeys(cat, 0)
		refPairs = 0

		for i in todo:

			genome = utils.myGenomes.Genome(arguments["genesFile"] % (i,phylTree.fileName[anc]))
			ancGenes = utils.myGenomes.Genome(arguments["ancGenesFile"] % (i,phylTree.fileName[anc]))
			lstDiags = utils.myFile.openFile(arguments["diagsFile"] % (i,phylTree.fileName[anc]), "r")

			refPairs += sum(len(x)-1 for x in genome.lstGenes.itervalues())
			do(genome, ancGenes, lstDiags)

			lstDiags.close()

		print "ALLDIAGS", anc, " ".join([str(scoresDiags[cat[x]]) for x in xrange(5)])
		print "ALLPAIRS", anc, " ".join([str(scoresPaires[cat[x]]) for x in xrange(5)])
		print "REFPAIRS", refPairs


