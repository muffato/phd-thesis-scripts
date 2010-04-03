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
arguments = utils.myTools.checkArgs( [("refGenome",str), ("reconstructedGenome",str)], [], __doc__)

cat = enum.Enum('Nb', 'Parfaites', 'Voisines', 'MemeChr', 'DiffChr')

def do(genome, lstDiags):
	for c in lstDiags.lstGenes:
		#t = l.split("\t")
		#diag = [int(x) for x in t[2].split()]
		#strand = [int(x) for x in t[3].split()]
		#score = [int(x) for x in t[4].split()]

		#if len(diag) == 1:
		if len(lstDiags.lstGenes[c]) == 1:
			continue

		scoreDiag = dict.fromkeys(cat, 0)

		for (gene1,gene2) in utils.myTools.myIterator.slidingTuple(lstDiags.lstGenes[c]):
			g1 = gene1.names[0]
			s1 = gene1.strand
			g2 = gene2.names[0]
			s2 = gene2.strand
		
		#for (((g1,s1),(g2,s2)),s) in itertools.izip(utils.myTools.myIterator.slidingTuple(zip(diag,strand)),score):
		
			#print gene1, gene2

			(c1,i1) = genome.dicGenes.get(g1, genome.dicGenes["0"])
			t1 = genome.lstGenes[c1][i1].strand
			
			#(c2,i2) = genome.dicGenes[g2]
			(c2,i2) = genome.dicGenes.get(g2, genome.dicGenes["0"])
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
				print "GAP", abs(i1-i2)
				res = cat.MemeChr
			scoreDiag[res] += 1
			print "PAIR", res.index

		for x in xrange(4, 0, -1):
			if scoreDiag[cat[x]] != 0:
				res = cat[x]
				break
		for x in cat:
			scoresPaires[x] += scoreDiag[x]
		
		scoresDiags[cat.Nb] += 1
		scoresDiags[res] += 1

		print "DIAG", " ".join([str(scoreDiag[cat[x]]) for x in xrange(5)])


scoresPaires = dict.fromkeys(cat, 0)
scoresDiags = dict.fromkeys(cat, 0)
refPairs = 0

for i in xrange(10):

	genome = utils.myGenomes.Genome(arguments["refGenome"] % i)
	if not utils.myFile.hasAccess(arguments["reconstructedGenome"] % i):
		continue
	lstDiags = utils.myGenomes.Genome(arguments["reconstructedGenome"] % i)

	refPairs += sum(len(x)-1 for x in genome.lstGenes.itervalues())
	do(genome, lstDiags)


print "ALLDIAGS", " ".join([str(scoresDiags[cat[x]]) for x in xrange(5)])
print "ALLPAIRS", " ".join([str(scoresPaires[cat[x]]) for x in xrange(5)])
print "REFPAIRS", refPairs


