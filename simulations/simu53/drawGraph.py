#! /users/ldog/muffato/python

__doc__ = """
Prend un genome ancestral et la liste des diagonales inferees a partir des especes modernes
Renvoie le pourcentage de qualite des diagonales.
"""

import sys
import enum
import itertools
import collections

import utils.myFile
import utils.myDiags
import utils.myTools
import utils.myGenomes
import utils.myPhylTree


arguments = utils.myTools.checkArgs( [("pairwiseDiagsFile",file), ("ancDiagsFile",file), ("genesFile",file), ("ancGenesFile",file)], [], __doc__ )

genome = utils.myGenomes.Genome(arguments["genesFile"])
ancGenes = utils.myGenomes.Genome(arguments["ancGenesFile"])
lstDiagsA = utils.myFile.openFile(arguments["ancDiagsFile"], "r")
lstDiagsPW = utils.myFile.openFile(arguments["pairwiseDiagsFile"], "r")

print >> sys.stderr, "Expected pairs ...",
expected = set()
for c in genome.lstChr:
	for (g1,g2) in utils.myTools.myIterator.slidingTuple(genome.lstGenes[c]):
		(_,i1) = ancGenes.dicGenes[g1.names[0]]
		(_,i2) = ancGenes.dicGenes[g2.names[0]]
		expected.add( ((i1,g1.strand),(i2,g2.strand)) )
		expected.add( ((i2,-g2.strand),(i1,-g1.strand)) )
print >> sys.stderr, "OK"

print >> sys.stderr, "Compared pairs ...",
compared = collections.defaultdict(int)
for l in lstDiagsPW:
	t = l.split("\t")
	diag = [int(x) for x in t[8].split()]
	strand = [int(x) for x in t[9].split()]
	for ((g1,s1),(g2,s2)) in utils.myTools.myIterator.slidingTuple(zip(diag,strand)):
		if g1 == g2:
			continue
		compared[ ((g1,s1),(g2,s2)) ] += 1
		compared[ ((g2,-s2),(g1,-s1)) ] += 1
lstDiagsPW.close()
print >> sys.stderr, "OK"

print >> sys.stderr, "Reconstructed pairs ...",
reconstructed = {}
for l in lstDiagsA:
	t = l.split("\t")
	diag = [int(x) for x in t[2].split()]
	strand = [int(x) for x in t[3].split()]
	score = [int(x) for x in t[4].split()]
	for (((g1,s1),(g2,s2)),s) in itertools.izip(utils.myTools.myIterator.slidingTuple(zip(diag,strand)),score):
		reconstructed[ ((g1,s1),(g2,s2)) ] = s
		reconstructed[ ((g2,-s2),(g1,-s1)) ] = s
lstDiagsA.close()
print >> sys.stderr, "OK"

print >> sys.stderr, "All edges ...",
allEdges = collections.defaultdict(list)
def addEdge((g1,s1), (g2,s2), s, color):
	allEdges[(g1,s1)].append( ((g2,s2), s, color) )

for ((g1,s1),(g2,s2)) in expected:
	addEdge((g1,s1), (g2,s2), "", "black")
for ((g1,s1),(g2,s2)) in reconstructed:
	addEdge((g1,s1), (g2,s2), reconstructed[((g1,s1),(g2,s2))], "red")
for ((g1,s1),(g2,s2)) in compared:
	addEdge((g1,s1), (g2,s2), compared[((g1,s1),(g2,s2))], "green")
print >> sys.stderr, "OK"

print >> sys.stderr, "Weird cases ...",
combin = utils.myTools.myCombinator()
union = expected.union(reconstructed, compared)
for ((g1,s1),(g2,s2)) in union.difference(expected).union(union.difference(reconstructed), union.difference(compared)):
	combin.addLink([(g1,s1), (g2,s2)])
print >> sys.stderr, "OK"

print >> sys.stderr, "Printing ...",
def drawEdge((g1,s1), (g2,s2), s, color):
	print '"%d/%d" -> "%d/%d" [color="%s", label="%s"]' % (g1,s1, g2,s2, color, s)

for grp in combin:
	print "digraph {"
	for (g1,s1) in grp:
		for ((g2,s2), s, color) in allEdges[(g1,s1)]:
			drawEdge((g1,s1), (g2,s2), s, color)
	print "}"
print >> sys.stderr, "OK"

