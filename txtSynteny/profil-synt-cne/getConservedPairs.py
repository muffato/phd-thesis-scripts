#!/usr/bin/env python2

__doc__ = """
	Lit les CNE et transforme les positions genomiques en positions indexees
"""

import sys
import bisect
import operator
import itertools
import collections

import utils.myTools
import utils.myGenomes

arguments = utils.myTools.checkArgs( [], [("ancDiags",str,""), ("genesFile",str,""), ("ancGenesFile",str,"")], __doc__)


# Chargement des diagonales
diags = collections.defaultdict(set)
anc = set()
print >> sys.stderr, "Loading synteny blocks ...",
ft = utils.myFile.myTSV.reader(sys.stdin)
for t in ft.csvobject:
	for (g1,g2) in utils.myTools.myIterator.slidingTuple(t[4].split()):
		diags[(t[2],t[5])].add( (g1,g2) )
		diags[(t[2],t[5])].add( (g2,g1) )
	for (g1,g2) in utils.myTools.myIterator.slidingTuple(t[7].split()):
		diags[(t[5],t[2])].add( (g1,g2) )
		diags[(t[5],t[2])].add( (g2,g1) )
	anc.add(t[0])
ft.file.close()
print >> sys.stderr, "OK"

esp = set(x[0] for x in diags)
for (e3,e1,e2) in itertools.permutations(esp, 3):
	nbTot = 0
	nbOK = 0
	for x in diags[(e1,e2)]:
		nbTot += 1
		if x in diags[(e1,e3)]:
			nbOK += 1
	print "Assuming %s/%s: %d/%d = %.2f%% conserved in %s" % (e1,e2,nbOK/2,nbTot/2,float(nbOK*100)/nbTot,e3)

genomes = {}
for e in esp:
	genomes[e] = utils.myGenomes.Genome(arguments["genesFile"] % e.replace(' ', '.'))

for a in anc:
	ancGenes = utils.myGenomes.Genome(arguments["ancGenesFile"] % a)
	ft = utils.myFile.myTSV.reader(arguments["ancDiags"] % a)
	nbTot = collections.defaultdict(int)
	nbOK = collections.defaultdict(int)
	for t in ft.csvobject:
		t = [int(x) for x in t[2].split()]
		for e in esp:
			for (g1,g2) in utils.myTools.myIterator.slidingTuple(t):
				g1 = genomes[e].getPosition(ancGenes.lstGenes[None][g1].names)
				g2 = genomes[e].getPosition(ancGenes.lstGenes[None][g2].names)
				if (len(g1) != 0) and (len(g2) != 0):
					nbTot[e] += 1
					for ((c1,i1),(c2,i2)) in itertools.product(g1, g2):
						if (c1==c2) and (abs(i1-i2)==1):
							nbOK[e] += 1
							break
	for e in nbTot:
		print "Assuming %s: %d/%d = %.2f%% conserved in %s" % (a, nbOK[e], nbTot[e], float(nbOK[e]*100)/nbTot[e], e)
	ft.file.close()

for (e1,e2) in itertools.permutations(esp, 2):
	nbTot = 0
	nbOK = 0
	for c1 in genomes[e1].lstGenes.itervalues():
		for (g1,g2) in utils.myTools.myIterator.slidingTuple(c1):
			nbTot += 1
			if (g1.names[0],g2.names[0]) in diags[(e1,e2)]:
				nbOK += 1
	print "Between %s and %s: %d/%d = %.2f%% conserved" % (e1,e2,nbOK,nbTot,float(nbOK*100)/nbTot)

for (e3,e1,e2) in itertools.permutations(esp, 3):
	nbTot = 0
	nbOK = 0
	for c1 in genomes[e1].lstGenes.itervalues():
		for (g1,g2) in utils.myTools.myIterator.slidingTuple(c1):
			if (g1.names[0],g2.names[0]) in diags[(e1,e2)]:
				nbTot += 1
				if (g1.names[0],g2.names[0]) in diags[(e1,e3)]:
					nbOK += 1
	#print "Assuming %s/%s: %d/%d = %.2f%% conserved in %s" % (e1,e2,nbOK,nbTot,float(nbOK*100)/nbTot,e3)

for (e3,e1,e2) in itertools.permutations(esp, 3):
	nbTot = 0
	nbOK = 0
	for c1 in genomes[e1].lstGenes.itervalues():
		for (g1,g2) in utils.myTools.myIterator.slidingTuple(c1):
			if (g1.names[0],g2.names[0]) not in diags[(e1,e2)]:
				nbTot += 1
				if (g1.names[0],g2.names[0]) in diags[(e1,e3)]:
					nbOK += 1
	print "Refusing %s/%s: %d/%d = %.2f%% conserved in %s" % (e1,e2,nbOK,nbTot,float(nbOK*100)/nbTot,e3)

