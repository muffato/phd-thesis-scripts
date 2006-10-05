#!/usr/bin/python2.4

__doc__ = """
Ce script cherche les suites de genes alignees entre deux genomes; des diagonales.
Il est capable de faire des diagonales de diagonales et ainsi de suite.
"""

# INITIALISATION #

# Librairies
import sys
import os

sys.path.append(os.environ['HOME'] + "/work/scripts/utils")
import myTools
import myOrthos
import myMaths


# Arguments
(noms_fichiers, options) = myTools.checkArgs( \
	["genome1", "genome2"], \
	[("orthologuesList",str,""), ("keepGaps",bool,True), ("fusionThreshold",int,-1)], \
	__doc__ \
)

genome1 = myOrthos.loadGenome(noms_fichiers[0])
genome2 = myOrthos.loadGenome(noms_fichiers[1])
if options["orthologuesList"] != "":
	genesAnc = myOrthos.AncestralGenome(options["orthologuesList"], False)
else:
	genesAnc = genome2

print >> sys.stderr, "Extraction des diagonales ",
nbOrthos = 0
nbDiag = 0.
nbDiag2 = 0.
lenDiag = 0.
lenDiag2 = 0.

for c1 in genome1.lstChr:
	for c2 in genome2.lstChr:
		tab1 = [-1 for x in range(len(genome1.lstGenes[c1]))]
		tab2 = [-2 for x in range(len(genome2.lstGenes[c2]))]
		for i in range(len(genome1.lstGenes[c1])):
			g1 = genome1.lstGenes[c1][i]
			for s in g1.names:
				if s not in genesAnc.dicGenes:
					continue
				(cA,iA) = genesAnc.dicGenes[s]
				gA = genesAnc.lstGenes[cA][iA]
				for ss in gA.names:
					if ss not in genome2.dicGenes:
						continue
					(cc2,i2) = genome2.dicGenes[ss]
					if cc2 != c2:
						continue
					tab1[i] = nbOrthos
					tab2[i2] = nbOrthos
					nbOrthos += 1
					break
				else:
					continue
				break
		if options["keepGaps"]:
			(t1,t2) = (tab1,tab2)
		else:
			(t1,t2) = [x for x in tab1 if x != -1], [x for x in tab2 if x != -2]
		diags = myMaths.extractDiags(t1, t2, options["fusionThreshold"])
		x = [len(d) for d in diags]
		for u in x:
			print u
		y = [l for l in x if l >= 2]
		nbDiag += len(x)
		nbDiag2 += len(y)
		lenDiag += sum(x)
		lenDiag2 += sum(y)

	sys.stderr.write(".")

print >> sys.stderr, " OK"
print >> sys.stderr, nbDiag, nbDiag2, lenDiag, lenDiag2
print >> sys.stderr, lenDiag/nbDiag
print >> sys.stderr, lenDiag2/nbDiag2

