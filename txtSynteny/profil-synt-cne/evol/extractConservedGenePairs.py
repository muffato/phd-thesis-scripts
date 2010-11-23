#!/usr/bin/env python2
# Copyright (c) 2010 IBENS/Dyogen : Matthieu MUFFATO, Alexandra LOUIS, Hugues ROEST CROLLIUS
# Mail : hrc@ens.fr or alouis@biologie.ens.fr
# Licences GLP v3 and CeCILL v2

__doc__ = """
	Renvoie les listes des devenirs de chaque gene le long des branches de l'arbre phylogenetique
"""

import sys
import collections
import itertools

import utils.myFile
import utils.myDiags
import utils.myMaths
import utils.myTools
import utils.myGenomes
import utils.myPhylTree


# Argument:
arguments = utils.myTools.checkArgs( \
	[("phylTree.conf",file)], \
	[("IN.genesFile",str,""), ("IN.ancGenesFile",str,""), ("IN.diagsFile",str,"")], \
	__doc__ \
)

# Chargement des tous les fichiers
###################################
phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])

todo = ["Homo-Pan-Gorilla_group", "Hominidae", "Catarrhini", "Euarchontoglires", "Boreoeutheria", "Theria", "Mammalia", "Amniota", "Tetrapoda", "Euteleostomi"]

lOK = set()
refGenes = utils.myGenomes.Genome(arguments["IN.genesFile"] % "Homo.sapiens")
for c in refGenes.lstGenes:
	for (g1,g2) in utils.myTools.myIterator.slidingTuple(refGenes.lstGenes[c]):
		lOK.add( (g1.names,g2.names, g1.strand,g2.strand, g1,g2) )
print >> sys.stderr, "ini", len(lOK)

lNC = set()
for anc in todo:

	# Chargement des diagonales
	dicDiags = {}
	for (i,l) in enumerate(utils.myFile.myTSV.readTabular(arguments["IN.diagsFile"] % anc, [str,int,(str,4)])):
		d = [int(x) for x in l[2].split()]
		s = [int(x) for x in l[3].split()]
		for (j,(g,s)) in enumerate(zip(d,s)):
			dicDiags[g] = (i, j, s)

	# Les noms des genes ancestraux
	ancGenes = utils.myGenomes.Genome(arguments["IN.ancGenesFile"] % anc)

	nbOK = 0
	nbNO = 0
	newOK = set()
	for (g1,g2, s1,s2, ref1,ref2) in lOK:
		g1 = ancGenes.getPosition(g1)
		g2 = ancGenes.getPosition(g2)
		data = (ref1.names[0],ref2.names[0], ref1.chromosome, ref1.end+1, ref2.beginning-1, ref1.strand,ref2.strand)

		if (len(g1) == 0) or (len(g2) == 0):
			print utils.myFile.myTSV.printLine((anc,"LOST", s1,s2) + data)
			print utils.myFile.myTSV.printLine((anc,"NON-CONSERVED", s1,s2) + data)
			lNC.add(data)
		else:
			assert len(g1) == len(g2) == 1
			g1 = g1.pop()[1]
			g2 = g2.pop()[1]
			(c1,i1,t1) = dicDiags[g1]
			(c2,i2,t2) = dicDiags[g2]

			assert ref1.beginning <= ref2.beginning
			assert ref1.chromosome == ref2.chromosome

			if (c1 == c2) and (abs(i1-i2) == 1):
				if i1 > i2:
					t1 = -t1
					t2 = -t2
				newOK.add( (ancGenes.lstGenes[None][g1].names,ancGenes.lstGenes[None][g2].names, t1,t2, ref1,ref2) )
				print utils.myFile.myTSV.printLine((anc,"CONSERVED", s1,s2, t1,t2) + data)
				nbOK += 1
			else:
				print utils.myFile.myTSV.printLine((anc,"BROKEN", s1,s2) + data)
				print utils.myFile.myTSV.printLine((anc,"NON-CONSERVED", s1,s2) + data)
				nbNO += 1
				lNC.add(data)

	for data in lNC:
		print utils.myFile.myTSV.printLine((anc,"CUMUL-NON-CONSERVED", s1,s2) + data)

	print >> sys.stderr, anc, nbOK, nbNO
	lOK = newOK

