#!/usr/bin/env python2
# Copyright (c) 2010 IBENS/Dyogen : Matthieu MUFFATO, Alexandra LOUIS, Hugues ROEST CROLLIUS
# Mail : hrc@ens.fr or alouis@biologie.ens.fr
# Licences GLP v3 and CeCILL v2

__doc__ = """
	Lit les exons d'une espece et cree le genome comme liste d'intervalles "exon" / ("intron" | "intergenic")
"""

import sys
import collections

import utils.myMaths
import utils.myTools
import utils.myGenomes

arguments = utils.myTools.checkArgs( [("exonsFile",file)], [], __doc__ )


# Lecture des exons
allExons = collections.defaultdict(set)
for line in utils.myFile.myTSV.readTabular(arguments["exonsFile"], [(str,4), (int,2), str]):
	if line[2] == "protein_coding":
		allExons[utils.myGenomes.commonChrName(line[6])].add( (line[4],line[5],frozenset([line[0]])) )


# Fusion des exons chevauchants
for c in allExons:
	found = True
	removed = 0
	print >> sys.stderr, c, len(allExons[c]), ">",
	while found:
		found = False
		toskip = False
		for (a,b) in utils.myTools.myIterator.slidingTuple(sorted(allExons[c])):
			if toskip:
				toskip = False
				continue
			if (b[0] <= a[1]+1):
				found = True
				allExons[c].remove(a)
				allExons[c].remove(b)
				allExons[c].add( (min(a[0],b[0]), max(a[1],b[1]), a[2]|b[2]) )
				toskip = True
				removed += ((a[1]-a[0]+1)+(b[1]-b[0]+1))-(max(a[1],b[1])-min(a[0],b[0])+1)
	print >> sys.stderr, len(allExons[c]), "=", -removed, "/", utils.myMaths.myStats.txtSummary([b[1]-b[0] for b in allExons[c]])

# Impression des intervalles
for c in allExons:
	lst = sorted(allExons[c])
	last = None
	for a in lst:
		if last != None:
			if len(last[2] & a[2]) > 0:
				print (c, last[1]+1, a[0]-1, "intron", last[2]|a[2])
			else:
				print (c, last[1]+1, a[0]-1, "intergenic", (last[2],a[2]))
		last = a
		print (c, a[0], a[1], "exon", a[2])

