#! /users/ldog/muffato/python

__doc__ = """
	Dessine le karyotype d'un genome face a l'autre
"""

import sys
import utils.myTools
import utils.myGenomes
import utils.myKaryoDrawer

arguments = utils.myTools.checkArgs( \
	[("studiedGenome",file), ("referenceGenome",file)], \
	[("orthologuesList",str,""), ("includeGaps",bool,False), ("includeScaffolds",bool,False), ("includeRandoms",bool,False), ("reverse",bool,False)], \
	__doc__
)


# Le genome avec les couleurs
genome2 = utils.myGenomes.Genome(arguments["referenceGenome"])

# Si on a utilise /dev/null, c'est que le caryotype est donne sous un autre format
if len(genome2.dicGenes) == 0:

	f = utils.myTools.myOpenFile(arguments["studiedGenome"], "r")
	table12 = {}
	for (i,l) in enumerate(f):
		nbChr = i+1
		chrom = []
		t = l.replace('\n', '').split(";")
		for s in t:
			if len(s) == 0:
				continue
			elif len(s) == 1:
				nbChr = s[0]
			c = s.split()
			chrom.extend( [[(c[0],None)]] *  int(10*float(c[1])) )
		table12[nbChr] = list(enumerate(chrom))
		table12[nbChr].reverse()

	chr1 = sorted(table12.keys())

# Sinon, procedure normale: on charge le genome avec les orthologues
else:
	genome1 = utils.myGenomes.Genome(arguments["studiedGenome"])
	if arguments["reverse"]:
		(genome1,genome2) = (genome2,genome1)

	genesAnc = utils.myGenomes.Genome(arguments["orthologuesList"]) if arguments["orthologuesList"] != "" else None

# Les chromosomes a etudier
chr1 = genome1.lstChr
chr2 = genome2.lstChr
if arguments["includeScaffolds"]:
	chr1.extend(genome1.lstScaff)
	chr2.extend(genome2.lstScaff)
if arguments["includeRandoms"]:
	chr1.extend(genome1.lstRand)
	chr2.extend(genome2.lstRand)

table12 = genome1.buildOrthosTable(chr1, genome2, chr2, arguments["includeGaps"], genesAnc)

data = {}
for (c,l) in table12.iteritems():
	newl = []
	for (_,val) in l:
		if len(val) == 0:
			newl.append( None )
		else:
			newl.append( val[0][0] )
	data[c] = newl
	
print >> sys.stderr, "Affichage ...",
utils.myKaryoDrawer.drawKaryo(data, arguments)
print >> sys.stderr, "OK"

