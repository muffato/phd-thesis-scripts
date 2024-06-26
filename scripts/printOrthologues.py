#!/usr/bin/env python2
# Copyright (c) 2010 IBENS/Dyogen : Matthieu MUFFATO, Alexandra LOUIS, Hugues ROEST CROLLIUS
# Mail : hrc@ens.fr or alouis@biologie.ens.fr
# Licences GLP v3 and CeCILL v2

__doc__ = """
	Compare deux genomes et renvoie la liste des couples de genes orthologues avec les details sur leurs positions
"""

##################
# INITIALISATION #
##################

# Librairies
import sys
import utils.myGenomes
import utils.myTools



########
# MAIN #
########

# Arguments
arguments = utils.myTools.checkArgs( \
	[("studiedGenome",file), ("referenceGenome",file)], \
	[("orthologuesList",str,""), ("includeGaps",bool,False), ("includeScaffolds",bool,False), ("includeRandoms",bool,False), ("includeNones",bool,False), ("reverse",bool,False)], \
	__doc__
)


# Chargement des fichiers
genome1 = utils.myGenomes.Genome(arguments["studiedGenome"])
genome2 = utils.myGenomes.Genome(arguments["referenceGenome"])
if arguments["reverse"]:
	x = genome1
	genome1 = genome2
	genome2 = x
if arguments["orthologuesList"] != "":
	genesAnc = utils.myGenomes.Genome(arguments["orthologuesList"])
else:
	genesAnc = genome2

# Les chromosomes a etudier
chr1 = []
chr2 = []
chr1.extend(genome1.chrList[utils.myGenomes.ContigType.Chromosome]
chr2.extend(genome2.chrList[utils.myGenomes.ContigType.Chromosome]
if arguments["includeScaffolds"]:
	chr1.extend(genome1.chrList[utils.myGenomes.ContigType.Scaffold]
	chr2.extend(genome2.chrList[utils.myGenomes.ContigType.Scaffold]
if arguments["includeRandoms"]:
	chr1.extend(genome1.chrList[utils.myGenomes.ContigType.Random]
	chr2.extend(genome2.chrList[utils.myGenomes.ContigType.Random]
if arguments["includeNones"]:
	chr1.extend(genome1.chrList[utils.myGenomes.ContigType.None]
	chr2.extend(genome2.chrList[utils.myGenomes.ContigType.None]

table12 = genome1.buildOrthosTable(chr1, genome2, chr2, arguments["includeGaps"], genesAnc)

# Fichier avec les noms des paires de genes orthologues et leurs coordonnees
for c1 in chr1:
	for (i1,t) in table12[c1]:
		g1 = genome1.lstGenes[c1][i1]
		for (c2,i2) in t:
			g2 = genome2.lstGenes[c2][i2]
			print utils.myFile.myTSV.printLine( [c1,g1.beginning,g1.end,g1.strand,"/".join(g1.names), c2,g2.beginning,g2.end,g2.strand,"/".join(g2.names)] )

