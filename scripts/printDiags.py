#! /users/ldog/muffato/python -OO

__doc__ = """
	Renvoie les diagonales de genes entre deux genomes
"""

##################
# INITIALISATION #
##################

# Librairies
import sys
import utils.myGenomes
import utils.myTools
import utils.myDiags



########
# MAIN #
########

# Arguments
arguments = utils.myTools.checkArgs( \
	[("studiedGenome",file), ("referenceGenome",file)], \
	[("orthologuesList",str,""), ("includeGaps",bool,False), ("includeScaffolds",bool,False), ("includeRandoms",bool,False), \
	("reverse",bool,False), ("fusionThreshold",int,-1), ("minimalLength",int,2), ("sameStrand",bool,True),], \
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
chr1 = genome1.lstChr
chr2 = genome2.lstChr
if arguments["includeScaffolds"]:
	chr1.extend(genome1.lstScaff)
	chr2.extend(genome2.lstScaff)
if arguments["includeRandoms"]:
	chr1.extend(genome1.lstRand)
	chr2.extend(genome2.lstRand)

print >> sys.stderr, "Extraction des diagonales ",
for ((c1,d1),(c2,d2),ds) in utils.myDiags.calcDiags(genome1, genome2, genesAnc, arguments["minimalLength"], arguments["fusionThreshold"], arguments["sameStrand"], not arguments["includeGaps"]):
	
	print utils.myTools.printLine([len(ds), c1," ".join([genome1.lstGenes[c1][i1].names[0] for i1 in d1]), c2," ".join([genome2.lstGenes[c2][i2].names[0] for i2 in d2]), ds])
	#res.append(str(len(ds)))
	#res.append(str(c1))
	#res.append(" ".join([genome1.lstGenes[c1][i1].names[0] for i1 in d1]))
	#res.append(str(c2))
	#res.append(" ".join([genome2.lstGenes[c2][i2].names[0] for i2 in d2]))
	#res.append(" ".join([str(x) for x in ds]))
	#print '\t'.join(res)

print >> sys.stderr, "OK"
