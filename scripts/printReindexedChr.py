#!/usr/bin/env python2
# Copyright (c) 2010 IBENS/Dyogen : Matthieu MUFFATO, Alexandra LOUIS, Hugues ROEST CROLLIUS
# Mail : hrc@ens.fr or alouis@biologie.ens.fr
# Licences GLP v3 and CeCILL v2

__doc__ = """
	Compare deux genomes et reordonne le genome 1 pour qu'il soit plus ressemblant au genome 2
"""

##################
# INITIALISATION #
##################

# Librairies
import sys
import utils.myGenomes
import utils.myTools


#############
# FONCTIONS #
#############

# Renvoie la liste des chromosomes orthologues
def getOrthosChr(table, chr):
	
	res = {}
	for c1 in chr:

		# On ne garde que les chromosomes orthologues majoritaires
		lst = [x for (x,_) in utils.myMaths.flatten(table[c1].itervalues())]
		count = [(lst.count(x),x) for x in set(lst)]
		count.sort()
		nb = (len(lst)*arguments["minHomology"])/100
		tmp = []
		for (n,c2) in count.__reversed__():
			tmp.append( (c2,n) )
			nb -= n
			if nb <= 0:
				break
		res[c1] = tmp
	return res




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

# On echange l'ordre des chromosomes pour que les deux genomes paraissent plus colineaires

# D'abord on fait la liste des paires de chromosomes homologues
res1 = getOrthosChr(table12, chr1)
besthits = []
for c1 in chr1:
	if len(res1[c1]) == 0:
		besthits.append( (100000,0,c1) )
	else:
		(c2,nb) = res[c1][0]
		besthits.append( (c2,-nb,c1) )

# On renvoie les chromosomes du genome 1 dans l'ordre des best hits avec le genome 2
besthits.sort()
#for i in xrange(len(besthits)):
#	(c2,nb,c1) = besthits[i]
for (i,(c2,nb,c1)) in enumerate(besthits):
	# Faut-il retourner le chromosome ?
	memeSens = 0
	# On restreint c1 a ses orthologues avec c2
	tmp = []
	for i1 in xrange(len(genome1.lstGenes[c1])):
		tg2 = [i2 for (cc2,i2) in table12[c1].get(i1,[]) if cc2 == c2]
		if len(tg2) > 0:
			tmp.append(tg2)
	# Le test de colinearite
	for j in xrange(len(tmp)-1):
		if max(tmp[j]) < min(tmp[j+1]):
			memeSens += 1
		elif min(tmp[j]) > max(tmp[j+1]):
			memeSens -= 1
	# Le resultat
	if memeSens > 0:
		res = genome1.lstGenes[c1].__iter__()
	else:
		res = genome1.lstGenes[c1].__reversed__()
	for g in res:
		print i+1, " ".join(g.names)

