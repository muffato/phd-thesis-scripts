#! /users/ldog/muffato/python

__doc__ = """
	Lit les CNE et transforme les positions genomiques en positions indexees
"""

import sys
import bisect
import operator
import itertools
import utils.myTools
import utils.myGenomes

arguments = utils.myTools.checkArgs( [("KX_HDMC_blocs.list",file), ("KX_HDMC_especes.list",file)], [("genesFile",str,"")], __doc__ )

dicSpeciesID = {'human': 60}
dicRealSpeciesNames = {60: 'Homo.sapiens'}


# Chargement des genomes des especes modernes
dicGenomes = {}
posCNE = utils.myTools.defaultdict(dict)


for (species_id,species_name) in dicRealSpeciesNames.iteritems():
	
	# Chargement du genome
	lst = utils.myTools.defaultdict(list)
	genome = utils.myGenomes.Genome(arguments["genesFile"] % species_name)
	for gene in genome:
		lst[gene.chromosome].append( (gene.beginning,gene.end,gene.names[0]) )
	
	# Tri des chromosomes
	for x in lst.itervalues():
		x.sort()

	dicGenomes[species_id] = lst


# Chargement des CNE
CNE = {}
for line in utils.myFile.myTSV.readTabular(arguments["KX_HDMC_blocs.list"], [str,int,int,int,str]):
	# Les infos generales du CNE
	CNE[line[0]] = (line[1:],[])


for line in utils.myFile.myTSV.readTabular(arguments["KX_HDMC_especes.list"], [str,str,int,int,str,str,int,str]):
	
	if (line[0] in dicSpeciesID) and (int(line[6]) >= 1):
		
		# Les infos de chaque CNE pour chaque espece
		species_id = dicSpeciesID[line[0]]
		c = utils.myGenomes.commonChrName(line[1])

		i = bisect.bisect_left(dicGenomes[species_id][c], (line[2],line[3]))
		posCNE[c][i] = posCNE[c].get(i,0) + (line[3]-line[2]+1)


posBlocks = utils.myTools.defaultdict(set)
for t in utils.myFile.myTSV.readTabular(sys.stdin, [str]*10):
	if "sapiens" in t[2]:
		do = t[4]
	elif "sapiens" in t[5]:
		do = t[7]
	else:
		continue

	for g in do.split():
		(c,i) = genome.dicGenes[g]
		posBlocks[c].update([i-1,i])

print >> sys.stderr, len(posBlocks), sum([len(x) for x in posBlocks.itervalues()])

for (c,l) in posCNE.iteritems():
	for i in xrange(min(l), max(l)+1):
		if (i >= 0) and (i+1 < len(genome.lstGenes[c])):
			d = genome.lstGenes[c][i+1].beginning - genome.lstGenes[c][i].beginning
		else:
			d = 0
		print c, i, l.get(i, 0), int(i in posBlocks[c]), d



