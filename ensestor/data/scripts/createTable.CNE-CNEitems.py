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
import utils.myPhylTree

arguments = utils.myTools.checkArgs( \
	[("phylTree.conf",file), ("KX_HDMC_blocs.list",file), ("KX_HDMC_especes.list",file)], \
	[("ancGenesFile",str,""), ("genesFile",str,""), ("syntenyCutoff",int,5), ("output",str,"")], \
	__doc__ \
)

phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])

dicScientificName = {
	'cow': 'Bos taurus', \
	'mouse': 'Mus musculus', \
	'human': 'Homo sapiens', \
	'dog': 'Canis familiaris', \
	'chicken': 'Gallus gallus', \
	'stickleback': 'Gasterosteus aculeatus' \
}


# Chargement des genomes des especes modernes si pas encore fait
dicGenomes = {}
posCNE1 = {}
posCNE2 = {}
def loadSpecies(species_name):
	if species_name in dicGenomes:
		return
	
	# Chargement du genome
	lst = collections.defaultdict(list)
	for gene in utils.myGenomes.Genome(arguments["genesFile"] % phylTree.fileName[species_name]):
		lst[str(gene.chromosome)].append( (gene.beginning,gene.end,gene.names[0]) )
	
	dicGenomes[species_name] = lst
	posCNE1[species_name] = collections.defaultdict(list)
	posCNE2[species_name] = collections.defaultdict(list)



# Chargement des CNE
CNE = {}
for line in utils.myFile.myTSV.readTabular(arguments["KX_HDMC_blocs.list"], [str,int,int,int,str]):
	# Les infos generales du CNE
	CNE[line[0]] = (line[1:],[])

levels = collections.defaultdict(set)
for line in utils.myFile.myTSV.readTabular(arguments["KX_HDMC_especes.list"], [str,str,int,int,str,str,int,str]):
	# Les infos de chaque CNE pour chaque espece
	species_name = dicScientificName[line[0]]
	loadSpecies(species_name)

	# Pas de vache !!
	if species_name == 52:
		continue

	# Il faut que le chromosome de l'alignement multiple soit dans le genome de reference
	if line[1] in dicGenomes[species_name]:
		# On recherche tout de suite l'index de position sur le chromosome
		i = bisect.bisect_left(dicGenomes[species_name][line[1]], (line[2],line[3]))
	else:
		i = 0
		#print >> sys.stderr, "no chromosome *%s* for *%s*" % (line[1],species_name)
	CNE[line[5]][1].append( (species_name,line[1],i,line) )
	levels[line[6]-1].add(species_name)

# Chargement des genes ancestraux
print >> sys.stderr, levels
ancNames = [(int(x),phylTree.lastCommonAncestor(list(levels[x]))) for x in levels]
ancNames.sort()
print >> sys.stderr, ancNames
ancGenes = [utils.myGenomes.Genome(arguments["ancGenesFile"] % x).dicGenes for (_,x) in ancNames]


# Filtre des CNE sur la syntenie
goodCNE = []
for (cne_id,(_,cne_items)) in CNE.iteritems():
	found = []
	ancGenesID = []
	# Parcours de chaque espece
	for (species_name,chrom,pos,line) in cne_items:
		lst = set()
		refGenes = ancGenes[int(line[6])-1]
		# On recupere les genes voisins
		for g in dicGenomes[species_name][chrom][max(0,pos-arguments["syntenyCutoff"]):pos+arguments["syntenyCutoff"]+1]:
			if g[-1] in refGenes:
				lst.add(refGenes[g[-1]])
		ancGenesID.append( (species_name,lst) )
	goodSpecies = set()
	# Comparaison des especes
	for ((e1,l1),(e2,l2)) in utils.myTools.myIterator.tupleOnStrictUpperList(ancGenesID):
		# Si deux especes ont un gene syntenique, on les garde
		if len(l1.intersection(l2)) > 0:
			goodSpecies.add(e1)
			goodSpecies.add(e2)
	goodCNE.append( (cne_id,goodSpecies) )


dataCNE = {}
f = utils.myFile.openFile(arguments["output"] % "CNE", "w")

# Parcours des CNE restants
for (realID,(cne_id,cne_spec)) in enumerate(goodCNE):

	# Parcours des especes
	for (species_name,chrom,pos,line) in CNE[cne_id][1]:
	
		level = line[6]
		syntenic = int(species_name in cne_spec)
		strand = int(line[4] + '1')
		intronic = -1

		# Si on a une info sur le chromosome
		if line[1] in dicGenomes[species_name]:
			i = bisect.bisect_left(dicGenomes[species_name][line[1]], (line[2],line[3]))
			if (i >= 1) and (dicGenomes[species_name][line[1]][i-1][1] > line[2]):
				intronic = i-1
		else:
			i = 0
		
		# On enregistre la position
		dataCNE[(realID,species_name)] = [realID,phylTree.indNames[species_name], line[1],None,None,  line[2],line[3],strand,intronic,syntenic,level,line[7]]
		posCNE1[species_name][line[1]].append( (i,line[2],realID) )
		if intronic >= 0:
			posCNE2[species_name][line[1]].append( (i-1,line[2],realID) )
		else:
			posCNE2[species_name][line[1]].append( (i,line[2],realID) )
	
	# On imprime les infos generales
	print >> f, utils.myFile.myTSV.printLine((realID,) + CNE[cne_id][0] + (level,))
f.close()


def writePositions(pos, index):

	# Parcours des CNE par espece
	for (species_name,x) in pos.iteritems():
		lasti = None
		# Puis par chromosome
		for (chrom,lst) in x.iteritems():
			lst.sort()
			# Puis par espace intergenique
			for (pos0,y) in itertools.groupby(lst, operator.itemgetter(0)):
				l = list(y)
				# On positionne les CNE dans le bon ordre et a intervalle regulier
				for (i,y) in enumerate(l):
					position = pos0 + (i+1.)/(len(l)+1.)
					dataCNE[(y[2],species_name)][index] = position

writePositions(posCNE1, 3)
writePositions(posCNE2, 4)

f = utils.myFile.openFile(arguments["output"] % "CNE_items", "w")
for x in dataCNE.itervalues():
	print >> f, utils.myFile.myTSV.printLine(x)
f.close()

