#!/usr/bin/env python2

__doc__ = """
	Ajoute aux blocs ancestraux les noms des chromosomes des especes modernes sur lesquels ils se retrouvent 
"""


import sys
import collections

import utils.myPhylTree
import utils.myGenomes
import utils.myFile
import utils.myTools
import utils.myMaths
import utils.myDiags


# Etablit la liste des chromosomes de chaque espece, avec les comptes
def getSpeciesCount(anc, diag):
	lstGenesAncAnc = genesAnc[anc].lstGenes[None]
	count = collections.defaultdict(int)
	countD = collections.defaultdict(int)
	for e in listSpecies:

		# L'ancetre commun
		a = phylTree.dicParents[anc][e]
		lstGenesAncA = genesAnc[a].lstGenes[None]
		genome = dicGenomes[e]

		for i in diag:
			# Les noms associes au gene ancestral
			names = lstGenesAncAnc[i].names
			if a != anc:
				names = set(genesAnc[a].getOtherNames(names[-1])) | set(names)
			
			lstByChrom = {}
			for t in utils.myGenomes.ContigType:
				lstByChrom[t] = [c for (c,_) in genome.getPosition(names) if c in genome.chrSet[t]]

			# Uniquement les chromosomes
			lstChrom = lstByChrom[utils.myGenomes.ContigType.Chromosome]
			# La liste des chromosomes + scaffolds
			lstChromScaff = lstChrom + lstByChrom[utils.myGenomes.ContigType.Scaffold]
			# On rajoute les chromosomes pointes par les randoms
			lstChromScaff.extend( utils.myGenomes.commonChrName(c.replace("_random","")) for c in lstByChrom[utils.myGenomes.ContigType.Random] if isinstance(c, str) )

			if len(lstChromScaff) == 0:
				# Aucun chrom
				pass
			elif len(lstChromScaff) == 1:
				# Un seul
				count[(e,lstChromScaff[0])] += 1
			elif len(lstChrom) == 1:
				# Plusieurs, mais parmi eux un seul vrai chromosome
				count[(e,lstChrom[0])] += 1
			else:
				for x in lstChromScaff:
					count[(e,x)] += 1
	s = 1 if len(diag) == 1 else 2
	return [(e,c,i) for ((e,c),i) in count.iteritems() if i >= s]



# Arguments
arguments = utils.myTools.checkArgs( \
	[("phylTree.conf",file), ("target",str), ("searchSpecies",str)], \
	[("IN.ancDiags",str,"proj/diags.%s.list.bz2"), \
	("OUT.ancDiags",str,"anc/diags.%s.list.bz2"), \
	("genesFile",str,"~/work/data/genes/genes.%s.list.bz2"), \
	("ancGenomesFile",str,"~/work/ancestralGenomes/Genome.%s.bz2"), \
	("ancGenesFile",str,"~/work/data/ancGenes/ancGenes.%s.list.bz2")], \
	__doc__ \
)


# L'arbre phylogenetique
phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])


# Les ancetres utilises pour trouver la localisation chromosomique
(listSpecies,listAnc) = utils.myDiags.getTargets(phylTree, arguments["searchSpecies"])
dicGenomes = {}
for e in listSpecies:
	if e in phylTree.listSpecies:
		dicGenomes[e] = utils.myGenomes.Genome(arguments["genesFile"] % phylTree.fileName[e])
	else:
		dicGenomes[e] = utils.myGenomes.Genome(arguments["ancGenomesFile"] % phylTree.fileName[e], withChr=True)
genesAnc = {}
for anc in listAnc:
	genesAnc[anc] = utils.myGenomes.Genome(arguments["ancGenesFile"] % phylTree.fileName[anc])


# Traitement final
for anc in utils.myDiags.getTargets(phylTree, arguments["target"])[1]:
	
	if anc not in genesAnc:
		genesAnc[anc] = utils.myGenomes.Genome(arguments["ancGenesFile"] % phylTree.fileName[anc])
	
	print >> sys.stderr, "Ajout des locations chromosomiques pour %s ..." % anc,
	n = 0
	fi = utils.myFile.myTSV.reader(arguments["IN.ancDiags"] % anc)
	fo = utils.myFile.myTSV.writer(arguments["OUT.ancDiags"] % anc)
	for t in fi.csvobject:
		d = [int(x) for x in t[2].split()]
		n += 1
		fo.csvobject.writerow(t + ["|".join("%s/%s/%d" % x for x in getSpeciesCount(anc, d))])
	fi.file.close()
	fo.file.close()
	print >> sys.stderr, "%d blocs OK" % n

