#!/usr/bin/env python2

__doc__ = """
	Extrait toutes les diagonales entre chaque paire d'especes.
"""


import sys

import utils.myPhylTree
import utils.myGenomes
import utils.myFile
import utils.myTools
import utils.myMaths
import utils.myDiags


# Comparaison de deux especes
def compare(e1, e2, toStudy):
	statsDiags = []
	print >> sys.stderr, "Extraction des diagonales entre %s et %s [%s] ..." % (e1,e2,phylTree.dicParents[e1][e2]),
	for ((c1,d1),(c2,d2),daa) in utils.myDiags.calcDiags(dicGenomes[e1], dicGenomes[e2], genesAnc[phylTree.dicParents[e1][e2]], \
		arguments["fusionThreshold"], arguments["sameStrand"], orthosFilter):

		l = len(daa)
		if l < arguments["minimalLength"]:
			continue
		statsDiags.append(l)
		
		dic1 = dicGenomes[e1].lstGenes[c1]
		dic2 = dicGenomes[e2].lstGenes[c2]
		
		# Impression de la diagonale pour chaque ancetre
		for anc in toStudy:
			
			# Les numeros des genes ancestraux
			dica = genesAnc[anc].dicGenes
			if phylTree.isChildOf(e1, anc):
				da = [dica[dic1[i1].names[-1]][1] for (i1,_) in d1]
			else:
				da = [dica[dic2[i2].names[-1]][1] for (i2,_) in d2]
			
			# Verification
			if anc == phylTree.dicParents[e1][e2]:
				assert da == daa, (da,daa)

			sizes[anc].append(l)
			res = [anc,l, \
				e1,c1," ".join(dicGenomes[e1].lstGenes[c1][i1].names[0] for (i1,_) in d1), \
				e2,c2," ".join(dicGenomes[e2].lstGenes[c2][i2].names[0] for (i2,_) in d2), \
				utils.myFile.myTSV.printLine(da, " ")]
			if arguments["sameStrand"]:
				# Un champ en plus pour l'orientation
				ds1 = [s1 for (_,s1) in d1]
				ds2 = [s2 for (_,s2) in d2]
				assert len(set(x/y for (x,y) in zip(ds1,ds2))) == 1, (ds1, ds2)
				res.append(utils.myFile.myTSV.printLine(ds1, " "))
			files[anc].csvobject.writerow(res)

	print >> sys.stderr, utils.myMaths.myStats.txtSummary(statsDiags), "OK"


# Arguments
modesOrthos = list(utils.myDiags.OrthosFilterType._keys)
arguments = utils.myTools.checkArgs( \
	[("phylTree.conf",file), ("target",str)], \
	[("fusionThreshold",int,-1), ("sameStrand",bool,True), ("orthosFilter",str,modesOrthos), ("minimalLength",int,2), ("withReverseCmp",bool,False), \
	("OUT.projDiags",str,"proj/diags.%s.list.bz2"), \
	("genesFile",str,"~/work/data/genes/genes.%s.list.bz2"), \
	("ancGenomesFile",str,"~/work/ancestralGenomes/Genome.%s.bz2"), \
	("ancGenesFile",str,"~/work/data/ancGenes/ancGenes.%s.list.bz2")], \
	__doc__ \
)


# L'arbre phylogenetique
phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])

orthosFilter = utils.myDiags.OrthosFilterType[modesOrthos.index(arguments["orthosFilter"])]

# Les especes a utiliser
listSpecies = []
for x in arguments["target"].split(','):
	if x[0] != '.':
		listSpecies.extend(phylTree.species[x])
	else:
		listSpecies.append(x[1:])

# La liste des ancetres edites
dicLinks = [(e1,e2,set(phylTree.dicLinks[e1][e2][1:-1] + [phylTree.dicParents[e1][e2]])) for (e1,e2) in utils.myTools.myIterator.tupleOnStrictUpperList(listSpecies)]
tmp = set()
for (_,_,s) in dicLinks:
	tmp.update(s)

# Fichiers relatifs aux ancetres
genesAnc = {}
files = {}
sizes = {}
for anc in tmp:
	genesAnc[anc] = utils.myGenomes.Genome(arguments["ancGenesFile"] % phylTree.fileName[anc])
	files[anc] = utils.myFile.myTSV.writer(arguments["OUT.projDiags"] % phylTree.fileName[anc])
	sizes[anc] = []

# Fichiers relatifs aux especes comparees
dicGenomes = {}
for esp in listSpecies:
	if esp in phylTree.listSpecies:
		dicGenomes[esp] = utils.myGenomes.Genome(arguments["genesFile"] % phylTree.fileName[esp])
	else:
		dicGenomes[esp] = utils.myGenomes.Genome(arguments["ancGenomesFile"] % phylTree.fileName[esp], ancGenes=genesAnc[esp])

# On compare toutes les especes entre elles
for (e1,e2,toStudy) in dicLinks:
	compare(e1, e2, toStudy)
	if arguments["withReverseCmp"]:
		compare(e2, e1, toStudy)

# Fin
for anc in tmp:
	files[anc].file.close()
	print >> sys.stderr, "Statistiques de", anc, utils.myMaths.myStats.txtSummary(sizes[anc])

