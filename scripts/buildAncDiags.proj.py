#! /users/ldog/muffato/python

__doc__ = """
	Extrait toutes les diagonales entre chaque paire d'especes.
	Les diagonales apportent les genes qui etaient sur un meme chromosome depuis leur ancetre commun dans les deux lignees.
"""

import sys

import utils.myPhylTree
import utils.myGenomes
import utils.myTools
import utils.myMaths
import utils.myDiags2

def do(e1, e2, toStudy):
	for ((c1,d1),(c2,d2),s) in utils.myDiags2.calcDiags(dicGenomes[e1], dicGenomes[e2], genesAnc[phylTree.dicParents[e1][e2]], arguments["minimalLength"], \
		arguments["fusionThreshold"], arguments["sameStrand"] and (e1 not in genesAnc) and (e2 not in genesAnc), arguments["keepOnlyOrthos"]):

		statsDiags.append(len(d1))
		
		dic1 = dicGenomes[e1].lstGenes[c1]
		dic2 = dicGenomes[e2].lstGenes[c2]
		d1 = tuple(dic1[i1].names[0] for i1 in d1)
		d2 = tuple(dic2[i2].names[0] for i2 in d2)
		s = utils.myTools.printLine(s, " ")
		sd1 = " ".join(d1)
		sd2 = " ".join(d2)

		for anc in toStudy:
			dic = genesAnc[anc].dicGenes
			if d1[0] in dic:
				d = d1
			else:
				d = d2
			stats[anc].append(len(d))
			fout[anc].csvobject.writerow( [anc,len(d), e1,c1,sd1, e2,c2,sd2, " ".join(str(dic[name][1]) for name in d), s] )
	print >> sys.stderr, utils.myMaths.myStats.txtSummary(statsDiags),


# Arguments
arguments = utils.myTools.checkArgs( \
	[("phylTree.conf",file), ("target",str)], \
	[("fusionThreshold",int,-1), ("minimalLength",int,2), ("sameStrand",bool,True), ("keepOnlyOrthos",bool,False), \
	("OUT.projDiags",str,"proj/diags.%s.list.bz2"), \
	("genesFile",str,"~/work/data/genes/genes.%s.list.bz2"), \
	("ancGenomesFile",str,"~/work/ancestralGenomes/Genome.%s.bz2"), \
	("ancGenesFile",str,"~/work/data/ancGenes/ancGenes.%s.list.bz2")], \
	__doc__ \
)


# L'arbre phylogenetique
phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])

# Les especes a utiliser
dicGenomes = {}
tmp = arguments["target"].split(',')
if len(arguments["target"]) == 0:
	print >> sys.stderr, "Aucune cible indiquee pour l'extraction des diagonales"
	sys.exit(1)
listSpecies = []
for x in tmp:
	if x[0] != '.':
		listSpecies.extend(phylTree.species[x])
		for e in phylTree.species[x]:
			dicGenomes[e] = utils.myGenomes.Genome(arguments["genesFile"] % phylTree.fileName[e])
	else:
		listSpecies.append(x[1:])
		dicGenomes[x[1:]] = utils.myGenomes.Genome(arguments["ancGenomesFile"] % phylTree.fileName[x[1:]], withChr=True)

# La liste des ancetres edites
dicLinks = [(e1,e2,set(phylTree.dicLinks[e1][e2][1:-1] + [phylTree.dicParents[e1][e2]])) for (e1,e2) in utils.myTools.myIterator.tupleOnStrictUpperList(listSpecies)]
tmp = set()
for (_,_,s) in dicLinks:
	tmp.update(s)
genesAnc = {}
fout = {}
stats = {}
for anc in tmp:
	genesAnc[anc] = utils.myGenomes.Genome(arguments["ancGenesFile"] % phylTree.fileName[anc])
	fout[anc] = utils.myTools.tsvWriter(arguments["OUT.projDiags"] % phylTree.fileName[anc])
	stats[anc] = []

# On compare toutes les especes entre elles
for (e1,e2,toStudy) in dicLinks:
	statsDiags = []
	print >> sys.stderr, "Extraction des diagonales entre %s et %s ..." % (e1,e2),
	do(e1, e2, toStudy)
	print >> sys.stderr, "OK"

# Traitement final
for (anc,lst) in stats.iteritems():
	print >> sys.stderr, "Statistiques de %s :" % anc, utils.myMaths.myStats.txtSummary(lst)
	# Fermeture du fichier de sortie
	fout[anc].file.close()

