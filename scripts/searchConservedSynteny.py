#! /users/ldog/muffato/python -OO

__doc__ = """
Extrait les paires de genes synteniques parmi toutes les especes
"""

# Librairies
import sys
import utils.myGenomes
import utils.myTools
import utils.myPhylTree


# Initialisation & Chargement des fichiers
(noms_fichiers, options) = utils.myTools.checkArgs( \
	["phylTree.conf"], \
	[("species",str,""), \
	("genesFile",str,"~/work/data/genes/genes.%s.list.bz2"), \
	("ancGenesFile",str,"~/work/data/ancGenes/ancGenes.%s.list.bz2")], \
	__doc__ \
)

phylTree = utils.myPhylTree.PhylogeneticTree(noms_fichiers["phylTree.conf"])
species = [phylTree.officialName[x] for x in options["species"].split(",")]
phylTree.loadSpeciesFromList(species, options["genesFile"])

espRef = species.pop()
genomeRef = phylTree.dicGenomes[espRef]
ancestr = set([phylTree.dicParents[espRef][e] for e in species])
ancFamilies = {}
for anc in ancestr:
	ancFamilies[anc] = utils.myGenomes.Genome(options["ancGenesFile"] % anc)
ancLst = [ancFamilies[phylTree.dicParents[espRef][e]] for e in species]

lastOK = False
for c in genomeRef.lstChr:
	genesChrRef = genomeRef.lstGenes[c]
	for i in xrange(len(genesChrRef)-1):

		# Pour ne pas reprendre le dernier gene (eviter les duplications locales)
		if lastOK:
			lastOK = False
			continue
		
		g1 = genesChrRef[i]
		g2 = genesChrRef[i+1]

		othersGenes = []
		for (i,e) in enumerate(species):
			anc = ancLst[i]
			genome = phylTree.dicGenomes[e]
			
			pos1 = genome.getPosition(anc.getOtherNames(g1.names[0]))
			pos2 = genome.getPosition(anc.getOtherNames(g2.names[0]))

			for ((c1,i1),(c2,i2)) in utils.myTools.myIterator.tupleOnTwoLists(pos1, pos2):
				if (c1 == c2) and (abs(i1-i2) == 1):
					othersGenes.append( (genome.lstGenes[c1][i1],genome.lstGenes[c2][i2]) )
					break
			else:
				break

		else:
			for (j,(gg1,gg2)) in enumerate(othersGenes):
				if gg1.beginning > gg2.beginning:
					x = gg1
					gg1 = gg2
					gg2 = x
				print species[j], gg1.chromosome, gg1.names[0],gg1.strand, gg2.names[0],gg2.strand
			print espRef, c, g1.names[0],g1.strand, g2.names[0],g2.strand
			print
			lastOK = True

	sys.stderr.write(".")

print >> sys.stderr, "OK"
