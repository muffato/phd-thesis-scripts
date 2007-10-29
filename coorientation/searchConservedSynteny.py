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

# Le genome de la 1ere espece est arbitrairement celui de reference
espRef = species.pop()
genomeRef = phylTree.dicGenomes[espRef]
# On charge tous les genes ancestraux necessaires
ancestr = set([phylTree.dicParents[espRef][e] for e in species])
ancFamilies = {}
for anc in ancestr:
	ancFamilies[anc] = utils.myGenomes.Genome(options["ancGenesFile"] % phylTree.fileName[anc])
ancLst = [ancFamilies[phylTree.dicParents[espRef][e]] for e in species]

# Le parcours du genome de reference
for c in genomeRef.lstChr:
	genesChrRef = genomeRef.lstGenes[c]
	g = 0
	# i et i+1 doivent rester sur le chromosome
	while g < len(genesChrRef)-1:

		# Les deux genes
		g1 = genesChrRef[g]
		g2 = genesChrRef[g+1]

		# On cherche a remplir othersGenes, la liste des paires de genes en syntenie des autres especes
		othersGenes = []
		# Parcours des autres especes
		for (i,e) in enumerate(species):
			
			anc = ancLst[i]
			genome = phylTree.dicGenomes[e]
			
			# Les genes orthologues dans l'autre espece
			pos1 = genome.getPosition(anc.getOtherNames(g1.names[0]))
			pos2 = genome.getPosition(anc.getOtherNames(g2.names[0]))

			# On cherche la meilleure combinaison (si il y a des paralogues)
			for ((c1,i1),(c2,i2)) in utils.myTools.myIterator.tupleOnTwoLists(pos1, pos2):
				if (c1 == c2) and (abs(i1-i2) == 1):
					# OK, on passe au suivant
					othersGenes.append( (genome.lstGenes[c1][i1],genome.lstGenes[c2][i2]) )
					break
			else:
				# NON, on sort de la boucle
				break

		else:
			# Ici, toutes les especes ont une paire de genes en syntenie
			for (j,(gg1,gg2)) in enumerate(othersGenes):
				# On echange les deux genes pour les mettre dans le bon sens
				if gg1.beginning > gg2.beginning:
					x = gg1
					gg1 = gg2
					gg2 = x
				# Impression
				print species[j], gg1.chromosome, gg1.names[0],gg1.strand, gg2.names[0],gg2.strand, gg1.beginning, gg1.end, gg2.beginning, gg2.end
			print espRef, c, g1.names[0],g1.strand, g2.names[0],g2.strand, g1.beginning, g1.end, g2.beginning, g2.end
			print
			# On saute le gene d'apres pour ne pas le compter 2 fois
			g += 1
		
		# On passe au suivant
		g += 1

	sys.stderr.write(".")

print >> sys.stderr, "OK"
