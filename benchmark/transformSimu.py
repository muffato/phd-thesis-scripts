#! /users/ldog/muffato/python -OO

__doc__ = """
Transforme mes fichiers de genomes / genes ancestraux
"""

# Librairies
import sys
import math
import random
import utils.myMaths
import utils.myTools
import utils.myGenomes
import utils.myPhylTree


# Arguments
(noms_fichiers, options) = utils.myTools.checkArgs( \
	["phylTree.conf"], \
	[("root",str,""), ("species",str,""), \
	("genomeFile",str,"simu/genes/genes.%s.list.bz2"), \
	("ancGenesFile",str,"simu/ancGenes/ancGenes.%s.list.bz2")], \
	__doc__ \
)

# L'arbre phylogenetique
phylTree = utils.myPhylTree.PhylogeneticTree(noms_fichiers["phylTree.conf"])
if options["root"] not in phylTree.listAncestr:
	print >> sys.stderr, "Unknown root '%s'" % options["root"]
	sys.exit(1)

if len(options["species"]) == 0:
	species = phylTree.species[options["root"]]
else:
	species = options["species"].split(",")
phylTree.loadSpeciesFromList(species, options["genomeFile"], True)

# Selection des marqueurs presents dans toutes les especes
ancGenes = utils.myGenomes.Genome(options["ancGenesFile"] %  options["root"])
genesOK = {}
nbOK = 0
for gene in ancGenes:
	speciesOK = utils.myTools.defaultdict(int)
	for n in gene.names:
		if n in phylTree.dicGenes:
			speciesOK[phylTree.dicGenes[n][0]] += 1
	# Unique et universel ...
	if len([e for (e,n) in speciesOK.iteritems() if n == 1]) == len(species):
		nbOK += 1
		for n in gene.names:
			genesOK[n] = nbOK
print >> sys.stderr, nbOK, "genes OK"

for esp in species:
	print ">", esp, nbOK
	genome = phylTree.dicGenomes[esp]
	for (c,lst) in genome.lstGenes.iteritems():
		print "# chr" + str(c)
		seq = ""
		for g in lst:
			for n in g.names:
				if n in genesOK:
					seq += "%d " % genesOK[n]
					break
		print seq + "$"

