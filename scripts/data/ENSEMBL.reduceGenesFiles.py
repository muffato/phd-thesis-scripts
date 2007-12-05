#! /users/ldog/muffato/python -OO

__doc__ = """
	Enleve de chaque genome les genes qui n'ont aucun homologue avec les autres genomes
"""

# Librairies
import os
import sys
import utils.myTools
import utils.myGenomes
import utils.myPhylTree


# Arguments
(noms_fichiers, options) = utils.myTools.checkArgs( ["phylTree.conf"], [("IN.ancGenesFile",str,""), ("IN.genesFile",str,""), ("OUT.genesFile",str,"")], __doc__)

# L'arbre phylogenetique
phylTree = utils.myPhylTree.PhylogeneticTree(noms_fichiers["phylTree.conf"])

utils.myTools.mkDir(options["OUT.genesFile"])


for esp in phylTree.listSpecies:
	genesAnc = utils.myGenomes.Genome(options["IN.ancGenesFile"] % phylTree.fileName[esp]).dicGenes
	genomeIni = utils.myGenomes.Genome(options["IN.genesFile"] % phylTree.fileName[esp])
	print >> sys.stderr, "Reducing %s ..." % esp,
	f = utils.myTools.myOpenFile(options["OUT.genesFile"] % phylTree.fileName[esp], "w")
	nb = 0
	for gene in genomeIni:
		if gene.names[0] in genesAnc:
			res = [gene.chromosome, gene.beginning, gene.end, gene.strand, gene.names[0]]
			print >> f, '\t'.join( [str(x) for x in res] )
			nb += 1
	f.close()
	print >> sys.stderr, "%d/%d OK" % (nb, len(genomeIni.dicGenes))

