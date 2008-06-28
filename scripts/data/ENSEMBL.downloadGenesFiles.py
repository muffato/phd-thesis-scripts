#! /users/ldog/muffato/python -OO

__doc__ = """
	Telecharge depuis le site d'Ensembl (>=42) les fichiers des listes de genes
"""

import os
import sys
import utils.myTools
import utils.myPhylTree


arguments = utils.myTools.checkArgs( \
	[("phylTree.conf",file)], \
	[("IN.EnsemblURL",str,"ftp://ftp.ensembl.org/pub/release-43/mart_43/data/mysql/ensembl_mart_43/%s_gene_ensembl__gene__main.txt.table.gz"), \
	("OUT.genesFile",str,"genes.%s.list.bz2")], \
	__doc__ \
)


# L'arbre phylogenetique
phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])

for esp in sorted(phylTree.listSpecies):

	# Les noms utilises dans les fichiers "Homo Sapiens" -> "hsap"
	tmp = esp.lower().split()
	tmp = tmp[0][0] + tmp[1]

	# Les fichiers de genes
	print >> sys.stderr, "Telechargement de la liste des genes de %s ..." % esp,
	fo = utils.myTools.myOpenFile(arguments["OUT.genesFile"] % phylTree.fileName[esp], 'w')
	nb = 0
	fi = utils.myTools.myOpenFile(arguments["IN.EnsemblURL"] % tmp,'r')
	for ligne in utils.myTools.MySQLFileLoader(fi):
		c = ligne.split('\t')
		print >> fo, "\t".join( [c[11],c[7],c[8],c[9],c[1],c[3]] )
		nb += 1
	fi.close()
	fo.close()
	print >> sys.stderr, "%d genes OK" % nb

