#! /users/ldog/muffato/python

__doc__ = """
	Telecharge depuis le site d'Ensembl (>=42) les fichiers des listes de genes
"""

import os
import sys

import utils.myFile
import utils.myTools
import utils.myPhylTree


arguments = utils.myTools.checkArgs( \
	[("phylTree.conf",file)], \
	[("IN.EnsemblURL",str,"ftp://ftp.ensembl.org/pub/release-43/mart_43/data/mysql/ensembl_mart_43/%s_gene_ensembl__gene__main.txt.table.gz"), \
	("OUT.genesFile",str,"genes.%s.list.bz2"), ("fields",str,"11,7,8,9,1,3")], \
	__doc__ \
)

# 19,1,5,6,17,12 pour Ensembl50

# L'arbre phylogenetique
phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])
fields = [int(x) for x in arguments["fields"].split(",")]

for esp in sorted(phylTree.listSpecies):

	# Les noms utilises dans les fichiers "Homo Sapiens" -> "hsap"
	tmp = esp.lower().split()
	tmp = tmp[0][0] + tmp[1]

	# Les fichiers de genes
	print >> sys.stderr, "Telechargement de la liste des genes de %s ..." % esp,
	fo = utils.myFile.openFile(arguments["OUT.genesFile"] % phylTree.fileName[esp], 'w')
	nb = 0
	fi = utils.myFile.openFile(arguments["IN.EnsemblURL"] % tmp,'r')
	for ligne in utils.myFile.MySQLFileLoader(fi):
		c = ligne.split('\t')
		print >> fo, "\t".join( [c[x] for x in fields] )
		nb += 1
	fi.close()
	fo.close()
	print >> sys.stderr, "%d genes OK" % nb

