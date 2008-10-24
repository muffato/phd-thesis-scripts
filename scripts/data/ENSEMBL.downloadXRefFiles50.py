#! /users/ldog/muffato/python

__doc__ = """
	Telecharge depuis le site d'Ensembl (>=42) les fichiers d'annotations xref
"""

import os
import sys
import utils.myTools
import utils.myPhylTree


arguments = utils.myTools.checkArgs( \
	[("phylTree.conf",file)], \
	[("IN.EnsemblURL",str,"ftp://ftp.ensembl.org/pub/release-50/mysql/ensembl_mart_50/"), \
	("transcriptList",str,"%s_gene_ensembl__translation__main.txt.gz"), \
	("xrefFiles",str,"%s_gene_ensembl__ox_*__dm.txt.gz"), \
	("OUT.transcriptFile",str,"transcripts/transcripts.%s.list.bz2"), \
	("OUT.xrefFile",str,"transcripts/xref.%s.list.bz2")], \
	__doc__ \
)


# L'arbre phylogenetique
phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])

for esp in sorted(phylTree.listSpecies):

	# Les noms utilises dans les fichiers "Homo Sapiens" -> "hsap"
	tmp = esp.lower().split()
	tmp = tmp[0][0] + tmp[1]

	# Les references dans xref
	print >> sys.stderr, "Telechargement des transcrits de %s ..." % esp,
	fi = utils.myTools.myOpenFile(arguments["IN.EnsemblURL"] + (arguments["transcriptList"] % tmp), 'r')
	fo = utils.myTools.myOpenFile(arguments["OUT.transcriptFile"] % phylTree.fileName[esp], 'w')
	n = 0
	for ligne in utils.myFile.MySQLFileLoader(fi):
		c = ligne.split('\t')
		print c
		print len(c)
		print >> fo, utils.myFile.myTSV.printLine([c[x-1] for x in [137,66,72,143,84,22,9]])
		n += 1
	fi.close()
	fo.close()
	print >> sys.stderr, "%d transcrits xref" % n
	
	print >> sys.stderr, "Telechargement des annotations xref de %s ..." % esp,
	fi = utils.myTools.myOpenFile(arguments["IN.EnsemblURL"] + (arguments["xrefFiles"] % tmp), 'r')
	fo = utils.myTools.myOpenFile(arguments["OUT.xrefFile"] % phylTree.fileName[esp], 'w')
	n = 0
	for l in fi:
		l = l.replace("\N", "")
		# Filtre pour les lignes vides (\t et \n)
		if len(set(l)) > 2:
			print >> fo, l.replace("\N", ""),
		n += 1
	fi.close()
	fo.close()
	print >> sys.stderr, "%d annotations xref" % n

