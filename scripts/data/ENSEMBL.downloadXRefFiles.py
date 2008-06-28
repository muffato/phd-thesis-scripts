#! /users/ldog/muffato/python -OO

__doc__ = """
	Telecharge depuis le site d'Ensembl (>=42) les fichiers d'annotations xref
"""

import os
import sys
import utils.myTools
import utils.myPhylTree


arguments = utils.myTools.checkArgs( \
	[("phylTree.conf",file)], \
	[("IN.EnsemblURL",str,"ftp://ftp.ensembl.org/pub/release-43/mart_43/data/mysql/ensembl_mart_43/%s_gene_ensembl__xref_{refseq_dna,pdb,unigene,refseq_peptide,mirbase,rfam,uniprot_swissprot,embl,protein_id,uniprot_accession,uniprot_id,uniprot_sptrembl,hugo,xref_ipi}__dm.txt.table.gz"), \
	("OUT.xrefFile",str,"xref.%s.list.bz2")], \
	__doc__ \
)


# L'arbre phylogenetique
phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])

for esp in sorted(phylTree.listSpecies):

	# Les noms utilises dans les fichiers "Homo Sapiens" -> "hsap"
	tmp = esp.lower().split()
	tmp = tmp[0][0] + tmp[1]

	# Les references dans xref
	print >> sys.stderr, "Telechargement des annotations xref de %s ..." % esp,
	dic = utils.myTools.defaultdict(set)
	fi = utils.myTools.myOpenFile(arguments["IN.EnsemblURL"] % tmp,'r')
	for ligne in utils.myTools.MySQLFileLoader(fi):
		c = ligne.split('\t')
		dic[(c[1],c[3],c[5])].update( [x for x in c[6:] if (x != "\\N") and (len(x) > 0)] )
	fi.close()
	
	fo = utils.myTools.myOpenFile(arguments["OUT.xrefFile"] % phylTree.fileName[esp], 'w')
	for ((gg,gt,gp),xref) in dic.iteritems():
		l = [gg,gt,gp]
		l.extend(xref)
		print >> fo, utils.myTools.printLine(l)
	fo.close()
	print >> sys.stderr, "%d annotations xref" % len(dic)

