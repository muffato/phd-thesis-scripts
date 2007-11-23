#! /users/ldog/muffato/python -OO

__doc__ = """
	Telecharge depuis le site d'Ensembl les fichiers des listes de genes et d'annotations xref
"""


##################
# INITIALISATION #
##################

# Librairies
import os
import sys
import utils.myTools
import utils.myPhylTree


#############
# FONCTIONS #
#############

def getEnsemblFile(nom):
	return utils.myTools.myOpenFile( ("%s/%s" % (options["IN.EnsemblURL"],nom)).replace("XXX", str(options["releaseID"])) , "r")


########
# MAIN #
########

# Arguments
(noms_fichiers, options) = utils.myTools.checkArgs( \
	["phylTree.conf"], \
	[("releaseID",int,[42,43,44,45,46,47]), ("OUT.directory",str,""), \
	("IN.EnsemblURL",str,"ftp://ftp.ensembl.org/pub/release-XXX/mart_XXX/data/mysql/"), \
	("IN.genesFile",str,"ensembl_mart_XXX/%s_gene_ensembl__gene__main.txt.table.gz"), \
	("IN.xrefFile",str,"ensembl_mart_XXX/%s_gene_ensembl__xref_{refseq_dna,pdb,unigene,refseq_peptide,mirbase,rfam,uniprot_swissprot,embl,protein_id,uniprot_accession,uniprot_id,uniprot_sptrembl,hugo,xref_ipi}__dm.txt.table.gz"), \
	("OUT.genesFile",str,"genes/genes.%s.list.bz2"), \
	("OUT.fullGenesFile",str,"genes/full/genes.%s.list.bz2"), \
	("OUT.xrefFile",str,"xref/xref.%s.list.bz2"), \
	], \
	__doc__ \
)



# L'arbre phylogenetique
phylTree = utils.myPhylTree.PhylogeneticTree(noms_fichiers["phylTree.conf"])

OUTgenesFile = os.path.join(options["OUT.directory"], options["OUT.genesFile"])
OUTfullGenesFile = os.path.join(options["OUT.directory"], options["OUT.fullGenesFile"])
OUTxrefFile = os.path.join(options["OUT.directory"], options["OUT.xrefFile"])
for dir in [OUTgenesFile, OUTfullGenesFile, OUTxrefFile]:
	try:
		os.makedirs(os.path.dirname(dir))
	except OSError:
		pass


# Les noms utilises dans les fichiers "Homo Sapiens" -> "hsap"
nomReel = []
for esp in sorted(phylTree.listSpecies):
	tmp = esp.lower().split()
	tmp = tmp[0][0] + tmp[1]
	nomReel.append( (tmp,esp) )

# Les fichiers de genes
for (tmp,esp) in nomReel:
	
	print >> sys.stderr, "Telechargement de la liste des genes de %s ..." % esp,
	fo1 = utils.myTools.myOpenFile(OUTgenesFile % phylTree.fileName[esp], 'w')
	fo2 = utils.myTools.myOpenFile(OUTfullGenesFile % phylTree.fileName[esp], 'w')
	nb1 = 0
	nb2 = 0
	fi = getEnsemblFile(options["IN.genesFile"] % tmp)
	for ligne in MySQLFileLoader(fi):
		c = ligne.split('\t')
		s = "\t".join( [c[11],c[7],c[8],c[9],c[1]] )
		if ("RNA" not in c[3]) and ("pseudogene" not in c[3]):
			print >> fo1, s
			nb1 += 1
		print >> fo2, s
		nb2 += 1
	fi.close()
	fo1.close()
	fo2.close()
	print >> sys.stderr, "%d/%d genes OK" % (nb1,nb2)
	
	# Les references dans xref
	print >> sys.stderr, "Telechargement des annotations xref de %s ..." % esp,
	dic = utils.myTools.defaultdict(set)
	fi = getEnsemblFile(options["IN.genesFile"] % tmp)
	for ligne in MySQLFileLoader(fi):
		c = ligne.split("\t")
		dic[(c[1],c[3],c[5])].update( [x for x in c[6:] if (x != "\\N") and (len(x) > 0)] )
	fi.close()
	
	fo = utils.myTools.myOpenFile(OUTxrefFile % phylTree.fileName[esp], 'w')
	for ((gg,gt,gp),xref) in dic.iteritems():
		print >> fo, "\t".join([gg,gt,gp] + list(xref))
	fo.close()
	print >> sys.stderr, "%d annotations xref" % len(dic)

