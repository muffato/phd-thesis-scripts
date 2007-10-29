#! /users/ldog/muffato/python -OO

__doc__ = """
	Telecharge depuis le site d'Ensembl les CDS des genes
"""


##################
# INITIALISATION #
##################

# Librairies
import os
import sys
import utils.myTools
import utils.myPhylTree


request = """<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query  virtualSchemaName = "default" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
				
	<Dataset name = "%s_gene_ensembl" interface = "default" >
		<Attribute name = "gene_stable_id" />
		<Attribute name = "coding" />
	</Dataset>
</Query>"""


########
# MAIN #
########

# Arguments
(noms_fichiers, options) = utils.myTools.checkArgs( \
	["phylTree.conf"], \
	[("OUT.file",str,"")], \
	__doc__ \
)


# L'arbre phylogenetique
phylTree = utils.myPhylTree.PhylogeneticTree(noms_fichiers["phylTree.conf"])

# Le repertoire
try:
	os.makedirs(os.path.dirname(options["OUT.file"]))
except OSError:
	pass


# Les noms utilises dans les fichiers "Homo Sapiens" -> "hsapiens"
nomReel = []
for esp in phylTree.listSpecies:
	tmp = esp.lower().split()
	tmp = tmp[0][0] + tmp[1]
	
	print >> sys.stderr, "Downloading %s (%s) ..." % (esp,tmp),
	nb = 0
	(stdin,stdout,stderr) = os.popen3( 'wget -O - http://www.biomart.org/biomart/martservice --post-data="query=`cat /dev/stdin`" ')
	stderr.close()
	print >> stdin, request % tmp
	stdin.close()
	dic = {}
	for l in stdout:
		try:
			(seq,gene) = l[:-1].split("\t")
			if seq == "Sequence unavailable":
				continue
			nb += 1
			if (gene in dic) and (len(seq) < len(dic[gene])):
				continue
			dic[gene] = seq
		except ValueError:
			print l,
	f = utils.myTools.myOpenFile( options["OUT.file"] % phylTree.fileName[esp], "w")
	for (gene,seq) in dic.iteritems():
		print >> f, "%s\t%s" % (gene,seq)
	f.close()
	print >> sys.stderr, "%d/%d OK" % (len(dic),nb)


