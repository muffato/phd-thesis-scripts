#! /users/ldog/muffato/python

import os
import sys
import urllib

import utils.myTools
import utils.myPhylTree


request = """<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query  virtualSchemaName = "default" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
				
	<Dataset name = "%s_gene_ensembl" interface = "default" >
		<Attribute name = "ensembl_gene_id" />
		<Attribute name = "coding" />
	</Dataset>
</Query>"""



# Arguments
arguments = utils.myTools.checkArgs( [("phylTree.conf",file)], [("OUT.file",str,"")], "Telecharge depuis le site d'Ensembl les CDS des genes")

# L'arbre phylogenetique
phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])

for esp in phylTree.listSpecies:
	
	# Les noms utilises dans les fichiers "Homo Sapiens" -> "hsapiens"
	tmp = esp.lower().split()
	tmp = tmp[0][0] + tmp[1]
	
	print >> sys.stderr, "Downloading %s (%s) ..." % (esp,tmp),
	dic = {}
	nb = 0
	for l in urllib.urlopen("http://www.biomart.org/biomart/martservice", data="query=" + (request % tmp)):
		try:
			(seq,gene) = l.replace('\n', '').split("\t")
			if seq == "Sequence unavailable":
				continue
			nb += 1
			if (gene in dic) and (len(seq) < len(dic[gene])):
				continue
			dic[gene] = seq
		except ValueError:
			print l,
	f = utils.myTools.myOpenFile( arguments["OUT.file"] % phylTree.fileName[esp], "w")
	for (gene,seq) in dic.iteritems():
		print >> f, "%s\t%s" % (gene,seq)
	f.close()
	print >> sys.stderr, "%d/%d OK" % (len(dic),nb)


