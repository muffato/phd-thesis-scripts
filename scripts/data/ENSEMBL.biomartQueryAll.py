#! /users/ldog/muffato/python


__doc__ = """
	Execute la requete biomart pour chaque espece et enregistre le resultat dans un fichier
"""

import os
import sys
import urllib

import utils.myFile
import utils.myTools
import utils.myPhylTree

# Arguments
arguments = utils.myTools.checkArgs( \
	[("phylTree.conf",file), ("xmlRequest",file)], \
	[("biomartServer",str,"http://www.biomart.org/biomart/martservice"), ("outputFileName",str,"output.%s.txt")], \
	__doc__ \
)

# L'arbre phylogenetique
phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])

# La requete
f = utils.myFile.openFile(arguments["xmlRequest"], "r")
request = f.read()
f.close()

for esp in phylTree.listSpecies:

	# Les noms utilises dans les fichiers "Homo Sapiens" -> "hsapiens"
	tmp = esp.lower().split()
	tmp = tmp[0][0] + tmp[1]
	
	print >> sys.stderr, "Downloading %s (%s) ..." % (esp,tmp),
	urllib.urlretrieve(arguments["biomartServer"], filename=arguments["outputFileName"] % phylTree.fileName[esp], data=urllib.urlencode({"query": request % tmp}))
	print >> sys.stderr, "OK"


