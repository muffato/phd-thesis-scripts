#! /users/ldog/muffato/python

__doc__ = """
	Cree la base Genomicus (table Species)
"""

import utils.myFile
import utils.myTools
import utils.myPhylTree


arguments = utils.myTools.checkArgs( [("phylTree.conf",file), ("genome_db.txt",file)], [], __doc__)

phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])

# Infos supplementaires d'assemblage et de genebuild
dicEsp = {}
for t in utils.myFile.myTSV.readTabular(arguments["genome_db.txt"], (int,int,str,str,int,str,str)):
	dicEsp[t[2]] = (t[1],t[3],t[5])


# Les groupes pour le menu
allgrp = []
todo = set(phylTree.indNames)


# Ajoute un nouveau groupe avec toutes les especes descedantes de node, qu'on n'a pas encore vues
def addModernSpecies(node, name):
	grp = todo.intersection(phylTree.species[node])
	todo.difference_update(grp)
	allgrp.append( (node,name,grp) )

addModernSpecies(phylTree.root, "Living yeasts")


# Ajoute un nouveau groupe avec tous les ancetres descendants de node, qui n'ont pas encore ete vus
def addAncestors(node, name):
	grp = todo.intersection(phylTree.allDescendants[node])
	todo.difference_update(grp)
	allgrp.append( (node,name,grp) )

addAncestors(phylTree.root, "Ancestral yeasts")


# Impression finale des groupes
maxid = len(phylTree.allNames)
for (i,(node,name,lst)) in enumerate(allgrp):

	# Entete
	res = [maxid+i, None, None, None, name, None,None, None, 100*i, 0]
	print utils.myFile.myTSV.MySQLFileWriter(res)

	# Contenu
	for (j,esp) in enumerate(lst):
		vers = "/".join(dicEsp[esp][1:]) if esp in dicEsp else ""
		othernames = [x for x in phylTree.commonNames[esp] if (x != esp) and (type(x) != int)]
		res = [phylTree.indNames[esp], esp, phylTree.ages[esp], vers, othernames[0] if len(othernames) > 0 else esp, None, None, None if esp in phylTree.items else esp.replace(' ', '_'), 100*i+j+1, int(esp in phylTree.lstEsp2X)]
		print utils.myFile.myTSV.MySQLFileWriter(res)


