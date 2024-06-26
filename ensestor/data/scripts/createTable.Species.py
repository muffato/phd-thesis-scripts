#!/usr/bin/env python2
# Copyright (c) 2010 IBENS/Dyogen : Matthieu MUFFATO, Alexandra LOUIS, Hugues ROEST CROLLIUS
# Mail : hrc@ens.fr or alouis@biologie.ens.fr
# Licences GLP v3 and CeCILL v2

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

addModernSpecies("Primates", "Primates")
addModernSpecies("Euarchontoglires", "Rodents etc.")
addModernSpecies("Laurasiatheria", "Laurasiatherias")
addModernSpecies("Xenarthra", "Xenarthras etc.")
addModernSpecies("Afrotheria", "Afrotherias etc.")
addModernSpecies("Mammalia", "Marsupials &amp; Monotremes")
addModernSpecies("Amniota", "Birds &amp; Reptiles")
addModernSpecies("Tetrapoda", "Amphibians")
addModernSpecies("Clupeocephala", "Fish")
# Avant ou apres l'insertion d'amphioxus
if "Deuterostomia" in phylTree.items:
	addModernSpecies("Deuterostomia", "Other chordates")
else:
	addModernSpecies("Chordata", "Other chordates")
addModernSpecies(phylTree.root, "Other eukaryotes")


# Ajoute un nouveau groupe avec tous les ancetres descendants de node, qui n'ont pas encore ete vus
def addAncestors(node, name):
	grp = todo.intersection(phylTree.allDescendants[node])
	todo.difference_update(grp)
	allgrp.append( (node,name,grp) )

addAncestors("Primates", "Ancestors in Primates")
addAncestors("Glires", "Ancestors in Rodents etc.")
addAncestors("Laurasiatheria", "Ancestors in Laurasiatherias")
addAncestors("Mammalia", "Ancestors in mammals")
addAncestors("Clupeocephala", "Ancestors in fish")
addAncestors("Euteleostomi", "Ancestors in vertebrates")
addAncestors(phylTree.root, "Other ancestors")


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
		if len(othernames) == 0:
			scientific_name = None
			common_name = esp
		else:
			scientific_name = esp
			common_name = othernames[0]
		res = [phylTree.indNames[esp], scientific_name, phylTree.ages[esp], vers, common_name, None, None, None if esp in phylTree.items else esp.replace(' ', '_'), 100*i+j+1, int(esp in phylTree.lstEsp2X)]
		print utils.myFile.myTSV.MySQLFileWriter(res)


