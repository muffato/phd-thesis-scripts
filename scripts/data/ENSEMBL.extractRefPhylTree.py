#!/usr/bin/env python2
# Copyright (c) 2010 IBENS/Dyogen : Matthieu MUFFATO, Alexandra LOUIS, Hugues ROEST CROLLIUS
# Mail : hrc@ens.fr or alouis@biologie.ens.fr
# Licences GLP v3 and CeCILL v2

__doc__ = """
	Extrait l'arbre phylogenetique des especes utilise par Ensembl
"""

import sys
import cStringIO

import utils.myFile
import utils.myTools
import utils.myPhylTree

arguments = utils.myTools.checkArgs( [("IN.protein_tree_tag",file)], [], __doc__)

dicTaxonName = {}
dicTaxonID = {}
dicTaxonAlias = {}

# Chargement des donnees
print >> sys.stderr, "Chargement des tags ...",
f = utils.myFile.openFile(arguments["IN.protein_tree_tag"], "r")
for ligne in f:
	t = ligne[:-1].split("\t")
	if t[1] == "taxon_name":
		dicTaxonName[t[0]] = t[2]
	elif t[1] == "taxon_id":
		dicTaxonID[t[0]] = t[2]
	elif (t[1] == "taxon_alias") and (t[0] not in dicTaxonAlias):
		dicTaxonAlias[t[0]] = t[2]
	elif t[1] == "taxon_alias_mya":
		dicTaxonAlias[t[0]] = t[2]
	elif t[1] == "species_tree_string":
		tree = t[2]
print >> sys.stderr, "OK (lengths:", len(dicTaxonName), len(dicTaxonID), len(dicTaxonAlias), ")"

# Recoupement des infos
resTaxon = {}
for x in dicTaxonName:
	i = dicTaxonID[x]
	s = dicTaxonName[x]
	a = dicTaxonAlias.get(x, s)
	if x in resTaxon:
		# qui sont censees dire la meme chose
		assert resTaxon[i] == (s,a)
	else:
		resTaxon[i] = (s,a)
print >> sys.stderr, len(resTaxon), "taxa named"

phylTree = utils.myPhylTree.PhylogeneticTree(cStringIO.StringIO(tree))

# Impression sous mon format, avec des indentations
def do(node, indent):
	node = node.replace("*", "")
	print ("\t" * indent) + "|".join(resTaxon[node])
	if node in phylTree.items:
		for (f,_) in phylTree.items[node]:
			do(f, indent+1)

do(phylTree.root, 0)

