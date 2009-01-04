#! /users/ldog/muffato/python

__doc__ = """
	Telecharge depuis le site d'Ensembl les fichiers des arbres de proteines
"""

import os
import sys
import collections

import utils.myFile
import utils.myTools
import utils.myPhylTree
import utils.myProteinTree

arguments = utils.myTools.checkArgs( \
	[("phylTree.conf",file)], \
	[("IN.EnsemblURL",str,"ftp://ftp.ensembl.org/pub/release-XXX/mysql/ensembl_compara_XXX"), \
	("IN.member",str,"member.txt.gz"), \
	("IN.genome_db",str,"genome_db.txt.gz"), \
	("IN.protein_tree_node",str,"protein_tree_node.txt.gz"), \
	("IN.protein_tree_member",str,"protein_tree_member.txt.gz"), \
	("IN.protein_tree_tag",str,"protein_tree_tag.txt.gz"), \
	], \
	__doc__ \
)


phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])


##############################################
# Chargement de la base de donnees d'Ensembl #
##############################################


# On charge les liens taxon_id -> species name
###############################################
print >> sys.stderr, "Chargement des liens taxon_id -> species_name ...",
taxonName = {}
f = utils.myFile.openFile(os.path.join(arguments["IN.EnsemblURL"], arguments["IN.genome_db"]), "r")
for ligne in utils.myFile.MySQLFileLoader(f):
	t = ligne.split("\t")
	taxonName[t[1]] = t[2]
f.close()
print >> sys.stderr, len(taxonName), "especes OK"


# On charge les liens member_id -> protein name
################################################
print >> sys.stderr, "Chargement des liens member_id -> protein_name ...",
tmpLinks = {}
f = utils.myFile.openFile(os.path.join(arguments["IN.EnsemblURL"], arguments["IN.member"]), "r")
for ligne in utils.myFile.MySQLFileLoader(f):
	t = ligne.split("\t")
	# A un numero member_id, on associe les noms (gene/transcrit/proteine) et l'espece
	if t[7] != "\\N":
		x = t[8].split()
		tmpLinks[t[0]] = ((x[1].split(':')[1], x[0].split(':')[1], t[1]), taxonName[t[4]])
f.close()
print >> sys.stderr, len(tmpLinks), "membres OK"


# On charge les liens node_id -> member_id
###########################################
print >> sys.stderr, "Chargement des liens node_id -> member_id ...",
x = 0
info = collections.defaultdict(dict)
f = utils.myFile.openFile(os.path.join(arguments["IN.EnsemblURL"], arguments["IN.protein_tree_member"]), "r")
for ligne in utils.myFile.MySQLFileLoader(f):
	t = ligne.split("\t")
	data = tmpLinks[t[1]]
	info[int(t[0])] = {'gene_name': data[0][0], 'transcript_name': data[0][1], 'protein_name': data[0][2], 'taxon_name': data[1]}
	x += 1
f.close()
print >> sys.stderr, x, "proteines OK"
del tmpLinks
del taxonName


# On charge les liens node_id -> infos 
#######################################
print >> sys.stderr, "Chargement des liens node_id -> infos ...",
f = utils.myFile.openFile(os.path.join(arguments["IN.EnsemblURL"], arguments["IN.protein_tree_tag"]), "r")
for ligne in utils.myFile.MySQLFileLoader(f):
	t = ligne.split("\t")

	# On extrait l'info sous forme numerique si possible
	try:
		t[2] = int(t[2])
	except ValueError:
		try:
			t[2] = float(t[2])
		except ValueError:
			pass
	
	info[ int(t[0]) ][t[1]] = t[2]

f.close()
print >> sys.stderr, len(info), "infos OK"

# On charge les liens node_id_father -> node_id_son
####################################################
print >> sys.stderr, "Chargement des arbres (node_id_father -> node_id_son) ...",
data = collections.defaultdict(list)
f = utils.myFile.openFile(os.path.join(arguments["IN.EnsemblURL"], arguments["IN.protein_tree_node"]), "r")
for ligne in utils.myFile.MySQLFileLoader(f):
	t = ligne.split("\t")
	# On rajoute le fils a son pere (avec une distance)
	node = int(t[1])
	data[ node ].append( (int(t[0]), float(t[5])) )
f.close()
print >> sys.stderr, len(data), "branches OK"


# On regle les parametres de info
##################################
for inf in info.itervalues():
	if 'taxon_name' in inf:
		# Pour passer des '/' et ' ' a '-' et '_', et enlever le 'Silurana'
		inf['taxon_name'] = phylTree.officialName[inf['taxon_name']]


# On a besoin des genomes modernes pour reconnaitre les genes
print >> sys.stderr, "Mise en forme des arbres ...",
for (root,_) in data[1]:
	utils.myProteinTree.printTree(sys.stdout, data, info, root)
print >> sys.stderr, len(data[1]), "arbres OK"

