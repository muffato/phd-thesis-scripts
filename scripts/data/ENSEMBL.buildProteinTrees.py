#! /users/ldog/muffato/python

__doc__ = """
	Telecharge depuis le site d'Ensembl les fichiers des arbres de proteines
"""

import os
import sys
import collections

import utils.myTools
import utils.myPhylTree
import utils.myProteinTree

arguments = utils.myTools.checkArgs( \
	[("phylTree.conf",file)], \
	[("minDuplicationScore",float,0), \
	("IN.EnsemblURL",str,"ftp://ftp.ensembl.org/pub/release-XXX/mysql/ensembl_compara_XXX"), \
	("IN.member",str,"member.txt.gz"), \
	("IN.genome_db",str,"genome_db.txt.gz"), \
	("IN.protein_tree_node",str,"protein_tree_node.txt.gz"), \
	("IN.protein_tree_member",str,"protein_tree_member.txt.gz"), \
	("IN.protein_tree_tag",str,"protein_tree_tag.txt.gz"), \
	("OUT.tree",str,"tree.%s.bz2"), \
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
f = utils.myTools.myOpenFile(os.path.join(arguments["IN.EnsemblURL"], arguments["IN.genome_db"]), "r")
for ligne in utils.myFile.MySQLFileLoader(f):
	t = ligne.split("\t")
	taxonName[t[1]] = t[2]
f.close()
print >> sys.stderr, len(taxonName), "especes OK"


# On charge les liens member_id -> protein name
################################################
print >> sys.stderr, "Chargement des liens member_id -> protein_name ...",
tmpLinks = {}
f = utils.myTools.myOpenFile(os.path.join(arguments["IN.EnsemblURL"], arguments["IN.member"]), "r")
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
f = utils.myTools.myOpenFile(os.path.join(arguments["IN.EnsemblURL"], arguments["IN.protein_tree_member"]), "r")
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
f = utils.myTools.myOpenFile(os.path.join(arguments["IN.EnsemblURL"], arguments["IN.protein_tree_tag"]), "r")
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
	
	# Le dictionnaire
	if (t[1] not in ['lost_taxon_id', 'taxon_id', 'taxon_alias', 'original_cluster_id', 'Sitewise_dNdS_subroot_id']) and ('SIS' not in t[1]):
		info[ int(t[0]) ][t[1]] = t[2]

f.close()
print >> sys.stderr, len(info), "infos OK"

# On charge les liens node_id_father -> node_id_son
####################################################
print >> sys.stderr, "Chargement des arbres (node_id_father -> node_id_son) ...",
data = collections.defaultdict(list)
nextNodeID = max(info)
f = utils.myTools.myOpenFile(os.path.join(arguments["IN.EnsemblURL"], arguments["IN.protein_tree_node"]), "r")
for ligne in utils.myFile.MySQLFileLoader(f):
	t = ligne.split("\t")
	# On rajoute le fils de son pere (avec une distance)
	node = int(t[1])
	nextNodeID = max(nextNodeID, node)
	data[ node ].append( (int(t[0]), float(t[5])) )
f.close()
print >> sys.stderr, len(data), "branches OK"


# On regle les parametres de info
##################################
for (node,inf) in info.iteritems():
	# On considere que les duplications 'dubious' ne sont pas valables pour une duplication
	if 'dubious_duplication' in inf:
		del inf['dubious_duplication']
		inf['Duplication'] = 1
	elif inf['Duplication'] != 0:
		if (inf.get('duplication_confidence_score', -1) >= arguments["minDuplicationScore"]) or phylTree.isChildOf(inf['taxon_name'], "Clupeocephala"):
			inf['Duplication'] = 2
		else:
			inf['Duplication'] = 1
	# Pour passer des '/' et ' ' a '-' et '_'
	inf['taxon_name'] = phylTree.officialName[inf['taxon_name']]


###################################################################################################################
# Aplatit un noeud et des descendants directs si ils representent le meme taxon et qu'il n'y a pas de duplication #
###################################################################################################################
def flattenTree(node, rec):
	
	# Fin du process sur une feuille
	if node not in data:
		return

	# Appels recursifs
	if rec:
		for (g,_) in data[node]:
			flattenTree(g, rec)

	# Si c'est une vraie duplication, on n'a plus rien a faire
	if info[node]['Duplication'] >= 2:
		return

	newData = []
	taxonName = info[node]['taxon_name']
	for (g,d) in data[node]:
		inf = info[g]
		# 2x le meme taxon et pas de duplication
		if (g in data) and (inf['taxon_name'] == taxonName) and (inf['Duplication'] < 2):
			newData.extend([(g2,d+d2) for (g2,d2) in data[g]])
		else:
			newData.append( (g,d) )
	data[node] = newData


##########################################################################
# Simple etape de contraction de l'arbre: flatten & renommage successifs #
##########################################################################
def contractTree(node):

	if node not in data:
		return

	while True:
		
		# On s'assure que le nom de l'ancetre est le bon
		l = [info[g]['taxon_name'] for (g,_) in data[node]]
		info[node]['taxon_name'] = phylTree.lastCommonAncestor(l)
		
		# On s'assure que le noeud est bien aplati
		#  - si on vient de changer le nom du taxon, celui-ci peut correspondre avec le fils, d'ou la fusion
		#  - si on vient de creer le noeud, des dubious_duplication successifs ont pu faire aterrir ensemble des branches de meme taxon (ex Euarchontoglires)
		backup = data[node]
		flattenTree(node, False)
		if backup == data[node]:
			break


########################################################
# Redonne la topologie attendue a l'arbre              #
#   Rassemble les noeuds equivalents sous le meme fils #
########################################################
def rebuildTree(node):

	# Fin du process sur une feuille
	if node not in data:
		return

	contractTree(node)

	# On ne change les fils que si ce n'est pas une vraie duplication
	if info[node]['Duplication'] < 2:

		# On redefinit les fils pour -notre- arbre phylogenetique en triant les enfants en paquets
		fils = collections.defaultdict(list)
		anc = info[node]['taxon_name']
		lfils = phylTree.items.get(anc, [])
		for (g,d) in data[node]:
			gname = info[g]['taxon_name']
			for (a,_) in lfils:
				if phylTree.isChildOf(gname, a):
					fils[a].append( (g,d) )
					break
			else:
				# Apparait si g est le meme ancetre que node et que g est une duplication, ce qui l'a empeche d'etre aplati par flatten
				assert (gname == anc), "ERREUR: name!=anc [%s / %s / %s]" % (node, anc, gname)
				# Le noeud courant sera donc un noeud de duplication
				info[node]['Duplication'] = 3
				fils[anc].append( (g,d) )
	
		# len(fils):
		#  1 -> uniquement anc
		#  2 ou 3 parmi F1/F2/anc
		assert (len(fils) != 1) or (anc in fils), "ERREUR: 1=anc [%s / %s]" % (node, fils)
		assert (len(fils) <= (1+len(lfils))), "ERREUR: len>(1+nbFils) [%s / %s]" % (node, fils)

		if len(fils) > 1:
			if anc in fils:
				lst1 = fils.pop(anc)
				lst2 = []
				for tmp in fils.itervalues():
					lst2.extend(tmp)
				items = [(anc,lst1), (anc,lst2)]
			else:
				items = fils.items()

			newData = []
			for (anc,lst) in items:
				if len(lst) == 1:
					newData.append( lst[0] )
				elif len(lst) > 1:
					global nextNodeID
					nextNodeID += 1
					length = min([d for (_,d) in lst]) / 2
					data[nextNodeID] = [(g,d-length) for (g,d) in lst]
					info[nextNodeID] = {'taxon_name':anc, 'Duplication':0, 'Bootstrap':-1}
					newData.append( (nextNodeID,length) )
			data[node] = newData

	# Appels recursifs
	for (g,_) in data[node]:
		rebuildTree(g)

	contractTree(node)


#################################################
# Retrouve les vraies racines dans les familles #
#################################################
def getRoots(node, previousAnc, lastWrittenAnc):

	newAnc = info[node]['taxon_name']
	(toWrite,newLastWritten,isroot) = utils.myProteinTree.getIntermediateAnc(phylTree, previousAnc, lastWrittenAnc, newAnc, info[node]['Duplication'] >=2)

	if isroot:
		return [node]

	# Les genes des descendants
	subRoots = []
	for (g,_) in data.get(node,[]):
		subRoots.extend( getRoots(g, newAnc, newLastWritten) )
	return subRoots


# On a besoin des genomes modernes pour reconnaitre les genes
ft1 = utils.myTools.myOpenFile(arguments["OUT.tree"] % "1.ensembl", "w")
ft2 = utils.myTools.myOpenFile(arguments["OUT.tree"] % "2.flatten", "w")
ft3 = utils.myTools.myOpenFile(arguments["OUT.tree"] % "3.rebuilt", "w")
ft4 = utils.myTools.myOpenFile(arguments["OUT.tree"] % "4.cut", "w")

print >> sys.stderr, "Mise en forme des arbres ...",
nb = 0
nbA = 0
for (root,_) in data[1]:
	# Permet d'eviter quelques noeuds artefacts
	if 'taxon_name' in info[root]:
		utils.myProteinTree.printTree(ft1, data, info, root)
		flattenTree(root, True)
		utils.myProteinTree.printTree(ft2, data, info, root)
		rebuildTree(root)
		utils.myProteinTree.printTree(ft3, data, info, root)
		for x in getRoots(root, None, None):
			utils.myProteinTree.printTree(ft4, data, info, x)
			nbA += 1
		nb += 1
print >> sys.stderr, "%d (%d) arbres OK" % (nbA,nb)
ft1.close()
ft2.close()
ft3.close()
ft4.close()

