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
	[("phylTree.conf",file), ("ensemblTree",file)], \
	[("minDuplicationScore",float,0), \
	("OUT.tree",str,"tree.%s.bz2"), \
	], \
	__doc__ \
)


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
	(_,newLastWritten,isroot) = utils.myProteinTree.getIntermediateAnc(phylTree, previousAnc, lastWrittenAnc, newAnc, info[node]['Duplication'] >=2)

	if isroot:
		return [node]

	# Les genes des descendants
	subRoots = []
	for (g,_) in data.get(node,[]):
		subRoots.extend( getRoots(g, newAnc, newLastWritten) )
	return subRoots



phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])

nextNodeID = 10000000
ft2 = utils.myFile.openFile(arguments["OUT.tree"] % "2.flatten", "w")
ft3 = utils.myFile.openFile(arguments["OUT.tree"] % "3.rebuilt", "w")
ft4 = utils.myFile.openFile(arguments["OUT.tree"] % "4.cut", "w")

nb = 0
nbA = 0
for (root,data,info) in utils.myProteinTree.loadTree(arguments["ensemblTree"]):

	# On regle les parametres de info
	##################################
	for inf in info.itervalues():
		for tagname in ['lost_taxon_id', 'taxon_id', 'taxon_alias', 'original_cluster_id', 'Sitewise_dNdS_subroot_id', 'SIS1', 'SIS2', 'SISi', 'SISu']:
			inf.pop(tagname, None)

		if 'Duplication' in inf:
			# On considere que les duplications 'dubious' ne sont pas valables pour une duplication
			assert (inf['Duplication'] in [0,1,2])
			if 'dubious_duplication' in inf:
				assert inf['Duplication'] == 1
				del inf['dubious_duplication']
			elif inf['Duplication'] != 0:
				if (inf.get('duplication_confidence_score', -1) < arguments["minDuplicationScore"]):
					inf['Duplication'] = 1
				else:
					inf['Duplication'] = 2

	flattenTree(root, True)
	utils.myProteinTree.printTree(ft2, data, info, root)
	rebuildTree(root)
	utils.myProteinTree.printTree(ft3, data, info, root)
	for x in getRoots(root, None, None):
		utils.myProteinTree.printTree(ft4, data, info, x)
		nbA += 1
	nb += 1

print >> sys.stderr, "Mise en forme des arbres : %d (%d) arbres OK" % (nbA,nb)

ft2.close()
ft3.close()
ft4.close()

