#! /users/ldog/muffato/python

__doc__ = """
	Corrige les arbres d'Ensembl en fonction du seuil minimal de duplication_score et de l'arbre des especes desire
"""

import sys
import itertools
import collections

import utils.myFile
import utils.myMaths
import utils.myTools
import utils.myPhylTree
import utils.myProteinTree

arguments = utils.myTools.checkArgs( \
	[("phylTree.conf",file), ("ensemblTree",file)], \
	[("minDuplicationScore",float,-1), ("minSyntenyScore",float,-1),
	("OUT.tree",str,"tree.%s.bz2"), ("genesFile",str,""),
	("defaultFamName",str,"GCUSGT%08d")], \
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
					info[nextNodeID] = {'taxon_name':anc, 'Duplication':0}
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


###########################################
# Sauvegarde toutes les familles de genes #
###########################################
def extractGeneFamilies(node, data, info, previousAnc, lastWrittenAnc):

	newAnc = info[node]['taxon_name']
	(toWrite,newLastWritten,isroot) = utils.myProteinTree.getIntermediateAnc(phylTree, previousAnc, lastWrittenAnc, newAnc, info[node]['Duplication'] >= 2)

	# Les genes des descendants
	if node in data:
		allGenes = []
		for (g,_) in data[node]:
			allGenes.extend( extractGeneFamilies(g, data, info, newAnc, newLastWritten) )
	else:
		allGenes = [ info[node]["gene_name"] ]

	for a in toWrite:
		i = len(geneFamilies[a])
		for x in allGenes:
			geneFamilies[a][x] = i

	return allGenes

@utils.myTools.memoize
def nextGene(g, anc):
	(e,c,i) = phylTree.dicGenes[g]
	c = phylTree.dicGenomes[e].lstGenes[c]
	l = len(c)
	res = []
	for step in [-1, 1]:
		j = i + step
		while (j >= 0) and (j < l):
			s = c[j].names[0]
			if s in geneFamilies[anc]:
				res.append(geneFamilies[anc][s])
				break
			j += step
		else:
			res.append(None)
	return (c[i].strand, res)

def isSynt(g1, g2):
	e1 = phylTree.dicGenes[g1].species
	e2 = phylTree.dicGenes[g2].species
	anc = phylTree.lastCommonAncestor([e1, e2])
	(s1,neigh1) = nextGene(g1, anc)
	(s2,neigh2) = nextGene(g2, anc)
	if s1 != s2:
		neigh2 = [neigh2[1], neigh2[0]]
	synt = 0
	all = 0
	for (j1,j2) in zip(neigh1, neigh2):
		if (j1 is None) or (j2 is None):
			print "notest", g1, g2
		else:
			all += 1
			synt += int(j1 == j2)
			if j1 == j2:
				print "synt", g1, g2
			else:
				print "nosynt", g1, g2
	return (synt,all)

baddup = []
def testDupSynt(node, data, info):
	if node in data:
		allgenes = []
		for (g,_) in data[node]:
			allgenes.append(testDupSynt(g, data, info))
		if (info[node]['Duplication'] >= 2) and (info[node]['taxon_name'] in phylTree.listAncestr):
			nbTot = 0
			nbSynt = 0
			for (l1,l2) in utils.myTools.myIterator.tupleOnStrictUpperList(allgenes):
				for x in itertools.product(l1, l2):
					y = isSynt(*x)
					nbTot += y[1]
					nbSynt += y[0]
			print node, info[node].get('duplication_confidence_score'), info[node].get('Bootstrap'), nbTot, nbSynt
			if nbSynt > arguments["minSyntenyScore"]*nbTot:
				baddup.append(node)
		return utils.myMaths.flatten(allgenes)
	else:
		return [info[node]['gene_name']]


phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])
phylTree.loadAllSpeciesBefore(None, arguments["genesFile"])

nextNodeID = 10000000
ft2 = utils.myFile.openFile(arguments["OUT.tree"] % "2.flatten", "w")
ft3 = utils.myFile.openFile(arguments["OUT.tree"] % "3.rebuilt", "w")
ft4 = utils.myFile.openFile(arguments["OUT.tree"] % "4.cut", "w")

geneFamilies = collections.defaultdict(dict)

nb = 0
nbA = 0
all = []
for (root,data,info) in utils.myProteinTree.loadTree(arguments["ensemblTree"]):
	extractGeneFamilies(root, data, info, None, None)
	all.append( (root,data,info) )

for (root,data,info) in [all[1]]:
	nextGene.cache.clear()
	testDupSynt(root, data, info)

	# On trie les bonnes duplications des mauvaises
	################################################
	for node in baddup.intersection(info):
		info[node]['Duplication'] = 1

	flattenTree(root, True)
	utils.myProteinTree.printTree(ft2, data, info, root)
	rebuildTree(root)
	utils.myProteinTree.printTree(ft3, data, info, root)
	for (i,x) in enumerate(getRoots(root, None, None)):
		if x != root:
			info[x]["tree_name"] = info[root]["tree_name"] + utils.myProteinTree.getDupSuffix(i+1, True)
		utils.myProteinTree.printTree(ft4, data, info, x)
		nbA += 1
	nb += 1

print >> sys.stderr, "Mise en forme des arbres : %d (%d) arbres OK" % (nbA,nb)

ft2.close()
ft3.close()
ft4.close()

