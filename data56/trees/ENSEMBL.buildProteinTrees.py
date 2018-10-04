#!/usr/bin/env python2

__doc__ = """
	Corrige les arbres d'Ensembl en fonction du seuil minimal de duplication_score et de l'arbre des especes desire
		1: score par defaut (0 -> 1)
		2: coef multiplicateur d'un duplication_confidence_score 6X_species / all_species
		3: duplication_confidence_score uniquement sur 6X_species et comparaison
"""

import sys
import collections

import utils.myFile
import utils.myTools
import utils.myPhylTree
import utils.myProteinTree

arguments = utils.myTools.checkArgs( \
	[("phylTree.conf",file), ("ensemblTree",file)], \
	[("minDuplicationScore",str,"-1"), ("OUT.tree",str,"tree.%s.bz2"), ("defaultFamName",str,"GCUSGT%08d"), ("scoreMethod",int,[1,2,3])], \
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



phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])

# Limites automatiques de score de duplication
if arguments["scoreMethod"] in [1, 3]:
	def calc(anc, val):
		return val
elif arguments["scoreMethod"] == 2:
	def calc(anc, val):
		nesp = len(phylTree.species[anc])
		n2X = len(phylTree.lstEsp2X.intersection(phylTree.species[anc]))
		# La moitie des especes non 2X a vu la duplication (au minimum 1 espece)
		return round(max(1., val*(nesp-n2X)) / nesp, 3) - 2e-3
else:
	assert False, "unknown scoreMethod"

minDuplicationScore = {}
try:
	# Une limite pour tout le monde
	val = float(arguments["minDuplicationScore"])
	for anc in phylTree.listAncestr:
		minDuplicationScore[anc] = calc(anc, val)
except ValueError:
	f = utils.myFile.openFile(arguments["cutoffList"], "r")
	for l in f:
		t = l.split()
		anc = phylTree.officialName[t[0]]
		minDuplicationScore[anc] = calc(anc, float(t[1]))
	f.close()
print >> sys.stderr, minDuplicationScore

# Les scores dans l'abre pour les especes modernes valent toujours 1, on doit toujours les accepter
for esp in phylTree.listSpecies:
	minDuplicationScore[esp] = 0

nextNodeID = 10000000
ft2 = utils.myFile.openFile(arguments["OUT.tree"] % "2.flatten", "w")
ft3 = utils.myFile.openFile(arguments["OUT.tree"] % "3.rebuilt", "w")
ft4 = utils.myFile.openFile(arguments["OUT.tree"] % "4.cut", "w")

@utils.myTools.memoize
def goodSpecies(anc):
	return phylTree.species[anc].difference(phylTree.lstEsp2X)


nb = 0
nbA = 0
for (root,data,info) in utils.myProteinTree.loadTree(arguments["ensemblTree"]):

	# Les arbres qui ne sont pas des arbres -> Plus necessaire
	#if 'taxon_name' not in info[root]:
	#	continue
	
	@utils.myTools.memoize
	def getSpeciesSets(node):
		if node in data:
			return set().union(*(getSpeciesSets(x) for (x,_) in data[node]))
		else:
			assert info[node]["taxon_name"] in phylTree.listSpecies
			return set([info[node]["taxon_name"]])

	if "tree_name" not in info[root]:
		info[root]["tree_name"] = arguments["defaultFamName"] % nb

	# On trie les bonnes duplications des mauvaises
	################################################
	for (node,inf) in info.iteritems():

		if 'Duplication' in inf:
			if 'dubious_duplication' in inf:
				# On considere que les duplications 'dubious' ne sont pas valables pour une duplication
				assert inf['Duplication'] == 1
				print "dubious\t%s\t0" % inf["taxon_name"]

			elif inf['Duplication'] != 0:
				#print node, inf
				speciessets = [getSpeciesSets(x) for (x,_) in data[node]]
				inters = speciessets[0].intersection(speciessets[1])
				#union = speciessets[0].union(*speciessets[1:])
				#print speciessets
				#print inters
				#print inf["taxon_name"], goodSpecies(inf["taxon_name"])
				#print len(inters.intersection(goodSpecies(inf["taxon_name"])))
				if	((arguments["scoreMethod"] < 3) and
						(inf['duplication_confidence_score'] < minDuplicationScore[inf["taxon_name"]])) or \
					((arguments["scoreMethod"] == 3) and
						(len(inters.intersection(goodSpecies(inf["taxon_name"]))) < (arguments["minHighCoverageGenomesProportion"]*len(goodSpecies(inf["taxon_name"]))))):

				# Attention: pour les arbres d'Ensembl dont la racine est une duplication, celle-ci n'ayant pas d'outgroup vaut 1
				#   Il faut la passer a 2 si le score est suffisant
					inf['Duplication'] = 1
					print "toolow\t%s\t%f" % (inf["taxon_name"], inf['duplication_confidence_score'])
					#if inf['duplication_confidence_score'] >= 0.75:
					#	print inf["taxon_name"], goodSpecies(inf["taxon_name"])
					#	print speciessets
					#	print inters
					#	print len(inters.intersection(goodSpecies(inf["taxon_name"])))
					#	print
				else:
					inf['Duplication'] = 2
					print "good\t%s\t%f" % (inf["taxon_name"], inf['duplication_confidence_score'])
				#print

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

