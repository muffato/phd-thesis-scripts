#! /users/ldog/muffato/python

__doc__ = """
	Lit les arbres de proteines et rajoute une nouvelle espece a partir d'orthologues
	speciesList contient la liste des especes utilisees pour calculer les orthologues
"""

import sys

import utils.myFile
import utils.myTools
import utils.myGenomes
import utils.myPhylTree
import utils.myProteinTree

# Arguments
############
arguments = utils.myTools.checkArgs( \
	[("phylTree.conf",file), ("proteinTree",file), ("orthologues",file), ("newSpecies",str), ("speciesList",str)], \
	[("genesFile",str,""), ("defaultFamName", str, "FAM%d")], \
	__doc__ \
)

# Chargement des fichiers
##########################
phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])
phylTree.loadSpeciesFromList(arguments["speciesList"].split(",") + [arguments["newSpecies"]], arguments["genesFile"], storeGenomes = False)
orthologues = utils.myGenomes.Genome(arguments["orthologues"])

# Chargement de tous les arbres
################################
roots = set()
data = {}
info = {}
for (r,d,i) in utils.myProteinTree.loadTree(arguments["proteinTree"]):
	roots.add(r)
	data.update(d)
	info.update(i)


# Creation d'un nouveau noeud dans l'arbre
###########################################
def mkNode(taxon_name, x, isdup):
	global maxNodeID
	maxNodeID += 1
	dup = 3 if isdup else 0
	if type(x) == str:
		info[maxNodeID] = {'taxon_name':taxon_name, 'Duplication':dup, 'gene_name':x}
	else:
		info[maxNodeID] = {'taxon_name':taxon_name, 'Duplication':dup}
		data[maxNodeID] = x
	return maxNodeID


# Creation d'un arbre pour une liste de genes d'une meme espece
################################################################
def mkSpeciesNode(taxon_name, lst):
	# Selection des genes de l'espece
	esp = [g for g in lst if (g in phylTree.dicGenes) and (phylTree.dicGenes[g][0] == taxon_name)]
	if len(esp) == 0:
		print >> sys.stderr, "NONE", taxon_name, lst
		return None
	
	# Creer les feuilles de la nouvelle espece
	esp = [mkNode(taxon_name, g, False) for g in esp]

	if len(esp) == 1:
		# Le noeud de l'unique gene
		return esp[0]
	else:
		# On englobe par un noeud de duplication
		return mkNode(taxon_name, [(x,0) for x in esp], True)


# Construction du dictionnaire: nom_gene -> noeud/arbre
########################################################
def storeTreeLinks(node, root):
	if node in data:
		for (fils,_) in data[node]:
			assert fils is not None, (node,data[node],info[node])
			storeTreeLinks(fils, root)
	else:
		assert node is not None
		assert node in info, node
		assert 'gene_name' in info[node], (node,info[node])
		assert info[node]['gene_name'] is not None, (node,info[node])
		dicGeneName2Tree[info[node]['gene_name']] = (root,node)


# Initialisation
#################
print >> sys.stderr, "Preparation des arbres ...",
maxNodeID = max(info)
newSpecies = phylTree.officialName[arguments["newSpecies"]]
# Association nom de gene -> (id de la racine de l'arbre, id du noeud de ce gene)
dicGeneName2Tree = {}
for r in roots:
	storeTreeLinks(r, r)
print >> sys.stderr, "OK"

# On definit les endroits ou rajouter les familles
###################################################
print >> sys.stderr, "Insertion des familles ...",
nb0 = 0
nb1 = 0
nb2 = 0
nbcastors = 0
for famille in orthologues:

	unknownGenes = [g for g in famille.names if g not in phylTree.dicGenes]
	if len(unknownGenes) > 0:
		print >> sys.stderr, "Unknown gene %s " % ("/".join(unknownGenes))

	# On verifie qu'il y a bien un gene de la nouvelle espece
	newGenes = [g for g in famille.names if (g in phylTree.dicGenes) and (phylTree.dicGenes[g].species == arguments["newSpecies"])]
	if len(newGenes) == 0:
		print >> sys.stderr, "No new gene in this family %s "% ("/".join(famille.names))
		continue

	# Les genes deja presents dans l'arbre
	lstAnchorGenes = [dicGeneName2Tree[g][1] for g in famille.names if g in dicGeneName2Tree]
	lstAnchorTrees = set(dicGeneName2Tree[g][0] for g in famille.names if g in dicGeneName2Tree)
	ntrees = len(lstAnchorTrees)
	print >> sys.stderr, ntrees
	
	if ntrees == 0:

		# On construit un nouvel arbre en copiant la phylogenie des especes
		####################################################################
		def mkGeneTreeFromSpeciesTree(node):
			if node in phylTree.items:
				# On calque les genes sur l'arbre phylogenetique de reference
				lst = [(mkGeneTreeFromSpeciesTree(f),l) for (f,l) in phylTree.items.get(node)]
				lst = [(x,l+l0) for ((x,l),l0) in lst if x != None]
				if len(lst) == 0:
					return (None,0)
				if len(lst) == 1:
					return lst[0]
				return (mkNode(node, lst, False), 0)
			else:
				return (mkSpeciesNode(node, famille.names),0)

		newNodeID = mkGeneTreeFromSpeciesTree(phylTree.root)[0]
		roots.add(newNodeID)
		storeTreeLinks(newNodeID, newNodeID)
		info[newNodeID]['tree_name'] = arguments["defaultFamName"] % newNodeID

		nb0 += 1


	elif ntrees > 1:
		
		# La nouvelle espece sert d'outgroup a une duplication
		########################################################

		# Les racines des arbres avec les noms de taxon
		lst = set((x,info[x]['taxon_name']) for x in lstAnchorTrees)
		
		# Les racines doivent etre en dessous de la cible
		commonAnc = phylTree.lastCommonAncestor([root_name for (_,root_name) in lst])
		if phylTree.isChildOf(newSpecies, commonAnc):
			continue

		# Les novueaux noeuds
		# TODO Il pourrait etre bien de chercher si ce ne sont pas des speciations (ex: on rassemble un arbre de primates et un arbre de rongeurs)
		n = mkNode(commonAnc, [(tree_node,0) for tree_node in lstAnchorTrees], True)
		par = phylTree.dicParents[commonAnc][newSpecies]
		newNodeID = mkNode(par, [(n,0),(mkSpeciesNode(newSpecies, famille.names),0)], False)
		roots.difference_update(lstAnchorTrees)
		roots.add(newNodeID)
		storeTreeLinks(newNodeID, newNodeID)
		info[newNodeID]['tree_name'] = arguments["defaultFamName"] % newNodeID

		nb2 += 1


	else:
		
		# Il suffit de trouver le bon point d'insertion
		#################################################
		def countNbHits(node):
			if node in data:
				s = 0
				for (f,_) in data[node]:
					dicNode2ParentNode[f] = node
					s += countNbHits(f)
			else:
				s = 1 if node in lstAnchorGenes else 0
			nbHits[node] = s
			return s
	
		root = lstAnchorTrees.pop()

		nbHits = {}
		# Association id du fils -> id du pere, pour pouvoir remonter l'arbre
		dicNode2ParentNode = {}
		countNbHits(root)

		print >> sys.stderr, "NEW TRIP"
		node = root
		# Parcourt l'arbre jusqu'a arriver sur le point d'insertion
		while node in data:
			
			currentNodeName = info[node]['taxon_name']
			print >> sys.stderr, "ON", node, currentNodeName

			# Duplication
			if info[node]['Duplication'] >= 2:
				# On cherche la sous-branche avec le plus de hits
				lst = [(nbHits[f],f) for (f,_) in data[node]]
				best = max(lst)
				print >> sys.stderr, "DUPLICATION", sorted(lst)
				
				if phylTree.isChildOf(newSpecies, currentNodeName):
					# On est au toujours au dessus du point d'insertion, on peut continuer a descendre
					node = best[1]
					continue
				else:
					# On est passe en dessous
					if (best[0] == nbHits[node]) and (nbHits[node] != 0):
						# On peut avoir un point d'insertion plus precis
						node = best[1]
						continue
					else:
						# Par defaut, on s'arrete
						break
			
			elif phylTree.isChildOf(newSpecies, currentNodeName):
				# On est toujours au dessus
				try:
					nextInterm = phylTree.dicLinks.get(newSpecies).get(currentNodeName)[-2]
				except IndexError:
					# TODO est-ce que ca arrive ?
					assert newSpecies == currentNodeName
					print >> sys.stderr, newSpecies, currentNodeName, phylTree.dicLinks.get(newSpecies).get(currentNodeName)
				lst = [f for (f,_) in data[node] if phylTree.isChildOf(info[f]['taxon_name'], nextInterm)]
				print >> sys.stderr, "SPECIATION", lst, nextInterm, [x+(info[x[0]]['taxon_name'],) for x in data[node]], currentNodeName
				if len(lst) == 0:
					print >> sys.stderr, famille
					(node,_) = data[node][0]
					nbcastors += 1
					break
				assert len(lst) == 1
				node = lst[0]
				continue
			else:
				print >> sys.stderr, "END", node
				break

		nb1 += 1
		print >> sys.stderr, "INSERTION", node

		# Insertion des nouvelles donnees
		##################################
		espNode = mkSpeciesNode(newSpecies, famille.names)
		# Il faut creer le noeud juste au dessus de node
		newName = phylTree.dicParents[info[node]['taxon_name']][newSpecies]

		if node in dicNode2ParentNode:
			# Mettre a jour le pere pour qu'il pointe vers newNodeID au lieu de node
			p = dicNode2ParentNode[node]
			if (info[p]['taxon_name'] == newName) and (info[p]['Duplication'] < 2):
				data[p].append( (espNode, 0) )
			else:
				newNodeID = mkNode(newName, [(node,0),(espNode,0)], newName in [info[node]['taxon_name'],newSpecies])
				data[p] = [(newNodeID,x[1]) if x[0] == node else x for x in data[p]]

				# Et ses grands-parents si il y a des duplications
				while phylTree.isChildOf(info[p]['taxon_name'], newName):
					info[p]['taxon_name'] = newName
					info[p]['Duplication'] = 3

					if p not in dicNode2ParentNode:
						break
					p = dicNode2ParentNode[p]
			storeTreeLinks(root, root)
		else:
			assert node == root
			newNodeID = mkNode(newName, [(node,0),(espNode,0)], newName in [info[node]['taxon_name'],newSpecies])
			storeTreeLinks(newNodeID, newNodeID)
			info[newNodeID]['tree_name'] = arguments["defaultFamName"] % newNodeID
			roots.remove(root)
			roots.add(newNodeID)

print >> sys.stderr, "%d-%d-%d/%d nouveaux genes OK" % (nb0,nb1,nb2,len(orthologues.lstGenes[None]))
print >> sys.stderr, nbcastors, "castors"

print >> sys.stderr, "Verification ...",
# Les racines doivent etre les noeuds sans parent
allnodesI = set(info)
allnodes1 = set(data)
allnodes2 = []
for l in data.itervalues():
	allnodes2.extend(x[0] for x in l)
assert len(set(allnodes2)) == len(allnodes2)
assert allnodesI == allnodes1.union(allnodes2).union(roots), (allnodesI.difference(allnodes1).difference(allnodes2), allnodes1.union(allnodes2).difference(allnodesI))
assert len(roots.intersection(allnodes2)) == 0, roots.intersection(allnodes2)
assert roots.issuperset(allnodes1.difference(allnodes2)), allnodes1.difference(allnodes2).difference(roots)
# Les racines doivent avoir un tree_name
for r in roots:
	assert 'tree_name' in info[r]
	def check(node):
		if node in data:
			for (fils,_) in data[node]:
				check(fils)
		else:
			assert dicGeneName2Tree[info[node]['gene_name']] == (r,node)
	check(r)
print >> sys.stderr, "OK"

print >> sys.stderr, "Mise en forme des arbres ...",
for r in roots:
	utils.myProteinTree.printTree(sys.stdout, data, info, r)
print >> sys.stderr, len(roots), "arbres OK"


