#! /users/ldog/muffato/python -OO

__doc__ = """
	Lit les arbres de proteines et rajoute une nouvelle espece a partir d'orthologues
"""

import sys
import utils.myTools
import utils.myGenomes
import utils.myPhylTree
import utils.myProteinTree

# Arguments
############
(noms_fichiers, options) = utils.myTools.checkArgs( \
	["phylTree.conf", "proteinTree", "orthologues"], \
	[("newSpecies",str,""), ("genesFile",str,"")], \
	__doc__ \
)

# Chargement des fichiers
##########################
phylTree = utils.myPhylTree.PhylogeneticTree(noms_fichiers["phylTree.conf"])
phylTree.loadAllSpeciesSince(None, options["genesFile"], storeGenomes = False)
orthologues = utils.myGenomes.Genome(noms_fichiers["orthologues"])
(data,info,roots) = utils.myProteinTree.loadTree(noms_fichiers["proteinTree"])


# Indique s'il s'agit d'un noeud de duplication (les duplications 'dubious' ne sont pas valables
#################################################################################################
def isDuplicatedNode(inf):
	return (inf['Duplication'] != 0) and ('dubious_duplication' not in inf)


# Creation d'un nouveau noeud dans l'arbre
###########################################
def mkNode(taxon_name, x, dup):
	global maxNodeID
	maxNodeID += 1
	dup = int(dup)
	if type(x) == str:
		info[maxNodeID] = {'taxon_name':taxon_name, 'Duplication':dup, 'gene_name':x}
	else:
		info[maxNodeID] = {'taxon_name':taxon_name, 'Duplication':dup, 'Bootstrap':-1}
		data[maxNodeID] = x
	return maxNodeID

# Creation d'un arbre pour une liste de genes d'une meme espece
################################################################
def mkSpeciesNode(taxon_name, lst):
	# Creer les feuilles de la nouvelle espece
	esp = []
	for g in lst:
		if phylTree.dicGenes[g][0] == taxon_name:
			esp.append(mkNode(taxon_name, g, False))
	if len(esp) == 0:
		return None
	# Fusions pour que l'arbre soit binaire
	while len(esp) >= 2:
		esp.append(mkNode(taxon_name, [(esp.pop(),0), (esp.pop(),0)], True))
	return esp[0]


# Initialisation
#################
print >> sys.stderr, "Preparation des arbres ...",
maxNodeID = max(info)
(newNodeName,age0) = phylTree.parent[options["newSpecies"]]

dicNames = {}
dicTrees = {}
for (i,r) in enumerate(roots):
	dicTrees[i] = r
	# Construction du dictionnaire: nom_gene -> noeud/arbre
	########################################################
	def store(node):
		if node in data:
			for (fils,_) in data[node]:
				store(fils)
		else:
			dicNames[info[node]['gene_name']] = (i,node)
	store(r)

# On definit les endroits ou rajouter les familles
###################################################
print >> sys.stderr, "Insertion des familles ...",
nb0 = 0
nb1 = 0
nb2 = 0
for famille in orthologues:

	# Les genes deja presents dans l'arbre
	lstPar = [dicNames[g] for g in famille.names if g in dicNames]
	dicPar = {}

	ntrees = len(set([x[0] for x in lstPar]))
	if ntrees == 0:

		# On construit un nouvel arbre avec toutes les especes
		#######################################################
		def mkExistingSpeciesNodes(node):
			if node in phylTree.items:
				# On calque les genes sur l'arbre phylogenetique de reference
				lst = [(mkExistingSpeciesNodes(f),l) for (f,l) in phylTree.items.get(node)]
				lst = [(x,l+l0) for ((x,l),l0) in lst if x != None]
				if len(lst) == 0:
					return (None,0)
				# Fusions pour que l'arbre soit binaire
				while len(lst) >= 2:
					lst.append( (mkNode(node, [lst.pop(), lst.pop()], False), 0) )
				return lst[0]
			else:
				return (mkSpeciesNode(node, famille.names),0)

		dicTrees[max(dicTrees)+1] = mkExistingSpeciesNodes(phylTree.root)[0]
		nb0 += 1
		continue

	elif ntrees > 1:
		
		# Les arbres doivent tous correspondre a des especes descendantes
		##################################################################

		# Les racines des arbres a relier
		lst = [dicTrees[x[0]] for x in lstPar]
		# Les taxa respectifs
		lst = set([(x,info[x]['taxon_name']) for x in lst])

		# Les racines doivent etre en dessous de la cible
		if len([x for (_,x) in lst if not phylTree.isChildOf(x, newNodeName)]) > 0:
			continue

		# On regroupe ces arbres
		while len(lst) >= 2:
			(e1,n1) = lst.pop()
			(e2,n2) = lst.pop()
			par = phylTree.dicParents.get(n1).get(n2)
			lst.add( (mkNode(par, [(e1,0), (e2,0)], (par == n1) or (par == n2)),par) )
		# Le nouveau noeud outgroup	
		newNodeID = mkNode(newNodeName, [(lst.pop()[0],1),(mkSpeciesNode(options["newSpecies"], famille.names),1)], False)

		# Mise a jour des racines
		for x in lstPar:
			dicTrees[x[0]] = newNodeID
		nb2 += 1
		continue

	else:
		# Les genes pointent sur un seul arbre, on doit retrouver leur ancetre commun
		##############################################################################
		tree = dicTrees[lstPar[0][0]]

		# La table des parents
		# Construction du dictionnaire: fils -> pere
		def mkPar(node):
			if node in data:
				for (fils,_) in data[node]:
					dicPar[fils] = node
					mkPar(fils)
		mkPar(tree)
		names = [x[1] for x in lstPar]
		# On remonte du gene de reference a la racine
		n0 = n = names.pop()
		refPath = {}
		while n != tree:
			n = dicPar[n]
			refPath[n] = (len(refPath),n)
		# On remonte les autres genes en conservant celui qui va le plus haut
		best = (-1,n0)
		while len(names) > 0:
			n = names.pop()
			while n not in refPath:
				n = dicPar[n]
			best = max(best, refPath[n])
		# -> ancestor = Le noeud de l'ancetre commun
		ancestor = best[1]

		# 2. Identifier le point d'insertion qui corresponde a l'age du nouveau noeud
		path = []
		# On va remonter l'arbre tant qu'on ne depasse pas la cible
		while (ancestor in dicPar) and (phylTree.ages.get(info[dicPar[ancestor]]['taxon_name']) < age0):
			path.append(ancestor)
			ancestor = dicPar[ancestor]
		# Puis on redescend pour eviter les duplications
		while (isDuplicatedNode(info[ancestor])) and (len(path) > 0):
			ancestor = path.pop()

		# L'ancetre commun trouve est trop vieux
		if phylTree.ages.get(info[ancestor]['taxon_name']) > age0:
				continue

		# Insertion des nouvelles donnees
		##################################
		espNode = mkSpeciesNode(options["newSpecies"], famille.names)
		# Il faut creer le noeud juste au dessus de ancestor
		newNodeID = mkNode(newNodeName, [(ancestor,0),(espNode,0)], False)

		if ancestor in dicPar:
			# Mettre a jour le pere
			p = dicPar[ancestor]
			newData = []
			for x in data[p]:
				if x[0] == ancestor:
					newData.append( (newNodeID,x[1]) )
				else:
					newData.append( x )
			data[p] = newData
			# Et ses grands-parents si il y a des duplications
			while isDuplicatedNode(info[p]) and (phylTree.ages.get(info[p]['taxon_name']) < age0):
				info[p]['taxon_name'] = newNodeName
				if p not in dicPar:
					break
				p = dicPar[p]
		else:
			dicTrees[lstPar[0][0]] = newNodeID
		nb1 += 1

print >> sys.stderr, "%d-%d-%d/%d nouveaux genes OK" % (nb0,nb1,nb2,len(orthologues.lstGenes[None]))

print >> sys.stderr, "Mise en forme des arbres ...",
nb = 0
for r in set(dicTrees.itervalues()):
	nb += 1
	utils.myProteinTree.printTree(sys.stdout, data, info, r)
print >> sys.stderr, "%d arbres OK" % nb


