#! /users/ldog/muffato/python

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
arguments = utils.myTools.checkArgs( \
	[("phylTree.conf",file), ("proteinTree",file), ("orthologues",file), ("newSpecies",str)], \
	[("genesFile",str,"")], \
	__doc__ \
)

# Chargement des fichiers
##########################
phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])
phylTree.loadAllSpeciesSince(None, arguments["genesFile"], storeGenomes = False)
orthologues = utils.myGenomes.Genome(arguments["orthologues"])

roots = []
data = {}
info = {}
for (r,d,i) in utils.myProteinTree.loadTree(arguments["proteinTree"]):
	roots.append(r)
	data.update(d)
	info.update(i)

# Creation d'un nouveau noeud dans l'arbre
###########################################
def mkNode(taxon_name, x, dup):
	global maxNodeID
	maxNodeID += 1
	dup = 3*int(dup)
	if type(x) == str:
		info[maxNodeID] = {'taxon_name':taxon_name, 'Duplication':dup, 'gene_name':x}
	else:
		info[maxNodeID] = {'taxon_name':taxon_name, 'Duplication':dup, 'Bootstrap':-1}
		data[maxNodeID] = x
	return maxNodeID

# Creation d'un arbre pour une liste de genes d'une meme espece
################################################################
def mkSpeciesNode(taxon_name, lst):
	# Selection des genes de l'espece
	esp = [g for g in lst if (g in phylTree.dicGenes) and (phylTree.dicGenes[g][0] == taxon_name)]
	
	if len(esp) == 0:
		return None
	
	# Creer les feuilles de la nouvelle espece
	esp = [mkNode(taxon_name, g, False) for g in esp]

	if len(esp) == 1:
		return esp[0]
	else:
		return mkNode(taxon_name, [(x,0) for x in esp], True)


# Initialisation
#################
print >> sys.stderr, "Preparation des arbres ...",
maxNodeID = max(info)
for lst in data.itervalues():
	maxNodeID = max(maxNodeID, max([x[0] for x in lst]))
newSpecies = phylTree.officialName[arguments["newSpecies"]]

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
				return (mkNode(node, lst, False), 0)
			else:
				return (mkSpeciesNode(node, famille.names),0)

		dicTrees[max(dicTrees)+1] = mkExistingSpeciesNodes(phylTree.root)[0]
		nb0 += 1

	elif ntrees > 1:
		
		# Les arbres doivent tous correspondre a des especes descendantes
		##################################################################

		# Les racines des arbres a relier
		lst = [dicTrees[x[0]] for x in lstPar]
		# Les taxa respectifs
		lst = set([(x,info[x]['taxon_name']) for x in lst])

		# Les racines doivent etre en dessous de la cible
		if phylTree.isChildOf(newSpecies, phylTree.lastCommonAncestor([x[1] for x in lst])):
			continue

		# On regroupe ces arbres
		while len(lst) >= 2:
			(e1,n1) = lst.pop()
			(e2,n2) = lst.pop()
			par = phylTree.dicParents.get(n1).get(n2)
			lst.add( (mkNode(par, [(e1,0), (e2,0)], (par == n1) or (par == n2)),par) )

		# Le nouveau noeud outgroup	
		(e,n) = lst.pop()
		par = phylTree.dicParents.get(n).get(newSpecies)
		newNodeID = mkNode(par, [(e,0),(mkSpeciesNode(newSpecies, famille.names),0)], False)

		# Mise a jour des racines
		for x in lstPar:
			dicTrees[x[0]] = newNodeID
		nb2 += 1

	else:
		
		def prepareTrip(node):
			if node in data:
				s = 0
				for (f,_) in data[node]:
					dicPar[f] = node
					s += prepareTrip(f)
			else:
				s = int(node in target)
			descDic[node] = s
			return s
		
		target = set([x[1] for x in lstPar])
		node = dicTrees[lstPar[0][0]]
		descDic = {}
		dicPar = {}
		prepareTrip(node)
		# Parcourt l'arbre jusqu'a arriver sur le point d'insertion
		while node in data:
			
			ref = info[node]['taxon_name']
			#print >> sys.stderr, "ON", node, ref

			# Duplication
			if info[node]['Duplication'] >= 2:
				# Le nombre de hits par branche
				lst = [(descDic[f],f) for (f,_) in data[node]]
				best = max(lst)
				#print >> sys.stderr, "DUPLICATION", sorted(lst)
				
				# On est au toujours au dessus
				if phylTree.isChildOf(newSpecies, ref):
					# On fonce dans la branche avec le plus de hits
					node = best[1]
					continue
				# On est passe en dessous
				else:
					# On prefere un point d'insertion plus precis
					if (best[0] == descDic[node]) and (descDic[node] != 0):
						node = best[1]
						continue
					# Par defaut, on s'arrete
					else:
						break
			
			# On est toujours au dessus
			if phylTree.isChildOf(newSpecies, ref):
				try:
					nextInterm = phylTree.dicLinks.get(newSpecies).get(ref)[-2]
				except IndexError:
					print newSpecies, ref, phylTree.dicLinks.get(newSpecies).get(ref)
				lst = [f for (f,_) in data[node] if phylTree.isChildOf(info[f]['taxon_name'], nextInterm)]
				#print >> sys.stderr, "SPECIATION", lst
				node = lst[0]
				continue
			else:
				#print >> sys.stderr, "END", node
				break

		# Insertion des nouvelles donnees
		##################################
		espNode = mkSpeciesNode(newSpecies, famille.names)
		# Il faut creer le noeud juste au dessus de node
		newName = phylTree.dicParents.get(info[node]['taxon_name']).get(newSpecies)
		newNodeID = mkNode(newName, [(node,0),(espNode,0)], newName in [info[node]['taxon_name'],newSpecies])

		if node in dicPar:
			# Mettre a jour le pere
			p = dicPar[node]
			newData = []
			for x in data[p]:
				if x[0] == node:
					newData.append( (newNodeID,x[1]) )
				else:
					newData.append( x )
			data[p] = newData

			# Et ses grands-parents si il y a des duplications
			while phylTree.isChildOf(info[p]['taxon_name'], newName):
				info[p]['taxon_name'] = newName
				info[p]['Duplication'] = 3

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


