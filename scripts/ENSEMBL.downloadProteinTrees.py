#! /users/ldog/muffato/python -OO

__doc__ = """
	Telecharge depuis le site d'Ensembl les fichiers des arbres de proteines
"""


# Librairies
import os
import sys
import utils.myTools
import utils.myPhylTree


# Arguments
(noms_fichiers, options) = utils.myTools.checkArgs( \
	["phylTree.conf"], \
	[("releaseID",int,[47]), \
	("IN.EnsemblURL",str,"ftp://ftp.ensembl.org/pub/release-XXX/mysql/ensembl_compara_XXX"), \
	("IN.member",str,"member.txt.gz"), \
	("IN.genome_db",str,"genome_db.txt.gz"), \
	("IN.protein_tree_node",str,"protein_tree_node.txt.gz"), \
	("IN.protein_tree_member",str,"protein_tree_member.txt.gz"), \
	("IN.protein_tree_tag",str,"protein_tree_tag.txt.gz"), \
	("OUT.ancGenesFile",str,""), \
	("OUT.tree",str,""), \
	], \
	__doc__ \
)


phylTree = utils.myPhylTree.PhylogeneticTree(noms_fichiers["phylTree.conf"])
ensemblURL = options["IN.EnsemblURL"].replace("XXX", str(options["releaseID"]))

# On charge les liens taxon_id -> species name
###############################################
print >> sys.stderr, "Chargement des liens taxon_id -> species name ...",
taxonName = {}
f = utils.myTools.myOpenFile(os.path.join(ensemblURL, options["IN.genome_db"]), "r")
for ligne in utils.myTools.MySQLFileLoader(f):
	t = ligne.split("\t")
	taxonName[t[1]] = t[2]
f.close()
print >> sys.stderr, len(taxonName), "OK"


# On charge les liens member_id -> protein name
################################################
print >> sys.stderr, "Chargement des liens member_id -> protein name ...",
tmpLinks = {}
f = utils.myTools.myOpenFile(os.path.join(ensemblURL, options["IN.member"]), "r")
for ligne in utils.myTools.MySQLFileLoader(f):
	t = ligne.split("\t")
	# A un numero member_id, on associe un nom et eventuellement un autre numero (pour passer de la proteine au gene)
	tmpLinks[t[0]] = (t[1], t[7], taxonName.get(t[4], None))
f.close()
print >> sys.stderr, len(tmpLinks), "OK"


# On charge les liens node_id -> member_id
###########################################
print >> sys.stderr, "Chargement des liens node_id -> member_id ...",
links = {}
f = utils.myTools.myOpenFile(os.path.join(ensemblURL, options["IN.protein_tree_member"]), "r")
for ligne in utils.myTools.MySQLFileLoader(f):
	t = ligne.split("\t")
	# Le nom de la proteine et l'indice du nom du gene
	(p,ig,esp) = tmpLinks[ t[1] ]
	# A un numero node_id, on associe une proteine et un gene
	links[ int(t[0]) ] = (p, tmpLinks[ig][0], esp)
f.close()
print >> sys.stderr, len(links), "liens OK"
del tmpLinks
del taxonName


# On charge les infos supplementaires
######################################
print >> sys.stderr, "Chargement des infos des noeuds ...",
info = utils.myTools.defaultdict(dict)
f = utils.myTools.myOpenFile(os.path.join(ensemblURL, options["IN.protein_tree_tag"]), "r")
for ligne in utils.myTools.MySQLFileLoader(f):
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
	info[ int(t[0]) ][t[1]] = t[2]
f.close()
print >> sys.stderr, len(info), "infos OK"


# On charge l'arbre
####################
print >> sys.stderr, "Chargement des arbres ...",
data = utils.myTools.defaultdict(list)
f = utils.myTools.myOpenFile(os.path.join(ensemblURL, options["IN.protein_tree_node"]), "r")
for ligne in utils.myTools.MySQLFileLoader(f):
	t = ligne.split("\t")
	# On rajoute le fils de son pere (avec une distance)
	data[ int(t[1]) ].append( (int(t[0]), float(t[5])) )
f.close()
print >> sys.stderr, len(data), "branches OK"


######################################################
# Imprime l'arbre sous forme indentee avec les infos #
######################################################
def printTree(f, n, node):
	indent = "\t" * n
	# L'id du noeud
	print >> f, "%sid\t%d" % (indent, node)
	if node in data:
		# Ses infos
		print >> f, "%sinfo\t%s" % (indent, info[node])
		for (g,d) in data[node]:
			# Chaque sous-branche avec sa longueur
			print >> f, "%s\tlen\t%f" % (indent, d)
			printTree(f, n+1, g)
	else:
		# Une feuille: un gene
		(p,g,e) = links[node]
		print >> f, "%snames\t%s\t%s" % (indent, p, g)
		print >> f, "%sspecies\t%s" % (indent, e)


##################################
# Nettoie l'arbre                #
# Renomme certains taxons:       #
#    Smegmamorpha -> Percomorpha #
#    Eutheria ~> Boreoeutheria   #
# Supprime certains infos:       #
#    - taxon_id                  #
#    - lost_taxon_id             #
##################################
def cleanTree(node):
	
	# Si on est sur une feuille (un gene)
	if node not in data:
		info[node]['taxon_name'] = links[node][2]
		return
	
	inf = info[node]
	# Suppression des infos non necessaires
	if 'lost_taxon_id' in inf:
		del inf['lost_taxon_id']
	if 'taxon_id' in inf:
		del inf['taxon_id']

	# Remplacement de Smegmamorpha par Percomorpha
	#if inf['taxon_name'] == 'Smegmamorpha':
	#	inf['taxon_name'] = 'Percomorpha'
	#if inf['taxon_name'] == 'Coelomata':
	#	inf['taxon_name'] = 'Bilateria'

	# Les appels recursifs
	for (g,d) in data[node]:
		cleanTree(g)


########################################################
# Redonne la topologie attendue a l'arbre              #
#   Rassemble les noeuds equivalents sous le meme fils #
########################################################
def rebuildTree(node):

	# Fin du process sur une feuille
	if node not in data:
		return

	inf = info[node]
	# Si c'est une vraie duplication, on n'a plus rien a faire
	if (inf['Duplication'] == 0) or ('dubious_duplication' in inf):
		# A. On trie les enfants en paquets selon l'arbre phylogenetique
		anc = inf['taxon_name']
		fils = utils.myTools.defaultdict(list)
		for (i,(g,d)) in enumerate(data[node]):
			gname = info[g]['taxon_name']
			for a in phylTree.branches[anc]:
				if phylTree.isChildOf(gname, a):
					fils[a].append( (g,d) )
					break
			else:
				fils[g].append( (g,d) )
		
		# B. Si un paquet a plus de 2 representants, on fonce dedans
		newData = []
		for (anc,lst) in fils.iteritems():
			if len(lst) == 1:
				newData.extend( lst )
			elif len(lst) > 1:
				global nextNodeID
				nextNodeID += 1
				newData.append( (nextNodeID,0) )
				data[nextNodeID] = lst
				info[nextNodeID] = {'taxon_name':anc, 'Duplication':0, 'Bootstrap':-1}
				flattenTree(nextNodeID, False)
	
		data[node] = newData
	
	for (g,_) in data[node]:
		rebuildTree(g)


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
	inf = info[node]
	# Si c'est une vraie duplication, on n'a plus rien a faire
	if (inf['Duplication'] != 0) and ('dubious_duplication' not in inf):
		return

	newData = []
	taxonName = inf['taxon_name']
	for (g,d) in data[node]:
		inf = info[g]
		# 2x le meme taxon et pas de duplication
		if (g in data) and (inf['taxon_name'] == taxonName) and ((inf['Duplication'] == 0) or ('dubious_duplication' in inf)):
			newData.extend([(g2,d+d2) for (g2,d2) in data[g]])
		else:
			newData.append( (g,d) )
	data[node] = newData

		
	

###########################################
# Sauvegarde toutes les familles de genes #
###########################################
#def extractGeneFamilies(node, previousAnc, previousDup):
def extractGeneFamilies(node, lastAnc):

	newAnc = info[node]['taxon_name']
	newDup = ((info[node]['Duplication'] != 0) and ('dubious_duplication' not in info[node]))

	# Les noeuds ou ecrire les familles
	if lastAnc == None:
		toWrite = [newAnc]
	else:
		toWrite = phylTree.dicLinks[lastAnc][newAnc][1:]
	if newDup:
		toWrite = toWrite[:-1]
	if len(toWrite) == 0:
		newAnc = lastAnc
	else:
		newAnc = toWrite[-1]

	#toWrite = phylTree.dicLinks[previousAnc][newAnc]
	#if newDup:
	#	toWrite = toWrite[:-1]
	#if not previousDup:
	#	toWrite = toWrite[1:]

	# Les genes des descendants
	if node in data:
		allGenes = []
		for (g,d) in data[node]:
			#allGenes.extend( extractGeneFamilies(g, newAnc, newDup) )
			allGenes.extend( extractGeneFamilies(g, newAnc) )
	else:
		allGenes = [links[node][1]]

	for a in toWrite:
		geneFamilies[a].append( allGenes )

	return allGenes


# On a besoin des genomes modernes pour reconnaitre les genes
geneFamilies = utils.myTools.defaultdict(list)
ft1 = utils.myTools.myOpenFile(options["OUT.tree"] + ".1", "w")
ft2 = utils.myTools.myOpenFile(options["OUT.tree"] + ".2", "w")
ft3 = utils.myTools.myOpenFile(options["OUT.tree"] + ".3", "w")
ft4 = utils.myTools.myOpenFile(options["OUT.tree"] + ".4", "w")
print >> sys.stderr, "Mise en forme des arbres ...",
nb = 0
nextNodeID = max(data) + 1
for (root,_) in data[1]:
	if 'taxon_name' in info[root]:
		printTree(ft1, 0, root)
		cleanTree(root)
		printTree(ft2, 0, root)
		flattenTree(root, True)
		printTree(ft3, 0, root)
		rebuildTree(root)
		#extractGeneFamilies(root, None, False)
		#extractGeneFamilies(root, info[root]['taxon_name'], True)
		extractGeneFamilies(root, None)
		printTree(ft4, 0, root)
		nb += 1
print >> sys.stderr, nb, "arbres OK"
ft1.close()
ft2.close()
ft3.close()
ft4.close()

for (anc,lst) in geneFamilies.iteritems():
	#if anc in phylTree.listSpecies:
	#	continue
	print >> sys.stderr, "Ecriture des familles de %s ..." % anc,
	f = utils.myTools.myOpenFile(options["OUT.ancGenesFile"] % phylTree.fileName[anc], "w")
	for gg in lst:
		print >> f, " ".join(gg)
	f.close()
	print >> sys.stderr, len(lst), "OK"

