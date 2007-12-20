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
	# A un numero member_id, on associe les noms (gene/transcrit/proteine) et l'espece
	if t[7] != "\\N":
		x = t[8].split()
		tmpLinks[t[0]] = (x[1].split(':')[1], x[0].split(':')[1], t[1], taxonName[t[4]])
f.close()
print >> sys.stderr, len(tmpLinks), "OK"


# On charge les liens node_id -> member_id
###########################################
print >> sys.stderr, "Chargement des liens node_id -> member_id ...",
links = {}
f = utils.myTools.myOpenFile(os.path.join(ensemblURL, options["IN.protein_tree_member"]), "r")
for ligne in utils.myTools.MySQLFileLoader(f):
	t = ligne.split("\t")
	links[ int(t[0]) ] = tmpLinks[t[1]]
f.close()
print >> sys.stderr, len(links), "liens OK"
del tmpLinks
del taxonName


# On charge les liens node_id -> infos 
#######################################
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


# On charge les liens node_id_father -> node_id_son
####################################################
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
		(g,t,p,e) = links[node]
		print >> f, "%snames\t%s\t%s\t%s" % (indent, g, t, p)
		print >> f, "%sspecies\t%s" % (indent, e)

####################################
# Imprime l'arbre au format Newick #
####################################
def printNewickTree(f, node):
	genes = []
	def rec(node):
		if node not in data:
			genes.append(links[node][0])
			return links[node][0]
		else:
			return "(" + ",".join([rec(x) + ":" + str(l) for (x,l)  in data[node]]) + ")"

	tr = rec(node)
	print >> f, " ".join(genes)
	print >> f, tr, ";"



##################################
# Nettoie l'arbre                #
#    - taxon_id                  #
#    - lost_taxon_id             #
##################################
def cleanTree(node):
	
	# Si on est sur une feuille (un gene), on rajoute le taxon_name
	if node not in data:
		info[node]['taxon_name'] = links[node][-1]
		return
	
	inf = info[node]
	# Suppression des infos non necessaires
	if 'lost_taxon_id' in inf:
		del inf['lost_taxon_id']
	if 'taxon_id' in inf:
		del inf['taxon_id']

	# Les appels recursifs
	for (g,d) in data[node]:
		cleanTree(g)


####################################################################
# On considere que les duplications 'dubious' ne sont pas valables #
####################################################################
def isDuplicatedNode(inf):
	return (inf['Duplication'] != 0) and ('dubious_duplication' not in inf)


########################################################
# Redonne la topologie attendue a l'arbre              #
#   Rassemble les noeuds equivalents sous le meme fils #
########################################################
def rebuildTree(node):

	# Fin du process sur une feuille
	if node not in data:
		return

	inf = info[node]

	# On ne change les fils que si ce n'est pas une vraie duplication
	if not isDuplicatedNode(inf):
		
		# Permet d'eviter les noeuds superflus en ayant juste le bon nom de taxon
		while True:

			# On s'assure que le noeud est bien aplati
			#  - si on vient de changer le nom du taxon, celui-ci peut correspondre avec le fils, d'ou la fusion
			#  - si on vient de creer le noeud, des dubious_duplication successifs ont pu faire aterrir ensemble des branches de meme taxon (ex Euarchontoglires)
			flattenTree(node, False)

			# A. On trie les enfants en paquets selon l'arbre phylogenetique
			fils = utils.myTools.defaultdict(list)
			anc = inf['taxon_name']
			for (i,(g,d)) in enumerate(data[node]):
				gname = info[g]['taxon_name']
				for a in phylTree.branches[anc]:
					if phylTree.isChildOf(gname, a):
						fils[a].append( (g,d) )
						break
				else:
					# Apparait si g est le meme ancetre que node et que g est une duplication, ce qui l'a empeche d'etre aplati par flatten
					#if gname != anc:
					#	print >> sys.stderr, "PB1", node, anc, gname
					fils[anc].append( (g,d) )
			
			# Si il y a un seul fils, qui correspond a un vrai enfant, on renomme et recommence
			if (len(fils) == 1) and (anc not in fils):
				info[node]['taxon_name'] = fils.popitem()[0]
			else:
				break

		if len(fils) == 1:
			# Ici, si il y a un seul fils, celui-ci a le meme taxon que son pere
			# -> On ne change rien
			pass
		else:
			if len(fils) == 2:
				# Fonctionnement normal
				items = fils.items()

			else:
				# Ici, len(fils) == 3, On a les deux sous-branches + le meme ancetre
				# On regroupe en deux
				lst1 = fils.pop(anc)
				lst2 = []
				for tmp in fils.itervalues():
					lst2.extend(tmp)
				items = [(anc,lst1), (anc,lst2)]

			# On construit la nouvelle structure data
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
	if isDuplicatedNode(inf):
		return

	newData = []
	taxonName = inf['taxon_name']
	for (g,d) in data[node]:
		inf = info[g]
		# 2x le meme taxon et pas de duplication
		if (g in data) and (inf['taxon_name'] == taxonName) and not isDuplicatedNode(inf):
			newData.extend([(g2,d+d2) for (g2,d2) in data[g]])
		else:
			#if (g in data) and (inf['taxon_name'] == taxonName) and isDuplicatedNode(inf):
			#	print >> sys.stderr, "SOURCE", node, info[node], (g,d), inf
			newData.append( (g,d) )
	data[node] = newData

		
	

###########################################
# Sauvegarde toutes les familles de genes #
###########################################
def extractGeneFamilies(node, previousAnc, lastWrittenAnc):

	newAnc = info[node]['taxon_name']

	# Les noeuds ou ecrire les familles
	if lastWrittenAnc == None:
		if previousAnc == None:
			toWrite = [newAnc]
		else:
			toWrite = phylTree.dicLinks[previousAnc][newAnc]
	else:
		toWrite = phylTree.dicLinks[lastWrittenAnc][newAnc][1:]
	if isDuplicatedNode(info[node]):
		toWrite = toWrite[:-1]
	
	if len(toWrite) == 0:
		newLastWritten = lastWrittenAnc
	else:
		newLastWritten = toWrite[-1]

	if (lastWrittenAnc == None) and (newLastWritten != None):
		trueRoots.append(node)
 
	# Les genes des descendants
	if node in data:
		allGenes = []
		for (g,d) in data[node]:
			allGenes.extend( extractGeneFamilies(g, newAnc, newLastWritten) )
	else:
		allGenes = [links[node][0]]

	for a in toWrite:
		geneFamilies[a].append( allGenes )

	return allGenes


# On a besoin des genomes modernes pour reconnaitre les genes
geneFamilies = utils.myTools.defaultdict(list)
ft1 = utils.myTools.myOpenFile(options["OUT.tree"] + ".1", "w")
ft2 = utils.myTools.myOpenFile(options["OUT.tree"] + ".2", "w")
ft3 = utils.myTools.myOpenFile(options["OUT.tree"] + ".3", "w")
ft4 = utils.myTools.myOpenFile(options["OUT.tree"] + ".4", "w")
ft5 = utils.myTools.myOpenFile(options["OUT.tree"] + ".5", "w")
print >> sys.stderr, "Mise en forme des arbres ...",
nb = 0
nextNodeID = max(data) + 1
for (root,_) in data[1]:
#for (root,_) in [(200984,0)]:
#for (root,_) in [(753697,0)]:
	if 'taxon_name' in info[root]:
		trueRoots = []
		printTree(ft1, 0, root)
		cleanTree(root)
		printTree(ft2, 0, root)
		flattenTree(root, True)
		printTree(ft3, 0, root)
		rebuildTree(root)
		trueRoots = []
		extractGeneFamilies(root, None, None)
		printTree(ft4, 0, root)
		for r in trueRoots:
			printNewickTree(ft5, r)
		nb += 1
print >> sys.stderr, nb, "arbres OK"
ft1.close()
ft2.close()
ft3.close()
ft4.close()
ft5.close()

for (anc,lst) in geneFamilies.iteritems():
	print >> sys.stderr, "Ecriture des familles de %s ..." % anc,
	f = utils.myTools.myOpenFile(options["OUT.ancGenesFile"] % phylTree.fileName[anc], "w")
	for gg in lst:
		print >> f, " ".join(gg)
	f.close()
	print >> sys.stderr, len(lst), "OK"

