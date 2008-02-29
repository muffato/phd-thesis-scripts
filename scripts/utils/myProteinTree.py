
import sys
import myTools


##########################################################
# Gere les arbres de proteines sous la forme (data,info) #
##########################################################

# Indique les noms des ancetres depuis le dernier lu
#####################################################
def getIntermediateAnc(phylTree, previousAnc, lastWrittenAnc, newAnc, isDuplication):

	# Les noeuds ou ecrire les familles
	if lastWrittenAnc == None:
		if previousAnc == None:
			toWrite = [newAnc]
		else:
			toWrite = phylTree.dicLinks[previousAnc][newAnc]
	else:
		toWrite = phylTree.dicLinks[lastWrittenAnc][newAnc][1:]
	if isDuplication:
		toWrite = toWrite[:-1]
	
	if len(toWrite) == 0:
		newLastWritten = lastWrittenAnc
	else:
		newLastWritten = toWrite[-1]

	root = ((lastWrittenAnc == None) and (newLastWritten != None))
	return (toWrite,newLastWritten,root)


# Imprime l'arbre sous forme indentee avec les infos
######################################################
def printTree(f, data, info, root):

	def rec(f, n, node):
		indent = "\t" * n
		# L'id du noeud
		print >> f, "%sid\t%d" % (indent, node)
		# Les infos
		print >> f, "%sinfo\t%s" % (indent, info[node])
		# Les enfants
		for (g,d) in data.get(node,[]):
			print >> f, "%s\tlen\t%f" % (indent, d)
			rec(f, n+1, g)

	rec(f, 0, root)


# Imprime l'arbre au format Newick
####################################
def printNewickTree(f, data, info, root):
	genes = []
	def rec(node):
		if node not in data:
			genes.append(info[node]['gene_name'])
			return links[node][0]
		else:
			return "(" + ",".join([rec(x) + ":" + str(l) for (x,l)  in data[node]]) + ") " + info[node]['family_name']

	tr = rec(root)
	print >> f, " ".join(genes)
	print >> f, tr, ";"


# Charge l'arbre depuis un fichier
###################################
def loadTree(name):

	# Chargement de toutes les lignes et mise en forme sommaire
	def loadFile(name):
		lignes = []
		f = myTools.myOpenFile(name, "r")
		for ligne in f:
			# On enleve le \n final et on coupe suivant les \t
			l = ligne.replace('\n', '').split('\t')
			# On stocke le triplet (indentation,key,value)
			lignes.append( (len(l)-2,l[-2],l[-1]) )
		f.close()
		return lignes
		
	# La procedure d'analyse des lignes du fichier
	def recLoad(indent):
		
		# l'ID du point
		currID = int(lignes.pop()[2])
		# Les infos associees
		info[currID] = eval(lignes.pop()[2])

		# Des fils ?
		child = []
		while (len(lignes) > 0) and (lignes[-1][0] == indent+1):
			length = float(lignes.pop()[2])
			child.append( (recLoad(indent+1), length) )
		if len(child) > 0:
			data[currID] = child
			
		return currID

	print >> sys.stderr, "Chargement du fichier d'arbres %s ..." % name,
	lignes = loadFile(name)
	lignes.reverse()

	roots = []
	info = {}
	data = {}
	
	print >> sys.stderr, "%d lignes, Analyse ..." % len(lignes),
	while len(lignes) > 0:
		roots.append(recLoad(0))
	print >> sys.stderr, "%d racines, %d branches, %d noeuds OK" % (len(roots),len(data),len(info)-len(data))

	return (data, info, roots)

