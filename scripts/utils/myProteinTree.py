# Copyright (c) 2010 IBENS/Dyogen : Matthieu MUFFATO, Alexandra LOUIS, Hugues ROEST CROLLIUS
# Mail : hrc@ens.fr or alouis@biologie.ens.fr
# Licences GLP v3 and CeCILL v2

import sys
import myFile


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


# Renvoie le suffixe associe a un numero de duplication ('a' -> 'z', 'aa' -> 'az', 'ba' ...)
def getDupSuffix(n, upper):
	base = 64 if upper else 96
	assert 1 <= n <= 26*27
	if n <= 26:
		return "." + chr(base + n)
	else:
		return "." + chr(base + (n-1)/26) + chr(base + 1+(n-1)%26)


# Charge l'arbre depuis un fichier
###################################
def loadTree(name):

	global curr

	# Lit la prochaine ligne du fichier (et bufferise la prochaine)
	def nextLine():
		global curr
		old = curr
		try:
			l = ""
			while (l == "") or l.startswith("#"):
				# On enleve le \n final et on coupe suivant les \t
				l = f.next().replace('\n', '')
			l = l.split('\t')
			# On stocke le triplet (indentation,key,value)
			curr = (len(l)-2, l[-2], l[-1])
			eval(l[-1])
		except StopIteration:
			curr = None
		return old

		
	# La procedure d'analyse des lignes du fichier
	def recLoad(indent):
		
		# l'ID du point
		currID = int(nextLine()[2])
		# Les infos associees
		info[currID] = eval(nextLine()[2])
		
		# Des fils ?
		child = []
		while (curr != None) and (curr[0] == indent+1):
			length = float(nextLine()[2])
			child.append( (recLoad(indent+1), length) )
		if len(child) > 0:
			data[currID] = child
			
		return currID

	print >> sys.stderr, "Chargement du fichier d'arbres %s ..." % name,
	f = myFile.openFile(name, "r")
	curr = None
	nextLine()
	n = (0,0,0)
	while True:
		info = {}
		data = {}
		root = recLoad(0)
		yield (root,data,info)
		n = (n[0]+1, n[1]+len(data), n[2]+len(info)-len(data))
		if curr == None:
			break
	print >> sys.stderr, "%d racines, %d branches, %d noeuds OK" % n

	f.close()

