#! /users/ldog/muffato/python -OO

#
# Fonctions communes de traitement des communautes
# Formats d'entree des donnees:
#  - Fonction de score
#  - Dictionnaire d'aretes
#  - Fichier "node1 node2 weight"
# Format de sortie des donnees:
#  - Liste pour chaque composante connexe
#    - Liste pour chaque score pertinent
#      - alpha
#      - relevance
#      - liste des clusters
#      - liste des noeuds en dehors des clusters


import os
import sys
import myMaths
import myGenomes
import myTools


#
# Creation du graphe (associations des noeuds avec poids des aretes)
# Separation en composantes connexes
#
def makeConnectedComponents(file = None, items = None, scoreFunc = None, edgesDict = {}):

	# Chargement depuis un fichier
	if file != None:
		for l in file:
			c = l.split()
			# Format = "node1 node2 weight"
			if c[0] not in edgesDict:
				edgesDict[c[0]] = {}
			edgesDict[c[0]][c[1]] = float(c[2])
	
	# Chargement avec une fonction
	if (items != None) and (scoreFunc != None):
		# Iteration des couples (i1,i2)
		for (x1,x2) in myTools.myMatrixIterator(items, None, myTools.myMatrixIterator.StrictUpperMatrix):
			# Le poids de l'arete
			score = scoreFunc(x1, x2)
			if score > 0:
				if x1 not in edgesDict:
					edgesDict[x1] = {}
				edgesDict[x1][x2] = score

	# Creation des composantes connexes
	combin = myTools.myCombinator([])
	for x in items:
		s = set(edgesDict.get(x,{}).iterkeys()).intersection(items)
		combin.addLink([x] + list(s))

	return (edgesDict, combin)


#
# Division en communautes
#
def launchCommunitiesBuild(**args):

	print >> sys.stderr, "Computing weights ...",

	# Creation des composantes connexes
	(edgesDict, combin) = makeConnectedComponents(**args)
	
	print >> sys.stderr, "Launching walktrap ",
	
	totalRes = []
	for nodes in combin:
		
		# Reindexation des noeuds
		indNodes = {}
		for i in xrange(len(nodes)):
			indNodes[nodes[i]] = i
		
		(stdin,stdout,stderr) = os.popen3('/users/ldog/muffato/work/scripts/walktrap/walktrap -t10')
		stderr.close()
		
		# Envoi des donnees du graphe
		for x in nodes:
			if x not in edgesDict:
				continue
			for y in nodes:
				if y in edgesDict[x]:
					print >> stdin, indNodes[x], indNodes[y], edgesDict[x][y]
		stdin.close()

		# Mise en forme des clusters (reindexation)
		res = []
		for (alpha,relevance,clusters) in loadCommunitiesFile(stdout):
			newC = [[nodes[i] for i in c] for c in clusters]
			lonelyNodes = set(nodes).difference(myMaths.flatten(newC))
			res.append( (alpha, relevance, newC, list(lonelyNodes)) )
		if len(res) == 0:
			res.append( (1,0,[nodes],[]) )
		
		totalRes.append(res)
		sys.stderr.write(".")
	
	print >> sys.stderr, " OK"
	return totalRes



#
# Chargement d'un fichier de resultat de walktrap
#
def loadCommunitiesFile(file):

	# Renvoie l'ensemble des fils d'un noeud, recursivement
	def getAllChildren(father):
		alreadySeen.add(father)
		if father in lstFils:
			res = []
			for child in lstFils[father]:
				res.extend( getAllChildren(child) )
			return res
		else:
			return [father]
	
	# On charge les fusions
	allMerges = []
	lstFils = {}
	for line in file:
		if line == "\n":
			break
		
		l = line.split(':')
		scale = float(l[0])
		l = l[1].split('-->')
		father = int(l[1])
		
		lstFils[father] = [int(x) for x in l[0].split('+')]
		allMerges.append( (scale,father) )

	allMerges.sort( reverse = True )

	res = []
	for line in file:
		try:
			# On cherche les lignes "alpha relevance)
			c = line.split()
			scale = float(c[0])
			relevance = float(c[1])

			# On extrait les communautes correspondantes
			alreadySeen = set()
			lstClusters = []
			for (s,pere) in allMerges:
				if (s < scale) and (pere not in alreadySeen):
					lstClusters.append( getAllChildren(pere) )
			res.append( (scale, relevance, lstClusters) )
			sys.stderr.write('+')

		except ValueError:
			pass
	
	return res

