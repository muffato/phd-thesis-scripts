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
import pywalktrap
import myMaths
import myTools


class WalktrapWrapper:

	def __init__(self):
		
		self.edges = {}


	def __internalAddEdge(self, x, y, weight):
		if x not in self.edges:
			self.edges[x] = {}
		self.edges[x][y] = weight

	def addEdge(self, x, y, weight):
		try:
			x = int(x)
		except ValueError:
			pass
		try:
			y = int(y)
		except ValueError:
			pass
		weight = float(weight)
		self.__internalAddEdge(x, y, weight)
		self.__internalAddEdge(y, x, weight)

	def updateFromFile(self, f):
		for l in f:
			c = l.split()
			self.addEdge(c[0], c[1], c[2])

	def updateFromFunc(self, items, func):
		for i1 in xrange(len(items)):
			x1 = items[i1]
			for i2 in xrange(i1+1,len(items)):
				x2 = items[i2]
				score = func(x1, x2)
				if score > 0:
					self.addEdge(x1, x2, score)

	def updateFromDict(self, d):
		for x1 in d:
			for (x2,v) in self.edges[x1].iteritems():
				self.addEdge(x1, x2, v)



	def doWalktrap(self):

		
		# Les composantes connexes
		combin = myTools.myCombinator([])
		for x in self.edges:
			combin.addLink([x] + list(self.edges[x]))

		self.res = []
		
		for nodes in combin:
		
			# Reindexation des noeuds
			indNodes = {}
			for i in xrange(len(nodes)):
				indNodes[nodes[i]] = i
		
			(scores,dend) = pywalktrap.doWalktrap(indNodes, self.edges)
			self.res.append( (nodes, scores, WalktrapDendrogram(dend, nodes)) )


class WalktrapDendrogram:

	def __init__(self, lstMerges, lstNodes):

		self.allMerges = lstMerges
		self.allMerges.sort(reverse = True)
		
		self.lstFils = {}
		for (_,fils,pere) in self.allMerges:
			self.lstFils[pere] = fils
		
		self.lstAll = lstNodes

	def cut(self, scale):

		# Renvoie l'ensemble des fils d'un noeud, recursivement
		def getAllChildren(father):
			if father in self.lstFils:
				fathersAlreadySeen.add(father)
				res = []
				for child in self.lstFils[father]:
					res.extend( getAllChildren(child) )
				return res
			else:
				nodesNotSeen.remove(father)
				return [father]

	
		# On extrait les communautes correspondantes
		lstClusters = []
		fathersAlreadySeen = set()
		nodesNotSeen = set(self.lstAll)
		for (s,_,pere) in self.allMerges:
			if (s < scale) and (pere not in fathersAlreadySeen):
				lstClusters.append( getAllChildren(pere) )
		return (lstClusters, list(nodesNotSeen))


