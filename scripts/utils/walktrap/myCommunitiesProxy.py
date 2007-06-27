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
import _walktrap
import utils.myTools


class WalktrapProxy:

	def __init__(self):
		self.edges = {}

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
		self.edges.setdefault(x,{})[y] = weight
		self.edges.setdefault(y,{})[x] = weight

	def updateFromFile(self, f):
		for l in f:
			c = l.split()
			self.addEdge(c[0], c[1], c[2])

	def updateFromFunc(self, items, func):
		for (x1,x2) in utils.myTools.myMatrixIterator(items, None, utils.myTools.myMatrixIterator.StrictUpperMatrix):
			score = func(x1, x2)
			if score > 0:
				self.addEdge(x1, x2, score)

	def updateFromDict(self, d):
		for x1 in d:
			for (x2,v) in self.edges[x1].iteritems():
				self.addEdge(x1, x2, v)


	# TODO Amener les options de walktrap
	def doWalktrap(self, internal = True):

		
		# Les composantes connexes
		combin = utils.myTools.myCombinator([])
		for x in self.edges:
			combin.addLink([x] + self.edges[x].keys())

		self.res = []
		
		for nodes in combin:
		
			# Reindexation des noeuds
			indNodes = {}
			for i in xrange(len(nodes)):
				indNodes[nodes[i]] = i
		
			if internal:
				(scores,dend) = _walktrap.doWalktrap(indNodes, self.edges)
			else:
				(stdin,stdout,stderr) = os.popen3('/users/ldog/muffato/work/scripts/utils/walktrap/walktrap -t10')
				stderr.close()
				
				# Envoi des donnees du graphe
				for x in nodes:
					for y in self.edges[x]:
						print >> stdin, indNodes[x], indNodes[y], self.edges[x][y]
				stdin.close()
				
				(scores,dend) = self.loadWalktrapOutput(stdout)
			
			self.res.append( (nodes, scores, WalktrapDendrogram(dend, nodes)) )


	#
	# Chargement d'un fichier de resultat de walktrap
	#
	def loadWalktrapOutput(self, file):

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
			
			allMerges.append( (scale,tuple([int(x) for x in l[0].split('+')]),int(l[1])) )

		allMerges.sort( reverse = True )

		lstCoup = []
		for line in file:
			try:
				# On extrait les lignes "alpha relevance"
				c = line.split()
				lstCoup.append( (float(c[0]),float(c[1])) )
			except ValueError:
				pass
		
		return (lstCoup, allMerges)



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


