#! /users/ldog/muffato/python -OO

#
# Fonctions communes de traitement des graphes
#

import sys
import myMaths
import myTools


class Graph:

	def __init__(self):

		self.lstAretes = dict()
		self.lstSommets = set([])
		self.lstAretesSortantes = dict()

	def addSommet(self, x):
		if x not in self.lstSommets:
			self.lstSommets.add(x)
			self.lstAretesSortantes[x] = set([])


	def addArete(self, x, y, poids=1):
		self.lstAretes[(x,y)] = poids
		self.lstAretesSortantes[x].add( (y,poids) )
		self.lstAretes[(y,x)] = poids
		self.lstAretesSortantes[y].add( (x,poids) )

	
	def buildReducedGraph(self):

		newGr = Graph()
		comb = myTools.myCombinator()
		for x in self.lstSommets:
			if len(self.lstAretesSortantes[x]) == 2:
				r = []
				for (y,_) in self.lstAretesSortantes[x]:
					r.extend( [(x,y),(y,x)] )
				comb.addLink(r)
				

