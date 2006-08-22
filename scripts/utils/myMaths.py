#! /usr/bin/python

#
# Fonctions communes de chargement et de traitement des donnees
#

import sys
import math


#
# Renvoie la moyenne d'une liste
#
def moyenne(lst):
	if len(lst) == 0:
		return 0.
	return float(sum(lst))/float(len(lst))

#
# Renvoie la moyenne ponderee d'une liste
#
def moyennePonderee(lst):
	sV = 0.
	sP = 0.
	for (val,poids) in lst:
		sV += val*poids
		sP += poids
	if sP == 0:
		return 0.
	return sV/sP




#
# Ecart type
def ecartType(lst):
	if len(lst) == 0:
		return 0
	m = moyenne(lst)
	return math.sqrt(moyenne([(x-m)*(x-m) for x in lst]))


#
# Min et max
#
def getMinMax(lst):
	mn = lst[0]
	mx = mn
	for x in lst:
		if x > mx:
			mx = x
		elif x < mn:
			mn = x
	return (mn, mx)


#
# Renvoie l'intersection de deux fenetres
#
def intersectionCouples(c1, c2):
	if c1[1] < c2[0] or c2[1] < c1[0]:
		return []
	else:
		return [ max(c1[0], c2[0]), min(c1[1], c2[1]) ]

#
# Cette classe permet de garder la valeur minimale qu'on lui a presente
# Elle utilise une fonction pour comparer les objets
#
class MinKeeper:

	def __init__(self, f):
		self.func = f
		self.empty = True
	
	def setBegin(self, r):
		self.res = r
		self.s = self.func(r)
		self.empty = False
	
	def check(self, r):
		if self.empty:
			self.setBegin(r)
			return
		s = self.func(r)
		if s < self.s:
			self.res = r
			self.s = s
			self.empty = False

