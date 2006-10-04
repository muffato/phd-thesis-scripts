#! /usr/bin/python2.4


import sys
import os
from bz2 import BZ2File


##################################################################
# Cette classe ouvre le fichier en le decompressant s'il le faut #
#   Retourne l'objet FILE et le nom complet du fichier           #
##################################################################
def myOpenFile(nom):
	nom = nom.replace("~", os.environ['HOME'])
	if nom.endswith("bz2"):
		f = BZ2File(nom, 'r')
	else:
		f = open(nom, 'r')
	return f

###################################################################
# Une classe pour avoir un iterateur a deux dimensions rapidement #
# Cela equivaut donc a un parcours de matrice                     #
###################################################################
class myMatrixIterator:

	WholeMatrix = 1
	UpperMatrix = 2
	StrictUpperMatrix = 3
	OnlyDiag = 4
	WholeWithoutDiag = 5

	def __init__(self, n, p, mode):
		self.n = n
		self.p = p
		self.mode = mode

	def __iter__(self):
		class __InternalIterator:
			def __init__(self, mat):
				self.mat = mat
				self.i = -1
				if mat.mode == myMatrixIterator.StrictUpperMatrix or mat.mode == myMatrixIterator.WholeWithoutDiag:
					self.i = 0
				self.j = mat.p
			def next(self):
				if self.mat.mode != myMatrixIterator.OnlyDiag:
					self.j += 1
					if self.mat.mode == myMatrixIterator.WholeWithoutDiag and self.i == self.j:
						self.j += 1
						
					if self.mat.mode == myMatrixIterator.UpperMatrix:
						borne = min(self.i + 1, self.mat.p)
					elif self.mat.mode == myMatrixIterator.StrictUpperMatrix:
						borne = min(self.i, self.mat.p)
					else:
						borne = self.mat.p
					
					if self.j >= borne:
						self.j = 0
						self.i += 1
						if self.i >= self.mat.n:
							raise StopIteration
				else:
					self.i += 1
					self.j = self.i
					if self.j == self.mat.p or self.i == self.mat.n:
						raise StopIteration
				return (self.i, self.j)
		return __InternalIterator(self)
			
########################################################################
# Cette classe permet de regrouper une liste d'elements                #
# Partant d'une liste initiale, on ajoute des liens entre des elements #
#  de la liste et la classe les regroupe                               #
########################################################################
class myCombinator:

	def __init__(self, ini):
		self.grp = ini
		self.dic = {}
		for i in range(len(ini)):
			for x in ini[i]:
				self.dic[x] = i
	
	def addLink(self, obj):
		if len(obj) == 0:
			return
		
		# Le 1er objet de la liste
		a = obj[0]
		if a in self.dic:
			i = self.dic[a]
		else:
			i = len(self.grp)
			self.grp.append(set([a]))
			self.dic[a] = i

		for x in obj[1:]:
			if x in self.dic:
				j = self.dic[x]
				if i == j:
					continue
				for b in self.grp[j]:
					self.dic[b] = i
				self.grp[i].update(self.grp[j])
				self.grp[j] = set([])
			else:
				self.grp[i].add(x)
				self.dic[x] = i

	def getGrp(self):
		return [g for g in self.grp if len(g) > 0]


#################################################################################
# Lit la ligne de commande et parse les arguments                               #
#  1. les arguments obligatoires (des noms de fichiers)                         #
#  2. des options sous la forme -opt=val                                        #
# En cas d'erreur, affiche la syntaxe demandee et une courte description (info) #
# -options- est un tableau de triplets (nom, constructeur, val_defaut)          #
#################################################################################
def checkArgs(args, options, info):

	#
	# Affiche le message d'erreur de mauvais arguments
	#
	def error_usage():
		s = "- ERREUR - Usage : " + sys.argv[0]
		for t in args:
			s += " " + t
		print >> sys.stderr, s
		for t in options:
			if t[1] == bool:
				invite = "+/-"
			else:
				invite = "-"
			print >> sys.stderr, "\t", invite + "%s [%s] (%s)" % t
		if info != "":
			print >> sys.stderr, "\n", info
		sys.exit(1)

	types = dict([ (x[0], x[1]) for x in options ])
	valOpt = dict([ (x[0], x[2]) for x in options ])
	valArg = []
	
	# On scanne les argumetns pour les compter et recuperer les valeurs
	for t in sys.argv[1:]:

		# Un petit peu d'aide
		if t == '-h' or t == '--help' or t == '-help':
			error_usage()
			
		# Un argument optionnel
		if t[0] in ['-', '+']:
		
			# Un parametre non bool
			if '=' in t:
				s = t[1:t.index('=')]
				v = t[t.index('=')+1:]
				if not s in valOpt:
					error_usage()
				if types[s] == bool:
					error_usage()
				valOpt[s] = types[s](v)
				
			else:
				s = t[1:]
				if not s in valOpt:
					error_usage()
				if types[s] != bool:
					error_usage()
				# Ici, on affecte False
				valOpt[s] = (t[0] == '+')
		elif os.access(t, os.R_OK):
			valArg.append(t)
		else:
			print >> sys.stderr, "No access to", t
			error_usage()

	# Il n'y a pas le nombre d'arguments minimal
	if len(valArg) != len(args):
		error_usage()
	
	return (valArg, valOpt)


