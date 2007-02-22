#! /users/ldog/muffato/python -OO


import sys
import os
from bz2 import BZ2File
from gzip import GzipFile


null = open('/dev/null', 'w')
stdin = sys.stdin
stdout = sys.stdout
stderr = sys.stderr

##################################################################
# Cette classe ouvre le fichier en le decompressant s'il le faut #
#   Retourne l'objet FILE et le nom complet du fichier           #
##################################################################
def myOpenFile(nom, mode):
	if "~" in nom:
		if "_CONDOR_SCRATCH_DIR" in os.environ:
			nom = nom.replace("~", "/users/ldog/muffato")
		else:
			nom = nom.replace("~", os.environ['HOME'])
	if nom.endswith(".bz2"):
		f = BZ2File(nom, mode)
	elif nom.endswith(".gz"):
		f = GzipFile(nom, mode)
	else:
		f = open(nom, mode)
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
	
		for i in xrange(0, self.n):
			for j in xrange(0, self.p):
			
				if self.mode == myMatrixIterator.OnlyDiag:
					if i != j:
						continue
				elif self.mode == myMatrixIterator.WholeWithoutDiag:
					if i == j:
						continue
				elif self.mode == myMatrixIterator.WholeMatrix:
					continue
				elif self.mode == myMatrixIterator.UpperMatrix:
					if j < i:
						continue
				elif self.mode == myMatrixIterator.StrictUpperMatrix:
					if j <= i:
						continue
				else:
					continue
				yield (i,j)
		

########################################################################
# Cette classe permet de regrouper une liste d'elements                #
# Partant d'une liste initiale, on ajoute des liens entre des elements #
#  de la liste et la classe les regroupe                               #
########################################################################
class myCombinator:

	def __init__(self, ini):
		self.grp = ini
		self.dic = {}
		for i in xrange(len(ini)):
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
			self.grp.append([a])
			self.dic[a] = i

		for x in obj[1:]:
			if x in self.dic:
				j = self.dic[x]
				if i == j:
					continue
				for b in self.grp[j]:
					self.dic[b] = i
				self.grp[i].extend(self.grp[j])
				self.grp[j] = []
			else:
				self.grp[i].append(x)
				self.dic[x] = i

	def getGrp(self):
		return [g for g in self.grp if len(g) > 0]

	def getNbGrp(self):
		nb = 0
		for g in self.grp:
			if len(g) > 0:
				nb += 1
		return nb

	def __iter__(self):
		for g in self.grp:
			if len(g) > 0:
				yield g



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
			print >> sys.stderr, "\t", invite + "%s %s (%s)" % t
		if info != "":
			print >> sys.stderr, "\n", info
		sys.exit(1)

	#types = {}
	valOpt = {}
	valArg = []
	opt = {}
	for (name,typ,val) in options:
		opt[name] = (typ,val)
		if type(val) == list:
			valOpt[name] = val[0]
		else:
			valOpt[name] = val
	
	# On scanne les argumetns pour les compter et recuperer les valeurs
	for tt in sys.argv[1:]:

		t = tt.replace('^', ' ')

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
				if (type(opt[s][1]) == list) and (v not in opt[s][1]):
					error_usage()
				if opt[s][0] == bool:
					error_usage()
				valOpt[s] = opt[s][0](v)
				
			else:
				s = t[1:]
				if not s in valOpt:
					error_usage()
				if opt[s][0] != bool:
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
	
	return (dict([(args[i],valArg[i]) for i in xrange(len(args))]), valOpt)


