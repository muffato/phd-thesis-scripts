#! /users/ldog/muffato/python -OO


import os
import sys
import bz2
import gzip


null = open('/dev/null', 'w')
stdin = sys.stdin
stdout = sys.stdout
stderr = sys.stderr

###########################################################################
# Consomme les elements d'un iterateur et renvoie la longueur de la liste #
###########################################################################
def leniter(it):
	nb = 0
	for _ in it:
		nb += 1
	return nb


##################################################################
# Cette classe ouvre le fichier en le decompressant s'il le faut #
#   Retourne l'objet FILE et le nom complet du fichier           #
##################################################################
def myOpenFile(nom, mode):
	nom = nom.replace("~", os.environ['HOME'])
	if nom.endswith(".bz2"):
		f = bz2.BZ2File(nom, mode)
	elif nom.endswith(".gz"):
		f = gzip.GzipFile(nom, mode)
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
			
				# Les valeurs a eviter pour chaque mode
				if self.mode == myMatrixIterator.OnlyDiag:
					if i != j:
						continue
				elif self.mode == myMatrixIterator.WholeWithoutDiag:
					if i == j:
						continue
				elif self.mode == myMatrixIterator.WholeMatrix:
					pass
				elif self.mode == myMatrixIterator.UpperMatrix:
					if j < i:
						continue
				elif self.mode == myMatrixIterator.StrictUpperMatrix:
					if j <= i:
						continue
				
				# Mode inconnu
				else:
					continue
				yield (i,j)
		

########################################################################
# Cette classe permet de regrouper une liste d'elements                #
# Partant d'une liste initiale, on ajoute des liens entre des elements #
#  de la liste et la classe les regroupe                               #
########################################################################
class myCombinator:

	#
	# Constructeur
	#
	def __init__(self, ini = []):
		self.grp = list(ini)
		self.dic = {}
		for i in xrange(len(self.grp)):
			for x in self.grp[i]:
				self.dic[x] = i
	
	#
	# Definit un lien entre tous les elements de obj
	#
	def addLink(self, obj):
	
		if len(obj) == 0:
			return
	
		d = set([self.dic[x] for x in obj if x in self.dic])
		
		if len(d) == 0:
			i = len(self.grp)
			self.grp.append(obj)
			for x in obj:
				self.dic[x] = i
		else:
			i = d.pop()
			for x in d:
				self.grp[i].extend(self.grp[x])
				for y in self.grp[x]:
					self.dic[y] = i
				self.grp[x] = []
			dd = [x for x in obj if x not in self.dic]
			for x in dd:
				self.dic[x] = i
			self.grp[i].extend(dd)
	

	#
	# Renvoie un iterateur sur les donnees
	#  Les ensembles vides sont donc elimines
	#
	def __iter__(self):
		for g in self.grp:
			if len(g) > 0:
				yield g
	#
	# Le nombre de groupes
	#
	def getNbGrp(self):
		return leniter(self)


	#
	# Enleve les ensembles vides
	#
	def reduce(self):
		newGrp = list(self)
		self.__init__(newGrp)

				

#######################################################################################
# Classe qui permet d'utiliser les cles d'un dictionnaire directement comme attributs #
#######################################################################################
class dicContainer(dict):

	def __getattr__(self, name):
		return dict.__getitem__(self, name)


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

	valOpt = dicContainer()
	valArg = []
	opt = {}
	for (name,typ,val) in options:
		opt[name] = (typ,val)
		if type(val) == list:
			valOpt[name] = val[0]
		else:
			valOpt[name] = val
	
	# On scanne les arguments pour les compter et recuperer les valeurs
	for tt in sys.argv[1:]:

		t = tt.replace('^', ' ')

		# Un petit peu d'aide
		if t == '-h' or t == '--help' or t == '-help':
			error_usage()
			
		# Un argument optionnel
		if t[0] in '-+':
		
			# Un parametre non bool
			try:
				i = t.index('=')
				s = t[1:i]
				v = t[i+1:]
				# Parametre non connu
				if not s in valOpt:
					error_usage()
				# Mauvaise syntaxe pour les bool
				if opt[s][0] == bool:
					error_usage()
				# Valeur de parametre non autorisee
				valOpt[s] = opt[s][0](v)
				if (type(opt[s][1]) == list) and (valOpt[s] not in opt[s][1]):
					error_usage()
				
			except ValueError:
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


