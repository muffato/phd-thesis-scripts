
import os
import sys
import bz2
import gzip

from collections import defaultdict

null = open('/dev/null', 'w')
stdin = sys.stdin
stdout = sys.stdout
stderr = sys.stderr
stdinInput = os.isatty(sys.stdin.fileno())

#########################################################################################################################
# Le but est de pouvoir acceder au fichier et lire la premiere ligne sans devoir le fermer pour le reouvrir juste apres #
#########################################################################################################################
class firstLineBuffer:
	def __init__(self, f):
		self.f = f
		try:
			self.firstLine = self.f.next()
		except StopIteration:
			self.firstLine = ""

	def __iter__(self):
		yield self.firstLine
		for l in self.f:
			if not l.startswith("#"):
				yield l

	def close(self):
		return self.f.close()


#############################################################################################
# Lit un fichier tabulaire, en convertissant les colonnes separees de delim selon type_list #
#############################################################################################
def readTabular(filename, type_list, delim = '\t'):
	f = myOpenFile(filename, 'r')
	for (i,line) in enumerate(f):
		current_line = line.replace('\n','').split(delim)
		if len(current_line) != len(type_list):
			raise IndexError,"Erreur nombre de colonne. Ligne:%d" % (i+1)
		yield [t(x) for (x,t) in zip(current_line,type_list)]
	f.close()


####################
# Fichier existant #
####################
def fileAccess(s):
	return os.access(s, os.R_OK)


###############################################
# Renvoie la ligne preparee pour l'impression "
###############################################
def printLine(line, delim = "\t", func = str):
	return delim.join([func(x) for x in line])


####################################################################
# Cette fonction ouvre le fichier en le decompressant s'il le faut #
#   Retourne l'objet FILE et le nom complet du fichier             #
####################################################################
def myOpenFile(nom, mode):
	if type(nom) != str:
		return nom
	elif nom.startswith("http://") or nom.startswith("ftp://"):
		comm = "wget %s -O -"
		if nom.endswith(".bz2"):
			comm += " | bunzip2"
		elif nom.endswith(".gz"):
			comm += " | gunzip"
		(stdin,f,stderr) = os.popen3( comm % nom )
		stdin.close()
		stderr.close()
	else:
		if "w" in mode:
			mkDir(nom)
		nom = nom.replace("~", os.environ['HOME'])
		i = nom.find(".zip/")
		if (mode == "r") and (i >= 0):
			import zipfile
			import cStringIO
			f = zipfile.ZipFile(nom[:i+4], "r")
			f = cStringIO.StringIO(f.read(nom[i+5:]))
		elif nom.endswith(".bz2"):
			f = bz2.BZ2File(nom, mode)
		elif nom.endswith(".gz"):
			f = gzip.GzipFile(nom, mode)
		else:
			f = open(nom, mode)
	return f


#####################################################
# Cree le repertoire pour les sorties dans fichiers #
#####################################################
def mkDir(fichier):
	try:
		os.makedirs(os.path.dirname(fichier))
	except OSError:
		pass

########################################
# Permet de charger les dumps de MySQL #
#    Raboute les lignes tronquees      #
########################################
def MySQLFileLoader(f):
	tmp = ""
	for ligne in f:
		ligne = ligne.replace('\n', '')
		if ligne[-1] == '\\':
			# Signe que la ligne n'est pas terminee
			tmp = ligne[:-1]
		else:
			yield tmp + ligne
			tmp = ""
	# Normalement, ici, tmp == ""
	if tmp != "":
		print >> sys.stderr, "File consistency error !"


#####################################################################
# Une classe pour avoir un iterateur a deux positions sur une liste #
#####################################################################
class myIterator:
	
	# Parcourt tous les (x,y)
	##########################
	@staticmethod
	def tupleOnWholeList(lst):
		for x in lst:
			for y in lst:
				yield (x,y)

	# Parcourt les (xi,yj) pour j>=i
	#################################
	@staticmethod
	def tupleOnUpperList(lst):
		for i in xrange(len(lst)):
			x = lst[i]
			for y in lst[i:]:
				yield (x,y)

	# Parcourt les (xi,yj) pour j>i
	################################
	@staticmethod
	def tupleOnStrictUpperList(lst):
		for i in xrange(len(lst)):
			x = lst[i]
			for y in lst[i+1:]:
				yield (x,y)
	
	# Parcourt tous les (x,y)
	##########################
	@staticmethod
	def tupleOnTwoLists(lstX, lstY):
		for x in lstX:
			for y in lstY:
				yield (x,y)

	# Fonction de Charles pour renvoyer toutes les combinaisons des elements des differentes listes passees en argument
	####################################################################################################################
	@staticmethod
	def tupleOnManyLists(*args):
		""" This generator combine all versus all sequences elements as follow:
		>>> args = [['A','C'],['A','C'],['A','C']]
		>>> [k for k in combination(args)]
		['AAA', 'AAC', 'ACA', 'ACC', 'CAA', 'CAC', 'CCA', 'CCC']
		"""
		import operator
		lengths = [len(seq) for seq in args]
		_tmp = lengths + [1] # append multiplicative identity
		range_len_args = range(len(args))
		dividers = [reduce(operator.mul, _tmp[-x-1:]) for x in range_len_args][::-1]
		for n in xrange(reduce(operator.mul, lengths)):
			yield tuple( args[r][(n/dividers[r])%lengths[r]] for r in range_len_args )
	
	# Couple (x,y) glissant
	########################
	@staticmethod
	def slidingTuple(lst):
		x = lst[0]
		for i in xrange(1, len(lst)):
			y = lst[i]
			yield (x,y)
			x = y
	
	# Consomme les elements d'un iterateur et renvoie la longueur de la liste
	##########################################################################
	@staticmethod
	def leniter(it):
		nb = 0
		for _ in it:
			nb += 1
		return nb

	
	# Construit les sous-ensembles de taille n de la liste
	#######################################################
	@staticmethod
	def buildSubsets(lst, s):
		l = len(lst)

		# Renvoie les sous-ensembles de taille n, avec des elements d'indices >= i
		@memoize
		def rec(i, n):
			if i >= l-n:
				return [lst[i:]]
			elif n <= 1:
				return [[x] for x in lst[i:]]
			else:
				ref = [lst[i]]
				res = list(rec(i+1, n))
				for x in rec(i+1, n-1):
					res.append(ref + x)
				return res

		if s <= 0:
			res = [lst]
			for i in xrange(1, l):
				res.extend(rec(0, i))
			return res
		else:
			return rec(0, s)


##########################################################
# Gestion du lancement multiple sur une plage de valeurs #
##########################################################
def getRange(s):
	if fileAccess(s):
		f = myOpenFile(s, "r")
		lst = []
		for l in f:
			lst.extend( [int(x) for x in l.replace('\n', '').split()] )
		f.close()
		return lst
	else:
		(start,_,end) = s.partition(':')
		return range(int(start), int(end)+1)


###########################################################
# Classe dict hashable, utile pour s'en servir comme clef #
###########################################################
class hashabledict(dict):
	def __hash__(self):
		return hash(tuple(sorted(self.items())))


################################################################################
# Enregistre les resultats d'une fonction pour chaque valeur de ses parametres #
################################################################################
class memoize:
	"""Decorator that caches a function's return value each time it is called.
	If called later with the same arguments, the cached value is returned, and
	not re-evaluated.
	"""
	def __init__(self, func):
		self.func = func
		self.nbcall = 0
		self.cache = {}

	def __repr__(self):
		return "[%s: %d values cached, %d calls]" % (self.func.__name__, len(self.cache), self.nbcall)

	def __call__(self, *args):
		self.nbcall += 1
		try:
			return self.cache[args]
		except KeyError:
			self.cache[args] = value = self.func(*args)
			return value
		except TypeError:
			# uncachable -- for instance, passing a list as an argument.
			# Better to not cache than to blow up entirely.
			return self.func(*args)
	
	def __repr__(self):
		"""Return the function's docstring."""
		return self.func.__doc__



########################################################################
# Cette classe permet de regrouper une liste d'elements                #
# Partant d'une liste initiale, on ajoute des liens entre des elements #
#  de la liste et la classe les regroupe                               #
########################################################################
class myCombinator:

	def __init__(self, ini = []):
		self.grp = list(ini)
		self.dic = {}
		for i in xrange(len(self.grp)):
			self.grp[i] = list(set(self.grp[i]))
			for x in self.grp[i]:
				self.dic[x] = i
	
	# Definit un lien entre tous les elements de obj
	# Met a jour les ensembles deja construits
	#################################################
	def addLink(self, obj):
	
		if len(obj) == 0:
			return []
	
		obj = set(obj)
		grp = self.grp
		dic = self.dic

		# Les elements de obj deja presents dans le combinateur
		d = set( dic[x] for x in obj if x in dic )
		
		if len(d) == 0:
			# Aucun, on rajoute tel quel l'objet alors
			i = len(grp)
			grp.append(list(obj))
			for x in obj:
				dic[x] = i
			return grp
		else:
			i = d.pop()
			grpiextend = grp[i].extend
			for x in d:
				grpx = grp[x]
				grpiextend(grpx)
				for y in grpx:
					dic[y] = i
				grp[x] = []
			dd = [x for x in obj if x not in dic]
			for x in dd:
				dic[x] = i
			grpiextend(dd)
			return grp[i]
	

	# Renvoie un iterateur sur les donnees
	# Les ensembles vides sont donc elimines
	#########################################
	def __iter__(self):
		for g in self.grp:
			if len(g) > 0:
				yield g

	# Enleve les ensembles vides
	# ###########################
	def reduce(self):
		self.__init__(self)

		

#######################################################################################
# Classe qui permet d'utiliser les cles d'un dictionnaire directement comme attributs #
#######################################################################################
class dicContainer(dict):

	def __getattr__(self, name):
		return dict.__getitem__(self, name)


#################################################################################
# Lit la ligne de commande et parse les arguments                               #
#  1. des arguments obligatoires (nom,constructeur)                             #
#  2. des options sous la forme -opt=val (nom, constructeur, val_defaut)        #
# En cas d'erreur, affiche la syntaxe demandee et une courte description (info) #
#################################################################################
def checkArgs(args, options, info):

	#
	# Affiche le message d'erreur de mauvais arguments
	#
	def error_usage(reason):
		print >> sys.stderr, "- ERREUR -", reason
		print >> sys.stderr, " Usage :", sys.argv[0]
		for (i,t) in enumerate(args):
			print >> sys.stderr, "\t", "%d:" % (i+1), "%s %s" % t
		for t in options:
			if t[1] == bool:
				invite = "+/-"
				print >> sys.stderr, "\t", "+/-%s (%s)" % (t[0],t[2])
			else:
				print >> sys.stderr, "\t", "  -%s %s (%s)" % t
		if info != "":
			print >> sys.stderr, "\n", info
		sys.exit(1)

	def putValue(typ, val, v):

		# Creation de la valeur en fonction du type
		if typ == bool:
			# Type booleen
			res = {"false": False, "true":True}[v.lower()]
		elif typ == file:
			# Type 'fichier': test de presence
			if not fileAccess(v):
				error_usage("Fichier '%s' non accessible" % v)
			else:
				res = v
		else:
			# Sinon, on utilise le constructeur
			res = typ(v)
		
		if (type(val) == list) and (res not in val):
			# Valeur de parametre non autorisee
			error_usage("'%s' n'est pas parmi %s" % (res,printLine(val, '/')))
		return res
		
	valOpt = dicContainer()
	valArg = dicContainer()
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

		# Un argument optionnel
		if t[0] in '-+':
		
			# Un parametre non bool
			try:
				i = t.index('=')
				s = t[1:i]
				
				# Le nom du parametre doit etre connnu
				if not s in valOpt:
					error_usage("Option '%s' non reconnue" % s)

				valOpt[s] = putValue(opt[s][0], opt[s][1], t[i+1:])
			
			# Si on ne trouve pas de '=', c'est une type bool
			except ValueError:
				s = t[1:]
				# Nom de parametre non attendu
				if s not in valOpt:

					# Valeurs predefinies
					if s.startswith("psyco"):
						if t[0] == '+':
							try:
								import utils.psyco
								from utils.psyco.classes import __metaclass__
								if s.startswith("psycoL"):
									utils.psyco.log()
								else:
									utils.psyco.full()
							except ImportError:
								print >> sys.stderr, "Unable to load psyco !"
					elif s == "bz2":
						if t[0] == '+':
							sys.stdout = bz2.BZ2File("/dev/stdout", "w")
					elif s == "gz":
						if t[0] == '+':
							sys.stdout = gzip.GzipFile("/dev/stdout", "w")
					else:
						error_usage("Option '%s' non reconnue" % s)
				elif opt[s][0] != bool:
					error_usage("Utiliser -%s=valeur" % s)
				else:
					# Ici, on affecte False
					valOpt[s] = (t[0] == '+')
		else:
			if len(valArg) < len(args):
				(s,typ) = args[len(valArg)]
				valArg[s] = putValue(typ, None, t)
			else:
				error_usage("Trop d'arguments sur '%s'" % t)

	# Il n'y a pas le nombre d'arguments minimal
	if len(valArg) < len(args):
		error_usage("Pas assez d'arguments")
	
	valArg.update(valOpt)
	return valArg


