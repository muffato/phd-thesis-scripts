
import os
import sys
import bz2
import gzip
import operator
import collections

defaultdict = collections.defaultdict

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
		nom = nom.replace("~", os.environ['HOME'])
		if nom.endswith(".bz2"):
			f = bz2.BZ2File(nom, mode)
		elif nom.endswith(".gz"):
			f = gzip.GzipFile(nom, mode)
		else:
			f = open(nom, mode)
	return f


#####################################################
# Cree le repertoire pour les sorties dans fichiers #
#####################################################
def mkDir(dir):
	try:
		os.makedirs(os.path.dirname(dir))
	except OSError:
		pass

########################################
# Permet de charger les dumps de MySQL #
#    Raboute les lignes tronquees      #
########################################
def MySQLFileLoader(f):
	tmp = ""
	for ligne in f:
		if ligne[-2] == '\\':
			# Signe que la ligne n'est pas terminee
			tmp = ligne[:-2]
		else:
			yield tmp + ligne[:-1]
			tmp = ""
	# Normalement, ici, tmp == ""
	if tmp != "":
		print >> sys.stderr, "File consistency error !"


#####################################################################
# Une classe pour avoir un iterateur a deux positions sur une liste #
#####################################################################
class myIterator:
	
	@staticmethod
	def tupleOnWholeList(lst):
		for x in lst:
			for y in lst:
				yield (x,y)

	@staticmethod
	def tupleOnUpperList(lst):
		for i in xrange(len(lst)):
			x = lst[i]
			for y in lst[i:]:
				yield (x,y)

	@staticmethod
	def tupleOnStrictUpperList(lst):
		for i in xrange(len(lst)):
			x = lst[i]
			for y in lst[i+1:]:
				yield (x,y)
	
	@staticmethod
	def tupleOnTwoLists(lstX, lstY):
		for x in lstX:
			for y in lstY:
				yield (x,y)

	# Fonction de Charles pour renvoyer toutes les combinaisons des elements des differentes listes passees en argument
	@staticmethod
	def tupleOnManyLists(*args):
		""" This generator combine all versus all sequences elements as follow:
		>>> args = [['A','C'],['A','C'],['A','C']]
		>>> [k for k in combination(args)]
		['AAA', 'AAC', 'ACA', 'ACC', 'CAA', 'CAC', 'CCA', 'CCC']
		"""
		lengths = [len(seq) for seq in args]
		_tmp = lengths + [1] # append multiplicative identity
		range_len_args = range(len(args))
		dividers = [reduce(operator.mul, _tmp[-x-1:]) for x in range_len_args][::-1]
		for n in xrange(reduce(operator.mul, lengths)):
			yield tuple( args[r][(n/dividers[r])%lengths[r]] for r in range_len_args )
	
	@staticmethod
	def slidingTuple(lst):
		x = lst[0]
		for i in xrange(1, len(lst)):
			y = lst[i]
			yield (x,y)
			x = y
	
	# Consomme les elements d'un iterateur et renvoie la longueur de la liste
	@staticmethod
	def leniter(it):
		nb = 0
		for _ in it:
			nb += 1
		return nb

	
	@staticmethod
	def buildSubsets(lst, n):
		l = len(lst)
		mem = {}

		def rec(i, n):
			if (i,n) in mem:
				return mem[(i,n)]
			
			if i >= l-n:
				res = [lst[i:]]
			elif n <= 1:
				res = [[x] for x in lst[i:]]
			else:
				ref = [lst[i]]
				res = rec(i+1, n)
				for x in rec(i+1, n-1):
					res.append(ref + x)

			mem[(i,n)] = res
			return res

		return rec(0, n)


class memoize:
	"""Decorator that caches a function's return value each time it is called.
	If called later with the same arguments, the cached value is returned, and
	not re-evaluated.
	"""
	def __init__(self, func):
		self.func = func
		self.cache = {}

	def __call__(self, *args):
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

	#
	# Constructeur
	#
	def __init__(self, ini = []):
		self.grp = list(ini)
		self.dic = {}
		for i in xrange(len(self.grp)):
			self.grp[i] = list(set(self.grp[i]))
			for x in self.grp[i]:
				self.dic[x] = i
	
	#
	# Definit un lien entre tous les elements de obj
	#
	def addLink(self, obj):
	
		if len(obj) == 0:
			return
	
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
		return myIterator.leniter(self)

	#
	# Enleve les ensembles vides
	#
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

		# Un argument optionnel
		if t[0] in '-+':
		
			# Un parametre non bool
			try:
				i = t.index('=')
				s = t[1:i]
				v = t[i+1:]
				# Parametre non connu
				if not s in valOpt:
					print >> sys.stderr, "Option non reconnue:", s
					error_usage()
				# Mauvaise syntaxe pour les bool
				if opt[s][0] == bool:
					print >> sys.stderr, "Utiliser +%s ou -%s, sans '='" % (s,s)
					error_usage()
				# Valeur de parametre non autorisee
				valOpt[s] = opt[s][0](v)
				if (type(opt[s][1]) == list) and (valOpt[s] not in opt[s][1]):
					print >> sys.stderr, "Les valeurs autorisees pour '%s' sont %s (ligne de commande: '%s')" \
						% (s,printLine(opt[s][1], '/'),valOpt[s])
					error_usage()
				
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
							#class bz2Proxy():
							#	def __init__(self):
							#		self.comp = bz2.BZ2Compressor()
							#		self.stdout = sys.stdout
							#	def __del__(self):
							#		self.stdout.write(self.comp.flush())
							#		self.stdout.flush()
							#	def write(self, data):
							#		self.stdout.write(self.comp.compress(data))
							#sys.stdout = bz2Proxy()
							sys.stdout = bz2.BZ2File("/dev/stdout", "w")
					elif s == "gz":
						if t[0] == '+':
							sys.stdout = gzip.GzipFile("/dev/stdout", "w")
					else:
						print >> sys.stderr, "Option non reconnue:", s
						error_usage()
				elif opt[s][0] != bool:
					print >> sys.stderr, "Utiliser -%s=valeur" % s
					error_usage()
				else:
					# Ici, on affecte False
					valOpt[s] = (t[0] == '+')
		elif os.access(t, os.R_OK):
			valArg.append(t)
		else:
			print >> sys.stderr, "Fichier %s non accessible" % t
			error_usage()

	# Il n'y a pas le nombre d'arguments minimal
	if len(valArg) != len(args):
		error_usage()
	
	return (dict([(arg,valArg[i]) for (i,arg) in enumerate(args)]), valOpt)


