
import sys
import bz2
import gzip
import itertools

import myFile

null = open('/dev/null', 'w')

def applyFunctions(fun, data):
	for (f,x) in itertools.izip(fun, data):
		yield f(x)

def funcFilter(fun):
	return lambda data: (f(x) for (f,x) in itertools.izip(fun, data))


def deprecated(func):
	"""This is a decorator which can be used to mark functions
	as deprecated. It will result in a warning being emmitted
	when the function is used."""
	import warnings
	def newFunc(*args, **kwargs):
		warnings.warn("Call to deprecated function %s." % func.__name__, category=DeprecationWarning, stacklevel=2)
		return func(*args, **kwargs)
	newFunc.__name__ = func.__name__
	newFunc.__doc__ = func.__doc__
	newFunc.__dict__.update(func.__dict__)
	return newFunc


#####################################################################
# Une classe pour avoir un iterateur a deux positions sur une liste #
#####################################################################
class myIterator:
	
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
	
	# Couple (x,y) glissant
	########################
	@staticmethod
	def slidingTuple(lst):
		if len(lst) > 0:
			x = lst[0]
			for i in xrange(1, len(lst)):
				y = lst[i]
				yield (x,y)
				x = y
	


##########################################################
# Gestion du lancement multiple sur une plage de valeurs #
##########################################################
def getRange(s):
	if myFile.hasAccess(s):
		f = myFile.openFile(s, "r")
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

########################
# Classe list hashable #
########################
class hashablelist(list):
	def __hash__(self):
		return hash(tuple(self))

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
		#self.nbcall += 1
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

def memoizeC(cache={}):
	def dowrap(f):
		def newf(*args):
			try:
				return cache[args]
			except KeyError:
				cache[args] = value = f(*args)
				return value
			except TypeError:
				# uncachable -- for instance, passing a list as an argument.
				# Better to not cache than to blow up entirely.
				return f(*args)
		return newf
	return dowrap


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

		
#####################################################
# Rajoute des options pour un module en particulier #
#####################################################
__moduleoptions = []
def addModuleOptions(namespace, options):
	for (name,typ,val) in options:
		__moduleoptions.append( (namespace+":"+name,typ,val) )


########################################################
# Permet de demander une liste de fichiers en argument #
########################################################
class FileList:
	def __init__(self, value):
		self.minNbFiles = value

	def __repr__(self):
		return '<FileList(%d)>' % self.minNbFiles


#################################################################################
# Lit la ligne de commande et parse les arguments                               #
#  1. des arguments obligatoires (nom,constructeur)                             #
#  2. des options sous la forme -opt=val (nom, constructeur, val_defaut)        #
# En cas d'erreur, affiche la syntaxe demandee et une courte description (info) #
#################################################################################
def checkArgs(args, options, info, showArgs=True):

	options = options + __moduleoptions
	#
	# Affiche le message d'erreur de mauvais arguments
	#
	def error_usage(reason):
		print >> sys.stderr, "- ERREUR -", reason
		print >> sys.stderr, " Usage :", sys.argv[0]
		for (i,t) in enumerate(args):
			print >> sys.stderr, "\t", "%d:" % (i+1), t[0], t[1]
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
			if not myFile.hasAccess(v):
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
		
	valOpt = {}
	valArg = {}
	opt = {}
	for (name,typ,val) in options:
		opt[name] = (typ,val)
		valOpt[name] = val[0] if type(val) == list else val
	
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
								import psyco
								from psyco.classes import __metaclass__
								psyco.full()
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
				if isinstance(typ, FileList):
					valArg[s] = list()
					assert len(valArg) == len(args)
					valArg[s].append(putValue(file, None, t))
				else:
					valArg[s] = putValue(typ, None, t)
			elif isinstance(args[-1][1], FileList):
				valArg[args[-1][0]].append(putValue(file, None, t))
			else:
				error_usage("Trop d'arguments sur '%s'" % t)

	if isinstance(args[-1][1], FileList):
		if args[-1][0] not in valArg:
			valArg[args[-1][0]] = []
		if len(valArg[args[-1][0]]) < args[-1][1].minNbFiles:
			error_usage("Pas assez de fichiers pour '%s'" % args[-1][0])

	# Il n'y a pas le nombre d'arguments minimal
	if len(valArg) < len(args):
		error_usage("Pas assez d'arguments")
	
	valArg.update(valOpt)
	if showArgs:
		print >> sys.stderr, "Arguments:", valArg
	return valArg


