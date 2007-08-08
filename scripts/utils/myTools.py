
import os
import sys
import bz2
import gzip
import operator
from collections import defaultdict

null = open('/dev/null', 'w')
stdin = sys.stdin
stdout = sys.stdout
stderr = sys.stderr
stdinInput = os.isatty(sys.stdin.fileno())

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
	if nom.startswith("http://") or nom.startswith("ftp://"):
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

#####################################################################
# Une classe pour avoir un iterateur a deux positions sur une liste #
#####################################################################
class myIterator:
	
	def _tupleOnWholeList(lst):
		for x in lst:
			for y in lst:
				yield (x,y)

	def _tupleOnUpperList(lst):
		for (i,x) in enumerate(lst):
			for y in lst[i:]:
				yield (x,y)

	def _tupleOnStrictUpperList(lst):
		for (i,x) in enumerate(lst):
			for y in lst[i+1:]:
				yield (x,y)
	
	def _tupleOnTwoLists(lstX, lstY):
		for x in lstX:
			for y in lstY:
				yield (x,y)
	#
	# Fonction de Charles pour renvoyer toutes les combinaisons des elements des differentes listes passees en argument
	#
	def _tupleOnManyLists(*args):
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
	
	
	
	def _buildSubsets(lst, n):
		
		# Cas special
		if n < 2:
			return [set([x]) for x in lst]
		
		ens = buildSubsets(lst, n-1)
		res = []
		for s in ens:
			m = max(s)
			for x in lst:
				if x <= m:
					continue
				ss = s.union([x])
				if len(ss) == n:
					res.append(ss)
		return res
	
	def _buildSubsets2(lst, n):
		if type(lst) != set:
			lst = frozenset(lst)
		
		# Cas special
		if n < 2:
			return [set([x]) for x in lst]
			#for x in lst:
			#	yield set([x])
		
		ens = _buildSubsets2(lst, n-1)
		res = []
		for s in ens:
			m = max(s)
			for x in lst:
				if x <= m:
					continue
				ss = s.union([x])
				if len(ss) == n:
					res.append(ss)
		return res

	buildSubsets = staticmethod(_buildSubsets)
	buildSubsets2 = staticmethod(_buildSubsets2)


	tupleOnWholeList = staticmethod(_tupleOnWholeList)
	tupleOnUpperList = staticmethod(_tupleOnUpperList)
	tupleOnStrictUpperList = staticmethod(_tupleOnStrictUpperList)
	tupleOnTwoLists = staticmethod(_tupleOnTwoLists)
	tupleOnManyLists = staticmethod(_tupleOnManyLists)



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
	
		obj = list(set(obj))
		# Les elements de obj deja presents dans le combinateur
		d = set( (self.dic[x] for x in obj if x in self.dic) )
		
		if len(d) == 0:
			# Aucun, on rajoute tel quel l'objet alors
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
						% (s,'/'.join([str(x) for x in opt[s][1]]),valOpt[s])
					error_usage()
				
			# Si on ne trouve pas de '=', c'est une type bool
			except ValueError:
				s = t[1:]
				if s == "psyco":
					try:
						import utils.psyco
						from utils.psyco.classes import __metaclass__
						#utils.psyco.log()
						utils.psyco.full()
						#utils.psyco.full(memory=100)
						#utils.psyco.profile(0.1, memory=100)
					except ImportError:
						print >> sys.stderr, "Unable to load psyco !"
				elif not s in valOpt:
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


