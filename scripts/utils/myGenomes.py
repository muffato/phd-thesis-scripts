
import sys
import myTools


##############
# Un gene :) #
##############
class Gene:

	def __init__(self, names, chromosome, beg, end, strand):

		self.names = tuple(intern(s) for s in names)
		self.chromosome = commonChrName(chromosome)
		self.beginning = beg
		self.end = end
		self.strand = strand

	def __repr__(self):
		return "Gene %s on chr %s from %d to %d on strand %d" % ("/".join(self.names), self.chromosome, self.beginning, self.end, self.strand)


###########################################################
# On convertit les noms de chromosomes en int si possible #
###########################################################
def commonChrName(x):
	try:
		return int(x)
	except TypeError:
		return None
	except ValueError:
		return intern(x)


##############################
# Charge un alignement FASTA #
##############################
def loadFastaFile(name):
	seq = {}
	f = myTools.myOpenFile(name, "r")
	name = None
	for ligne in f:
		ligne = ligne.replace('\n', '').strip()
		if len(ligne) == 0:
			continue
		if ligne[0] == ">":
			name = ligne[1:].strip()
			seq[name] = ""
		elif name != None:
			seq[name] += ligne
	f.close()
	return seq



##########################################
# Classe generale pour stocker un genome #
##########################################
class Genome:

	#
	# Constructeur
	#
	def __init__(self, fichier, **args):

		# le nom et l'instance de file
		f = myTools.myOpenFile(fichier, 'r')
		self.nom = f.name
		self.f = myTools.firstLineBuffer(f)
		
		# la liste des genes par chromosome
		self.lstGenes = myTools.defaultdict(list)
		# Associer un nom de gene a sa position sur le genome
		self.dicGenes = {}
		# La liste triee des noms des chromosomes, des scaffolds et des randoms
		self.lstChr = []
		self.lstScaff = []
		self.lstRand = []

		# Le choix de la fonction de chargement
		c = self.f.firstLine.split()
		
		# Fichier de genome d'Ensembl: les trois champs du milieu sont des nombres
		try:
			x = int(c[1]) + int(c[2]) + int(c[3])
			self.__loadFromEnsemblGenome__(**args)
		except (ValueError, IndexError):
				
			# Fichier d'orthologues classique, 4 champs indiquent les scores de similarite
			try:
				x = int(c[7]) + int(c[8]) + int(c[9]) + int(c[10])
				self.__loadFromEnsemblOrthologs__(**args)
			except (ValueError, IndexError):
				
				# 1. Noms de chromosomes ?
				if "withChr" not in args:
					firstElt = ""
					try:
						firstElt = c[0]
						int(firstElt[0])
						withChr = True
					except (ValueError, IndexError):
						withChr = (len(firstElt) < 4)    # Les noms de genes font au minimum 4 caracteres
						# Les quelques exceptions
						withChr |= firstElt.upper().strip() in ["ALPHA", "BETA", "DELTA", "EPSILON", "GAMMA", "PHI"]
				else:
					withChr = args.pop("withChr")

				# 2. Facteur de qualite de concorde
				if "concordeQualityFactor" not in args:
					try:
						x = int(c[1])
						conc = True
					except (ValueError, IndexError):
						conc = False
				else:
					conc = args.pop("concordeQualityFactor")
					
				self.__loadFromAncestralGenome__(withChr=withChr, concordeQualityFactor=conc, **args)
		self.f.close()
		self.sortGenome()


	#
	# Rajoute un gene au genome
	#
	def addGene(self, gene):
	
		self.lstGenes[gene.chromosome].append(gene)
	
	#
	# Trie les chromsomoses
	#
	def sortGenome(self):
		
		import operator

		self.lstChr.sort()
		self.lstScaff.sort()
		self.lstRand.sort()
		
		self.dicGenes = {}
		for c in self.lstGenes:
			self.lstGenes[c].sort(key = operator.attrgetter('beginning'))
			for (i,g) in enumerate(self.lstGenes[c]):
				g.chromosome = c
				for s in g.names:
					self.dicGenes[s] = (c,i)
		
	#
	# Renvoie les genes presents sur le chromosome donne a certaines positions
	#	
	def getGenesAt(self, chr, beg = 0, end = sys.maxint):
		if chr not in self.lstGenes:
			return
		lst = self.lstGenes[chr]
		
		# Recherche dichotomique d'un gene dans la fenetre
		def dichotFind(a, b):
			i = (a+b)/2
			if a == i:
				return a
			if lst[i].end < beg:
				return dichotFind(i,b)
			elif lst[i].beginning > end:
				return dichotFind(a,i)
			else:
				return i
		index = dichotFind(0, len(lst)-1)
		
		for i in xrange(index, len(lst)):
			g = lst[i]
			if g.beginning > end:
				break
			yield g
		
		for i in xrange(index-1, -1, -1):
			g = lst[i]
			if g.end < beg:
				break
			yield g
		
	#
	# Renvoie les noms des genes presents aux alentours d'un gene donne (fenetre l en nombre de genes)
	#	
	def getGenesNear(self, chr, index, l):
		if chr not in self.lstGenes:
			return []
		
		return self.lstGenes[chr][max(0, index-l):index+l+1]

	#
	# Renvoie tous les genes
	#
	def __iter__(self):
		for t in self.lstGenes.itervalues():
			for g in t:
				yield g

	#
	# Cherche la position d'un gene donne par ses noms)
	#
	def getPosition(self, names):
		return set( (self.dicGenes[s] for s in names if s in self.dicGenes) )

	#
	# Renvoie les autres noms d'un gene
	#
	def getOtherNames(self, name):
		if name not in self.dicGenes:
			return []
		(c,i) = self.dicGenes[name]
		return [x for x in self.lstGenes[c][i].names if x != name]

	#
	# Fabrique la liste des orthologues entre les deux genomes
	#
	def buildOrthosTable(self, chr1, genome2, chr2, includeGaps, genesAnc = None):

		# Tous les orthologues entre pour les chromosomes OK du genome 1
		res = {}
		for c1 in chr1:
			res[c1] = []
			for (i1,g1) in enumerate(self.lstGenes[c1]):
				tmp = genome2.getPosition(g1.names)
				if genesAnc != None:
					for (c,i) in genesAnc.getPosition(g1.names):
						tmp.update(genome2.getPosition(genesAnc.lstGenes[c][i].names))
				
				# On ne garde que les chromsomes OK du genome 2
				tmp = [(c,i) for (c,i) in tmp if c in chr2]
				# +/- includeGaps
				if includeGaps or (len(tmp) > 0):
					res[c1].append( (i1,tmp) )
			res[c1].reverse()

		return res



	################################################################
	# Cette fonction gere un fichier de liste de genes d'Ensembl   #
	#   "Chr Debut Fin Brin ENSFAMxxxx Nom"                        #
	################################################################
	def __loadFromEnsemblGenome__(self):
		
		print >> sys.stderr, "Chargement du genome de", self.nom, "...",
		
		# On lit chaque ligne
		for ligne in self.f:
			champs = ligne.split()
			self.addGene( Gene(champs[4:], champs[0], int(champs[1]), int(champs[2]), int(champs[3])) )

		dicConvRomain = {"I":1, "II":2, "III":3, "IV":4, "V":5, "VI":6, "VII":7, "VIII":8, "IX":9, "X":10, "XI":11, "XII":12, "XIII":13, "XIV":14, "XV":15, "XVI":16, "XVII":17, "XVIII":18, "XIX":19, "XX":20, "XXI":21, "XXII":22, "XXIII":23, "XXIV":24, "XXV":25}

		# Les veritables chromosomes sont des entiers < 50 ou des entiers suivis de p/q/L/R/a/b ou W/X/Y/Z
		# Les chromosomes '*random*' ou 'UNKN' sont des scaffold mis bout a bout -> pas d'ordre utilisable
		# Le reste correspond aux scaffolds
		for c in self.lstGenes:
			if (type(c) == int) or (type(c) == long):
				if c < 50:
					self.lstChr.append(c)
				else:
					self.lstScaff.append(c)
			elif ((c[-1] in "pqLRab") and (c != "c6_QBL")) or ((c[0] in "WXYZ") and len(c) <= 2) or (c in dicConvRomain) or (c[:5] == "group"):
				self.lstChr.append(c)
			elif ("andom" in c) or (c == "UNKN"):
				self.lstRand.append(c)
			else:
				self.lstScaff.append(c)
		self.sortGenome()
		
		print >> sys.stderr, "OK"


	##################################################################################
	# Cette fonction gere un fichier d'orthologues a partir duquel on cree un genome #
	# Les noms des genes sont en positions 1 et 4                                    #
	##################################################################################
	def __loadFromEnsemblOrthologs__(self, homologyFilter=[], ancFilter=[]):
		
		print >> sys.stderr, "Chargement des orthologues de", self.nom, "...",
		
		combin = myTools.myCombinator([])
		fH = (len(homologyFilter) != 0)
		fA = (len(ancFilter) != 0)
		
		# On lit chaque ligne
		for ligne in self.f:
			champs = ligne.replace('\n', '').split('\t')
			if len(champs) == 12:
				# Nouvelle version des fichiers
				if (champs[-1] not in homologyFilter) and fH:
					continue
				if (champs[6] not in ancFilter) and fA:
					continue
			else:
				# Ancienne version des fichiers
				if (champs[6] not in homologyFilter) and fH:
					continue
			combin.addLink([champs[0], champs[3]])
		
		for (nb,g) in enumerate(combin):
			self.addGene( Gene(g, None, nb, nb, 0) )

		print >> sys.stderr, "OK"


	####################################################################################
	# Cette fonction gere un fichier de genome ancestral                               #
	# il s'agit d'une liste de lignes de la forme "CHROMOSOME GENE_ESP1 GENE_ESP2 ..." #
	####################################################################################
	def __loadFromAncestralGenome__(self, withChr=False, concordeQualityFactor=False):

		if withChr:
			print >> sys.stderr, "Chargement du genome ancestral de", self.nom, "...",
		else:
			print >> sys.stderr, "Chargement des genes ancestraux de", self.nom, "...",
		
		# On initialise tout
		strand = 0
		c = None
		i = 0
		for ligne in self.f:
			champs = ligne.split()

			if len(champs) == 0:
				continue

			# Le chromosome du gene lu
			if withChr:
				c = champs.pop(0)
			
			if concordeQualityFactor:
				# Fichiers de genomes ancestraux avec score de Concorde
				strand = int(champs.pop(0))
				
			# On ajoute le gene
			self.addGene( Gene(champs, c, i, i, strand) )
			i += 1
		
		self.lstChr = self.lstGenes.keys()
		
		print >> sys.stderr, "OK"


