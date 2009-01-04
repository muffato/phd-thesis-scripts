
import sys
import collections

import myFile
import myTools


##############
# Un gene :) #
##############

Gene = collections.namedtuple("gene", ['chromosome', 'beginning', 'end', 'strand', 'names'])
GenePosition = collections.namedtuple("geneposition", ['chromosome', 'index'])


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
	f = myFile.openFile(name, "r")
	name = None
	for ligne in f:
		ligne = ligne.replace('\n', '').strip()
		if len(ligne) == 0:
			continue
		# Les chevrons indiquent le debut d'une nouvelle sequence
		if ligne[0] == ">":
			name = ligne[1:].strip()
			seq[name] = ""
		# Les lignes doivent etre concatenees
		elif name != None:
			seq[name] += ligne
	f.close()
	return seq



##########################################
# Classe generale pour stocker un genome #
##########################################
class Genome:

	# Constructeur
	###############
	def __init__(self, fichier, **args):

		# le nom et l'instance de file
		f = myFile.openFile(fichier, 'r')
		self.nom = f.name
		self.f = myFile.firstLineBuffer(f)
		
		# la liste des genes par chromosome
		self.lstGenes = collections.defaultdict(list)
		# Associer un nom de gene a sa position sur le genome
		self.dicGenes = {}

		# Le choix de la fonction de chargement
		c = self.f.firstLine.split()
		
		# Fichier de genome d'Ensembl: les trois champs du milieu sont des nombres
		try:
			x = int(c[1]) + int(c[2]) + int(c[3])
			mode = 1
		except (ValueError, IndexError):
				
			# Fichier d'orthologues classique, 4 champs indiquent les scores de similarite
			try:
				x = int(c[7]) + int(c[8]) + int(c[9]) + int(c[10])
				mode = 2
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
				mode = 3	
		if mode == 1:
			self.__loadFromEnsemblGenome__(**args)
		elif mode == 2:
			self.__loadFromEnsemblOrthologs__(**args)
		elif mode == 3:
			self.__loadFromAncestralGenome__(withChr, conc, **args)
		self.f.close()

		# Classification en chromosomes/contigs/scaffolds
		self.lstChr = []
		self.lstScaff = []
		self.lstRand = []
		self.lstNone = []
		for chrom in self.lstGenes.keys():
			if chrom in [None, "Un_random", "UNKN", "Un"]:
				self.lstNone.append(chrom)
				continue
			try:
				x = int(chrom)
				if x < 100:
					self.lstChr.append(chrom)
				else:
					self.lstScaff.append(chrom)
			except:
				c = chrom.lower()
				if "rand" in c:
					self.lstRand.append(chrom)
				else:
					keys = ["cont", "scaff", "ultra", "reftig", "_", "un", "mt"]
					for x in keys:
						if x in c:
							self.lstScaff.append(chrom)
							break
					else:
						if (chrom in ["U", "E64", "2-micron"]) or chrom.endswith("Het"):
							self.lstScaff.append(chrom)
						else:
							self.lstChr.append(chrom)

		# Trie les chromsomoses et cree l'index
		self.lstChr.sort()
		self.lstChrS = frozenset(self.lstChr)
		self.lstScaff.sort()
		self.lstScaffS = frozenset(self.lstScaff)
		self.lstRand.sort()
		self.lstRandS = frozenset(self.lstRand)
		self.lstNone.sort()
		self.lstNoneS = frozenset(self.lstNone)
		
		self.dicGenes = {}
		for c in self.lstGenes:
			self.lstGenes[c].sort()
			for (i,g) in enumerate(self.lstGenes[c]):
				assert g.chromosome == c
				for s in g.names:
					self.dicGenes[s] = GenePosition(c,i)
		

	# Rajoute un gene au genome
	############################
	def addGene(self, names, chromosome, beg, end, strand):

		assert 0 <= beg <= end
		assert strand in [-1,0,1]
		
		names = [intern(s) for s in names]
		chromosome = commonChrName(chromosome)
		self.lstGenes[chromosome].append( Gene(chromosome, beg, end, strand, names) )

	
	# Renvoie les genes presents sur le chromosome donne a certaines positions
	###########################################################################
	def getGenesAt(self, chr, beg, end):
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

		
	# Renvoie les genes presents aux alentours d'un gene donne (fenetre l en nombre de genes)
	##########################################################################################
	def getGenesNear(self, chr, index, l):
		if chr not in self.lstGenes:
			return []
		
		return self.lstGenes[chr][max(0, index-l):index+l+1]


	# Renvoie tous les genes
	#########################
	def __iter__(self):
		for t in self.lstGenes.itervalues():
			for g in t:
				yield g


	# Cherche la position d'un gene donne par ses noms)
	####################################################
	def getPosition(self, names):
		return set( (self.dicGenes[s] for s in names if s in self.dicGenes) )


	# Renvoie les autres noms d'un gene
	####################################
	def getOtherNames(self, name):
		if name not in self.dicGenes:
			return []
		(c,i) = self.dicGenes[name]
		return [x for x in self.lstGenes[c][i].names if x != name]


	# Fabrique la liste des orthologues entre les deux genomes
	###########################################################
	def buildOrthosTable(self, chr1, genome2, chr2, includeGaps, genesAnc):

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
		for l in self.f:
			c = l.replace('\n', '').split('\t')
			self.addGene(c[4].split(), c[0], int(c[1]), int(c[2]), int(c[3]))
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
			self.addGene(g, None, nb, nb, 0)

		print >> sys.stderr, "OK"


	####################################################################################
	# Cette fonction gere un fichier de genome ancestral                               #
	# il s'agit d'une liste de lignes de la forme "CHROMOSOME GENE_ESP1 GENE_ESP2 ..." #
	####################################################################################
	def __loadFromAncestralGenome__(self, withChr, concordeQualityFactor):

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
			self.addGene(champs, c, i, i+1, strand)
			i += 1
		
		print >> sys.stderr, "OK"


