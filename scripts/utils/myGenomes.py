
import sys
import operator
import myMaths
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


# On convertit les noms de chromosomes en int si possible
def commonChrName(x):
	try:
		return int(x)
	except Exception:
		return intern(x)



#######################################################################################
# Cette fonction determine le type du fichier de genome a partir de la premiere ligne #
#######################################################################################
def loadGenome(nom):

	# Classe qui simule un fichier
	# Le but est de pouvoir acceder au fichier et lire la premiere ligne sans devoir le fermer pour le reouvrir juste apres
	class myFileChecker:
		def __init__(self, nom):
			self.nom = nom
			self.f = myTools.myOpenFile(nom, 'r')
			try:
				self.firstLine = self.f.next()
			except StopIteration:
				self.firstLine = ""

		def __iter__(self):
			return self

		def next(self):
			self.next = self.f.next
			return self.firstLine

		def close(self):
			return self.f.close()
	
	f = myFileChecker(nom)
	c = f.firstLine.split()
	
	# Fichier de genome d'Ensembl: les trois champs du milieu sont des nombres
	try:
		x = int(c[1]) + int(c[2]) + int(c[3])
		return EnsemblGenome(f)
	except (ValueError, IndexError):
		# On a un genome ancestral
		pass
		
	# Fichier d'orthologues classique, 4 champs indiquent les scores de similarite
	try:
		x = int(c[7]) + int(c[8]) + int(c[9]) + int(c[10])
		return EnsemblOrthosListGenome(f)
	except (ValueError, IndexError):
		# On a un genome ancestral
		pass
	
	# 1. Noms de chromosomes ?
	firstElt = ""
	try:
		firstElt = c[0]
		int(firstElt[0])
		withChr = True
	except (ValueError, IndexError):
		withChr = (len(firstElt) < 4)    # Les noms de genes font au minimum 4 caracteres
		withChr |= firstElt.upper().strip() in ["ALPHA", "BETA", "DELTA", "EPSILON", "GAMMA", "PHI"]

	# 2. Facteur de qualite de concorde
	try:
		x = int(c[1])
		conc = True
	except (ValueError, IndexError):
		conc = False
		
	return AncestralGenome(f, chromPresents=withChr, concordeQualityFactor=conc)


##########################################
# Classe generale pour stocker un genome #
##########################################
class Genome:

	# Le chromosome par defaut (liste de genes ancestraux sans chromosomes)
	defaultChr = "-"
	
	#
	# Constructeur
	#
	def __init__(self, nom):

		# le nom et l'instance de file
		if type(nom) == str:
			self.nom = nom
			self.f = myTools.myOpenFile(nom, 'r')
		else:
			self.nom = nom.nom
			self.f = nom
		
		# la liste des genes par chromosome
		self.lstGenes = myTools.defaultdict(list)
		# Associer un nom de gene a sa position sur le genome
		self.dicGenes = {}
		# La liste triee des noms des chromosomes, des scaffolds et des randoms
		self.lstChr = []
		self.lstScaff = []
		self.lstRand = []
		
	#
	# Rajoute un gene au genome
	#
	def addGene(self, gene):
	
		self.lstGenes[gene.chromosome].append(gene)
	
	#
	# Trie les chromsomoses
	#
	def sortGenome(self):
		
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

	def __iter__(self):
		for c in self.lstGenes:
			for g in self.lstGenes[c]:
				yield g

	def getPosition(self, gene):
		return set( (self.dicGenes[s] for s in gene.names if s in self.dicGenes) )

	def getOtherNames(self, name):
		if name not in self.dicGenes:
			return []
		(c,i) = self.dicGenes[name]
		return [x for x in self.lstGenes[c][i].names if x != name]

	# Fabrique la liste des orthologues entre les deux genomes
	def buildOrthosTable(self, chr1, genome2, chr2, includeGaps, genesAnc = None):

		# Tous les orthologues entre pour les chromosomes OK du genome 1
		res = {}
		for c1 in chr1:
			res[c1] = []
			for (i1,g1) in enumerate(self.lstGenes[c1]):
				tmp = genome2.getPosition(g1)
				if genesAnc != None:
					for (c,i) in genesAnc.getPosition(g1):
						tmp.update(genome2.getPosition(genesAnc.lstGenes[c][i]))
				
				# On ne garde que les chromsomes OK du genome 2
				tmp = [(c,i) for (c,i) in tmp if c in chr2]
				# +/- includeGaps
				if includeGaps or (len(tmp) > 0):
					res[c1].append( (i1,tmp) )
			res[c1].reverse()

		return res



##############################################################
# Cette classe gere un fichier de liste de genes d'Ensembl   #
#   "Chr Debut Fin Brin ENSFAMxxxx Nom"                      #
# NE Convertit PAS les nombres romains en arabe              #
##############################################################
class EnsemblGenome(Genome):

	dicConvRomain = {"I":1, "II":2, "III":3, "IV":4, "V":5, "VI":6, "VII":7, "VIII":8, "IX":9, "X":10, "XI":11, "XII":12, "XIII":13, "XIV":14, "XV":15, "XVI":16, "XVII":17, "XVIII":18, "XIX":19, "XX":20, "XXI":21, "XXII":22, "XXIII":23, "XXIV":24, "XXV":25}

	#
	# Constructeur
	#
	def __init__(self, nom):
		
		Genome.__init__(self, nom)
		
		print >> sys.stderr, "Chargement de", self.nom, "...",
		
		# On lit chaque ligne
		for ligne in self.f:
			champs = ligne.split()
			self.addGene( Gene(champs[4:], champs[0], int(champs[1]), int(champs[2]), int(champs[3])) )
			
		self.f.close()
		
		# Les veritables chromosomes sont des entiers < 50 ou des entiers suivis de p/q/L/R/a/b ou W/X/Y/Z
		# Les chromosomes '*random*' ou 'UNKN' sont des scaffold mis bout a bout -> pas d'ordre utilisable
		# Le reste correspond aux scaffolds
		for c in self.lstGenes:
			if (type(c) == int) or (type(c) == long):
				if c < 50:
					self.lstChr.append(c)
				else:
					self.lstScaff.append(c)
			elif ((c[-1] in "pqLRab") and (c != "c6_QBL")) or ((c[0] in "WXYZ") and len(c) <= 2) or (c in self.dicConvRomain) or (c[:5] == "group"):
				self.lstChr.append(c)
			elif ("andom" in c) or (c == "UNKN"):
				self.lstRand.append(c)
			else:
				self.lstScaff.append(c)
		self.sortGenome()
		
		print >> sys.stderr, "OK"


################################################################################
# Cette classe gere un fichier d'orthologues a partir duquel on cree un genome #
# Les noms des genes sont en positions 1 et 4                                  #
################################################################################
class EnsemblOrthosListGenome(Genome):

	#
	# Constructeur
	#
	def __init__(self, nom, homologyFilter=[], ancFilter=[]):
		
		Genome.__init__(self, nom)
		
		print >> sys.stderr, "Chargement de", self.nom, "...",
		
		combin = myTools.myCombinator([])
		
		# On lit chaque ligne
		for ligne in self.f:
			champs = ligne[:-1].split('\t')
			if len(champs) == 12:
				# Nouvelle version des fichiers
				if (champs[-1] not in homologyFilter) and (len(homologyFilter) != 0):
					continue
				if (champs[6] not in ancFilter) and (len(ancFilter) != 0):
					continue
			else:
				# Ancienne version des fichiers
				if (champs[6] not in homologyFilter) and (len(homologyFilter) != 0):
					continue
			combin.addLink([champs[0], champs[3]])
		self.f.close()
		
		nb = 0
		for g in combin:
			self.addGene( Gene(g, Genome.defaultChr, nb, nb, 0) )
			nb += 1

		print >> sys.stderr, "OK"


####################################################################################
# Cette classe gere un fichier de genome ancestral                                 #
# il s'agit d'une liste de lignes de la forme "CHROMOSOME GENE_ESP1 GENE_ESP2 ..." #
####################################################################################
class AncestralGenome(Genome):

	
	def __init__(self, nom, chromPresents=False, concordeQualityFactor=False):

		Genome.__init__(self, nom)
		
		if chromPresents:
			print >> sys.stderr, "Chargement du genome ancestral de", self.nom, "...",
		else:
			print >> sys.stderr, "Chargement des genes ancestraux de", self.nom, "...",
		
		# On initialise tout
		strand = 0
		c = Genome.defaultChr
		i = 0
		for ligne in self.f:
			champs = ligne.split()

			if len(champs) == 0:
				continue

			# Le chromosome du gene lu
			if chromPresents:
				c = champs.pop(0)
			
			if concordeQualityFactor:
				# Fichiers de genomes ancestraux avec score de Concorde
				strand = int(champs.pop(0))
				
			# On ajoute le gene
			self.addGene( Gene(champs, c, i, i, strand) )
			i += 1
		
		self.f.close()
		self.lstChr = self.lstGenes.keys()
		self.sortGenome()
		
		print >> sys.stderr, "OK"


