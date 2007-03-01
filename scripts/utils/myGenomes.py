#! /users/ldog/muffato/python -OO

import sys
import myMaths
import myTools
import myBioObjects



#######################################################################################
# Cette fonction determine le type du fichier de genome a partir de la premiere ligne #
#######################################################################################
def loadGenome(nom):
	
	f = myTools.myOpenFile(nom, 'r')
	c = f.readline().split()
	f.close()
	
	# Fichier d'Ensembl: les trois champs du milieu sont des nombres
	try:
		x = int(c[1]) + int(c[2]) + int(c[3])
		return EnsemblGenome(nom)
	except Exception:
		# On a un genome ancestral
		pass
	
	# 1. Noms de chromosomes ?
	withChr = (len(c[0]) < 4)    # Les noms de genes font au minimum 4 caracteres
	withChr |= c[0] in ["ALPHA", "BETA", "DELTA", "EPSILON", "GAMMA", "PHI"]

	# 2. Facteur de qualite de concorde
	try:
		x = int(c[1])
		conc = True
	except Exception:
		conc = False
		
	return AncestralGenome(nom, chromPresents=withChr, concordeQualityFactor=conc)


##########################################
# Classe generale pour stocker un genome #
##########################################
class Genome:

	defaultChr = "-"

	#
	# Constructeur
	#
	def __init__(self, nom):

		# la liste des genes par chromosome
		self.lstGenes = {}
		# le nom	
		self.nom = nom
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
	
		c = gene.chromosome
		if c not in self.lstGenes:
			self.lstGenes[c] = []

		n = len(self.lstGenes[c])
		self.lstGenes[c].append( gene )
		
		for s in gene.names:
			self.dicGenes[s] = (c, n)
	
	#
	# Trie les chromsomoses
	#
	def sortGenome(self):
		
		self.lstChr.sort()
		self.lstScaff.sort()
		self.lstRand.sort()
		#print >> sys.stderr, self.lstChr
		#print >> sys.stderr, self.lstScaff
		#print >> sys.stderr, self.lstRand
		
		self.dicGenes = {}
		for c in self.lstGenes:
			self.lstGenes[c].sort(lambda g1, g2: cmp(g1.beginning, g2.beginning))
			for i in xrange(len(self.lstGenes[c])):
				for s in self.lstGenes[c][i].names:
					self.dicGenes[s] = (c, i)
		
	#
	# Renvoie les noms des genes presents sur le chromosome donne a certaines positions
	#	
	def getGenesAt(self, chr, beg = 0, end = sys.maxint):
		for g in self.lstGenes[chr]:
			if g.end >= beg and g.beginning <= end:
				yield g
	#
	# Renvoie les noms des genes presents aux alentours d'un gene donne
	#	
	def getGenesNearB(self, chr, index, l):
		if chr not in self.lstGenes:
			return
		g = self.lstGenes[chr][index]
		x1 = g.beginning-l
		x2 = g.end+l
		for i in xrange(index+1, len(self.lstGenes[chr])):
			g = self.lstGenes[chr][i]
			if g.beginning > x2:
				break
			yield g
		
		for i in xrange(index-1, -1, -1):
			g = self.lstGenes[chr][i]
			if g.end < x1:
				break
			yield g
	#
	# Renvoie les noms des genes presents aux alentours d'un gene donne
	#	
	def getGenesNearN(self, chr, index, l):
		if chr not in self.lstGenes:
			return
		
		for i in xrange(index+1, min(len(self.lstGenes[chr]), index+l+1)):
			yield self.lstGenes[chr][i]
		
		for i in xrange(index-1, max(-1, index-1-l), -1):
			yield self.lstGenes[chr][i]


	def iterOnChromosome(self, c):
		for g in self.lstGenes[c]:
			yield g
	
	def __iter__(self):
		for c in self.lstGenes:
			for g in self.lstGenes[c]:
				yield g

	def getPosition(self, gene):
		return myMaths.unique([self.dicGenes[s] for s in gene.names if s in self.dicGenes])



##############################################################
# Cette classe gere un fichier de liste de genes d'Ensembl   #
#   "Chr Debut Fin Brin ENSFAMxxxx Nom"                      #
# Convertit automatiquement les nombres romains en arabe     #
##############################################################
class EnsemblGenome(Genome):

	dicConvRomain = {"I":1, "II":2, "III":3, "IV":4, "V":5, "VI":6, "VII":7, "VIII":8, "IX":9, "X":10, "XI":11, "XII":12, "XIII":13, "XIV":14, "XV":15, "XVI":16, "XVII":17, "XVIII":18, "XIX":19, "XX":20, "XXI":21, "XXII":22, "XXIII":23, "XXIV":24, "XXV":25}

	#
	# Constructeur
	#
	def __init__(self, nom):
		
		Genome.__init__(self, nom)
		
		print >> sys.stderr, "Chargement de", nom, "...",
		
		f = myTools.myOpenFile(nom, 'r')
		
		# On lit chaque ligne
		for ligne in f:
			champs = ligne.split()
			
			# On convertit en entier le nom du chromosome si possible
			try:
				champs[0] = int(champs[0])
			except ValueError:
				pass
	
			self.addGene( myBioObjects.Gene([champs[-1]], champs[0], int(champs[1]), int(champs[2]), int(champs[3])) )
			
		f.close()
		
		# Les veritables chromosomes sont des entiers < 50 ou des entiers suivis de p/q/L/R/a/b ou W/X/Y/Z
		# Les chromosomes '*random*' ou 'UNKN' sont des scaffold mis bout a bout -> pas d'ordre utilisable
		# Le reste correspond aux scaffolds
		for c in self.lstGenes:
			if (type(c) == int) or (type(c) == long):
				if c < 50:
					self.lstChr.append(c)
				else:
					self.lstScaff.append(c)
			elif ((c[-1] in "pqLRab") and (c != "c6_QBL")) or (c in "WXYZ"):
				self.lstChr.append(c)
			elif ("ando" in c) or (c == "UNKN"):
				self.lstRand.append(c)
			else:
				self.lstScaff.append(c)
		self.sortGenome()
		
		print >> sys.stderr, "OK"


################################################################################
# Cette classe gere un fichier d'orthologues a partir duquel on cree un genome #
# Les noms des genes sont en positions 1 et 4                                  #
################################################################################
class GenomeFromOrthosList(Genome):

	#
	# Constructeur
	#
	def __init__(self, nom, filter=[]):
		
		Genome.__init__(self, nom)
		
		print >> sys.stderr, "Chargement de", nom, "...",
		
		f = myTools.myOpenFile(nom, 'r')
		combin = myTools.myCombinator([])
		
		# On lit chaque ligne
		for ligne in f:
			champs = ligne.split()
			if (champs[6] in filter) or (len(filter) == 0):
				combin.addLink([champs[0], champs[3]])
		f.close()
		
		nb = 0
		for g in combin:
			self.addGene( myBioObjects.Gene(tuple(g), Genome.defaultChr, nb, nb, 0) )
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
			print >> sys.stderr, "Chargement du genome ancestral de", nom, "... ",
		else:
			print >> sys.stderr, "Chargement des genes ancestraux de", nom, "... ",
		
		# On initialise tout
		f = myTools.myOpenFile(nom, 'r')
		strand = 0
		for ligne in f:
			champs = ligne.split()

			# Le chromosome du gene lu
			if not chromPresents:
				c = Genome.defaultChr
			else:
				try:
					c = int(champs[0])
				except ValueError:
					c = champs[0]
				del champs[0]
			
			if concordeQualityFactor:
				# Fichiers de genomes ancestraux avec score de Concorde
				strand = int(champs[0])
				del champs[0]
				
			
			# La position du gene lu
			if c in self.lstGenes:
				i = len(self.lstGenes[c])
			else:
				i = 0
			
			# On ajoute le gene
			self.addGene( myBioObjects.Gene(champs, c, i, i, strand) )
		
		f.close()
		self.lstChr = sorted(self.lstGenes)
		
		print >> sys.stderr, "OK"

	#
	# Cette fonction separe un genome ancestral en blocs dans chaque
	# chromosome de chaque espece
	#
	def splitChr(self, geneBank, chrAnc):
		
		# 1ere etape: separer pour chaque espece en chromosomes
		blocsAnc = {}
		j = 0
		for g in self.lstGenes[chrAnc]:
			for s in g.names:
				if s not in geneBank.dicGenes:
					continue
				(e,c,i) = geneBank.dicGenes[s]
				if e not in blocsAnc:
					blocsAnc[e] = {}
				if c not in blocsAnc[e]:
					blocsAnc[e][c] = []
				blocsAnc[e][c].append( (i,s,j) )
			j += 1

		# 2eme etape: trier chacun de ces ensembles et creer un
		# dictionnaire qui pemet de retrouver la position de chaque gene
		dico = {}
		for e in blocsAnc:
			for c in blocsAnc[e]:
				blocsAnc[e][c].sort()
				for i in xrange(len(blocsAnc[e][c])):
					(_,s,j) = blocsAnc[e][c][i]
					dico[s] = (e,c,i,j)
		
		return (blocsAnc, dico)


