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
		# La liste triee des chromosomes
		self.lstChr = []
		
	#
	# Rajoute un gene au genome
	#
	def addGene(self, gene):
	
		c = gene.chromosome
		if c not in self.lstGenes:
			self.lstChr.append(c)
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
	def getGenesNear(self, chr, index, l):
		g = self.lstGenes[chr][index]
		x1 = g.beginning
		x2 = g.end
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

	def iterOnChromosome(self, c):
		for g in self.lstGenes[c]:
			yield g
	
	def __iter__(self):
		for c in self.lstChr:
			for g in self.lstGenes[c]:
				yield g

	def getPosition(self, gene):
		return myMaths.unique([self.dicGenes[s] for s in gene.names if s in self.dicGenes])

##############################################################
# Cette classe gere un fichier de liste de genes d'Ensembl   #
#   "Chr Debut Fin Brin Nom"                                 #
##############################################################
class EnsemblGenome(Genome):

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
	
			self.addGene( myBioObjects.Gene(champs[4:], champs[0], int(champs[1]), int(champs[2]), int(champs[3])) )
		f.close()

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
	def __init__(self, nom):
		
		Genome.__init__(self, nom)
		
		print >> sys.stderr, "Chargement de", nom, "...",
		
		f = myTools.myOpenFile(nom, 'r')
		combin = myTools.myCombinator([])
		
		# On lit chaque ligne
		for ligne in f:
			champs = ligne.split()
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
			if c in self.lstChr:
				i = len(self.lstGenes[c])
			else:
				i = 0
			
			# On ajoute le gene
			self.addGene( myBioObjects.Gene(champs, c, i, i, strand) )
		
		f.close()
		self.lstChr.sort()
		
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


###########################################################
# Cette classe gere une liste de fichiers de genes.       #
# Le format des lignes est "Nom_Espece nom_fichier_genes" #
# Un symbole indique le type de l'espece                  #
#   Rien -> Vertebre non duplique                         #
#    . -> Vertebre duplique                               #
#    * -> Non vertebre = outgroup                         #
###########################################################
class GeneBank:

	def __init__(self, nom, only = []):

		self.dicGenes = {}
		self.dicEspeces = {}
		self.lstEspecesNonDup = []
		self.lstEspecesDup = []
		self.lstEspecesOutgroup = []
		f = myTools.myOpenFile(nom, 'r')
		
		for ligne in f:

			champs = ligne.split()
			
			if len(champs) != 2:
				continue

			# Un commentaire
			if ligne[0] == '#':
				continue

			# Le nom de l'espece
			if champs[0][0] in ['.', '-']:
				s = champs[0][1:]
				dup = champs[0][0]
			else:
				s = champs[0]
				dup = ""
			
			# Est-ce que l'espece a ete filtree
			if s not in only and len(only) > 0:
				continue

			g = EnsemblGenome(champs[1])
			
			if dup == "":
				self.lstEspecesNonDup.append(s)
			elif dup == '.':
				self.lstEspecesDup.append(s)
			elif dup == '-':
				self.lstEspecesOutgroup.append(s)
			self.dicEspeces[s] = g
			
			for x in g.dicGenes:
				self.dicGenes[x] = (s, g.dicGenes[x][0], g.dicGenes[x][1])

		f.close()


