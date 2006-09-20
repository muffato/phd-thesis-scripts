#! /usr/bin/python2.4

import sys
import string
import myTools
import myMaths


##############
# Un gene :) #
##############
class Gene:

	def __init__(self, names, chromosome, beg, end, strand):
	
		self.names = names
		self.chromosome = chromosome
		self.beg = beg
		self.end = end
		self.strand = stand

#######################################################################################
# Cette fonction determine le type du fichier de genome a partir de la premiere ligne #
#######################################################################################
def loadGenome(nom):
	
	f = myTools.myOpenFile(nom)
	c = f.readline().split()
	f.close()
	try:
		x = int(c[1]) + int(c[2]) + int(c[3])
		return EnsemblGenome(nom)
	except Exception:
		c = nom.split('/')
		return AncestralGenome(nom, c[-1][0] in string.uppercase)


##########################################
# Classe generale pour stocker un genome #
##########################################
class Genome:

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
	def addGene(self, names, chromosome, beg, end, strand):
	
		if chromosome not in self.lstGenes:
			self.lstChr.append(chromosome)
			self.lstChr.sort()
			self.lstGenes[chromosome] = []

		n = len(self.lstGenes[chromosome])
		self.lstGenes[chromosome].append( (beg, end, strand, tuple(names)) )
		
		for s in names:
			self.dicGenes[s] = (chromosome, n)
	
	#
	# Trie les chromsomoses
	#
	def sortGenome(self):
		
		self.dicGenes = {}
		for c in self.lstGenes:
			self.lstGenes[c].sort()
			for i in range(len(self.lstGenes[c])):
				#s = self.lstGenes[c][i][3]
				for s in self.lstGenes[c][i][3]:
					self.dicGenes[s] = (c, i)
		
		# La liste triee des chromosomes
		self.lstChr = self.lstGenes.keys()
		self.lstChr.sort()

	#
	# Renvoie les noms des genes presents sur le chromosome donne a certaines positions
	#	
	def getGenesAt(self, chr, beg = 0, end = sys.maxint):
		return [s for (x1,x2,_,s) in self.lstGenes[chr] if x2 >= beg and x1 <= end]



##############################################################
# Cette classe gere un fichier de liste de genes d'Ensembl   #
#   "Chr Debut Fin Brin Nom"                                 #
# Le filtre correspond a une liste de chromosomes a eliminer #
##############################################################
class EnsemblGenome(Genome):

	#
	# Constructeur
	#
	def __init__(self, nom):
		
		Genome.__init__(self, nom)
		
		print >> sys.stderr, "Chargement de", nom, "...",
		
		f = myTools.myOpenFile(nom)
		
		# On lit chaque ligne
		for ligne in f:
			champs = ligne.split()
			
			# On convertit en entier le nom du chromosome si possible
			try:
				champs[0] = int(champs[0])
			except ValueError:
				pass
	
			self.addGene(champs[4:], champs[0], int(champs[1]), int(champs[2]), int(champs[3]))
		
		f.close()

		self.sortGenome()
		
		print >> sys.stderr, "OK"



####################################################################################
# Cette classe gere un fichier de genome ancestral                                 #
# il s'agit d'une liste de lignes de la forme "CHROMOSOME GENE_ESP1 GENE_ESP2 ..." #
####################################################################################
class AncestralGenome(Genome):

	defaultChr = "-"
	
	def __init__(self, nom, chromPresents):

		Genome.__init__(self, nom)
		
		if chromPresents:
			print >> sys.stderr, "Chargement du genome ancestral de", nom, "... ",
		else:
			print >> sys.stderr, "Chargement des genes ancestraux de", nom, "... ",
		
		# On initialise tout
		f = myTools.myOpenFile(nom)

		for ligne in f:
			champs = ligne.split()

			# Le chromosome du gene lu
			if not chromPresents:
				c = AncestralGenome.defaultChr
			else:
				try:
					c = int(champs[0])
				except ValueError:
					c = champs[0]
				del champs[0]
			
			# On efface un entier s'il y en a un
			# C'est le cas des fichiers de genomes ancestraux avec score de Concorde
			try:
				int(champs[0])
				del champs[0]
			except ValueError:
				pass
			
			# La position du gene lu
			if c in self.lstChr:
				i = len(self.lstGenes[c])
			else:
				i = 0
			
			# On ajoute le gene
			self.addGene(champs, c, i, i, 0)
		
		f.close()
		
		print >> sys.stderr, "OK"

	#
	# Cette fonction separe un genome ancestral en blocs dans chaque
	# chromosome de chaque espece
	#
	def splitChr(self, geneBank, chrAnc):
		
		# 1ere etape: separer pour chaque espece en chromosomes
		blocsAnc = {}
		j = 0
		for (_,_,_,x) in self.lstGenes[chrAnc]:
			for s in x:
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
				for i in range(len(blocsAnc[e][c])):
					(_,s,j) = blocsAnc[e][c][i]
					dico[s] = (e,c,i,j)
		
		return (blocsAnc, dico)


#########################################################################
# Cette classe gere une liste de fichiers de genes.                     #
# Le format des lignes est "Nom_Espece nom_fichier_genes chr1 chr2 ..." #
#    avec une liste de chromosomes a ne pas prendre en compte           #
#########################################################################
class GeneBank:

	def __init__(self, nom, only = []):

		self.dicGenes = {}
		self.dicEspeces = {}
		self.lstEspecesNonDup = []
		self.lstEspecesDup = []
		self.lstEspecesOutgroup = []
		f = open(nom, 'r')
		
		for ligne in f:

			champs = ligne.split()
			
			# Un commentaire
			if ligne[0] == '#':
				continue

			# Le nom de l'espece
			if champs[0][-1] in ['.', '*']:
				s = champs[0][:-1]
				dup = champs[0][-1]
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
			elif dup == '*':
				self.lstEspecesOutgroup.append(s)
			self.dicEspeces[s] = g
			
			for x in g.dicGenes:
				self.dicGenes[x] = (s, g.dicGenes[x][0], g.dicGenes[x][1])

		f.close()




#############################################
# Cette classe gere un arbre phylogenetique #
#  "A = [('*2',140), ('T',450)]"            #
#############################################
class PhylogeneticTree:

	def __init__(self, nom):
		
		print >> sys.stderr, "Chargement de l'arbre phylogenique de %s ..." % nom,
		
		self.items = {}
		f = open(nom, 'r')
		for ligne in f.readlines():
			# On separe la ligne
			i = ligne.index('=')
			# A gauche, le nom du resultat
			nom = ligne[:i].strip()
			# A droite, les descendants
			self.items[nom] = eval(ligne[i+1:])
		
		print >> sys.stderr, "OK"
		
	def getDescendants(self, anc):
	
		if anc not in self.items:
			return []
			
		return [myMaths.flatten(self.getDescendants(e)) for (e,_) in self.items[anc]]
	
	def getFils(self, anc):
		if anc not in self.items:
			return []
			
		return [e for (e,_) in self.items[anc]]
		

##################################################################################################
# Cette classe gere un fichier d'orthologues d'Ensembl                                           #
#    "Khumain start stop brin gene_HS gene_Tn KTn start stop ortho_type %ID %cov"                #
#    "1 914833  920104  1  ENSG00000187634 GSTENG00028217001.1  9  2679651 2687663 UBRH  86  58" #
##################################################################################################
class EnsemblTwoSpeciesOrthos:

	def __init__(self, nom):

		print >> sys.stderr, "Chargement de %s ... " % nom,
		
		# Les tableaux des genes 1 et des genes 2
		self.tabGenes1 = []
		self.tabGenes2 = []
		self.dictGenes1 = {}
		self.dictGenes2 = {}
		self.lstChr1 = {}
		self.lstChr2 = {}
		
		self.info = []
		# Le nombre de genes
		j = 0

		# On lit les lignes
		f = myTools.myOpenFile(nom)
		self.nom = nom

		for ligne in f:
			champs = ligne.split()

			# On convertit en nombre les chaines de caracteres
			for i in range(len(champs)):
				try:
					champs[i] = int(champs[i])
				except ValueError:
					pass
			# On enregistre des lignes "Chr - debut - fin - nom - indice"

			# Le gene 1
			nom1 = champs[4].split(".")[0]
			self.tabGenes1.append( (champs[0], champs[1], champs[2], nom1, j) )
			if champs[0] not in self.lstChr1:
				self.lstChr1[champs[0]] = []
			self.lstChr1[champs[0]].append(j)
			self.dictGenes1[nom1] = j
			
			# Le gene 2
			nom2 = champs[5].split(".")[0]
			self.tabGenes2.append( (champs[6], champs[7], champs[8], nom2, j) )
			if champs[6] not in self.lstChr2:
				self.lstChr2[champs[6]] = []
			self.lstChr2[champs[6]].append(j)
			self.dictGenes2[nom2] = j
			
			# Les infos supplementaires
			self.info.append( [j, champs[3]] + champs[9:] )
			
			# Et un gene en plus !
			j += 1
		
		f.close()
		self.nbGenes = j
		self.lstChr1Sorted = self.lstChr1.keys()
		self.lstChr1Sorted.sort()
		self.lstChr2Sorted = self.lstChr2.keys()
		self.lstChr2Sorted.sort()
		
		self.sort()
		print >> sys.stderr, "OK"

	#
	# Renvoie un iterateur pour parcourir tous les orthologues
	#
	def __iter__(self):
		class OrthoIterator:
			def __init__(self, ortho):
				self.ortho = ortho
				self.index = -1
			def next(self):
				self.index += 1
				if self.index >= self.ortho.nbGenes:
					raise StopIteration
				return (self.ortho.tabGenes1[self.ortho.indGenes1[self.index]], self.ortho.tabGenes2[self.ortho.indGenes2[self.index]], self.ortho.info[self.index], self.index)
		return OrthoIterator(self)
	
	#
	# Renvoie la liste des orthologues, tries sur le gene1
	#
	def iterGenes1(self):
		return [ (g, self.tabGenes2[self.indGenes2[g[4]]], self.info[g[4]], self.indGenes1[g[4]]) for g in self.tabGenes1]

	#
	# Renvoie la liste des orthologues, tries sur le gene1
	#
	def iterGenes2(self):
		return [ (self.tabGenes1[self.indGenes1[g[4]]], g, self.info[g[4]], self.indGenes2[g[4]]) for g in self.tabGenes2]

	
	
	#
	# Inverse les references gene1 et gene2
	#
	def swapGene1Gene2(self):
		t = self.tabGenes1
		self.tabGenes1 = self.tabGenes2
		self.tabGenes2 = t
		
		t = self.dictGenes1
		self.dictGenes1 = self.dictGenes2
		self.dictGenes2 = t

		t = self.indGenes1
		self.indGenes1 = self.indGenes2
		self.indGenes2 = t

		t = self.lstChr1
		self.lstChr1 = self.lstChr2
		self.lstChr2 = t

		t = self.lstChr1Sorted
		self.lstChr1Sorted = self.lstChr2Sorted
		self.lstChr2Sorted = t

		
	#
	# Trie les genes selon les coordonnes genomiques
	#
	def sort(self):
	
		# On trie les deux listes de genes suivant leurs coordonnees
		self.tabGenes1.sort()
		self.tabGenes2.sort()

		# Permet d'associer les indices reels aux indices tries
		self.indGenes1 = range(self.nbGenes)
		self.indGenes2 = range(self.nbGenes)
		for i in range(self.nbGenes):
			self.indGenes1[self.tabGenes1[i][4]] = i
			self.indGenes2[self.tabGenes2[i][4]] = i

	#
	# Rajoute une paire de genes
	#
	def addOrtho(self, g1, g2):

		# D'abord le gene 1
		self.tabGenes1.append( (g1[0], g1[1], g1[2], g1[3], self.nbGenes) )
		if g1[0] not in self.lstChr1:
			self.lstChr1[g1[0]] = []
		self.lstChr1[g1[0]].append(self.nbGenes)
		self.dictGenes1[g1[3]] = self.nbGenes
		
		# Puis le gene 2
		self.tabGenes2.append( (g2[0], g2[1], g2[2], g2[3], self.nbGenes) )
		if g2[0] not in self.lstChr2:
			self.lstChr2[g2[0]] = []
		self.lstChr2[g2[0]].append(self.nbGenes)
		self.dictGenes2[g2[3]] = self.nbGenes

		self.info.append( ["", 0, 0, 0, self.nbGenes] )
		self.nbGenes += 1
		self.sort()



	#
	# Filtre les genes 1 et ne garde que ceux qui sont sur un certain
	# chromosome, entre certaines coordonnees
	#
	def getGene1Filtered(self, chr, debut = -1, fin = sys.maxint):
		return [g for g in self.tabGenes1 if g[0] == chr and g[1] >= debut and g[2] <= fin]

	#
	# Filtre les genes 2 et ne garde que ceux qui sont sur un certain
	# chromosome, entre certaines coordonnees
	#
	def getGene2Filtered(self, chr, debut = -1, fin = sys.maxint):
		return [g for g in self.tabGenes2 if g[0] == chr and g[1] >= debut and g[2] <= fin]




	#
	# Renvoie les genes 1 correspondant aux genes 2
	#
	def getGenes1(self, tab2):
		return [self.tabGenes1[self.indGenes1[g[4]]] for g in tab2]

	#
	# Renvoie les genes 2 correspondant aux genes 1
	#
	def getGenes2(self, tab1):
		return [self.tabGenes2[self.indGenes2[g[4]]] for g in tab1]




#####################################################
# Cette classe gere un fichier resultat de Concorde #
#####################################################

class ConcordeFile:

	def __init__(self, nom):
		tmp = []
		f = open(nom)
		for ligne in f:
			for x in ligne.split():
				tmp.append(int(x))
		
		f.close()
		i = tmp.index(0)
		self.res = tmp[i+1:] + tmp[1:i]
		self.buildIndex()

	def buildIndex(self):
		self.dic = {}
		for i in range(len(self.res)):
			self.dic[self.res[i]] = i
	
	def isMemeSens(self, other):
		s = 0
		for i in range(len(self.res)-1):
			if other.dic[self.res[i]] < other.dic[self.res[i+1]]:
				s += 1
			else:
				s -= 1
		return s>0

	def reverse(self):
		self.res.reverse()
		self.buildIndex()
	
