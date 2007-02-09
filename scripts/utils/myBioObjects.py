#! /users/ldog/muffato/python -OO

import sys
import myTools
import myMaths
import myGenomes


##############
# Un gene :) #
##############
class Gene:

	def __init__(self, names, chromosome, beg, end, strand):
	
		self.names = tuple(names)
		self.chromosome = chromosome
		self.beginning = beg
		self.end = end
		self.strand = strand


#############################################
# Cette classe gere un arbre phylogenetique #
#############################################
class PhylogeneticTree:

	def __init__(self, nom):
		
		# Procedure de chargement du fichier
		def loadFile(s):
			f = myTools.myOpenFile(s, 'r')
			lignes = []
			for ligne in f:

				# On enleve le \n final et on coupe suivant les \t
				l = ligne[:-1].split('\t')
				
				# Le niveau d'indentation
				indent = 0
				for x in l:
					if x != '':
						break
					indent += 1

				# On stocke le triplet (indentation,noms,age)
				noms = l[indent].split('|')
				if len(l) > (indent+1):
					age = int(l[indent+1])
				else:
					age = 0
				lignes.append( (indent,noms,age) )
				
			f.close()
			lignes.reverse()
			return lignes
			

		# La procedure d'analyse des lignes du fichier
		def recLoad(indent):
			
			# Est-ce que la ligne suivante est une ligne fille ?
			if len(lignes) == 0  or lignes[-1][0] != indent:
				return None

			# On charge toutes les sous arbres-fils
			currLine = lignes.pop()
			fils = []
			while True:
				tmp = recLoad(indent + 1)
				if tmp == None:
					break
				fils.append( (tmp, currLine[2]-self.ages[tmp]) )
				self.parent[tmp] = currLine[1][0]
			
			if len(fils) != 0:
				self.items[currLine[1][0]] = fils
			self.commonNames[currLine[1][0]] = currLine[1][1:]
			for s in currLine[1]:
				self.ages[s] = currLine[2]
				self.officialName[s] = currLine[1][0]
				
			return currLine[1][0]
		
		
		# La procedure d'analyse de l'arbre
		def recInitialize(node):
			self.branchesSpecies[node] = []
			self.species[node] = []
			self.branches[node] = []
			self.outgroupNode[node] = None
			#if len(self.items[node]) != 0:
			if node in self.items:
				self.listAncestr.append(node)
				for (f,_) in self.items[node]:
					recInitialize(f)
					self.branches[node].append(f)
					self.branchesSpecies[node].append(self.species[f])
					self.species[node].extend(self.species[f])
				for (f,_) in self.items[node]:
					self.outgroupNode[f] = self.branches[node][:]
					self.outgroupNode[f].remove(f)
					
			else:
				self.listSpecies.append(node)
				self.branchesSpecies[node].append(node)
				self.species[node].append(node)
	
				
		self.items = {}
		self.commonNames = {}
		self.officialName = {}
		self.parent = {}
		self.ages = {}
		self.branches = {}
		self.branchesSpecies = {}
		self.species = {}
		self.listSpecies = []
		self.listAncestr = []
		self.outgroupNode = {}
		self.outgroupSpecies = {}
		self.dicGenes = {}
		self.dicGenomes = {}

		print >> sys.stderr, "Chargement de l'arbre phylogenique de %s ..." % nom,
		lignes = loadFile(nom)
		self.root = recLoad(0)
		print >> sys.stderr, "Analyse des donnees ...",
		recInitialize(self.root)
		
		for node in self.items:
			self.outgroupSpecies[node] = list(set(self.listSpecies).difference(self.species[node]))

		print >> sys.stderr, "OK"


		

	# Renvoie le dernier ancetre commun des deux entites
	def getFirstParent(self, anc1, anc2):
		if anc1 in self.species[anc2]:
			return anc2
		anc = anc1
		while anc2 not in self.species[anc]:
			anc = self.parent[anc]
		return anc


	# Renvoie l'arbre au format avec des parentheses
	def convertToFlatFile(self, anc):

		a = anc.replace(' ', '.')
		if anc in self.listSpecies:
			return a
		else:
			return "(" + ",".join(["%s:%d" % (self.convertToFlatFile(e),l) for (e,l) in self.items[anc]]) + ")%s|%d" % (a,self.ages[anc])



	# Charge toutes les especes qui descendent d'un ancetre
	def loadAllSpeciesSince(self, ancestr, template):
		self.loadSpeciesFromList(self.species.get(ancestr,self.listSpecies), template)
	
	# Charge toutes les especes outgroup d'un ancetre
	def loadAllSpeciesBefore(self, ancestr, template):
		self.loadSpeciesFromList(self.outgroupSpecies.get(ancestr,self.listSpecies), template)

	# Charge toutes les especes d'une liste
	def loadSpeciesFromList(self, lst, template):

		for esp in lst:
			g = myGenomes.EnsemblGenome(template % self.commonNames[esp][0])
			self.dicGenomes[esp] = g
			for x in g.dicGenes:
				self.dicGenes[x] = (esp, g.dicGenes[x][0], g.dicGenes[x][1])

	# Renvoie le nombre de genes dans chaque espece pour une famille donnee
	def findFamilyComposition(self, fam, breakWhenIncomplete=False):
		
		score = dict( [(e,[]) for e in self.dicGenomes] )
		
		for g in fam:
			if (g not in self.dicGenes) and breakWhenIncomplete:
				return None
			(e,_,_) = self.dicGenes[g]
			score[e].append(g)
		
		return score







#####################################################
# Cette classe gere un fichier resultat de Concorde #
#####################################################
class ConcordeFile:

	def __init__(self, nom):
		tmp = []
		f = myTools.myOpenFile(nom, 'r')
		for ligne in f:
			for x in ligne.split():
				tmp.append(int(x))
		
		f.close()
		i = tmp.index(0)
		self.res = tmp[i+1:] + tmp[1:i]
		self.buildIndex()

	def buildIndex(self):
		self.dic = {}
		for i in xrange(len(self.res)):
			self.dic[self.res[i]] = i
	
	def isMemeSens(self, other):
		s = 0
		for i in xrange(len(self.res)-1):
			if other.dic[self.res[i]] < other.dic[self.res[i+1]]:
				s += 1
			else:
				s -= 1
		return s>0

	def reverse(self):
		self.res.reverse()
		self.buildIndex()
	



	

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
		f = myTools.myOpenFile(nom, 'r')
		self.nom = nom

		for ligne in f:
			champs = ligne.split()

			# On convertit en nombre les chaines de caracteres
			for i in xrange(len(champs)):
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
		for i in xrange(self.nbGenes):
			yield (self.tabGenes1[self.indGenes1[i]], self.tabGenes2[self.indGenes2[i]], self.info[i], i)
	
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
		for i in xrange(self.nbGenes):
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



