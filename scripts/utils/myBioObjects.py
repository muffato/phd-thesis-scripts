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
	

