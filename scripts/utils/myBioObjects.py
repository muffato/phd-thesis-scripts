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
				noms = [x.strip() for x in l[indent].split('|')]
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
			
			# On ne continue que si la ligne suivante est indentee comme prevu
			if len(lignes) == 0  or lignes[-1][0] != indent:
				return None

			# On charge la ligne
			currLine = lignes.pop()
			
			# On charge toutes les sous arbres-fils tant que possible
			fils = []
			while True:
				tmp = recLoad(indent + 1)
				if tmp == None:
					break
				# On stocke (nom_du_fils, temps_d_evolution)
				fils.append( (tmp, currLine[2]-self.ages[tmp]) )
				self.parent[tmp] = currLine[1][0]
			
			# Un seul fils, on remonte le noeud
			if len(fils) == 1:
				return fils[0][0]
			# Plusieurs fils, on les enregistre
			elif len(fils) > 1:
				self.items[currLine[1][0]] = fils

			# Info standard
			self.commonNames[currLine[1][0]] = currLine[1][1:]
			self.ages[currLine[1][0]] = currLine[2]
			for s in currLine[1]:
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
	
				
		self.officialName = {}
		self.commonNames = self.newCommonNamesMapperInstance()
		self.items = self.newCommonNamesMapperInstance()
		self.parent = self.newCommonNamesMapperInstance()
		self.ages = self.newCommonNamesMapperInstance()
		self.branches = self.newCommonNamesMapperInstance()
		self.branchesSpecies = self.newCommonNamesMapperInstance()
		self.species = self.newCommonNamesMapperInstance()
		self.listSpecies = []
		self.listAncestr = []
		self.outgroupNode = self.newCommonNamesMapperInstance()
		self.dicGenes = self.newCommonNamesMapperInstance()
		self.dicGenomes = self.newCommonNamesMapperInstance()

		print >> sys.stderr, "Chargement de l'arbre phylogenique de %s ..." % nom,
		lignes = loadFile(nom)
		self.root = recLoad(0)
		
		print >> sys.stderr, "Analyse des donnees ...",
		recInitialize(self.root)
		
		# Les especes qui peuvent servir d'outgroup
		self.outgroupSpecies = self.newCommonNamesMapperInstance()
		for node in self.items:
			self.outgroupSpecies[node] = list(set(self.listSpecies).difference(self.species[node]))
		# Les noms a donner aux fichiers
		self.fileName = self.newCommonNamesMapperInstance()
		for n in self.commonNames:
			if n in self.items:
				self.fileName[n] = n.replace(' ', '_').replace('/', '-')
			else:
				self.fileName[n] = n.replace(' ', '.')
			# Ancienne version
			#self.fileName[n] = (self.commonNames[n]+[n])[0].replace(' ', '_').replace('/', '_')

		print >> sys.stderr, "OK"
		

	#
	# Renvoie un dictionnaire qui utilise en interne les noms officiels des taxons
	#  mais que l'on peut utiliser avec les noms communs
	#
	def newCommonNamesMapperInstance(self):
		class commonNamesMapper(dict):

			def __getitem__(d, name):
				if name in self.officialName:
					return dict.__getitem__(d, self.officialName[name])
				else:
					return dict.__getitem__(d, name)
			
			def __setitem__(d, name, value):
				if name in self.officialName:
					return dict.__setitem__(d, self.officialName[name], value)
				else:
					return dict.__setitem__(d, name, value)
			
		return commonNamesMapper()
		

	# Renvoie le dernier ancetre commun des deux entites
	def getFirstParent(self, anc1, anc2):
		if anc1 not in self.officialName:
			return None
		elif anc1 in self.listSpecies:
			e1 = set([anc1])
		else:
			e1 = set(self.species[anc1])
		
		if anc2 not in self.officialName:
			return None
		elif anc2 in self.listSpecies:
			e2 = set([anc2])
		else:
			e2 = set(self.species[anc2])
		
		if e1.issubset(e2):
			return anc2
		elif e2.issubset(e1):
			return anc1
			
		anc = self.parent[anc1]
		while not e2.issubset(self.species[anc]):
			anc = self.parent[anc]
			
		return anc

	# Renvoie tous les noeuds de l'arbre entre les deux especes (en remontant jusqu'a leur ancetre commun)
	def getNodesBetween(self, anc1, anc2):
		anc = self.getFirstParent(anc1, anc2)
		if anc == None:
			return []

		res = [anc]
		if anc != anc1:
			a = self.parent[anc1]
			while a != anc:
				res.append(a)
				a = self.parent[a]
		if anc != anc2:
			a = self.parent[anc2]
			while a != anc:
				res.append(a)
				a = self.parent[a]
		
		return res


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
			g = myGenomes.EnsemblGenome(template % self.fileName[esp])
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

	# A partir d'un noeud racine, liste les sous-arbres qui en partent dans toutes les directions
	# Appelle la fonction func sur chacun de ces sous-arbres avec le poids dependant de la position et de l'eloignement
	def travelFunc(self, node, func, useOutgroups):

		# Les fils a egalite avec un poids de 1
		for f in self.branches[node]:
			func(f, 1.)
		if not useOutgroups:
			return
		outgroup = []
		anc = node
		while anc in self.parent:
			par = self.parent[anc]
			outgroup.extend([(e,1./float(2*self.ages[par]-self.ages[node])) for (e,_) in self.items[par] if e != anc])
			anc = par
		s = sum([a for (_,a) in outgroup])
		for (e,a) in outgroup:
			func(e, a/s)






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
	

