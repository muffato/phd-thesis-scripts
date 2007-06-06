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

	def __repr__(self):
		return "Gene %s on chr %s from %d to %d on strand %d" % ("/".join(self.names), self.chromosome, self.beginning, self.end, self.strand)


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
				for f in self.branches[node]:
					self.outgroupNode[f] = self.branches[node][:]
					self.outgroupNode[f].remove(f)
					
			else:
				self.listSpecies.append(node)
				self.branchesSpecies[node].append(node)
				self.species[node].append(node)
		
		def buildPhylLinks():
			
			# Initialisation de la table de tous les liens entre les objets
			self.dicLinks = self.newCommonNamesMapperInstance()
			self.dicParents = self.newCommonNamesMapperInstance()
			for f1 in self.commonNames:
				self.dicLinks[f1] = self.newCommonNamesMapperInstance()
				self.dicParents[f1] = self.newCommonNamesMapperInstance()
				for f2 in self.commonNames:
					self.dicLinks[f1][f2] = []
				self.dicLinks[f1][f1] = [f1]
				self.dicParents[f1][f1] = f1
			
			# Remplissage de chaque objet avec tous ses parents et vice-versa
			for f1 in self.commonNames:
				parent = f1
				while parent != self.root:
					f2 = parent
					parent = self.parent[f2]
					self.dicLinks[f1][parent] = self.dicLinks[f1][f2] + [parent]
					self.dicLinks[parent][f1] = [parent] + self.dicLinks[f2][f1]
					self.dicParents[f1][parent] = self.dicParents[parent][f1] = parent
					
			# Liens entre objets de branches differentes
			tmp = set()
			for f1 in self.commonNames:
				for f2 in self.commonNames:
					if len(self.dicLinks[f1][f2]) != 0:
						continue
					for f3 in self.commonNames:
						if len(set(self.dicLinks[f1][f3]).intersection(self.dicLinks[f3][f2])) == 1:
							tmp.add( (f1, f2, f3) )
							break
			# On enregistre le tout
			for (f1, f2, f3) in tmp:
				self.dicLinks[f1][f2] = self.dicLinks[f1][f3][:-1] + self.dicLinks[f3][f2]
				self.dicParents[f1][f2] = f3
				
		self.officialName = {}
		self.commonNames = self.newCommonNamesMapperInstance()
		self.items = self.newCommonNamesMapperInstance()
		self.parent = self.newCommonNamesMapperInstance()
		self.ages = self.newCommonNamesMapperInstance()

		print >> sys.stderr, "Chargement de l'arbre phylogenique de %s ..." % nom,
		lignes = loadFile(nom)
		self.root = recLoad(0)
		
		print >> sys.stderr, "Analyse des donnees ...",
		
		self.branches = self.newCommonNamesMapperInstance()
		self.branchesSpecies = self.newCommonNamesMapperInstance()
		self.species = self.newCommonNamesMapperInstance()
		self.listSpecies = []
		self.listAncestr = []
		self.outgroupNode = self.newCommonNamesMapperInstance()
		self.dicGenes = self.newCommonNamesMapperInstance()
		self.dicGenomes = self.newCommonNamesMapperInstance()
		recInitialize(self.root)
		
		# Les especes qui peuvent servir d'outgroup
		self.outgroupSpecies = self.newCommonNamesMapperInstance()
		for node in self.listAncestr:
			self.outgroupSpecies[node] = list(set(self.listSpecies).difference(self.species[node]))
		# Les noms a donner aux fichiers
		self.fileName = self.newCommonNamesMapperInstance()
		for n in self.listAncestr:
			self.fileName[n] = n.replace(' ', '_').replace('/', '-')
		for n in self.listSpecies:
			self.fileName[n] = n.replace(' ', '.')

		buildPhylLinks()

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
			esp = self.officialName[esp]
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


	#
	# Cree une structure d'arbre secondaire pour le calcul d'une moyenne sur les especes
	# Si on utilise les outgroups, il faut faire basculer l'arbre pour que les outgroups
	#    soient des fils (de fils de fils ...) de plus en plus eloignes
	#
	def initCalcDist(self, node, useOutgroups):
		self.tmpItems = self.items.copy()
		if useOutgroups:
			anc = node
			while anc in self.parent:
				par = self.parent[anc]
				self.tmpItems[anc].append( (par,self.ages[par]-self.ages[anc]) )
				self.tmpItems[par] = [x for x in self.tmpItems[par] if x[0] != anc]
				anc = par
		self.tmpItems[0] = self.tmpItems[node]

	#
	# Lance le calcul de la moyenne etant donne les valeurs stockees dans values
	#
	def calcDist(self, values, init=0):

		# La partie recursive
		def recCalc(anc):
			d = 0.
			s = 0.
			nb = 0
		
			# Moyenne sur tous les fils du noeud
			for (e,p) in self.tmpItems[anc]:
				# Des valeurs d'ancetres / d'especes fixees
				if e in values:
					val = values[e]
				# Le fils est un nouveau noeud - Appel recursif
				elif e in self.tmpItems:
					(val,x) = recCalc(e)
					p += x
				# Aucune info
				else:
					continue
				
				# Si on a une vraie valeur de distance, on continue la moyenne
				if val != None:
					poids = 1./float(p)
					d += val*poids
					s += poids
					nb += 1

			# Test final
			if nb != 0:
				d /= s
				if nb == 1:
					return (d,1./s)
				return (d,0)
			return (None,0)

		return recCalc(init)[0]




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
	

