
import sys
import myTools
import myGenomes


dgi = dict.__getitem__
dsi = dict.__setitem__

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

				# Un commentaire
				if '#' in ligne:
					continue

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
				fils.append( (tmp, currLine[2]-dgi(self.ages,tmp)) )
				dsi(self.parent, tmp, currLine[1][0])
			
			n = currLine[1][0]
			
			# Un seul fils, on remonte le noeud
			if len(fils) == 1:
				return fils[0][0]
			# Plusieurs fils, on les enregistre
			elif len(fils) > 1:
				dsi(self.items, n, fils)

			# Info standard
			dsi(self.commonNames, n, currLine[1][1:])
			dsi(self.ages, n, currLine[2])
			for s in currLine[1]:
				self.officialName[s] = n
				
			return n
		
		# La procedure d'analyse de l'arbre
		def recInitialize(node):
			if node in self.items:
				dsi(self.outgroupNode, node, None)
				b = []
				bs = []
				s = []
				self.listAncestr.append(node)
				for (f,_) in dgi(self.items, node):
					recInitialize(f)
					b.append(f)
					x = dgi(self.species, f)
					bs.append(x)
					s.extend(x)
				dsi(self.branches, node, b)
				dsi(self.branchesSpecies, node, bs)
				dsi(self.species, node, s)
				for f in b:
					x = b[:]
					x.remove(f)
					dsi(self.outgroupNode, f, x)
			else:
				dsi(self.branchesSpecies, node, [node])
				dsi(self.species, node, [node])
				dsi(self.branches, node, [])
				dsi(self.outgroupNode, node, None)
				self.listSpecies.append(node)
		
		def buildPhylLinks():
			
			esp = self.commonNames.keys()
			nb = len(esp)
			tab = range(nb)
			dic = {}
			tmpPar = range(nb)
			for i in tab:
				dic[esp[i]] = i
			root = dic[self.root]
			for i in tab:
				if i != root:
					tmpPar[i] = dic[dgi(self.parent,esp[i])]

			# Initialisation de la table de tous les liens entre les objets
			dicLinks = [[[] for _ in tab] for _ in tab]
			dicParents = [[None for _ in tab] for _ in tab]
			for i in tab:
				dicLinks[i][i] = [i]
				dicParents[i][i] = i
			
			# Remplissage de chaque objet avec tous ses parents et vice-versa
			for f1 in tab:
				parent = f1
				while parent != root:
					f2 = parent
					parent = tmpPar[f2]
					dicLinks[f1][parent] = dicLinks[f1][f2] + [parent]
					dicLinks[parent][f1] = [parent] + dicLinks[f2][f1]
					dicParents[f1][parent] = dicParents[parent][f1] = parent

			# Liens entre objets de branches differentes
			todo = []
			for f1 in tab:
				for f2 in tab:
					if len(dicLinks[f1][f2]) != 0:
						continue
					for f3 in tab:
						f = False
						try:
							# Pour continuer il faut que les extremites concordent
							if (dicLinks[f1][f3][-1] == dicLinks[f3][f2][0]):
								f = True
								# Mais pas trop
								if (dicLinks[f1][f3][-2] == dicLinks[f3][f2][1]):
									f = False
						except Exception:
							pass
						if f:
							# Pour ne pas que les modifications se perturbent entre elles
							todo.append( (f1,f2,f3) )
							break
			
			for (f1,f2,f3) in todo:
				dicLinks[f1][f2] = dicLinks[f1][f3][:-1]  + dicLinks[f3][f2]
				dicParents[f1][f2] = f3
			
			# Le dictionnaire final
			self.dicLinks = self.newCommonNamesMapperInstance()
			self.dicParents = self.newCommonNamesMapperInstance()
			for i in tab:
				t1 = self.newCommonNamesMapperInstance()
				t2 = self.newCommonNamesMapperInstance()
				dsi(self.dicLinks, esp[i], t1)
				dsi(self.dicParents, esp[i], t2)
				for j in tab:
					dsi(t1, esp[j], [esp[x] for x in dicLinks[i][j]])
					dsi(t2, esp[j], esp[dicParents[i][j]])
				
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
			dsi(self.outgroupSpecies, node, list(set(self.listSpecies).difference(dgi(self.species,node))))
		# Les noms a donner aux fichiers
		self.fileName = self.newCommonNamesMapperInstance()
		for n in self.listAncestr:
			dsi(self.fileName, n, n.replace(' ', '_').replace('/', '-'))
		for n in self.listSpecies:
			dsi(self.fileName, n, n.replace(' ', '.'))

		buildPhylLinks()

		# Les numeros des noms d'ancetres/d'especes
		self.indNames = self.newCommonNamesMapperInstance()
		self.allNames = self.listAncestr + self.listSpecies
		for (i,e) in enumerate(self.allNames):
			self.indNames[e] = i

		print >> sys.stderr, "OK"
		

	#
	# Renvoie un dictionnaire qui utilise en interne les noms officiels des taxons
	#  mais que l'on peut utiliser avec les noms communs
	#
	def newCommonNamesMapperInstance(self):
		class commonNamesMapper(dict):

			def __getitem__(d, name):
				return dgi(d, self.officialName.get(name, name))
			
			def __setitem__(d, name, value):
				return dsi(d, self.officialName.get(name, name), value)
			
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
			g = myGenomes.Genome(template % self.fileName[esp])
			self.dicGenomes[esp] = g
			for (x,(c,i)) in g.dicGenes.iteritems():
				self.dicGenes[x] = (esp, c, i)

	# Renvoie le nombre de genes dans chaque espece pour une famille donnee
	def findFamilyComposition(self, fam):
		
		score = self.newCommonNamesMapperInstance()
		for e in self.officialName:
			score[e] = []
		
		for g in fam:
			if g in self.dicGenes:
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
			
			# Des valeurs d'ancetres / d'especes fixees
			if anc in values:
				return (values[anc],0)
			
			# Initialisation du parcours des fils
			d = 0.
			s = 0.
			nb = 0

			# Moyenne sur tous les fils du noeud
			for (e,p) in self.tmpItems.get(anc, []):
				
				(val,x) = recCalc(e)
				# Si on a une vraie valeur de distance, on continue la moyenne
				if val != None:
					poids = 1./float(p+x)
					d += val*poids
					s += poids
					nb += 1
		
			# Test final
			if nb != 0:
				d /= s
				if nb == 1:
					return (d,1./s)
				return (d,0)
			else:
				return (None,0)

		return recCalc(init)[0]


	#
	# Cree une structure d'arbre secondaire pour le calcul d'une moyenne sur les especes
	# Si on utilise les outgroups, il faut faire basculer l'arbre pour que les outgroups
	#    soient des fils (de fils de fils ...) de plus en plus eloignes
	#
	def initCalcDist2(self, node, useOutgroups):
		tmp = self.items.copy()
		if useOutgroups:
			anc = node
			while anc in self.parent:
				par = self.parent[anc]
				tmp[anc].append( (par,self.ages[par]-self.ages[anc]) )
				tmp[par] = [x for x in tmp[par] if x[0] != anc]
				anc = par
		self.tmpItems = [[(self.indNames[fils],age) for (fils,age) in tmp.get(e,[])] for e in self.allNames]
		self.tmpItems.append(self.tmpItems[self.indNames[node]])

	#
	# Lance le calcul de la moyenne etant donne les valeurs stockees dans values
	#
	def calcDist2(self, values, init=-1):

		# La partie recursive
		def recCalc(anc):
			d = 0.
			s = 0.
			nb = 0
		
			# Moyenne sur tous les fils du noeud
			for (e,p) in self.tmpItems[anc]:
				val = values[e]
				if val == None:
					(val,x) = recCalc(e)
					p += x
				
				# Si on a une vraie valeur de distance, on continue la moyenne
				if val != None:
					d += val/float(p)
					s += 1./p
					nb += 1

			# Test final
			if nb != 0:
				d /= s
				if nb == 1:
					return (d,1./s)
				return (d,0)
			else:
				return (None,0)

		return recCalc(init)[0]

