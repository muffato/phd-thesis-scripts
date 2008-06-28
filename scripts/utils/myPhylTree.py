
import sys
import myTools

dsi = dict.__setitem__
dgi = dict.__getitem__


##########################################
# Classe mere de tous les types d'arbres #
##########################################
class PhylogeneticTree:
	
	def __init__(self, fichier):
		
		if type(fichier) == tuple:
			print >> sys.stderr, "Creation de l'arbre phylogenetique ...",
			(self.items, self.root, self.officialName) = fichier
		else:
			print >> sys.stderr, "Chargement de l'arbre phylogenique de %s ..." % fichier,
			self.officialName = {}
			self.items = self.newCommonNamesMapperInstance()

			# le nom et l'instance de file
			f = myTools.myOpenFile(fichier, 'r')
			try:
				self.nom = f.name
			except AttributeError:
				self.nom = fichier
			
			f = myTools.firstLineBuffer(f)
			if ';' in f.firstLine:
				self.__loadFromNewick__(f.firstLine)
			else:
				self.__loadFromMyFormat__(f)
			f.close()

	
		# La procedure d'analyse de l'arbre
		def recInitialize(node):
			
			sys.stderr.write(".")

			self.dicLinks.setdefault(node, self.newCommonNamesMapperInstance())
			self.dicLinks.get(node).setdefault(node, [node])
			self.dicParents.setdefault(node, self.newCommonNamesMapperInstance())
			self.dicParents.get(node).setdefault(node, node)
			res = [node]

			if node in self.items:
				self.fileName.setdefault(node, str(node).replace(' ', '_').replace('/', '-'))
				s = []
				ld = []
				self.listAncestr.append(node)
				for (f,l) in self.items.get(node):
					self.parent.setdefault(f, (node,l))
					desc = recInitialize(f)
					s.extend(self.species.get(f))
					# On remonte a node
					res.extend(desc)
					ld.append(desc)
					for e in desc:
						self.dicParents.get(e).setdefault(node, node)
						self.dicParents.get(node).setdefault(e, node)
						self.dicLinks.get(e).setdefault(node, self.dicLinks.get(e).get(f) + [node])
						self.dicLinks.get(node).setdefault(e, [node] + self.dicLinks.get(f).get(e))
				self.species.setdefault(node, frozenset(s))
				# Liens de parente
				for (s1,s2) in myTools.myIterator.tupleOnStrictUpperList(ld):
					for e1 in s1:
						for e2 in s2:
							self.dicParents.get(e1).setdefault(e2, node)
							self.dicParents.get(e2).setdefault(e1, node)
							self.dicLinks.get(e1).setdefault(e2, self.dicLinks.get(e1).get(node) + self.dicLinks.get(node).get(e2)[1:])
							self.dicLinks.get(e2).setdefault(e1, self.dicLinks.get(e2).get(node) + self.dicLinks.get(node).get(e1)[1:])
			else:
				self.fileName.setdefault(node, str(node).replace(' ', '.'))
				self.species.setdefault(node, frozenset([node]))
				self.listSpecies.append(node)

			self.allDescendants.setdefault(node, frozenset(res))
			
			return res

		print >> sys.stderr, "Analyse des donnees ",
		
		# Initialisation des structures
		self.parent = self.newCommonNamesMapperInstance()
		self.species = self.newCommonNamesMapperInstance()
		self.fileName = self.newCommonNamesMapperInstance()
		self.dicLinks = self.newCommonNamesMapperInstance()
		self.dicParents = self.newCommonNamesMapperInstance()
		self.allDescendants = self.newCommonNamesMapperInstance()
		self.listSpecies = []
		self.listAncestr = []

		# Remplissage
		recInitialize(self.root)

		# Structures post-analyse
		self.allNames = self.listAncestr + self.listSpecies
		self.outgroupSpecies = self.newCommonNamesMapperInstance()
		self.indNames = self.newCommonNamesMapperInstance()
		for (i,e) in enumerate(self.allNames):
			self.outgroupSpecies.setdefault(e, frozenset(self.listSpecies).difference(self.species.get(e)))
			self.indNames.setdefault(e, i)
			self.officialName[i] = e
		self.numItems = [[]] * len(self.allNames)
		for a in self.listAncestr:
			self.numItems[self.indNames.get(a)] = [(self.indNames.get(f),l) for (f,l) in self.items.get(a)]
		self.numParent = [None] * len(self.allNames)
		for (e,(par,l)) in self.parent.iteritems():
			self.numParent[self.indNames.get(e)] = (self.indNames.get(par),l)
		# Les noms officiels / courants
		tmp = myTools.defaultdict(list)
		for (curr,off) in self.officialName.iteritems():
			tmp[off].append(curr)
		self.commonNames = self.newCommonNamesMapperInstance()
		self.commonNames.update(tmp)
		self.dicGenes = {}
		self.dicGenomes = self.newCommonNamesMapperInstance()

		print >> sys.stderr, " OK"


	# Renvoie le nom de l'ancetre commun de plusieurs especes
	##########################################################
	def lastCommonAncestor(self, species):
		anc = species[0]
		for e in species[1:]:
			anc = self.dicParents[anc][e]
			if anc == self.root:
				return self.root
		return anc

	# Teste si le fils est bien un fils de son pere
	################################################
	def isChildOf(self, fils, pere):
		return self.dicParents[fils][pere] == pere


	# Calcule les valeurs sur les noeuds de l'arbre
	#  - values represente des valeurs definies (noeuds ou feuilles)
	#  - notdefined est la valeur a renvoyer si pas de resultat
	#  - resultNode indique de quel ancetre on veut le resultat
	#########################################################################
	def calcWeightedValue(self, values, notdefined, resultNode = None):
		
		import numpy
		
		n = len(self.allNames)
		# Par defaut, les resultats seront "notdefined"
		matriceA = numpy.identity(n)
		matriceB = numpy.empty((n,))
		matriceB.fill(notdefined)

		# Construit recursivement la matrice
		def recBuild(ianc, father, length):
			anc = self.allNames[ianc]

			# Appels recursifs: on selectionne les fils OK (en supposant qu'on le soit nous meme)
			items = [(e,p) for (e,p) in self.numItems[ianc] if recBuild(e, ianc, p)]
			# Le pere si disponible
			if father != None:
				items.append( (father,length) )

			if anc in values:
				# Si on a une valeur, on met "x = val"
				matriceB[ianc] = values.get(anc)
				return True

			elif len(items) >= 2:
				# S'il a suffisament de voisins, on ecrit l'equation
				s = 0.
				for (e,p) in items:
					p = 1./max(p,0.00001)
					matriceA[ianc][e] = p
					s += p
				matriceA[ianc][ianc] = -s
				matriceB[ianc] = 0
				return True
			else:
				return False
		
		# Construction de la matrice
		if len(values) == 0:
			return matriceB
			
		rootNode = self.indNames[self.lastCommonAncestor(values.keys())]
		recBuild(rootNode, None, 0)
		# Resolution de l'equation
		try:
			res = numpy.linalg.solve(matriceA, matriceB)
		except numpy.linalg.linalg.LinAlgError:
			# Ne se produit jamais
			return None


		# Cherche la valeur la plus proche de node sachant que notre parcours nous fait venir de origin
		def searchBestValue(node, origin, tripLength):

			# Est-ce que node a une valeur definie
			val = res[node]
			if val != notdefined:
				return (tripLength,node,val)
			
			# Les chemins a parcourir
			test = self.numItems[node][:]
			par = self.numParent[node]
			if par != None:
				test.append(par)
			# Iteration
			best = (1e300,node,notdefined)
			for (e,l) in test:
				# Pour ne pas boucler
				if e == origin:
					continue
				# Coupure
				if best[0] < tripLength+l:
					continue
				tmp = searchBestValue(e, node, tripLength+l)
				if tmp[0] < best[0]:
					best = tmp
			return best

		
		if resultNode == None:
			return res
		resultNode = self.indNames[resultNode]
		best = searchBestValue(resultNode, None, 0)
		return (best[0],self.allNames[best[1]],best[2])


	# Renvoie un dictionnaire qui utilise en interne les noms officiels des taxons
	#  mais que l'on peut utiliser avec les noms communs
	# #############################################################################
	def newCommonNamesMapperInstance(self):
		class commonNamesMapper(dict):

			def __getitem__(d, name):
				if name in self.officialName:
					return dgi(d, self.officialName[name])
				else:
					return dgi(d, name)
			
			def __setitem__(d, name, value):
				if name in self.officialName:
					dsi(d, self.officialName[name], value)
				else:
					dsi(d, name, value)
			
		return commonNamesMapper()
		
	# Charge toutes les especes qui descendent d'un ancetre
	def loadAllSpeciesSince(self, ancestr, template, **args):
		if ancestr in self.species:
			l = self.species[ancestr]
		else:
			l = self.listSpecies
		self.loadSpeciesFromList(l, template, **args)
	
	# Charge toutes les especes outgroup d'un ancetre
	def loadAllSpeciesBefore(self, ancestr, template, **args):
		if ancestr in self.outgroupSpecies:
			l = self.outgroupSpecies[ancestr]
		else:
			l = self.listSpecies
		self.loadSpeciesFromList(l, template, **args)

	# Charge toutes les especes d'une liste
	def loadSpeciesFromList(self, lst, template, storeGenomes = True):

		import myGenomes

		for esp in lst:
			esp = self.officialName[esp]
			g = myGenomes.Genome(template % self.fileName[esp])
			if storeGenomes:
				self.dicGenomes[esp] = g
			for (x,(c,i)) in g.dicGenes.iteritems():
				self.dicGenes[x] = (esp, c, i)
	
	# Renvoie le nombre de genes dans chaque espece pour une famille donnee
	def findFamilyComposition(self, fam):
		
		score = self.newCommonNamesMapperInstance()
		score.update(dict.fromkeys(self.officialName,[]))
		
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
				(par,l) = self.parent[anc]
				self.tmpItems[anc].append( (par,l) )
				self.tmpItems[par].remove( (anc,l) )
				anc = par
		for esp in self.allNames:
			self.tmpItems.setdefault(esp, [])
		self.tmpItems[0] = self.tmpItems[node]

	#
	# Lance le calcul de la moyenne etant donne les valeurs stockees dans values
	#
	def calcDist(self, values, init=0):
		
		it = self.tmpItems
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
			#for (e,p) in self.tmpItems.get(anc, []):
			for (e,p) in it[anc]:
				
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




	# Charge un arbre phylogenetique dans mon format a moi
	#######################################################
	def __loadFromMyFormat__(self, f):
		
		# Procedure de chargement du fichier
		def loadFile():
			lignes = []
			for ligne in f:

				# Un commentaire
				if '#' in ligne:
					continue

				# On enleve le \n final et on coupe suivant les \t
				l = ligne.replace('\n', '').split('\t')
				
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
				fils.append( (tmp, currLine[2]-self.ages.get(tmp)) )
			
			n = currLine[1][0]
			if len(fils) == 0:
				if n[0] == ".":
					currLine[1][0] = n = n[1:]
					self.lstEsp6X.add(n)
				elif n[0] == "*":
					currLine[1][0] = n = n[1:]
					self.lstEsp2X.add(n)
				else:
					self.lstEspFull.add(n)
				
			# Un seul fils, on remonte le noeud
			if len(fils) == 1:
				return fils[0][0]
			# Plusieurs fils, on les enregistre
			elif len(fils) > 1:
				self.items.setdefault(n, fils)

			# Info standard
			self.ages.setdefault(n, currLine[2])
			for s in currLine[1]:
				self.officialName[s] = n
				
			return n
					
		self.ages = self.newCommonNamesMapperInstance()
		self.lstEsp2X = set()
		self.lstEsp6X = set()
		self.lstEspFull = set()
		lignes = loadFile()
		self.root = recLoad(0)
	

	# Arbre au format Newick (parentheses)
	#######################################
	def __loadFromNewick__(self, s):

		import string

		# Lit les nb prochains caracteres de l'arbre
		def readStr(nb):
			ret = s[self.pos:self.pos+nb]
			self.pos += nb
			return ret

		# Enleve les caracteres interdits
		def keepWhile(car):
			x = self.pos
			while s[x] in car:
				x += 1
			return readStr(x-self.pos)

		# Garde les caracteres autorises
		def keepUntil(car):
			x = self.pos
			while s[x] not in car:
				x += 1
			return readStr(x-self.pos)

		# Lit l'arbre sous forme de texte
		def readTree():
			keepWhile(' ')
			if s[self.pos] == '(':
				items = []
				# '(' la premiere fois, puis des ',' jusqu'au ')' final
				while readStr(1) != ')':
					items.append( readTree() )
					keepWhile(' ')
				keepWhile(' ')
				# Le resultat est la liste des fils + le nom
				elt = (items, keepUntil("),:; "))
			else:
				# Le resultat est un nom
				elt = keepUntil("),:; ")
		
			keepWhile(' ')
			# Eventuellement une longueur de branche non nulle
			if s[self.pos] == ':':
				readStr(1) # ':"
				length = float(keepWhile("0123456789.eE-"))
			else:
				length = 0
			keepWhile(' ')
			return (elt,length)

		# Remplit les variables d'un GenericTree
		def storeTree(data):
			(elt,length) = data

			if type(elt) == tuple:
				(fils,nom) = elt
			else:
				fils = []
				nom = elt

			if (nom == '') or (nom in self.officialName):
				nom = "NAME_%d" % self.pos
				self.pos += 1
			self.officialName[nom] = nom

			if len(fils) > 0:
				items = []
				for arbre in fils:
					items.append( storeTree(arbre) )
				self.items.setdefault(nom, items)
			
			self.root = nom
			return (nom,length)

		self.pos = 0
		data = readTree()
		self.pos = 0
		storeTree(data)

