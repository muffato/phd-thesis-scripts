
###################################################
# Fonctions communes de traitement des diagonales #
###################################################

import sys
import operator
import itertools
import collections

import myFile
import myMaths
import myTools

slidingTuple = myTools.myIterator.slidingTuple


#
# Extrait toutes les diagonales entre deux genomes (eventuellement des singletons)
# Pour optimiser, on demande 
#   genome1 qui est un dictionnaire qui associe a chaque chromosome 
#     la liste des numeros des genes ancestraux sur ce chromosome
#   dic2 qui associe a un numero de gene ancestral ses positions sur le genome 2
########################################################################################
def iterateDiags(genome1, dic2, sameStrand):

	diag = collections.deque()
	listI1 = []
	listI2 = []
	lastPos2 = []
	listStrand = []
	lastS1 = 0
	
	# Parcours du genome 1
	for (i1,(j1,s1)) in enumerate(genome1):
		if j1 < 0:
			presI2 = []
		else:
			presI2 = dic2[j1]
		
		# On regarde chaque orthologue du gene
		for ((c2,i2,s2), (lastC2,lastI2,lastS2)) in itertools.product(presI2, lastPos2):
			# Chromosomes differents -> indiscutable
			if c2 != lastC2:
				continue
			# Meme brin
			if sameStrand:
				# Les brins initiaux imposent le sens de parcours (+1 ou -1)
				if i2 != lastI2 + lastS1*lastS2:
					continue
				# Le nouveau brin doit etre coherent
				if lastS1*s1 != lastS2*s2:
					continue
			else:
				# On demande juste a ce que les deux genes soient cote a cote
				if abs(i2-lastI2) != 1:
					continue
			
			# On a passe les test, c'est OK
			# On ecrit l'orthologue que l'on a choisi pour le coup d'avant (aucun effet si one2one)
			listI2[-1] = lastI2
			listI2.append(i2)
			lastPos2 = [(c2,i2,s2)]
			break

		# On n'a pas trouve de i2 satisfaisant, c'est la fin de la diagonale
		else:
			# On l'enregistre si elle n'est pas vide
			if len(listI2) > 0:
				yield (listI1,listI2, lastPos2[0][0], listStrand, (deb1,fin1), (min(listI2),max(listI2)) )
			# On recommence a zero
			deb1 = i1
			lastPos2 = presI2
			listI1 = []
			listStrand = []
			# Pour que les diagonales de longueur 1 soient correctes
			listI2 = [i2 for (lastC2,i2,_) in presI2[:1]]
		
		listI1.append(i1)
		listStrand.append(s1)
		lastS1 = s1
		fin1 = i1
	
	if len(listI2) > 0:
		yield (listI1,listI2, lastC2, listStrand, (deb1,fin1), (min(listI2),max(listI2)))


class queueWithBackup:

	def __init__(self, gen):
		self.gen = gen
		self.backup = collections.deque()
		self.todofirst = collections.deque()
	
	def __iter__(self):
		return self
	
	def next(self):
		if len(self.todofirst) > 0:
			return self.todofirst.popleft()
		return self.gen.next()
	
	def putBack(self, x):
		self.backup.append(x)
	
	def rewind(self):
		self.todofirst = self.backup
		self.backup = collections.deque()



def diagMerger(diagGen, largeurTrou):
	# On rassemble des diagonales separees par une espace pas trop large
	for (d1,d2,c2,s,(deb1,fin1),(deb2,fin2)) in diagGen:
		for curr in diagGen:
			(dd1,dd2,cc2,ss,(debb1,finn1),(debb2,finn2)) = curr

			# Aucune chance de poursuivre la diagonale
			if debb1 > (fin1+largeurTrou+1):
				diagGen.putBack(curr)
				break
			elif (min(abs(deb2-finn2), abs(debb2-fin2)) <= (largeurTrou+1)) and (c2 == cc2):
				d1.extend(dd1)
				d2.extend(dd2)
				s.extend(ss)
				fin1 = finn1
				deb2 = min(deb2,debb2)
				fin2 = max(fin2,finn2)
			else:
				diagGen.putBack(curr)

		yield (c2,d1,d2,s)
		diagGen.rewind()

#
# Procedure complete de calculs des diagonales a partir de 2 genomes, des orthologues et de certains parametres
################################################################################################################
def calcDiags(g1, g2, orthos, minimalLength, fusionThreshold, sameStrand, keepOnlyOrthos):


	# Ecrire un genome en suite de genes ancestraux
	def translateGenome(genome, dicOrthos, onlyOrthos):
		newGenome = {}
		transNewOld = {}
		for c in genome.lstChr + genome.lstScaff:
			newGenome[c] = [(dicOrthos.get(g.names[0], (None,-1))[1],g.strand) for g in genome.lstGenes[c]]
			# Si on a garde uniquement les genes avec des orthologues
			# On doit construire un dictionnaire pour revenir aux positions originales
			if onlyOrthos:
				tmp = [x for x in newGenome[c] if x[0] != -1]
				tmpD = {}
				last = -1
				for (i,obj) in enumerate(tmp):
					last = newGenome[c].index(obj, last + 1)
					tmpD[i] = last
				newGenome[c] = tmp
				transNewOld[c] = tmpD
					
		return (newGenome,transNewOld)

	# Les dictionnaires pour accelerer la recherche de diagonales
	(newGen,trans1) = translateGenome(g1, orthos.dicGenes, keepOnlyOrthos)
	(tmp,trans2) = translateGenome(g2, orthos.dicGenes, keepOnlyOrthos)
	
	newLoc = [[] for x in xrange(len(orthos.lstGenes[None]))]
	for c in g2.lstChr + g2.lstScaff:
		for (i,(ianc,s)) in enumerate(tmp[c]):
			if ianc != -1:
				newLoc[ianc].append( (c,i,s) )

	for (c1,tmpGen1) in newGen.iteritems():
		for (c2,d1,d2,s) in diagMerger(queueWithBackup(iterateDiags(tmpGen1, newLoc, sameStrand)), fusionThreshold):

			if len(d1) < minimalLength:
				continue
			
			# Si on a garde uniquement les genes avec des orthologues, il faut revenir aux positions reelles dans le genome
			if keepOnlyOrthos:
				yield ((c1,[trans1[c1][i] for i in d1]), (c2,[trans2[c2][i] for i in d2]), s)
			else:
				yield ((c1,d1), (c2,d2), s)


#
# Un ensemble de diagonales que l'on represente comme un graphe ou les noeuds sont les genes
##############################################################################################
class WeightedDiagGraph:

	#
	# Initialisation du graphe a partir de la liste des diagonales
	# On construit la liste des genes voisins
	################################################################
	def __init__(self, lstDiags):
		
		def newDicInt():
			return collections.defaultdict(int)
		def newDicList():
			def newStrandCount():
				return collections.defaultdict(int)
			return collections.defaultdict(newStrandCount)

		# Les aretes du graphe et les orientations relatives des genes
		self.aretes = collections.defaultdict(newDicInt)
		self.strand = collections.defaultdict(newDicList)
		for d in lstDiags:
			for ((x,sx),(y,sy)) in slidingTuple(d):
				self.aretes[x]
				self.aretes[y]
				if x != y:
					# On compte pour chaque arete le nombre de fois qu'on l'a vue
					self.aretes[x][y] += 1
					self.aretes[y][x] += 1
					# On enregistre l'orientation
					self.strand[x][y][(sx,sy)] += 1
					self.strand[y][x][(-sy,-sx)] += 1
		
		# Pour "normaliser"
		for (x,sx) in self.strand.iteritems():
			for (y,l) in sx.iteritems():
				sx[y] = max( [(v,strand) for (strand,v) in l.iteritems()] )

	#
	# Construit un graphe dans lequel les suites d'aretes sans carrefours ont ete reduites
	#
	def getBestDiags(self):
	
		# Renvoie le chemin qui part de src en parcourant les aretes ar
		def followSommet(src):
			res = []
			pred = None
			curr = src
			# On part de src et on prend les successeurs jusqu'a la fin du chemin
			while True:
				# On marque notre passage
				res.append( curr )
				try:
					# Le prochain noeud a visiter
					(next,_) = self.aretes[curr].popitem()
					# On supprime l'arete inverse pour ne pas revenir en arriere
					del self.aretes[next][curr]
					# L'orientation relative des deux genes
					pred = curr
					curr = next
				# Jusqu'a ce qu'il n'y ait plus de successeur
				except KeyError:
					# On reconstruit l'orientation
					# La liste des co-orientations
					iniStrand = [self.strand[x][y] for (x,y) in slidingTuple(res)]
					# La premiere orientation
					strand = [iniStrand[0][1][0]]
					# On scanne les paires de co-orientations pour voir si elles sont consistantes
					for ((v1,(_,sy1)),(v2,(sx2,_))) in slidingTuple(iniStrand):
						if (sy1 == sx2) or (v1 > v2):
							strand.append( sy1 )
						else:
							strand.append( sx2 )
					# La derniere orientation
					strand.append( iniStrand[-1][1][1] )
					return (res,strand)

		nr = len(self.aretes)
		no = 0
		# 1. On supprime les bifurcations
		todo = []
		for (x,sx) in self.aretes.iteritems():
			# Pour les noeuds qui ont plus de 2 voisins ...
			if len(sx) > 2:
				vois = sx.items()
				# On trie selon la certitude de chaque lien
				vois.sort(key = operator.itemgetter(1))
				# On ne garde que les deux premiers et les autres seront supprimes
				todo.append( (vois[-3][1]-vois[-2][1], x, tuple(i for (i,_) in vois[:-2])) )
		
		# Les liens a supprimer en premier sont ceux qui ont un plus fort differentiel par rapport aux liens gardes
		todo.sort()
		for (_,x,s) in todo:
			# Au cours du processus, on a peut-etre rendu inutile la coupe
			# On verifie qu'il s'agit toujours d'une bifurcation pour x
			if len(self.aretes[x]) <= 2:
				continue
			for y in s:
				# Que y est toujours lie a x
				if y not in self.aretes[x]:
					continue
				del self.aretes[x][y]
				del self.aretes[y][x]

		
		# On a isole quelques noeuds
		for (x,sx) in self.aretes.iteritems():
			if len(sx) == 0:
				yield ([x], [1])
				no += 1

		# Maintenant, on n'a plus de bifurcations
		for (x,sx) in self.aretes.iteritems():
			if len(sx) == 0:
				pass
			elif len(sx) == 1:
				# On peut suivre les chemins non cycliques	
				t = followSommet(x)
				yield t
				no += len(t[0])
				#yield followSommet(x)
				assert len(sx) == 0
			else:
				# Les autres attendent
				assert len(sx) == 2

		# 2. On coupe les boucles en chemins
		todo = []
		for (x,sx) in self.aretes.iteritems():
			if len(sx) == 2:
				for (y,val) in sx.iteritems():
					todo.append( (val,x,y) )
			else:
				# Ici, on est deja passe par le noeud
				assert len(sx) == 0
		
		# Les liens a supprimer en premier sont les plus faibles
		todo.sort()
		for (_,x,y) in todo:
			try:
				del self.aretes[x][y]
				del self.aretes[y][x]
				#yield followSommet(x)
				t = followSommet(x)
				yield t
				no += len(t[0])
			except KeyError:
				pass
	
		assert no == nr
		assert (len(self.aretes) == 0) or (max(len(sx) for sx in self.aretes.itervalues()) == 0)


####################################
# Charge le fichier des diagonales #
####################################
def loadDiagsFile(nom, ancName, officialName):
	
	print >> sys.stderr, "Chargement des diagonales de %s ..." % nom,
	lst = collections.defaultdict(list)
	for (anc,_,diag,strands,e1,e2) in myFile.myTSV.readTabular(nom, [str,int,str,str,str,str]):
		if anc not in ancName:
			continue
		
		# La diagonale
		d = [int(x) for x in diag.split(' ')]
		
		# On joint les especes qui ont vu la diagonale et celles qui n'apportent que le chromosome
		tmp = [y.split("/")[:2] for y in (e1+'|'+e2).split("|") if len(y) > 0]
		# Les chromosomes de ces especes
		espChr = frozenset( [(officialName[e],c) for (e,c) in tmp if ('Un' not in c) and (e in officialName)] )
		esp = frozenset( [e for (e,_) in espChr] )
		
		lst[anc].append( (d,espChr,esp,diag,strands) )

	print >> sys.stderr, "OK (%d diagonales)" % sum([len(lst[x]) for x in ancName])
	return lst


