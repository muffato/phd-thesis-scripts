
###################################################
# Fonctions communes de traitement des diagonales #
###################################################

import operator
import itertools
import collections

import enum

import utils.myTools
import utils.myGenomes


OrthosFilterType = enum.Enum('None', 'InCommonAncestor', 'InBothSpecies')


#
# Forme la liste des especes a comparer
# Les ancetres doivent etre prefixes d'un '.' pour etre pris en tant que genomes
#   sinon, ce sont leurs especes descendantes qui sont utilisees
##################################################################################
def getTargets(phylTree, s):
	listSpecies = []
	for x in s.split(','):
		if x[0] != '.':
			listSpecies.extend(phylTree.species[x])
		else:
			listSpecies.append(x[1:])

	# La liste des ancetres edites
	dicLinks = [(e1,e2,set(phylTree.dicLinks[e1][e2][1:-1] + [phylTree.dicParents[e1][e2]])) for (e1,e2) in utils.myTools.myIterator.tupleOnStrictUpperList(listSpecies)]
	tmp = set()
	for (_,_,s) in dicLinks:
		tmp.update(s)
	
	return (listSpecies, tmp)


#
# Extrait toutes les diagonales entre deux genomes (eventuellement des singletons)
# Pour optimiser, on demande 
#   genome1 qui est un dictionnaire qui associe a chaque chromosome 
#     la liste des numeros des genes ancestraux sur ce chromosome
#   dic2 qui associe a un numero de gene ancestral ses positions sur le genome 2
########################################################################################
def iterateDiags(genome1, dic2, sameStrand):

	l1 = []
	l2 = []
	la = []
	lastPos2 = []
	lastS1 = 0
	
	# Parcours du genome 1
	for (i1,(j1,s1)) in enumerate(genome1):
		presI2 = dic2[j1] if j1 >= 0 else []
		
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
				# On demande juste que les deux genes soient cote a cote
				if abs(i2-lastI2) != 1:
					continue
			
			# On a passe tous les tests, c'est OK
			# On ecrit l'orthologue que l'on a choisi pour le coup d'avant (utile en cas de one2many, aucun effet si one2one)
			l2[-1] = (lastI2,lastS2)
			l2.append((i2,s2))
			lastPos2 = [(c2,i2,s2)]
			break

		# On n'a pas trouve de i2 satisfaisant, c'est la fin de la diagonale
		else:
			# On l'enregistre si elle n'est pas vide
			if len(l2) > 0:
				yield (l1,l2,la, lastPos2[0][0])
			# On recommence a zero
			lastPos2 = presI2
			l1 = []
			la = []
			# Pour que les diagonales de longueur 1 soient correctes
			l2 = [presI2[0][1:]] if len(presI2) > 0 else []
		
		l1.append((i1,s1))
		la.append(j1)
		lastS1 = s1
	
	if len(l2) > 0:
		yield (l1,l2,la, lastPos2[0][0])


#
# Proxy de generateur gerant la remise differee d'elements
############################################################
class queueWithBackup:

	def __init__(self, gen):
		self.gen = gen
		self.backup = collections.deque()
		self.todofirst = []
	
	def __iter__(self):
		return self
	
	# Le prochain element renvoye vient soit du buffer, soit du generateur principal
	def next(self):
		if len(self.todofirst) > 0:
			return self.todofirst.pop()
		return self.gen.next()
	
	# L'element reinsere est mis en attente
	def putBack(self, x):
		self.backup.appendleft(x)
	
	# Les elements sauvegardes sont mis dans un buffer de sortie et deviennent prioritaires
	def rewind(self):
		self.todofirst.extend(self.backup)
		self.backup = collections.deque()


#
# Lit les diagonales et les fusionne si elles sont separees par un trou pas trop grand
########################################################################################
def diagMerger(diagGen, sameStrand, largeurTrou):
	
	diagGen = queueWithBackup( (l1, l2, la, c2, l1[0][1]/l2[0][1], (min(l1)[0],max(l1)[0]), (min(l2)[0],max(l2)[0])) for (l1,l2,la,c2) in diagGen )

	# On rassemble des diagonales separees par une espace pas trop large
	for (la1,la2,laa,ca2,sa,(_,fina1),(deba2,fina2)) in diagGen:
		for curr in diagGen:
			(lb1,lb2,lba,cb2,sb,(debb1,finb1),(debb2,finb2)) = curr
			
			# Trou trop grand sur l'espece 1, aucune chance de le continuer
			if debb1 > (fina1+largeurTrou+1):
				diagGen.putBack(curr)
				break
			
			# Chromosomes differents de l'espece 2
			if ca2 != cb2:
				ok = False
			elif sameStrand:
				if sa != sb:
					ok = False
				elif sa > 0:
					ok = (fina2 < debb2 <= (fina2+largeurTrou+1))
				else:
					ok = (finb2 < deba2 <= (finb2+largeurTrou+1))
			else:
				ok = (min(abs(deba2-finb2), abs(debb2-fina2)) <= (largeurTrou+1))
			
			if ok:
				la1.extend(lb1)
				la2.extend(lb2)
				laa.extend(lba)
				fina1 = finb1
				deba2 = min(deba2,debb2)
				fina2 = max(fina2,finb2)
			else:
				diagGen.putBack(curr)

		yield (ca2,la1,la2,laa)
		diagGen.rewind()

#
# Procedure complete de calculs des diagonales a partir de 2 genomes, des orthologues et de certains parametres
################################################################################################################
def calcDiags(g1, g2, orthos, fusionThreshold, sameStrand, orthosFilter):

	# Ecrit les genomes comme suites de numeros de genes ancestraux
	def translateGenome(genome):
		newGenome = {}
		for c in genome.chrList[utils.myGenomes.ContigType.Chromosome] + genome.chrList[utils.myGenomes.ContigType.Scaffold]:
			newGenome[c] = [(orthos.dicGenes[g.names[-1]].index if g.names[-1] in orthos.dicGenes else -1, g.strand) for g in genome.lstGenes[c]]
		return newGenome
	newg1 = translateGenome(g1)
	newg2 = translateGenome(g2)

	# Ne conserve que les genes presents dans les deux genomes
	if orthosFilter == OrthosFilterType.InBothSpecies:
		def usedValues(genome):
			val = set()
			for x in genome.itervalues():
				val.update(i for (i,_) in x)
			return val
		def rewrite(genome, inters):
			for c in genome:
				genome[c] = [(i,s) if i in inters else (-1,s) for (i,s) in genome[c]]
		val1 = usedValues(newg1)
		val2 = usedValues(newg2)
		inters = val1.intersection(val2)
		rewrite(newg1, inters)
		rewrite(newg2, inters)

	# Enleve les genes sans lien
	if orthosFilter != OrthosFilterType.None:
		def filterGenome(genome):
			trans = {}
			for c in genome:
				tmp = [(i,x) for (i,x) in enumerate(genome[c]) if x[0] != -1]
				genome[c] = [x for (_,x) in tmp]
				trans[c] = dict((newi,oldi) for (newi,(oldi,_)) in enumerate(tmp))
			return trans
		trans1 = filterGenome(newg1)
		trans2 = filterGenome(newg2)

	# Pour chaque gene ancestral, ses positions dans le genome 2
	newLoc = [[] for x in xrange(len(orthos.lstGenes[None]))]
	for c in newg2:
		for (i,(ianc,s)) in enumerate(newg2[c]):
			if ianc != -1:
				newLoc[ianc].append( (c,i,s) )

	for c1 in newg1:
		for (c2,d1,d2,da) in diagMerger(iterateDiags(newg1[c1], newLoc, sameStrand), sameStrand, fusionThreshold):
			# Si on a garde uniquement les genes avec des orthologues, il faut revenir aux positions reelles dans le genome
			if orthosFilter != OrthosFilterType.None:
				yield ((c1,[(trans1[c1][i1],s1) for (i1,s1) in d1]), (c2,[(trans2[c2][i2],s2) for (i2,s2) in d2]), da)
			else:
				yield ((c1,d1), (c2,d2), da)

#
# Un ensemble de diagonales que l'on represente comme un graphe ou les noeuds sont les genes
##############################################################################################
class WeightedDiagGraph:

	#
	# Constructeur
	################
	def __init__(self):
		# Les aretes du graphe et les orientations relatives des genes
		def newDicInt():
			return collections.defaultdict(int)
		self.aretes = collections.defaultdict(newDicInt)
	
	
	#
	# Insere une diagonale avec un certain poids
	################################################################
	def addDiag(self, diag, weight=1):
		
		for ((x,sx),(y,sy)) in utils.myTools.myIterator.slidingTuple(diag):
			self.aretes[(x,sx)]
			self.aretes[(x,-sx)]
			self.aretes[(y,sy)]
			self.aretes[(y,-sy)]
			if x != y:
				# On compte pour chaque arete le nombre de fois qu'on l'a vue
				self.aretes[(x,sx)][(y,sy)] += weight
				self.aretes[(y,-sy)][(x,-sx)] += weight
		
	#
	# Affiche le graphe initial
	#############################
	def printIniGraph(self, ):
		print "INIGRAPH %d {" % len(self.aretes)
		for x in self.aretes:
			print "\t", x, "> {%d}" % len(self.aretes[x])
			for y in self.aretes[x]:
				print "\t\t", y, "[%s]" % self.aretes[x][y]
		print "}"


	#
	# Garde successivement les aretes de meilleur poids tant qu'elles n'introduisent pas de carrefour ou de cycle
	##############################################################################################################
	def cleanGraphTopDown(self, minimalWeight):
		allEdges = []
		res = {}
		allSucc = {}
		allPred = {}
		allNodes = set(self.aretes)
		
		for (xsx,l) in self.aretes.iteritems():
			allSucc[xsx] = []
			allPred[xsx] = []
			for (ysy,c) in l.iteritems():
				#if c >= minimalWeight:
				if (c >= minimalWeight) and (xsx < ysy):
					allEdges.append( (c,xsx,ysy) )
		
		def addEdge(c, xsx, ysy):
			# On ecrit l'arete
			res[xsx] = (ysy, c)
			# Tables de detection de cycles
			allSucc[xsx] = [ysy] + allSucc[ysy]
			assert len(allSucc[ysy]) == len(set(allSucc[ysy]))
			for t in allPred[xsx]:
				allSucc[t].extend(allSucc[xsx])
				assert len(allSucc[t]) == len(set(allSucc[t]))
			allPred[ysy] = [xsx] + allPred[xsx]
			assert len(allPred[ysy]) == len(set(allPred[ysy]))
			for t in allSucc[ysy]:
				allPred[t].extend(allPred[ysy])
				assert len(allPred[t]) == len(set(allPred[t]))

		allEdges.sort(reverse = True)
		for (c,(x,sx),(y,sy)) in allEdges:
			# 3 cas ambigus
			if (x,sx) in res:
				# > 1 successeur
				print "not used /successor", (x,sx), (y,sy), c
				continue
			if (y,-sy) in res:
				# > 1 predecesseur
				print "not used /predecessor", (x,sx), (y,sy), c
				continue
			if (x,sx) in allSucc[(y,sy)]:
				# Cycle
				print "not used /cycle", (x,sx), (y,sy), c
				continue
			assert (y,-sy) not in allSucc[(x,-sx)]

			addEdge(c, (x,sx), (y,sy))
			addEdge(c, (y,-sy), (x,-sx))
			
	
		allNodes.difference_update(res)
		allNodes.difference_update(ysy for (ysy,_) in res.itervalues())

		self.aretes = res
		self.singletons = allNodes

		# Affichage du graphe
		print "GRAPH %d {" % len(self.aretes)
		for (tx,(ty,c)) in self.aretes.iteritems():
			print "\t", tx, ">", ty, "[%s]" % c
		print "}"



	#
	# Construit un graphe dans lequel les suites d'aretes sans carrefours ont ete reduites
	#
	def getBestDiags(self):

		# Renvoie le chemin qui part de src en parcourant les aretes
		def followSommet(x, sx):
		
			res = []
			scores = []
			print "begin", (x,sx)
			
			# On part de src et on prend les successeurs jusqu'a la fin du chemin
			while (x,sx) in self.aretes:
				# On marque notre passage
				res.append( (x,sx) )
				
				# Le prochain noeud a visiter
				((y,sy),c) = self.aretes.pop( (x,sx) )
				print "pop edge", (x,sx), (y,sy), c
				scores.append(c)

				# Traitement de l'arete inverse
				assert self.aretes.pop( (y,-sy) ) == ((x,-sx),c)

				# Passage au suivant
				(x,sx) = (y,sy)

			res.append( (x,sx) )
			assert (x,-sx) not in self.aretes
			assert len(res) >= 2

			print "end", len(res), res

			return (res,scores)

		
		# On cherche les extremites pour lancer les blocs integres
		todo = [(x,sx) for (x,sx) in self.aretes if (x,-sx) not in self.aretes]
		for (x,sx) in todo:
			if (x,sx) in self.aretes:
				yield followSommet(x, sx)
			
		assert len(self.aretes) == 0

		# On a isole quelques noeuds
		print self.singletons
		r = set()
		for (x,sx) in self.singletons:
			if (x,sx) not in r:
				r.add( (x,-sx) )

		for (x,sx) in r:
			print "singleton", x
			yield ([(x,1)],[])
			self.singletons.remove( (x,sx) )
			self.singletons.remove( (x,-sx) )

		assert len(self.singletons) == 0



#
# Un ensemble de diagonales que l'on represente comme un graphe ou les noeuds sont les genes
##############################################################################################
class WeightedDiagGraphWithoutStrand:

	#
	# Constructeur
	################
	def __init__(self):
		# Les aretes du graphe et les orientations relatives des genes
		def newDicInt():
			return collections.defaultdict(int)
		self.aretes = collections.defaultdict(newDicInt)
	
	
	#
	# Insere une diagonale avec un certain poids
	################################################################
	def addDiag(self, diag, weight=1):
		
		for (x,y) in utils.myTools.myIterator.slidingTuple(diag):
			self.aretes[x]
			self.aretes[y]
			if x != y:
				# On compte pour chaque arete le nombre de fois qu'on l'a vue
				self.aretes[x][y] += weight
				self.aretes[y][x] += weight
		
	#
	# Affiche le graphe initial
	#############################
	def printIniGraph(self, ):
		print "INIGRAPH %d {" % len(self.aretes)
		for x in self.aretes:
			print "\t", x, "> {%d}" % len(self.aretes[x])
			for y in self.aretes[x]:
				print "\t\t", y, "[%s]" % self.aretes[x][y]
		print "}"


	#
	# Garde successivement les aretes de meilleur poids tant qu'elles n'introduisent pas de carrefour ou de cycle
	##############################################################################################################
	def cleanGraphTopDown(self, minimalWeight):
		allEdges = []
		res = collections.defaultdict(list)
		allSucc = {}
		
		for (x,l) in self.aretes.iteritems():
			allSucc[x] = []
			for (y,c) in l.iteritems():
				if (c >= minimalWeight) and (x < y):
					allEdges.append( (c,x,y) )
		
		def addEdge(c, x, y):
			# On ecrit l'arete
			res[x].append( (y, c) )
			res[y].append( (x, c) )
			
			# Tables de detection de cycles
			lx = allSucc[x] + [x]
			ly = allSucc[y] + [y]
			assert len(set(lx).intersection(ly)) == 0
			for t in lx:
				allSucc[t].extend(ly)
				assert len(allSucc[t]) == len(set(allSucc[t]))
			for t in ly:
				allSucc[t].extend(lx)
				assert len(allSucc[t]) == len(set(allSucc[t]))

		end0 = set(self.aretes)
		end1 = set()
		end2 = set()
		def upgrade(node):
			if node in end0:
				end0.discard(node)
				end1.add(node)
			else:
				assert node in end1
				end1.discard(node)
				end2.add(node)

		allEdges.sort(reverse = True)
		for (c,x,y) in allEdges:
			# 3 cas ambigus
			if x in end2:
				# > 1 successeur
				print "not used /x", x, y, c
				continue
			if y in end2:
				# > 1 predecesseur
				print "not used /y", x, y, c
				continue
			if x in allSucc[y]:
				# Cycle
				print "not used /cycle", x, y, c
				continue
			assert y not in allSucc[x]
			addEdge(c, x, y)
			upgrade(x)
			upgrade(y)
			
		self.aretes = res
		self.singletons = end0

		# Affichage du graphe
		print "GRAPH %d {" % len(self.aretes)
		for (tx,l) in self.aretes.iteritems():
			assert ((len(l) == 1) and (tx in end1)) or ((len(l) == 2) and (tx in end2))
			print "\t", tx, ">", " / ".join("%d [%s]" % ty for ty in l)
		print "}"



	#
	# Construit un graphe dans lequel les suites d'aretes sans carrefours ont ete reduites
	#
	def getBestDiags(self):

		# Renvoie le chemin qui part de src en parcourant les aretes
		def followSommet(x):
		
			res = []
			scores = []
			print "begin", x
			
			# On part de src et on prend les successeurs jusqu'a la fin du chemin
			while x in self.aretes:
				# On marque notre passage
				res.append(x)
				
				# Le prochain noeud a visiter
				assert len(self.aretes[x]) == 1, self.aretes[x]
				[(y,c)] = self.aretes.pop(x)
				print "pop edge", x, y, c
				scores.append(c)

				# Traitement de l'arete inverse
				self.aretes[y].remove( (x,c) )
				if len(self.aretes[y]) == 0:
					del self.aretes[y]

				# Passage au suivant
				x = y

			res.append(x)
			assert x not in self.aretes
			assert len(res) >= 2

			print "end", len(res), res

			return (res,scores)

		
		# On cherche les extremites pour lancer les blocs integres
		todo = [x for x in self.aretes if len(self.aretes[x]) == 1]
		for x in todo:
			if x in self.aretes:
				yield followSommet(x)
			
		assert len(self.aretes) == 0

		# On a isole quelques noeuds
		for x in self.singletons:
			print "singleton", x
			yield ([x],[])


