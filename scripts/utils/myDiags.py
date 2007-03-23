#! /users/ldog/muffato/python -OO

#
# Fonctions communes de traitement des diagonales
#

import sys
import myMaths
import myGenomes
import myTools

#
# Extrait toutes les diagonales entre deux genomes (eventuellement des singletons)
# Pour optimiser, on demande 
#   genome1 qui est un dictionnaire qui associe a chaque chromosome 
#     la liste des numeros des genes ancestraux sur ce chromosome
#   dic2 qui associe a un numero de gene ancestral ses positions sur le genome 2
#

def iterateDiags(genome1, dic2, largeurTrou, sameStrand, callBackFunc):

	def getMinMaxDiag(lst):
		a = lst[0]
		b = lst[-1]
		if abs(a-b) < len(lst)-1:
			return myMaths.getMinMax(lst)
		elif a < b:
			return (a,b)
		else:
			return (b,a)

	for c1 in genome1:
		
		diag = []
		lastI2 = []
		lastC2 = []
		lastS2 = []
		listI1 = []
		listI2 = []
		
		# Parcours du genome 1
		for i1 in xrange(len(genome1[c1])):
			(j1,s1) = genome1[c1][i1]
			if j1 < 0:
				presI2 = []
			else:
				presI2 = dic2[j1]

			# On regarde chaque orthologue du gene
			for (c2,i2,s2) in presI2:
				# Est-ce qu'on est dans le prolongement d'une diagonale
				if c2 in lastC2 and (((i2+1) in lastI2) or ((i2-1) in lastI2)):
					# Cas special: la diagonale commence avec un g1 qui a plusieurs orthologues
					# Il faut corriger listI2 et lastI2 avec les bons orthologues
					if len(listI1) == 1 and len(lastI2) > 1:
						i = lastC2.index(c2)
						newCS = currentStrand * lastS2[i]
						if sameStrand and (s1*s2 != newCS):
							continue
						listI2 = [lastI2[i]]
						currentStrand = newCS
					listI2.append(i2)
					lastI2 = [i2]
					lastC2 = [c2]
					break
			else:
				# On n'a pas trouve de i2 satisfaisant, c'est la fin de la diagonale
				if len(listI1) > 0 and len(listI2) > 0:
					diag.append( (listI1,listI2,lastC2[0],(deb1,fin1),getMinMaxDiag(listI2)) )
				deb1 = i1
				listI1 = []
				lastC2 = [c for (c,_,_) in presI2]
				lastI2 = [i for (_,i,_) in presI2]
				lastS2 = [s for (_,_,s) in presI2]
				listI2 = [i for (_,i,_) in presI2[:1]] # Pour que les diagonales de longueur 1 soient correctes
				currentStrand = s1
				if len(presI2) == 1:
					currentStrand *= presI2[0][2]
				
			listI1.append(i1)
			fin1 = i1
		
		if len(listI1) > 0 and len(listI2) > 0:
			diag.append( (listI1,listI2,lastC2[0],(deb1,fin1),getMinMaxDiag(listI2)) )

		# On rassemble des diagonales separees par une espace pas trop large
		while len(diag) > 0:
			#print >> sys.stderr, "on est sur", diag[0]
			(d1,d2,c2,(deb1,fin1),(deb2,fin2)) = diag.pop(0)
			i = 0
			while i < len(diag):
				(dd1,dd2,cc2,(debb1,finn1),(debb2,finn2)) = diag[i]
				#print >> sys.stderr, "test de", diag[i]

				# Aucune chance de poursuivre la diagonale
				if debb1 > (fin1+largeurTrou+1):
					#print >> sys.stderr, "meme pas la peine"
					break
				a = min(abs(deb2-finn2), abs(debb2-fin2))
				if a <= (largeurTrou+1) and c2==cc2:
					#print >> sys.stderr, "OK"
					d1.extend(dd1)
					d2.extend(dd2)
					fin1 = finn1
					deb2 = min(deb2,debb2)
					fin2 = max(fin2,finn2)
					del diag[i]
				else:
					#print >> sys.stderr, "non"
					i += 1
			#print >> sys.stderr, "envoi"
			#yield (d1,d2,c2,(deb1,fin1),(deb2,fin2))
			callBackFunc(c1, c2, d1, d2)


#
# Un ensemble de diagonales que l'on represente comme un graphe ou les noeuds sont les genes
#
class DiagGraph:


	#
	# Initialisation du graphe a partir de la liste des diagonales
	# On construit la liste des genes voisins
	#
	def __init__(self, lstDiags):
		
		# La liste des sommets
		self.sommets = set()
		for d in lstDiags:
			self.sommets.update(d)
		
		# Les aretes du graphe
		self.aretes = dict([(x,dict()) for x in self.sommets])
		for d in lstDiags:
			for i in xrange(len(d)-1):
				x = d[i]
				y = d[i+1]
				self.aretes[x][y] = self.aretes[x].get(y,[]) + [y]
				self.aretes[y][x] = self.aretes[y].get(x,[]) + [x]



	#
	# Construit un graphe dans lequel les suites d'aretes sans carrefours ont ete reduites
	#
	def reduceGraph(self):

		# Renvoie le chemin qui part de s (en venant de pred) tant qu'il ne croise pas de carrefours
		def followSommet(s, pred):
			if s in newSommets:
				return self.aretes[pred][s]
			next = [x for x in self.aretes[s] if x != pred][0]
			return self.aretes[pred][s] + followSommet(next, s)

		# Cree les aretes du nouveau graphe, avec les chemins les plus longs entre les nouveaux sommets
		newSommets = set([x for x in self.sommets if len(self.aretes[x]) != 2])
		newAretes = dict([(x,dict()) for x in newSommets])
		for x in newSommets:
			for v in self.aretes[x]:
				l = followSommet(v, x)
				if len(l) > len(newAretes[x].get(l[-1],[])):
					newAretes[x][l[-1]] = l
		
		self.sommets = newSommets
		self.aretes = newAretes

	
	#
	# Recherche le plus long chemin en les testant tous
	#
	def doSearchLongestPath(self):
		todo = [ [i] for i in self.newSommets ]
		res = []
		max = -1
		while len(todo) > 0:
			path = todo.pop()
			for j in self.aretes[path[-1]]:
				if (j not in path) and (j in self.newSommets):
					todo.append( path + self.aretes[path[-1]][j] )
			if len(path) > max:
				max = len(path)
				res = path
		return res

	
	def doFloydWarshall():
		# Tous les plus longs chemins
		vide = set()
		chemins = dict([(x,dict([(y,vide) for y in self.newSommets])) for x in self.newSommets])
		for x in self.newSommets:
			for y in self.aretes[x]:
				chemins[x][y] = set([x] + self.aretes[x][y])
		for z in self.newSommets:
			for x in self.newSommets:
				for y in self.newSommets:
					c1 = chemins[x][z]
					c2 = chemins[z][y]
					
					# Cycle absorbant
					if len(c1 & c2) > 1:
						continue
					
					new = c1 | c2
					if len(new) > len(chemins[x][y]):
						chemins[x][y] = new
		
		# Le plus long plus long chemin
		best = []
		bestS = 0
		for x in chemins:
			all2 = chemins[x]
			for y in all2:
				if len(all2[y]) > bestS:
					best = all2[y]
					bestS = len(best)
		return best


	def getBestDiags(self):
		backupSommets = self.newSommets.copy()
		while len(self.newSommets) > 0:
			best = self.doFloydWarshall()
			#best = doSearchLongestPath()
			
			if len(best) < 2:
				break
			yield best
			newSommets.difference_update(best)
			print >> sys.stderr, "{%s}" % best
		
		self.newSommets = backupSommets


	def printGraph(self, subset):
		print "graph {"
		for x in self.sommets:
			if x not in subset:
				continue
			for y in self.aretes[x]:
				if y not in subset or y >= x:
					continue
				print "%s -- %s [label=\"%d\"]" % (x,y,len(self.aretes[x][y]))
		print "}"



#
# Un ensemble de diagonales que l'on represente comme un graphe ou les noeuds sont les genes
#
class WeightedDiagGraph:


	#
	# Initialisation du graphe a partir de la liste des diagonales
	# On construit la liste des genes voisins
	#
	def __init__(self, lstDiags):
		
		# La liste des sommets
		self.sommets = set()
		for d in lstDiags:
			self.sommets.update(d)
		
		# Les aretes du graphe
		self.aretes = dict([(x,{}) for x in self.sommets])
		for d in lstDiags:
			for i in xrange(len(d)-1):
				x = d[i]
				y = d[i+1]
				if x != y:
					self.aretes[x][y] = self.aretes[x].get(y, 0) + 1
					self.aretes[y][x] = self.aretes[y].get(x, 0) + 1



	#
	# Construit un graphe dans lequel les suites d'aretes sans carrefours ont ete reduites
	#
	def getBestDiags(self):

		# Pour les noeuds qui ont plus de deux voisins, on ne considere que les meilleurs hits
		for x in self.sommets:
			if len(self.aretes[x]) <= 2:
				continue
			vois = self.aretes[x].keys()
			vois.sort(lambda a, b: cmp(self.aretes[x][b], self.aretes[x][a]))
			print >> sys.stderr, "Reduction", x, " ".join([str(self.aretes[x][y]) for y in vois])
			for y in vois[2:]:
				del self.aretes[x][y]
				del self.aretes[y][x]
			#print >> sys.stderr, "Resultat", self.aretes[x]
		
		# Renvoie le chemin qui part de s (en venant de pred) tant qu'il ne croise pas de carrefours
		def followSommet(s):
			res = []
			pred = None
			# On part de s et on prend les successeurs jusqu'a la fin du chemin
			#print >> sys.stderr, "DEBUT"
			while True:
				# On marque notre passage
				res.append(s)
				alreadySeen.add(s)
				# Les prochains noeuds a visiter
				# On doit verifier qu'on est pas deja passer par les noeuds sinon on boucle indefiniment
				next = [x for x in self.aretes[s] if (x != pred) and (x not in alreadySeen)]
				#print >> sys.stderr, pred, s, self.aretes[s], next
				if len(next) == 0:
					#print >> sys.stderr, "FIN"
					#print >> sys.stderr
					return res
				# N'arrive jamais
				#elif len(next) >= 2:
				#	print >> sys.stderr, "boudiou"
				pred = s
				s = next[0]
				# N'arrive jamais
				#if s in alreadySeen:
				#	print >> sys.stderr, "drame"
				

		alreadySeen = set()
		for x in self.sommets:
			if (len(self.aretes[x]) == 1) and (x not in alreadySeen):
				yield followSommet(x)

		#print >> sys.stderr, "BOUCLES-DEBUT"
		for x in self.sommets:
			if (len(self.aretes[x]) == 2) and (x not in alreadySeen):
				yield followSommet(x)
		#print >> sys.stderr, "BOUCLES-FIN"



