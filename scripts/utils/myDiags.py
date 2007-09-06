
#
# Fonctions communes de traitement des diagonales
#

import sys
import operator
import myMaths
import myGenomes
import myTools
from collections import defaultdict

#
# Extrait toutes les diagonales entre deux genomes (eventuellement des singletons)
# Pour optimiser, on demande 
#   genome1 qui est un dictionnaire qui associe a chaque chromosome 
#     la liste des numeros des genes ancestraux sur ce chromosome
#   dic2 qui associe a un numero de gene ancestral ses positions sur le genome 2
#

def iterateDiags(genome1, dic2, largeurTrou, sameStrand):

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
		listI1 = []
		listI2 = []
		lastPos2 = []
		listStrand = []
		lastS1 = 0
		
		# Parcours du genome 1
		for (i1,(j1,s1)) in enumerate(genome1[c1]):
			if j1 < 0:
				presI2 = []
			else:
				presI2 = dic2[j1]
			
			# On regarde chaque orthologue du gene
			for ((c2,i2,s2), (lastC2,lastI2,lastS2)) in myTools.myIterator.tupleOnTwoLists(presI2, lastPos2):
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
					diag.append( (listI1,listI2, lastPos2[0][0], listStrand, (deb1,fin1),getMinMaxDiag(listI2)) )
				# On recommence a zero
				deb1 = i1
				lastPos2 = presI2
				listI1 = []
				# Pour que les diagonales de longueur 1 soient correctes
				listI2 = [i2 for (lastC2,i2,_) in presI2[:1]]
			
			listI1.append(i1)
			listStrand.append(s1)
			lastS1 = s1
			fin1 = i1
		
		if len(listI2) > 0:
			diag.append( (listI1,listI2, lastC2, listStrand, (deb1,fin1),getMinMaxDiag(listI2)) )

		# On rassemble des diagonales separees par une espace pas trop large
		while len(diag) > 0:
			#print >> sys.stderr, "on est sur", diag[0]
			(d1,d2,c2,s,(deb1,fin1),(deb2,fin2)) = diag.pop(0)
			i = 0
			while i < len(diag):
				(dd1,dd2,cc2,ss(debb1,finn1),(debb2,finn2)) = diag[i]
				#print >> sys.stderr, "test de", diag[i]

				# Aucune chance de poursuivre la diagonale
				if debb1 > (fin1+largeurTrou+1):
					#print >> sys.stderr, "meme pas la peine"
					break
				if (min(abs(deb2-finn2), abs(debb2-fin2)) <= (largeurTrou+1)) and (c2 == cc2):
					#print >> sys.stderr, "OK"
					d1.extend(dd1)
					d2.extend(dd2)
					s.extend(ss)
					fin1 = finn1
					deb2 = min(deb2,debb2)
					fin2 = max(fin2,finn2)
					del diag[i]
				else:
					#print >> sys.stderr, "non"
					i += 1
			#print >> sys.stderr, "envoi", (d1,d2,c2,(deb1,fin1),(deb2,fin2))
			yield (c1, c2, d1, d2, s)



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
	sys.stderr.write(".")
	(tmp,trans2) = translateGenome(g2, orthos.dicGenes, keepOnlyOrthos)
	sys.stderr.write(".")
	
	newLoc = [[] for x in xrange(len(orthos.lstGenes[None]))]
	for c in g2.lstChr + g2.lstScaff:
		for (i,(ianc,s)) in enumerate(tmp[c]):
			if ianc != -1:
				newLoc[ianc].append( (c,i,s) )
	sys.stderr.write(". ")

	statsDiags = []
	for (c1,c2, d1,d2, s) in iterateDiags(newGen, newLoc, fusionThreshold, sameStrand):

		if len(d1) < minimalLength:
			continue
		
		statsDiags.append(len(d1))

		# Si on a garde uniquement les genes avec des orthologues, il faut revenir aux positions reelles dans le genome
		if keepOnlyOrthos:
			yield ((c1,[trans1[c1][i] for i in d1]), (c2,[trans2[c2][i] for i in d2]), s)
		else:
			yield ((c1,d1), (c2,d2), s)
	
	print >> sys.stderr, myMaths.myStats(statsDiags),





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
		self.aretes = dict([(x,defaultdict(int)) for x in self.sommets])
		for d in lstDiags:
			y = d[0]
			for x in d:
				if x != y:
					# On compte pour chaque arete le nombre de fois qu'on l'a vue
					self.aretes[x][y] += 1
					self.aretes[y][x] += 1
				y = x

	#
	# Construit un graphe dans lequel les suites d'aretes sans carrefours ont ete reduites
	#
	def getBestDiags(self):

		todo = []
		for x in self.sommets:
			# Pour les noeuds qui ont plus de 2 voisins ...
			if len(self.aretes[x]) <= 2:
				continue
			vois = self.aretes[x].items()
			# On trie selon la certitude de chaque lien
			vois.sort(key = operator.itemgetter(1), reverse=True)
			# On ne garde que les deux premiers et les autres seront supprimes
			todo.append( (vois[2][1]-vois[1][1], x, tuple(i for (i,_) in vois[2:])) )
		
		# Les liens a supprimer en premier sont ceux qui ont un plus fort differentiel par rapport aux liens gardes
		todo.sort()

		for (_,x,s) in todo:
			# Au cours du processus, on a peut-etre rendu inutile la coupe
			if len(self.aretes[x]) <= 2:
				continue
			for y in s:
				if y not in self.aretes[x]:
					continue
				del self.aretes[x][y]
				del self.aretes[y][x]
		
		# Renvoie le chemin qui part de s (en venant de pred) tant qu'il ne croise pas de carrefours
		def followSommet(ar, src):
			res = []
			pred = None
			curr = src
			# On part de src et on prend les successeurs jusqu'a la fin du chemin
			while True:
				# On marque notre passage
				res.append(curr)
				try:
					# Le prochain noeud a visiter
					(next,_) = ar[curr].popitem()
					del ar[next][curr]
					pred = curr
					curr = next
				except KeyError:
					return res

		for x in self.sommets:
			if (len(self.aretes[x]) == 1):
				yield followSommet(self.aretes, x)

		for x in self.sommets:
			if (len(self.aretes[x]) == 2):
				yield followSommet(self.aretes, x)


