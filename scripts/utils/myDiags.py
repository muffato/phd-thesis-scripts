#! /users/ldog/muffato/python

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
					diag.append( (listI1,listI2,lastC2[0],(deb1,fin1),myMaths.getMinMax(listI2)) )
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
			diag.append( (listI1,listI2,lastC2[0],(deb1,fin1),myMaths.getMinMax(listI2)) )

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
# A partir d'une liste de diagonales, contruit la liste des voisins de chaque gene
#
def buildVoisins(lstDiags):
	voisins = {}
	
	for c in lstDiags:
		if len(c) == 0:
			continue
		if len(c) == 1 and c[0] not in voisins:
			voisins[c[0]] = set([])
		else:
			for i in xrange(len(c)-1):
				x = c[i]
				y = c[i+1]
				if x not in voisins:
					voisins[x] = set([])
				if y not in voisins:
					voisins[y] = set([])
				voisins[x].add(y)
				voisins[y].add(x)
	return voisins


# lstTout est une liste de diagonales chevauchantes
# Renvoie les plus longs chemins de genes ancestraux issus de ces diagonales,
# puis les plus longs chemins avec les genes non encore utilises et ainsi de suite
def getLongestPath(lstTout):

	# prend une liste de liste
	# renvoie la liste des listes de longueur maximale
	def selectLongest(lst):
		if len(lst) == 0:
			return []
		tutu = max([len(x) for x in lst])
		return [x for x in lst if len(x) == tutu]
		
	# prend une liste (le chemin de depart)
	# renvoie la liste des chemins maximaux en partant de ce chemin de depart
	def recLongestPath(currPath):
		toto = [recLongestPath(currPath + [j]) for j in voisins[currPath[-1]] if j not in currPath]
		return selectLongest([currPath] + myMaths.flatten(toto))
		
	ens = set(myMaths.flatten(lstTout))
	voisins = buildVoisins(lstTout)

	res = []
	while len(ens) > 0:
		res.append(selectLongest(myMaths.flatten([recLongestPath([i]) for i in ens])))
		ens.difference_update(res[-1][0])
	return res




#
# Classe censee gerer un ensemble de diagonales en evitant la redondance ...
# encore utile ??
#
class DiagRepository:

	__vide = set([])

	def __init__(self, eliminateSubDiags):
		self.lstDiags = []
		self.lstApp = []
		self.lstDiagsSet = []
		self.genesToDiags = {}
		self.voisins = {}
		self.nbRealDiag = 0
		self.eliminateSubDiags = eliminateSubDiags
	
	def addDiag(self, diag, appar):

		if len(diag) == 0:
			return

		# On doit verifier si on est une sous-diagonale
		# Les diagonales potentielles
		lst = set(myMaths.flatten([self.genesToDiags[x] for x in diag if x in self.genesToDiags]))
		flag = False
		diagS = set(diag)
		if self.eliminateSubDiags:
			for j in lst:
				if self.lstDiagsSet[j].issubset(diagS):
					for x in self.lstDiagsSet[j]:
						self.genesToDiags[x] = [u for u in self.genesToDiags[x] if u != j]
					self.lstDiags[j] = []
					self.lstDiagsSet[j] = self.__vide
					self.lstApp[j] = self.__vide
					self.nbRealDiag -= 1
				elif diagS.issubset(self.lstDiagsSet[j]):
					self.lstApp[j].update(appar)
					flag = True
				
		if not flag:
			n = len(self.lstDiags)
			self.lstDiags.append( diag )
			self.lstDiagsSet.append( diagS )
			self.lstApp.append( set(appar) )
			for x in diagS:
				if x not in self.genesToDiags:
					self.genesToDiags[x] = [n]
				else:
					self.genesToDiags[x].append(n)
			self.nbRealDiag += 1
	
	def __iter__(self):
		for i in xrange(len(self.lstDiags)):
			d = self.lstDiags[i]
			if len(d) == 0:
				continue
			yield (d, self.lstDiagsSet[i], self.lstApp[i])

	def addRepository(self, init):
		for (d,_,a) in init:
			self.addDiag(s, a)

	
	
	def buildVoisins(self):
		self.voisins = buildVoisins(self.lstDiags)
		return
	
	def buildOverlapScores(self):
		
		nbDiags = len(self.lstDiags)
		self.overlapScores = [{} for i in xrange(nbDiags)]
		
		for i in xrange(nbDiags):
			if len(self.lstDiags[i]) == 0:
				continue
			lst = set(myMaths.flatten([self.genesToDiags[x] for x in self.lstDiags[i]]))
			lst.remove(i)
			for j in lst:
				s = self.lstDiagsSet[i].intersection(self.lstDiagsSet[j])
				nb = 0
				for x in s:
					nb += min(self.lstDiags[i].count(x), self.lstDiags[j].count(x))
				self.overlapScores[j][i] = nb
				self.overlapScores[i][j] = nb



	def checkInsert(self):
		
		self.buildVoisins()
		flag = False
		# On cherche des elements qui ont deux voisins eux memes voisins l'un de l'autre
		# L'element est alors a inserer entre les deux pour qu'ils ne soient plus voisins
		for x in self.voisins:
			s = set(self.voisins[x])
			if len(s) != 2:
				continue
			i = s.pop()
			j = s.pop()
			if j not in self.voisins[i]:
				continue
			diags = set(self.genesToDiags[i]).intersection(self.genesToDiags[j])
			
			for ind in diags:
				d = self.lstDiags[ind]
				for k in xrange(len(d)-1):
					if (d[k],d[k+1]) not in [(i,j), (j,i)]:
						continue
					d.insert(k+1, x)
					self.lstDiagsSet[ind].add(x)
					self.genesToDiags[x].append(ind)
					flag = True
					break
		return flag

	
	def __extendExtrem(self, lastPos, beforeLastPos, fun):
		i = 0
		toUpdate = set([])
		while i < len(self.lstDiags):
			curr = self.lstDiags[i]
			if len(curr) >= 2:
				last = curr[lastPos]
				last2 = curr[beforeLastPos]
				v = [y for y in self.voisins[last] if y != last2]
				if len(v) == 1:
					x = v[0]
					if x not in self.lstDiagsSet[i]:
						fun(curr, x)
						self.lstDiagsSet[i].add(x)
						self.genesToDiags[x].append(i)
						toUpdate.add(i)
						continue
			i += 1
		
		for i in toUpdate:
			self.addDiag( self.lstDiags[i], self.lstApp[i] )
		self.buildVoisins()
		return (len(toUpdate) > 0)

	def extendLeft(self):
		return self.__extendExtrem(0, 1, lambda l,x: l.insert(0, x))

	def extendRight(self):
		return self.__extendExtrem(-1, -2, lambda l,x: l.append(x))



	def combinOverlap(self, chev):
	
		self.buildOverlapScores()
		nbDiags = len(self.lstDiags)
		combin = myTools.myCombinator([])
		for i in xrange(nbDiags):
			if len(self.lstDiags[i]) == 0:
				continue
			lst = [i]
			scores = self.overlapScores[i]
			for j in scores:
				if (chev >= 1) and (scores[j] >= chev):
					lst.append(j)
				elif scores[j] >= (min(len(self.lstDiagsSet[i]),len(self.lstDiagsSet[j]))*chev):
					lst.append(j)
			
			if len(lst) == 1:
				continue
			combin.addLink(lst)
		for lst in combin:
			s = set(myMaths.flatten([self.lstDiags[j] for j in lst]))
			a = set(myMaths.flatten([self.lstApp[j] for j in lst]))
			self.addDiag(s, a)
	
	def combinDiags(self, fils, func):
		combin = myTools.myCombinator([])
		for (i,j) in myTools.myMatrixIterator(len(self.lstDiags), len(self.lstDiags), myTools.myMatrixIterator.StrictUpperMatrix):
			commun = set([e for (e,_) in self.lstApp[i].intersection(self.lstApp[j])])
			diff = set([e for (e,_) in self.lstApp[i].symmetric_difference(self.lstApp[j])])
			filsOK = [len(commun.intersection(x)) for x in fils]
			filsNO = [len(diff.intersection(x)) for x in fils]
			outgroupOK = len(commun) - sum(filsOK)
			outgroupNO = len(diff) - sum(filsNO)
			if func(filsOK, filsNO, outgroupOK, outgroupNO):
			#if outgroupOK >= 1 and min(filsOK) >= 1 and max(filsOK) >= 1:
				combin.addLink([i,j])
		res = []
		for g in combin:
			diag = myMaths.unique(myMaths.flatten([self.lstDiags[i] for i in g]))
			orig = myMaths.unique(myMaths.flatten([self.lstApp[i] for i in g]))
			self.addDiag(diag, orig)
			

	def buildCliques(self):

		def recBuild(last):
			newList = set([])
			newListEnd = set([])
			for cl in last:
				# On est dans une clique et on va tenter de l'etendre
				inter = set([])
				beginning = True
				for d in cl:
					s = self.overlapScores[d].keys()
					if beginning:
						beginning = False
						inter.update(s)
					else:
						inter.intersection_update(s)
						if len(inter) == 0:
							break
				if len(inter) == 0:
					continue
			
				# inter contient les voisins de tout le monde
				for newV in inter:
					newCL = list(cl)
					newCL.append(newV)
					newCL.sort()
					newCL = tuple(newCL)
					if len(inter) == 1:
						newListEnd.add(newCL)
					else:
						newList.add(newCL)
			newList.difference_update(newListEnd)
			return (newList,newListEnd)
		
		self.buildOverlapScores()
		
		initCliques = set([])
		for i in xrange(len(self.lstDiags)):
			for j in self.overlapScores[i]:
				if i < j:
					cl = (i,j)
				else:
					cl = (j,i)
				initCliques.add(cl)

		last = initCliques
		last2 = initCliques
		self.cliquesList = [[], [(i,) for i in xrange(len(self.lstDiags)) if len(self.lstDiags[i]) > 0], initCliques]
		while len(last2) > 0:
			(last,last2) = recBuild(last)
			last2.update(last)
			self.cliquesList.append(last2)
		del self.cliquesList[-1]


