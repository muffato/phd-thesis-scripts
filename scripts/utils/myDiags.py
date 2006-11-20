#! /users/ldog/muffato/python

#
# Fonctions communes de traitement des diagonales
#

import sys
import myMaths
import myGenomes
import myTools

#
# Fait la liste de tous les genes de tab1 en diagonale avec ceux de tab2
#  (eventuellement des singletons)
# Pour optimiser, on demande 
#   tab1 qui est la liste des numeros des genes ancestraux du chromosome 1 (cf translateGenome)
#   dic2 qui associe a un numero de gene ancestral ses positions sur le genome 2
#

def __extractDiags(tab1, dic2, largeurTrou, sameStrand):
	
	diag = []
	lastI2 = []
	lastC2 = []
	lastS2 = []
	listI1 = []
	listI2 = []
	
	# Parcours du genome 1
	for i1 in xrange(len(tab1)):
		(j1,s1) = tab1[i1]
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
		yield (d1,d2,c2,(deb1,fin1),(deb2,fin2))


#
# Convertit un genome en liste de genes ancestraux 
#
def translateGenome(genome, genesAnc):
	newGenome = {}
	for c in genome.lstChr:
		newGenome[c] = [(genesAnc.dicGenes.get(g.names[0], (0,-1))[1],g.strand) for g in genome.lstGenes[c]]

	return newGenome


def iterateDiags(genome1, dic2, threshold, sameStrand, callBackFunc):

	for c1 in genome1:
		for (d1,d2,c2,aa,bb) in __extractDiags(genome1[c1], dic2, threshold, sameStrand):
			callBackFunc(c1,c2,d1,d2)

class DiagRepository:

	__vide = set([])

	def __init__(self):
		self.lstDiags = []
		self.lstApp = []
		self.lstDiagsSet = []
		self.genesToDiags = {}
		self.voisins = {}
	
	def addDiag(self, diag, appar):

		# On doit verifier si on est une sous-diagonale
		# Les diagonales potentielles
		lst = set(myMaths.flatten([self.genesToDiags[x] for x in diag if x in self.genesToDiags]))
		flag = False
		diagS = set(diag)
		#diagR = diag[:]
		#diagR.reverse()
		for j in lst:
			#elif myMaths.issublist(diags[j][0], diag) or myMaths.issublist(diags[j][0], diagS):
			if self.lstDiagsSet[j].issubset(diagS):
				for x in self.lstDiagsSet[j]:
					self.genesToDiags[x] = [u for u in self.genesToDiags[x] if u != j]
				self.lstDiags[j] = []
				self.lstDiagsSet[j] = self.__vide
				self.lstApp[j] = self.__vide
			#if myMaths.issublist(diag, diags[j][0]) or myMaths.issublist(diagR, diags[j][0]):
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

	def nbRealDiags(self):
		nb = 0
		for d in self.lstDiags:
			if len(d) > 0:
				nb += 1
		return nb
	
	def __iter__(self):
		for i in xrange(len(self.lstDiags)):
			d = self.lstDiags[i]
			if len(d) == 0:
				continue
			yield (d, self.lstDiagSet[i], self.lstApp[i])

	def buildVoisins(self):
		self.voisins = {}
		
		for c in self.lstDiags:
			if len(c) == 0:
				continue
			if len(c) == 1 and c[0] not in self.voisins:
				self.voisins[c[0]] = set([])
			else:
				for i in xrange(len(c)-1):
					x = c[i]
					y = c[i+1]
					if x not in self.voisins:
						self.voisins[x] = set([])
					if y not in self.voisins:
						self.voisins[y] = set([])
					self.voisins[x].add(y)
					self.voisins[y].add(x)
	
	def checkInsert(self):
		
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
			#print "insertion de", x, "entre", i, "et", j, "(", self.voisins[x], self.voisins[i], diags, ")"
			
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


	def buildOverlap(self, chev):
		combin = myTools.myCombinator([])
		newDiags = DiagRepository()
		for (i,j) in myTools.myMatrixIterator(len(self.lstDiags), len(self.lstDiags), myTools.myMatrixIterator.StrictUpperMatrix):
			s = self.lstDiagsSet[i].intersection(self.lstDiagsSet[j])
			if len(s) < 0:
				continue
			nb = 0
			for x in s:
				nb += min(self.lstDiags[i].count(x), self.lstDiags[j].count(x))
			if nb >= chev:
				combin.addLink([i,j])
		for g in combin:
			s = set(myMaths.flatten([self.lstDiags[i] for i in g]))
			a = set(myMaths.flatten([self.lstApp[i] for i in g]))
			newDiags.addDiag(s, a)

		return newDiags
	
	def buildOverlap2(self, chev):
		nbDiags = len(self.lstDiags)
		for i in xrange(nbDiags):
			if len(self.lstDiags[i]) == 0:
				continue
			lst = [i]
			for j in xrange(i+1,nbDiags):
				s = self.lstDiagsSet[i].intersection(self.lstDiagsSet[j])
				if len(s) == 0:
					continue
				nb = 0
				for x in s:
					nb += min(self.lstDiags[i].count(x), self.lstDiags[j].count(x))
				if nb >= chev:
					lst.append(j)
			
			if len(lst) == 1:
				continue
			s = set(myMaths.flatten([self.lstDiags[j] for j in lst]))
			a = set(myMaths.flatten([self.lstApp[j] for j in lst]))
			self.addDiag(s, a)
	
def combinDiags(anc, diags):
	combin = utils.myTools.myCombinator([])
	fils = phylTree.getSpecies(anc)
	for (i,j) in utils.myTools.myMatrixIterator(len(diags), len(diags), utils.myTools.myMatrixIterator.StrictUpperMatrix):
		combin.addLink([i])
		(_,orig1) = diags[i]
		(_,orig2) = diags[j]
		commun = [x for x in orig1.intersection(orig2)]
		filsOK = 0
		outgroupOK = 0
		for (k,l) in utils.myTools.myMatrixIterator(len(commun), len(commun), utils.myTools.myMatrixIterator.StrictUpperMatrix):
			(e1,c1) = commun[k]
			(e2,c2) = commun[l]
			if phylTree.getFirstParent(e1,e2) == anc:
				filsOK += 1
			if ((e1 in fils) and (e2 not in fils)) or ((e2 in fils) and (e1 not in fils)):
				outgroupOK += 1
		if (outgroupOK > 0) and (filsOK > 0):
			combin.addLink([i,j])
	res = []
	for g in combin:
		diag = utils.myMaths.unique(utils.myMaths.flatten([diags[i][0] for i in g]))
		orig = utils.myMaths.unique(utils.myMaths.flatten([diags[i][1] for i in g]))
		res.append( (diag,orig) )
	return res



