#! /users/ldog/muffato/python

#
# Fonctions communes de traitement des diagonales
#

import sys
import myMaths
import myGenomes

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
	
	#for d in diag:
	#	#print >> sys.stderr, d
	#	yield d

	#return
	#sys.exit(0)

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

# Ajoute une diagonale a la liste
def addDiag(repos, diag, appar ):

	diags = repos[0]
	dic = repos[1]
	lst = set(myMaths.flatten([dic[x] for x in diag if x in dic]))
	flag = False
	dd = diag[:]
	dd.reverse()
	for j in lst:
		#if myMaths.issublist(diag, diags[j][0]) or myMaths.issublist(dd, diags[j][0]):
		if set(diag).issubset(set(diags[j][0])) or set(dd).issubset(set(diags[j][0])):
			diags[j][1].update(appar)
			flag = True
		#elif myMaths.issublist(diags[j][0], diag) or myMaths.issublist(diags[j][0], dd):
		elif set(diags[j][0]).issubset(set(diag)) or set(diags[j][0]).issubset(set(dd)):
			for x in diags[j][0]:
				dic[x] = [u for u in dic[x] if u != j]
			diags[j] = ([], set([]))
	if not flag:
		n = len(diags)
		diags.append( (diag,set(appar)) )
		for x in diag:
			if x not in dic:
				dic[x] = []
			dic[x].append(n)

def buildVoisins(lst):
	vois = {}
	
	for c in lst:
		if len(c) == 0:
			continue
		if len(c) == 1 and c[0] not in vois:
			vois[c[0]] = []
		else:
			for i in xrange(len(c)-1):
				x = c[i]
				y = c[i+1]
				if x not in vois:
					vois[x] = set([])
				if y not in vois:
					vois[y] = set([])
				vois[x].add(y)
				vois[y].add(x)
	return vois
