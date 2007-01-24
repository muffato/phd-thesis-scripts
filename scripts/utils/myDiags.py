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
# A partir d'une liste de diagonales, contruit la liste des voisins de chaque gene
#
def buildVoisins(lstDiags):
	voisins = {}
	
	for c in lstDiags:
		if len(c) == 0:
			continue
		y = c[0]
		if len(c) == 1 and y not in voisins:
			voisins[y] = set([])
		for i in xrange(1, len(c)):
			x = y
			y = c[i]
			if x not in voisins:
				voisins[x] = set([y])
			else:
				voisins[x].add(y)
			if y not in voisins:
				voisins[y] = set([x])
			else:
				voisins[y].add(x)
	return voisins


def buildReducedGraph(voisins):
	newSommets = [x for x in voisins if len(voisins[x]) != 2]
	for (i1,i2) in myTools.myMatrixIterator(len(newSommets), len(newSommets), myTools.myMaths.WholeWithoutDiag):
		x = newSommets[i1]
		y = newSommets[i2]
		#if x[-1] in voisins[y[0]] or x[-1]



def FloydWarshall(voisins):

	dist = {}
	pred = {}
	for i in voisins:
		for j in voisins:
			pred[(i,j)] = None
			if i in voisins[j]:
				dist[(i,j)] = 1
			else:
				dist[(i,j)] = 1000000

	for k in voisins:
		for i in voisins:
			for j in voisins:
				if dist[(i,j)] < dist[(i,k)]+dist[(k,j)]:
					pass
"""
function fw(int[1..n,1..n] graph) {
    // Initialization
    var int[1..n,1..n] dist := graph
    var int[1..n,1..n] pred
    for i from 1 to n
        for j from 1 to n
            pred[i,j] := nil
            if (dist[i,j] > 0) and (dist[i,j] < Infinity)
                pred[i,j] := i
    // Main loop of the algorithm
    for k from 1 to n
        for i from 1 to n
            for j from 1 to n
                if dist[i,j] > dist[i,k] + dist[k,j]
                    dist[i,j] = dist[i,k] + dist[k,j]
                    pred[i,j] = pred[k,j]
    return ( dist, pred ) // Tuple of the distance and predecessor matrices
}
"""




def extractLongestOverlappingDiags(oldDiags, genesAnc):

	# lstTout est une liste de diagonales chevauchantes
	# Renvoie les plus longs chemins de genes ancestraux issus de ces diagonales,
	# puis les plus longs chemins avec les genes non encore utilises et ainsi de suite
	def getLongestPath(lstTout):


		# prend une liste de liste
		# renvoie la liste des listes de longueur maximale
		def selectLongest(lst):
			m = -1
			r = []
			for x in lst:
				n = len(x)
				if n < m:
					continue
				if n == m:
					r.append(x)
				else:
					m = n
					r = [x]
			return r
			
		# prend une liste (le chemin de depart)
		# renvoie la liste des chemins maximaux en partant de ce chemin de depart
		def recLongestPath(currPath):
			toto = [currPath]
			for j in voisins[currPath[-1]]:
				if j not in currPath:
					toto = selectLongest(toto + recLongestPath(currPath + [j]))
					if len(toto[0]) == len(voisins):
						return [toto[0]]
			return toto
			
		ens = set(myMaths.flatten(lstTout))
		voisins = buildVoisins(lstTout)
		
		#print >> sys.stderr, "longestPath", len(lstTout), len(ens)

		res = []
		while len(ens) > 0:
			res.append( selectLongest(myMaths.flatten([recLongestPath([i]) for i in ens]))[0] )
			ens.difference_update(res[-1])
		return res


	dic = {}
	diags = []
	combin = myTools.myCombinator([])
	for i in xrange(len(oldDiags)):
		d1 = oldDiags[i][0][2]
		d2 = oldDiags[i][1][2]
		da1 = [genesAnc.dicGenes.get(s,("",""))[1] for s in d1]
		da2 = [genesAnc.dicGenes.get(s,("",""))[1] for s in d2]
		if "" in da1:
			diags.append(da2)
		else:
			diags.append(da1)
		for s in d1+d2:
			if s not in dic:
				dic[s] = []
			dic[s].append(i)
		combin.addLink([i])
	
	#combin = myTools.myCombinator([[x] for x in xrange(len(oldDiags))])
	for s in dic:
		combin.addLink(dic[s])

	newDiags = []
	for g in combin:
		for res in getLongestPath([diags[i] for i in g]):
			ok = set([])
			for i in g:
				d = diags[i]
				flag = False
				for j in xrange(len(d)-1):
					if (d[j] not in res[0]) or (d[j+1] not in res[0]):
						continue
					if abs(res[0].index(d[j])-res[0].index(d[j+1])) == 1:
						flag = True
						break
				if flag:
					ok.add( (oldDiags[i][0][0],oldDiags[i][0][1]) )
					ok.add( (oldDiags[i][1][0],oldDiags[i][1][1]) )
			#print "%s\t%d\t%s\t%s" % (anc, len(res[0]), " ".join([str(x) for x in res[0]]), " ".join(["%s.%s" % (e,c) for (e,c) in ok]))
			newDiags.append( (len(res[0]), res[0], list(ok)) )
	return newDiags
	

