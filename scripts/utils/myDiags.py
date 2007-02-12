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


# lstTout est une liste de diagonales chevauchantes
# Renvoie les plus longs chemins de genes ancestraux issus de ces diagonales,
# puis les plus longs chemins avec les genes non encore utilises et ainsi de suite
def getLongestPath(lstTout):

	# A partir d'une liste de diagonales, contruit la liste des voisins de chaque gene
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


	# Construit un graphe dans lequel les suites d'aretes sans carrefours ont ete reduites
	def buildReducedGraph():
		newSommets = set([x for x in voisins if len(voisins[x]) != 2])
		aretes = dict([(x,dict()) for x in newSommets])

		# Renvoie le chemin maximal qui part de s (en venant de pred)
		# qui ne croise pas de carrefours
		def followSommet(s, pred):
			if s in newSommets:
				return [s]
			next = [x for x in voisins[s] if x != pred][0]
			l = followSommet(next, s)
			return [s]+l

		for x in newSommets:
			for v in voisins[x]:
				l = followSommet(v, x)
				if (l[-1] not in aretes[x]) or (len(l) > len(aretes[x][l[-1]])):
					aretes[x][l[-1]] = l

		return (newSommets, aretes)


	# prend une liste de liste
	# renvoie la liste des listes de longueur maximale
	def selectLongest(lst):
		m = -1
		r = []
		for (x,n) in lst:
			if n < m:
				continue
			if n != m:
				m = n
				r = []
			r.append( (x,n) )
		return r
		
	# prend une liste (le chemin de depart)
	# renvoie la liste des chemins maximaux en partant de ce chemin de depart
	# Chaque liste est en fait un couple (liste ds noeuds, poids)
	def recLongestPath( (currPath,currLength) ):
		toto = [ (currPath,currLength) ]
		for j in aretes[currPath[-1]]:
			if j in currPath:
				continue
			tmp = recLongestPath( (currPath+[j], currLength+len(aretes[currPath[-1]][j])) )
			toto = selectLongest(toto + tmp)
			if toto[0][1] == len(voisins):
				return [toto[0]]
		return toto
		
	def doSearchLongestPath():
		todo = [ ([i],0) for i in newSommets ]
		res = []
		max = -1
		while len(todo) > 0:
			(path,weight) = todo.pop()
			for j in aretes[path[-1]]:
				if j not in path:
					todo.append( (path+[j], weight+len(aretes[path[-1]])) )
			if weight > max:
				max = weight
				res = path
		return (res,max)

	def doFloydWarshall():
		chemins = dict([(x,dict([(y,(None,0)) for y in newSommets])) for x in newSommets])
		for x in newSommets:
			for y in newSommets:
				if y in aretes[x]:
					chemins[x][y] = (set([x,y]), len(aretes[x][y]))
		for z in newSommets:
			for x in newSommets:
				for y in newSommets:
					(c1,l1) = chemins[x][z]
					(c2,l2) = chemins[z][y]
					newL = l1 + l2
					if (l1 > 0) and (l2 > 0) and (newL > chemins[x][y][1]) and (len(c1 & c2) <= 1):
						chemins[x][y] = (c1 | c2, newL)
		return chemins


	# 1. On construit le graphe original
	voisins = buildVoisins(lstTout)
	sys.stderr.write('.')

	# 2. On reduit le graphe
	(newSommets, aretes) = buildReducedGraph()
	#print "graph {"
	#for x in aretes:
	#	for y in aretes[x]:
	#		print x, " -- ", y
	#print "}"

	print >> sys.stderr, '%d->%d/%d' % (len(voisins), len(newSommets), len(aretes))
	sys.stderr.write('.')
	
	# 3. On extrait les chemins les plus longs
	res = []
	while len(newSommets) > 0:
		#sys.stderr.write('*')
		#allChemins = doFloydWarshall()
		#sys.stderr.write('*')
		#best = []
		#bestS = 0
		#for x in allChemins:
		#	all2 = allChemins[x]
		#	for y in all2:
		#		if all2[y][1] > bestS:
		#			(best,bestS) = all2[y]
		#new = best
		#sys.stderr.write('*')
		#(long,_) = selectLongest(myMaths.flatten([recLongestPath( ([i],0) ) for i in newSommets]))[0]
		(long,_) = doSearchLongestPath()
		new = []
		for x in long:
			if len(new) == 0:
				new.append(x)
			else:
				new.extend(aretes[new[-1]][x])
		if len(new) < 2:
			break
		res.append(new)
		newSommets.difference_update(long)
		print >> sys.stderr, "%d/%d" % (len(new),len(long))
	
	sys.stderr.write('.')
	
	return res



def extractLongestOverlappingDiags(oldDiags, genesAnc):


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
					if (d[j] not in res) or (d[j+1] not in res):
						continue
					if abs(res.index(d[j])-res.index(d[j+1])) == 1:
						flag = True
						break
				if flag:
					ok.add( (oldDiags[i][0][0],oldDiags[i][0][1]) )
					ok.add( (oldDiags[i][1][0],oldDiags[i][1][1]) )
			#print "%s\t%d\t%s\t%s" % (anc, len(res[0]), " ".join([str(x) for x in res[0]]), " ".join(["%s.%s" % (e,c) for (e,c) in ok]))
			newDiags.append( (len(res), res, list(ok)) )
	return newDiags
	

