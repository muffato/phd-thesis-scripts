#! /usr/bin/python2.4

#
# Fonctions communes de chargement et de traitement des donnees
#

import math
import sys

#
# Renvoie les cles d'un dictionnaire, triees suivant les valeurs associees
#
def sortDict(d):
	k = d.keys()
	k.sort(lambda x, y: cmp(d[y], d[x]))
	return k


#
# Renvoie la moyenne d'une liste
#
def moyenne(lst):
	if len(lst) == 0:
		return 0.
	return float(sum(lst))/float(len(lst))

#
# Renvoie la moyenne ponderee d'une liste
#
def moyennePonderee(lst):
	sV = 0.
	sP = 0.
	for (val,poids) in lst:
		sV += val*poids
		sP += poids
	if sP == 0:
		return 0.
	return sV/sP




#
# Ecart type
def ecartType(lst):
	if len(lst) == 0:
		return 0
	m = moyenne(lst)
	return math.sqrt(moyenne([(x-m)*(x-m) for x in lst]))


#
# Min et max
#
def getMinMax(lst):
	if len(lst) == 0:
		return (None, None)
	mn = lst[0]
	mx = mn
	for x in lst:
		if x > mx:
			mx = x
		elif x < mn:
			mn = x
	return (mn, mx)


#
# Renvoie l'intersection de deux fenetres
#
def intersectionCouples(c1, c2):
	if c1[1] < c2[0] or c2[1] < c1[0]:
		return []
	else:
		return [ max(c1[0], c2[0]), min(c1[1], c2[1]) ]

#
# Cette classe permet de garder la valeur minimale qu'on lui a presente
# Elle utilise une fonction pour comparer les objets
#
class MinKeeper:

	def __init__(self, f):
		self.func = f
		self.empty = True
	
	def setBegin(self, r):
		self.res = r
		self.s = self.func(r)
		self.empty = False
	
	def check(self, r):
		if self.empty:
			self.setBegin(r)
			return
		s = self.func(r)
		if s < self.s:
			self.res = r
			self.s = s
			self.empty = False



#
# Renvoie tous les ensembles de n elements de la liste
#
def buildSubsets(lst, n):
	if n < 2:
		return [set([x]) for x in lst]
	ens = buildSubsets(lst, n-1)
	res = []
	for s in ens:
		m = max(s)
		for x in lst:
			if x <= m:
				continue
			ss = s.union([x])
			if len(ss) == n:
				res.append(ss)
	return res


#
# Applatit une liste
#
def flatten(lst):
	res = []
	for x in lst:
		res.extend(x)
	return res

#
# Fait la liste de tous les genes de tab1 en diagonale avec ceux de tab2
#  (eventuellement des singletons)
#
def extractDiags(tab1, dic2, largeurTrou):
	
	#dic2 = {}
	#for i in range(len(tab2)):
	#	x = tab2[i]
	#	if x in dic2:
	#		dic2[x].append(i)
	#	else:
	#		dic2[x] = [i]

	diag = []
	cour = []
	lastI2 = []
	listI2 = []
	for i1 in xrange(len(tab1)):
		g1 = tab1[i1]
		presI2 = dic2.get(g1, [])
		
		for i2 in presI2:
			# Est-ce qu'on est dans le prolongement d'une diagonale
			if ((i2+1) in lastI2) or ((i2-1) in lastI2):
				listI2.append(i2)
				break
		else:
			# On n'a pas trouve de i2 satisfaisant, c'est la fin de la diagonale
			if len(cour) > 0 and len(listI2) > 0:
				diag.append( (cour,(deb1,fin1),getMinMax(listI2)) )
			deb1 = i1
			cour = []
			listI2 = presI2
		cour.append(g1)
		lastI2 = presI2
		fin1 = i1
	
	if len(cour) > 0 and len(listI2) > 0:
		diag.append( (cour,(deb1,fin1),getMinMax(listI2)) )
	elif len(diag) == 0:
		return []

	# On rassemble des diagonales separees par une espace pas trop large
	combiDiag = []
	(cour,(x1,x2),(y1,y2)) = diag.pop()
	while len(diag) > 0:
		(c2,(xx1,xx2),(yy1,yy2)) = diag.pop()
		
		a = max(min(abs(x1-xx2), abs(xx1-x2)), min(abs(y1-yy2), abs(yy1-y2)))
		if a <= (largeurTrou+1):
			cour.extend(c2)
			x1 = min(x1,xx1)
			y1 = max(y1,yy1)
			x2 = min(x2,xx2)
			y2 = max(y2,yy2)
		else:
			combiDiag.append(cour)
			cour = c2
			x1 = xx1
			y1 = yy1
			x2 = xx2
			y2 = yy2
	combiDiag.append(cour)
	return combiDiag


def translateGenome(genome, default, genesAnc):

	# On construit les couleurs
	newGen = {}
	for c in genome.lstChr:

		newChr = []
		for gene in genome.lstGenes[c]:
			
			for g in gene.names:
				if g in genesAnc.dicGenes:
					newChr.append( genesAnc.dicGenes[g][1] )
					break
			else:
				newChr.append(default)
		newGen[c] = newChr
	return newGen


def iterateDiags(genome1, genome2, keepGaps, threshold, callBackFunc):

	#genome1 = translateGenome(g1, "genome1")
	#genome2 = translateGenome(g2, "genome2")

	for c1 in genome1:
	
		tab1 = genome1[c1]
		dic1 = {}
		for i in range(len(tab1)):
			x = tab1[i]
			if x in dic1:
				dic1[x].append(i)
			else:
				dic1[x] = [i]

		for c2 in genome2:
			if not keepGaps:
				s1 = set(genome1[c1])
				s2 = set(genome2[c2])
				t1 = [x for x in genome1[c1] if x in s2]
				t2 = [x for x in genome2[c2] if x in s1]
			else:
				t1 = genome1[c1]
				t2 = genome2[c2]
			diags = extractDiags(t2, dic1, threshold)
			for d in diags:
				callBackFunc(c1, c2, len(diags), d)


def unique(s):
    """Return a list of the elements in s, but without duplicates.

    For example, unique([1,2,3,1,2,3]) is some permutation of [1,2,3],
    unique("abcabc") some permutation of ["a", "b", "c"], and
    unique(([1, 2], [2, 3], [1, 2])) some permutation of
    [[2, 3], [1, 2]].

    For best speed, all sequence elements should be hashable.  Then
    unique() will usually work in linear time.

    If not possible, the sequence elements should enjoy a total
    ordering, and if list(s).sort() doesn't raise TypeError it's
    assumed that they do enjoy a total ordering.  Then unique() will
    usually work in O(N*log2(N)) time.

    If that's not possible either, the sequence elements must support
    equality-testing.  Then unique() will usually work in quadratic
    time.
    """

    n = len(s)
    if n == 0:
        return []

    # Try using a dict first, as that's the fastest and will usually
    # work.  If it doesn't work, it will usually fail quickly, so it
    # usually doesn't cost much to *try* it.  It requires that all the
    # sequence elements be hashable, and support equality comparison.
    u = {}
    try:
        for x in s:
            u[x] = 1
    except TypeError:
        del u  # move on to the next method
    else:
        return u.keys()

    # We can't hash all the elements.  Second fastest is to sort,
    # which brings the equal elements together; then duplicates are
    # easy to weed out in a single pass.
    # NOTE:  Python's list.sort() was designed to be efficient in the
    # presence of many duplicate elements.  This isn't true of all
    # sort functions in all languages or libraries, so this approach
    # is more effective in Python than it may be elsewhere.
    try:
        t = list(s)
        t.sort()
    except TypeError:
        del t  # move on to the next method
    else:
        assert n > 0
        last = t[0]
        lasti = i = 1
        while i < n:
            if t[i] != last:
                t[lasti] = last = t[i]
                lasti += 1
            i += 1
        return t[:lasti]

    # Brute force is all that's left.
    u = []
    for x in s:
        if x not in u:
            u.append(x)
    return u

