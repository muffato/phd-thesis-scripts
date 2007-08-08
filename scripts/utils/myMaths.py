
#
# Fonctions communes de chargement et de traitement des donnees
#

import sys
import math

#
# Renvoie la moyenne d'une liste
#
def mean(lst):
	if len(lst) == 0:
		return 0.
	return float(sum(lst))/float(len(lst))


#
# Ecart type
#
def stddev(lst):
	if len(lst) == 0:
		return 0
	m = mean(lst)
	return math.sqrt(mean([(x-m)*(x-m) for x in lst]))

#
# Mediane
#
def median(lst):
	if len(lst) == 0:
		return 0
	l = lst[:]
	l.sort()
	return l[len(l)/2]


#
# Min et max
#
def getMinMax(lst):
	mn = mx = None
	for x in lst:
		if mn == None:
			mn = x
			mx = x
		elif x > mx:
			mx = x
		elif x < mn:
			mn = x
	return (mn, mx)


#
# Contient les stats classiques d'une serie de nombres
#
class myStats:

	def __init__(self, lst):

		self.data = list(lst)
		self.data.sort()
		self.update()
		

	# Met a jour la moyenne et l'ecart-type
	def update(self):

		if len(self.data) == 0:
			self.mean = self.stddev = 0
		else:
			s = 0
			for x in self.data:
				s += x
			self.mean = float(s) / float(len(self.data))
			
			s = 0
			for x in self.data:
				s += (x-self.mean)*(x-self.mean)
			self.stddev = math.sqrt(float(s) / float(len(self.data)))
	

	# Renvoie la valeur comprise a x% des donnees (x=50 -> mediane)
	def getValue(self, x):
		if len(self.data) == 0:
			return None
		return self.data[int((x*(len(self.data)-1))/100.)]

	# renvoie la valeur telle que x% soit au dessus (x=50 -> N50)
	def getValueNX(self, x):
		if len(self.data) == 0:
			return None
		tmp = (sum(self.data) * float(x)) / 100.
		for x in self.data.__reversed__():
			tmp -= x
			if tmp <= 0:
				return x


	def __repr__(self):
		# min quart1 median quart3 max mean stddev len
		return "%s [%s/%s/%s] [%s/%s/%s] %s [%.2f/%.2f-%d]" % \
		(self.getValue(0), self.getValue(25),self.getValue(50),self.getValue(75), self.getValueNX(75),self.getValueNX(50),self.getValueNX(25),
			self.getValue(100), self.mean,self.stddev, len(self.data))




def issublist(l1, l2):

	if len(l1) > len(l2):
		return False
	
	n = len(l1)
	for i in xrange(len(l2)+1-n):
		if l1 == l2[i:i+n]:
			return True
	return False


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

