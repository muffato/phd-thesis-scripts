
#
# Fonctions communes de chargement et de traitement des donnees
#

import sys
import math


# Renvoie la moyenne d'une liste
#################################
def mean(lst):
	if len(lst) == 0:
		return 0.
	return float(sum(lst))/float(len(lst))

# Renvoie la moyenne ponderee d'une liste
##############################
def weightedMean(lst):
	sV = 0.
	sP = 0.
	for (val,poids) in lst:
		sV += val*poids
		sP += poids
	if sP == 0:
		return 0.
	return sV/sP



# Min et max
#############
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

# Renvoie la correlation entre les deux variables
##################################################
def correlation(x, y):
	N = len(x)
	if N != len(y):
		N = min(N, len(y))
		x = x[:N]
		y = y[:N]

	sum_sq_x = 0.
	sum_sq_y = 0.
	sum_coproduct = 0.
	mean_x = x[0]
	mean_y = y[0]
	for i in range(1,N):
		sweep = i / (i + 1.0)
		delta_x = x[i] - mean_x
		delta_y = y[i] - mean_y
		sum_sq_x += delta_x * delta_x * sweep
		sum_sq_y += delta_y * delta_y * sweep
		sum_coproduct += delta_x * delta_y * sweep
		mean_x += delta_x / (i + 1.0)
		mean_y += delta_y / (i + 1.0)
	pop_sd_x = math.sqrt( sum_sq_x / N )
	pop_sd_y = math.sqrt( sum_sq_y / N )
	cov_x_y = sum_coproduct / N
	correlation = cov_x_y / (pop_sd_x * pop_sd_y)
	return correlation


# Calcule la proba dans le cas d'une distribution binomiale
#  On observe l parmi ll, alors qu'on attendait une proportion pi
############################################################################
def binom(pi, l, ll):
	p = pow(pi, l) * pow(1.-pi, ll-l)
	for i in xrange(l):
		p *= float(ll-i)/float(i+1)
	return p


# Calcule la p-value de l'observation de l/ll (attendu: pi)
############################################################
def binomPvalue(pi, l, ll, above):
	# La p-value pour > l
	s = 0
	lastS = 0
	for i in xrange(l,ll):
		s += binom(pi, i+1, ll)
		if s == lastS:
			break
		lastS = s
	# Selon above, on renvoit ...
	if above:
		# pvalue pour >= l
		return binom(pi, l, ll) + s
	else:
		# pvalue pour <= l
		return 1. - s

# Calcule le log de cette meme proba
#################################################
def binomLog(pi, l, ll):
	p = l*math.log10(pi) + (ll-l)*math.log10(1-pi)
	for i in xrange(l):
		p += math.log10(ll-i) - math.log10(i+1)
	return p


########################################################
# Contient les stats classiques d'une serie de nombres #
########################################################
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


####################################
# Generateur de valeurs aleatoires #
####################################
class randomValue:
	
	vonmisesmean = 0.5
	vonmiseskappa = 2.
	import random
	import bisect

	# Renvoie un indice au hasard dans l, compte tenu des valeurs de l, qui indiquent une frequence d'apparition
	#############################################################################################################
	def initBisect(self, l):
		self.l = [0] * (len(l)+1)
		for i in xrange(len(l)):
			self.l[i+1] = self.l[i] + l[i]
		self.max = self.l[-1]

	def getRandomBisectPos(self):
		return self.bisect.bisect_left(self.l, self.random.random() * self.max) - 1


	# Distribution VonMises restreinte a [0,1] et centree sur une moyenne
	######################################################################
	def myVonMises(self):
		# Un nombre entre 0 et 1
		r = self.random.vonmisesvariate(0, self.vonmiseskappa) / (2*math.pi) + .5
		# On decale la distribution vers la moyenne voulue
		x0 = 2*self.vonmisesmean - 1
		if x0 < 0:
			return abs(x0 + r*(1-x0))
		else:
			return (1 - abs(-x0 + r*(1+x0)))




def issublist(l1, l2):

	if len(l1) > len(l2):
		return False
	
	n = len(l1)
	for i in xrange(len(l2)+1-n):
		if l1 == l2[i:i+n]:
			return True
	return False




# Applatit une liste
#####################
def flatten(lst):
	res = []
	for x in lst:
		res.extend(x)
	return res


# Eneleve les elements redondants d'une liste
##############################################
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


