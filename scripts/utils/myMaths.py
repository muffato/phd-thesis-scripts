
#
# Fonctions communes de chargement et de traitement des donnees
#

import sys
import math

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


	# Renvoie la moyenne d'une liste
	#################################
	@staticmethod
	def mean(lst):
		if len(lst) == 0:
			return 0.
		return float(sum(lst))/float(len(lst))

	# Renvoie l'ecart type
	#######################
	@staticmethod
	def stddev(lst, m = None):
		if len(lst) == 0:
			return 0.
		if m == None:
			m = myStats.mean(lst)
		s = 0
		for x in lst:
			s += (x-m) ** 2
		return math.sqrt(float(s) / float(len(lst)))

	# Renvoie la valeur comprise a x% des donnees (x=50 -> mediane)
	################################################################
	@staticmethod
	def getValue(lst, x):
		if len(lst) == 0:
			return None
		return lst[int((x*(len(lst)-1))/100.)]

	# Renvoie la valeur telle que x% soit au dessus (x=50 -> N50)
	#############################################################
	@staticmethod
	def getValueNX(lst, x):
		if len(lst) == 0:
			return None
		tmp = (sum(lst) * float(x)) / 100.
		for x in lst.__reversed__():
			tmp -= x
			if tmp <= 0:
				return x

	
	# Renvoie (min quart1 median quart3 max mean stddev len)
	#########################################################
	@staticmethod
	def valSummary(lst):
		l = list(lst)
		l.sort()
		m = myStats.mean(l)
		return (myStats.getValue(l, 0), myStats.getValue(l, 25),myStats.getValue(l, 50),myStats.getValue(l, 75), \
			myStats.getValueNX(l, 75),myStats.getValueNX(l, 50),myStats.getValueNX(l, 25), myStats.getValue(l, 100), m, myStats.stddev(l, m), len(l))

	@staticmethod
	def txtSummary(lst):
		return "%s [%s/%s/%s] [%s/%s/%s] %s [%.2f/%.2f-%d]" % myStats.valSummary(lst)


	# Renvoie la moyenne ponderee d'une liste
	##############################
	@staticmethod
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
	@staticmethod
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
	@staticmethod
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


def gcd(a, b):
	'''Greatest common divisor function; Euclid's algorithm.
		[ a and b are integers ->
		return the greatest common divisor of a and b ]
	'''
	while b != 0:
		(a, b) = (b, a%b)
	return a

class Rational:

	'''rational.py:  Module to do rational arithmetic.

	  For full documentation, see http://www.nmt.edu/tcc/help/lang/python/examples/rational/.
	  Exports:
		Rational ( a, b ):
		  [ (a is a nonnegative integer) and (b is a positive integer)
			-> return a new Rational instance with numerator a and denominator b ]
		.n:    [ the numerator ]
		.d:    [ the denominator ]   
	'''

	def __init__ ( self, a, b ):
		"""Constructor for Rational.
		"""
		if  b == 0:
			raise ZeroDivisionError, ( "Denominator of a rational may not be zero." )
		else:
			g  =  gcd ( a, b )
			self.n  =  a / g
			self.d  =  b / g
	def __add__ ( self, other ):
		"""Add two rational numbers.
		"""
		return Rational ( self.n * other.d + other.n * self.d, self.d * other.d )
	def __sub__ ( self, other ):
		"""Return self minus other.
		"""
		return Rational ( self.n * other.d - other.n * self.d, self.d * other.d )
	def __mul__ ( self, other ):
		"""Implement multiplication.
		"""
		return  Rational ( self.n * other.n, self.d * other.d )
	def __div__ ( self, other ):
		"""Implement division.
		"""
		return  Rational ( self.n * other.d, self.d * other.n )

	def __str__ ( self ):
		''' Return a string representation of self '''
		return "%d/%d" % ( self.n, self.d )
	def __repr__ ( self ):
		''' Return a string representation of self '''
		return "%d/%d" % ( self.n, self.d )
	def mixed ( self ):
		""" Render self as a mixed fraction in string form. """
		whole, n2  =  divmod ( self.n, self.d )
		if  whole == 0:
			if  n2 == 0:  return "0"
			else:         return ("%s/%s" % (n2, self.d) )
		else:
			if  n2 == 0:  return str(whole)
			else:         return ("%s and %s/%s" % (whole, n2, self.d) )

	def __hash__ ( self):
		return hash( (self.n,self.d) )
	def __eq__ ( self, other):
		return (self.n == other.n) and (self.d == other.d)

	def __float__ ( self ):
		""" Implement the float() conversion function. """
		return  float ( self.n ) / float ( self.d )


def sqrti(n, part, nbdec):
	if (len(n) % 2) == 1:
		n = "0" + n
	backup = len(n) / 2
	iniL = len(n) + len(part)
	s = n + (part + ("00" * nbdec))[:2*nbdec]
	maxL = len(s)
	c = 0
	p = 0
	i = 0
	while ((c != 0) or (i < iniL)) and (i<maxL):
		c = 100*c + int(s[i:i+2])
		i += 2
		for x in [9,8,7,6,5,4,3,2,1,0]:
			y = (20*p+x)*x
			if y <= c:
				break
		p = 10*p+x
		c -= y
	sp = str(p)
	return (sp[:backup], sp[backup:])

