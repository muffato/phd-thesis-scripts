
import math

# Calcule la proba dans le cas d'une distribution binomiale
#  On observe l parmi ll, alors qu'on attendait une proportion pi
############################################################################
def binom(pi, l, ll):
	#print pi, l, ll
	p = pow(pi, l) * pow(1.-pi, ll-l)
	#print pow(pi, l), 1.-pi, ll-l, pow(1.-pi, ll-l)
	for i in xrange(l):
		#print p
		p *= float(ll-i)/float(i+1)
	return p


# Calcule la p-value de l'observation de l/ll (attendu: pi)
############################################################
def binomPvalue(pi, l, ll, above):
	# La p-value pour > l
	s = 0
	lastS = 0
	for i in xrange(l,ll):
		x = binom(pi, i+1, ll)
		s += x
		#print pi, i+1, ll, "=", x
		#print lastS, ">", s
		#s += binom(pi, i+1, ll)
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


# Calcule la proba dans le cas d'une distribution binomiale
#  On observe l parmi ll, alors qu'on attendait une proportion pi
############################################################################
def binomF(pi, l, ll):
	#p = pow(pi, l) * pow(1-pi, ll-l)
	num = pi.numerator**l * (pi.denominator-pi.numerator)**(ll-l)
	den = pi.denominator**l
	import fractions
	for i in xrange(l):
		#print i
		num *= (ll-i)
		den *= (i+1)
		#p *= fractions.Fraction(ll-i, i+1)
		#p *= (ll-i)
		#p /= (i+1)
	#print "frac", pi, l, ll
	return fractions.Fraction(num, den)


# Calcule la p-value de l'observation de l/ll (attendu: pi)
############################################################
def binomPvalueF(pi, l, ll, above):
	# La p-value pour > l
	ref = binomF(pi, l, ll)
	tmp = pi * (1 - pi)
	import fractions
	s = ref
	for i in xrange(l+1,ll):
		print "add", i
		ref *= (fractions.Fraction(ll-(i-1), i ) * tmp)
		#ref *= ((ll-(i-1))*pi)
		#ref /= (i*(1-pi))
		s += ref
		#s += binomF(pi, i+1, ll)
	# Selon above, on renvoit ...
	return s
	if above:
		# pvalue pour >= l
		return binomF(pi, l, ll) + s
	else:
		# pvalue pour <= l
		return 1 - s


# Calcule le log de cette meme proba
#################################################
def binomLogF(pi, l, ll):
	x = binomF(pi, l, ll)
	return math.log10(x.numerator) - math.log10(x.denominator)


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


	# Renvoie la quantite de chaque valeur dans la liste
	#####################################################
	@staticmethod
	def count(l):
		import collections
		d = collections.defaultdict(int)
		for x in l:
			d[x] += 1
		return d


######################################################
# Genere les permutations de [1..n]                  #
# Chaque permutation est referee par un index unique #
######################################################
class permutationGenerator:

	def __init__(self, l):

		self.l = l
		self.n = len(l)
		fact = [1] * self.n
		for i in xrange(self.n,1,-1):
			fact[i-2] = i * fact[i-1]
		self.fact = fact
		self.a = [None] * self.n

	def getNbPerm(self, p):
		return self.fact[self.n-p-1]

	def getPermutation(self, k, p):

		assert k >= 0
		assert k < self.fact[self.n-p-1]

		for i in xrange(1,self.n):
			(self.a[i],k) = divmod(k, self.fact[i])

		b = self.l[:]
		for i in xrange(1,self.n):
			(b[i],b[self.a[i]]) = (b[self.a[i]],b[i])

		return b[-p:]


#############################
# Fonctions d'interpolation #
#############################
class myInterpolator:

	########################################
	# Interpolation lineaire (1 dimension) #
	########################################
	@staticmethod
	def oneDimLinear(points, val):
		
		intervals = []
		tmp = zip(points, val)
		for ((u,FU), (v,FV)) in zip(tmp[:-1], tmp[1:]):
			A = (FU - FV) / (u - v)
			B = FU - A * u
			intervals.append((A, B))
		
		import bisect
		def tmp(x):
			i = bisect.bisect_right(points, x)
			if i >= len(points):
				i = len(points)-1
			(A, B) = intervals[i-1]
			return (A * x + B)
		
		return tmp


	#######################################
	# Interpolation cubique (1 dimension) #
	#######################################
	@staticmethod
	def oneDimCubic(points, val):
		
		intervals = []
		tmp = zip(points, val)
		# Generate a cubic spline for each interpolation interval.
		for ((u,FU), (v,FV)) in zip(tmp[:-1], tmp[1:]):
			denom = (u - v)**3
			A = (2 * FV - 2 * FU) / denom
			B = -(u * (3 * FV - 3 * FU) + 3 * FV * v - 3 * FU * v) / denom
			C = (u * (6 * FV * v - 6 * FU * v)) / denom
			D = -(u *(-3 * FU * v**2) + FU * v**3 + u**2 * (3 * FV * v)) / denom
			D = -(u *(-3 * FU * v**2) + FU * v**3 + u**2 * (3 * FV * v) - u**3 * FV) / denom
			intervals.append((A, B, C, D))

		import bisect
		def tmp(x):
			i = bisect.bisect_right(points, x)
			if i >= len(points):
				i = len(points)-1
			(A, B, C, D) = intervals[i-1]
			return ((A * x + B) * x + C) * x + D

		def c_code():
			def codeChoice(intervalList):
				n = len(intervalList)
				if n < 2:
					return ("A=%.10e;B=%.10e;C=%.10e;D=%.10e;" % intervalList[0])
				n2 = n / 2
				return ("if (x < %.10e) {%s} else {%s}" % (points[n2], codeChoice(intervalList[:n2]), codeChoice(intervalList[n2:])))
			return ("double interpolator(double x) {double A,B,C,D;" + codeChoice(intervals) + "return ((A * x + B) * x + C) * x + D;}")
		
		return tmp

	
	###############################################
	# Cree un interpolateur pour chaque dimension #
	###############################################
	@staticmethod
	def getMultDim(interpolator, points, val):
		linter = [interpolator(points, [x[i] for x in val]) for i in xrange(len(val[0]))]
		return lambda x: tuple(inter(x) for inter in linter)


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


##########################################################################
# Renvoie la racine carre d'un entier, sous forme de chaine de caractere #
##########################################################################
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

