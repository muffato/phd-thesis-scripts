#! /usr/bin/python2.4


import sys
import os


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
def extractDiags(t1, t2, byCouple = False, largeurTrou = 0):
	diag = []
	cour = []
	lastI2 = -2
	
	# Permet de ne considerer que la case seule, et donc de supprimer les
	# trous induits par de genes qui ont saute sur d'autre chromosomes
	if byCouple:
		s = [u for u in set(t1) & set(t2) if t1.count(u) == 1 and t2.count(u) == 1]
		(tab1,tab2) = ([x for x in t1 if x in s], [x for x in t2 if x in s])
	else:
		(tab1,tab2) = (t1, t2)
	
	# Recherche simple de diagonale
	for g1 in tab1:
		if g1 not in tab2:
			continue
		i2 = tab2.index(g1)
		# Le point d'apres doit se trouver a une distance de 1
		if abs(i2-lastI2) > 1:
			if len(cour) > 0:
				diag.append(cour)
			cour = []
		lastI2 = i2
		cour.append(g1)
	if len(cour) > 0:
		diag.append(cour)

	if len(diag) == 0:
		return []
		
	# On rassemble des diagonales separees par une espace pas trop large
	diag2 = []
	cour = diag.pop()
	while len(diag) > 0:
		c2 = diag.pop()
		if abs(tab1.index(cour[-1]) - tab1.index(c2[0])) <= (largeurTrou+1):
			cour.extend(c2)
		else:
			diag2.append(cour)
			cour = c2
	diag2.append(cour)
	return diag2



#
# Renvoie les cles d'un dictionnaire, triees suivant les valeurs associees
#
def sort_dict(d):
	k = d.keys()
	k.sort(lambda x, y: cmp(d[y], d[x]))
	return k


#
# Verifie que les arguments obligatoires de la ligne de commande sont presents
# Recupere les valeurs des parametres optionnels
# -options- est un tableau de triplets (nom, constructeur, val_defaut)
#
def checkArgs(args, options = [], info = ""):

	#
	# Affiche le message d'erreur de mauvais arguments
	#
	def error_usage():
		s = "- ERROR - Usage : " + sys.argv[0]
		for t in args:
			s += " " + t
		print >> sys.stderr, s
		for t in options:
			if t[1] == bool:
				invite = "+/-"
			else:
				invite = "-"
			print >> sys.stderr, "  ", invite + "%s [%s] (%s)" % t
		if info != "":
			print >> sys.stderr, info
		sys.exit(1)

	types = dict([ (x[0], x[1]) for x in options ])
	valOpt = dict([ (x[0], x[2]) for x in options ])
	valArg = []
	
	# On scanne les argumetns pour les compter et recuperer les valeurs
	for t in sys.argv[1:]:

		# Un petit peu d'aide
		if t == '-h' or t == '--help':
			error_usage()
			
		# Un argument optionnel
		if t[0] in ['-', '+']:
		
			# Un parametre non bool
			if '=' in t:
				s = t[1:t.index('=')]
				v = t[t.index('=')+1:]
				if not s in valOpt:
					error_usage()
				if types[s] == bool:
					error_usage()
				valOpt[s] = types[s](v)
				
			else:
				s = t[1:]
				if not s in valOpt:
					error_usage()
				if types[s] != bool:
					error_usage()
				# Ici, on affecte False
				valOpt[s] = (t[0] == '+')
		elif os.access(t, os.R_OK):
			valArg.append(t)
		else:
			print >> sys.stderr, "No access to", t
			error_usage()

	# Il n'y a pas le nombre d'arguments minimal
	if len(valArg) != len(args):
		error_usage()
	
	return (valArg, valOpt)



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

