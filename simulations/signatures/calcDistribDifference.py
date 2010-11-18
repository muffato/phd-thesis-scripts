#!/usr/bin/env python2

__doc__ = """
	Renvoie les listes des devenirs de chaque gene le long des branches de l'arbre phylogenetique
"""

import sys
import math
import random
import utils.myMaths
import utils.myTools


# Argument:
arguments = utils.myTools.checkArgs( [("scoresOK",file), ("scoresNO",file)], [], __doc__ )

def calc():
	scoresOK = utils.myTools.myOpenFile(arguments["scoresOK"], "r")
	listeOK = [float(x) for x in scoresOK]
	scoresOK.close()

	scoresNO = utils.myTools.myOpenFile(arguments["scoresNO"], "r")
	listeNO = [float(x) for x in scoresNO]
	scoresNO.close()

	def stats(liste):
		mean = utils.myMaths.myStats.mean(liste)
		var = sum( (x-mean)**2 for x in liste ) / float(len(liste))
		return (mean,var)

	def do(l1, l2):
		(m1,s21) = stats(l1)
		(m2,s22) = stats(l2)

		return ((m1,s21,len(l1)), (m2,s22,len(l2)), abs(m1-m2)/math.sqrt(s21/len(l1)+s22/len(l2)))
	
	(x,y,z) = do(listeOK,listeNO)
	print x
	print y
	print z
	print

	n = 0
	for _ in xrange(1000):
		l1 = random.sample(listeOK, 1000)
		l2 = random.sample(listeNO, 1000)
		(x,y,z) = do(l1, l2)
		print z
		if z >= 1.96:
			n += 1
	
	print
	print "*", n
	print


calc()

