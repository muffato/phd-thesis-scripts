#! /users/ldog/muffato/python -OO

__doc__ = """
Teste differents modeles de cassures/fusions de chromosomes
"""

import sys
import math
import random
import utils.myMaths
from utils.myTools import defaultdict


def countAltern(lst):
	
	# La liste des chromosomes de l'alternance
	lst.reverse()

	# Compte le nombre d'occurrences de c dans la liste courante
	def countChr(c):
		nb = 0
		for x in lst:
			if c not in x:
				break
			nb += 1
		return nb
	
	# Le compte final
	count = defaultdict(int)
	# La derniere ligne lue
	last = defaultdict(int)
	# On parcourt la liste
	while len(lst) > 0:
		curr = lst.pop()
		for x in curr:
			# Les alternances sont mesurees entre deux positions consecutives
			for y in last:
				if y == x:
					continue
				count[(x,y)] += (countChr(x)+1) * last[y]
				count[(y,x)] = count[(x,y)]
			# Et aussi entre les paralogues
			for y in curr:
				if y >= x:
					continue
				count[(x,y)] += 1
				count[(y,x)] = count[(x,y)]
		
		# On met a jour last
		for y in last:
			if y not in curr:
				last[y] = 0
		for x in curr:
			last[x] += 1

	return count


lst = [[9], [9], [11], [11], [9], [11], [11], [9,11], [9], [11], [9], [9,11], [9], [11], [9], [11]]
print countAltern(lst)

lst = [[9], [9], [11], [17], [9], [11], [11], [9,11], [9], [11], [9], [9,11], [9], [17], [9], [11]]
print countAltern(lst)

lst = [[9], [9], [11], [11], [9], [11]]
print countAltern(lst)


lst = [[9], [9], [11], [11], [9]] #, [11], [11]] #, [9], [9], [11], [9], [9], [9], [11], [9], [11]]
print countAltern(lst)



