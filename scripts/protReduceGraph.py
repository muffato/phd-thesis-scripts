#! /users/ldog/muffato/python -OO

__doc__ = """
	Lit un graphe (aretes + noms) et fusionne les noeuds relies par un score > seuil
	  Produit un nouveau graphe (aretes + noms)
"""

import sys
import math
import utils.myTools



(noms_fichiers, _) = utils.myTools.checkArgs( ["grapheInit", "nomsNoeuds", "lstFusions"], [], __doc__ )

# MAIN #

# On lit la liste des fusions
comb = utils.myTools.myCombinator([])
f = utils.myTools.myOpenFile(noms_fichiers["lstFusions"], 'r')
for l in f:
	t = l.split()
	comb.addLink(t)

comb.reduce()


def getInd(node):
	
	if node not in comb.dic:
		comb.addLink( [node])
		
	i = comb.dic[node]
	t = comb.grp[i]

	if len(t) == 1:
		return i
	return -1 - i

def revert(ind):
	if ind < 0:
		return -1-ind
	else:
		return ind

# Le nouveau graphe
todo = {}
fg = utils.myTools.myOpenFile(noms_fichiers["grapheInit"], 'r')
for l in fg:
	
	t = l.split()

	a = getInd(t[0])
	b = getInd(t[1])

	if a == b:
		continue
	
	if a < 0 or b < 0:
		todo[(a,b)] = todo.get( (a,b), 0) + float(t[2])
	else:
		print a, b, t[2]

fg.close()

# On rajoute les aretes qui touchent des nouveaux noeuds fusionnes
for (a,b) in todo:
	if (a < b) and ((b,a) in todo):
		continue
	print revert(a), revert(b), todo[(a,b)] + todo.get( (b,a), 0)

todo = None


# On lit les anciens noms des noeuds
f = utils.myTools.myOpenFile(noms_fichiers["nomsNoeuds"], 'r')
dicNodes = {}
for l in f:
	t = l.split()
	dicNodes[t[0]] = t[1:]
f.close()

# On ecrit les nouveaux noms de noeuds
n = 0
for g in comb:
	s = set()
	for x in g:
		s.update(dicNodes[x])
	print >> sys.stderr, n, " ".join(s)
	n += 1

