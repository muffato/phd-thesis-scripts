#! /users/ldog/muffato/python

__doc__ = """
	Lit un graphe (aretes + noms) et fusionne les noeuds relies par un score > seuil
	  Produit un nouveau graphe (aretes + noms)
"""

import sys
import math
import utils.myTools



arguments = utils.myTools.checkArgs( [("grapheInit",file), ("nomsNoeuds",file), ("lstFusions",file)], [], __doc__ )

# MAIN #

# On lit la liste des fusions
comb = utils.myTools.myCombinator([])
f = utils.myTools.myOpenFile(arguments["lstFusions"], 'r')
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
todo = collections.defaultdict(float)
fg = utils.myTools.myOpenFile(arguments["grapheInit"], 'r')
for l in fg:
	
	t = l.split()

	a = getInd(t[0])
	b = getInd(t[1])

	if a == b:
		continue
	
	if a < 0 or b < 0:
		if (b,a) in todo:
			todo[(b,a)] += float(t[2])
		else:
			todo[(a,b)] += float(t[2])
	elif (b,a) not in todo:
		todo[(a,b)] = float(t[2])


fg.close()

# On rajoute les aretes qui touchent des nouveaux noeuds fusionnes
for ((a,b),v) in todo.iteritems():
	print revert(a), revert(b), v

todo = None


# On lit les anciens noms des noeuds
f = utils.myTools.myOpenFile(arguments["nomsNoeuds"], 'r')
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

