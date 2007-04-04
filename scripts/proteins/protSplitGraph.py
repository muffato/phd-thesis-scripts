#! /users/ldog/muffato/python -OO

__doc__ = """
	Lit le graphe (aretes + noms) et les composantes connexes puis cree les sous-graphes (aretes + noms)
"""

import sys
import math
import utils.myTools




(noms_fichiers, _) = utils.myTools.checkArgs( ["grapheInit", "nomsNoeuds", "lstComposantesConnexes", "dossierGraphes"], [], __doc__ )



# On lit les noms des noeuds
f = utils.myTools.myOpenFile(noms_fichiers["nomsNoeuds"], 'r')
dicNodes = {}
for l in f:
	t = l.split()
	dicNodes[t[0]] = t[1:]
f.close()

print >> sys.stderr, "Noms des noeuds lus"

# On lit les composantes
f = utils.myTools.myOpenFile(noms_fichiers["lstComposantesConnexes"], 'r')
dic = {}
n = 0
for l in f:
	t = [x for x in l.split()]
	fout = utils.myTools.myOpenFile(noms_fichiers["dossierGraphes"] + "nodes.%d" % n, 'w')
	for i in xrange(len(t)):
		dic[t[i]] = (n,i)
		print >> fout, i, " ".join(dicNodes[t[i]])
	fout.close()
	n += 1
f.close()
dicNodes = None

print >> sys.stderr, "Composantes connexes chargees et noeuds mis a jour"

# On lit le graphe et cree les sous-graphes
lastG = -1
f = None
fg = utils.myTools.myOpenFile(noms_fichiers["grapheInit"], 'r')
openFiles = {}

for l in fg:
	
	t = l.split()

	(g,ia) = dic[t[0]]
	(_,ib) = dic[t[1]]

	if g not in openFiles:
		if len(openFiles) > 1000:
			(x,f) = openFiles.popitem()
			f.close()
		f = utils.myTools.myOpenFile(noms_fichiers["dossierGraphes"] + "graph.%d" % g, 'a')
		openFiles[g] = f
	else:
		f = openFiles[g]

	print >> f, ia, ib, t[2]

fg.close()

for x in openFiles:
	openFiles[x].close()

print >> sys.stderr, "Graphe divise"

