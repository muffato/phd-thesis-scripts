#!/usr/bin/env python2
# Copyright (c) 2010 IBENS/Dyogen : Matthieu MUFFATO, Alexandra LOUIS, Hugues ROEST CROLLIUS
# Mail : hrc@ens.fr or alouis@biologie.ens.fr
# Licences GLP v3 and CeCILL v2

__doc__ = """
	Charge le fichier de definition de GO et les associations genes Ensembl <-> GO
"""

import sys
import math
import collections

import utils.myFile
import utils.myTools
import utils.myMaths
import utils.myGenomes

# Arguments
arguments = utils.myTools.checkArgs( [("GO_db.obo",file), ("ensembl_export",file)], [], __doc__)

GO_down = collections.defaultdict(set)
GO_up = collections.defaultdict(set)
GO_name = {}
currid = None

f = utils.myFile.openFile(arguments["GO_db.obo"], "r")
for l in f:
	if l.startswith("id"):
		currid = l.split()[1]
	if l.startswith("name"):
		GO_name[currid] = l[6:-1]
	if l.startswith("is_a"):
		link = l.split()[1]
		if link.startswith("GO"):
			GO_up[currid].add(link)
			GO_down[link].add(currid)
f.close()

GO_roots = [x for x in GO_down if x not in GO_up]
print >> sys.stderr, GO_roots

GO_level = collections.defaultdict(int)
GO_allTerms = collections.defaultdict(set)
def setLevel(node, l):
	GO_level[node] = max(l, GO_level[node])
	GO_allTerms[node].add(node)
	for x in GO_up[node]:
		GO_allTerms[node].update(GO_allTerms[x])
	for x in GO_down[node]:
		setLevel(x, l+1)

for x in GO_roots:
	setLevel(x, 3)


assoc = collections.defaultdict(set)
f = utils.myFile.openFile(arguments["ensembl_export"], "r")
for l in f:
	t = l.split("\t")
	for x in [1]:
		if (t[x] in GO_up) or (t[x] in GO_down):
			assoc[t[0]].update(GO_allTerms[t[x]])
f.close()

for (x,l) in assoc.iteritems():
	print x, " ".join([y+":"+str(GO_level[y]) for y in l])

