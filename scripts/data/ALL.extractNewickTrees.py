#!/usr/bin/env python2
# Copyright (c) 2010 IBENS/Dyogen : Matthieu MUFFATO, Alexandra LOUIS, Hugues ROEST CROLLIUS
# Mail : hrc@ens.fr or alouis@biologie.ens.fr
# Licences GLP v3 and CeCILL v2

__doc__ = """
	Lit les arbres de proteines et cree les fichiers Newick
"""

# Librairies
import sys

import utils.myTools
import utils.myProteinTree

# Arguments
arguments = utils.myTools.checkArgs( [("proteinTree",file)], [("withDist",bool,False)], __doc__ )

def printNewickTree(node):
	if node in data:
		return "(" + ",".join(printNewickTree(g) + ((":%f" % l) if arguments["withDist"] else "") for (g,l) in data[node]) + ")"
	else:
		return info[node]['transcript_name'].split("/")[0]

print >> sys.stderr, "Mise en forme des arbres ...",
nb = 0
for (r,data,info) in utils.myProteinTree.loadTree(arguments["proteinTree"]):
	print printNewickTree(r)
	nb += 1
print >> sys.stderr, "%d arbres OK" % nb

