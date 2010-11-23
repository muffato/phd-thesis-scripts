#!/usr/bin/env python2
# Copyright (c) 2010 IBENS/Dyogen : Matthieu MUFFATO, Alexandra LOUIS, Hugues ROEST CROLLIUS
# Mail : hrc@ens.fr or alouis@biologie.ens.fr
# Licences GLP v3 and CeCILL v2

__doc__ = """
	Lit les arbres de proteines et verifie que la structure est compatible avec l'arbre des especes
	Le mode "strict" ne s'applique pas a l'arbre tree.1.ensembl
"""

import sys
import collections

import utils.myFile
import utils.myTools
import utils.myPhylTree
import utils.myProteinTree

arguments = utils.myTools.checkArgs( [("phylTree.conf",file), ("proteinTree",file)], [("strict",bool,True)], __doc__ )

phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])

dupVal = 2 if arguments["strict"] else 1

def check(node):
	txt = "Node %d %s" % (node, info[node])
	if node in data:
		assert (info[node]['taxon_name'] in phylTree.listAncestr) or (info[node]['Duplication'] >= dupVal), ("SPECIES INTERNAL NODE WITHOUT DUPLICATION", txt)
		anc = None
		for (f,l) in data[node]:
			newtxt = "%s / %d %s" % (txt, f, info[f])
			assert l >= 0, ("NEGATIVE BRANCH LENGTH", l, newtxt)
			assert phylTree.isChildOf(info[f]['taxon_name'], info[node]['taxon_name']), ("TAXON_NAME NOT A CHILD OF LAST NODE (%s)" % info[node]['taxon_name'], newtxt)
			if arguments["strict"]:
				assert (info[f]['taxon_name'] != info[node]['taxon_name']) or (info[node]['Duplication'] >= dupVal), ("SAME CONSECUTIVE ANC AND NO DUPLICATION", newtxt)
			if anc is None:
				anc = info[f]['taxon_name']
			else:
				anc = phylTree.dicParents[anc][info[f]['taxon_name']]
			check(f)
		if arguments["strict"]:
			assert anc == info[node]['taxon_name'], ("ANC != LAST_COMMON_ANCESTOR (%s)" % anc, txt)
		else:
			assert phylTree.isChildOf(anc, info[node]['taxon_name']), ("ANC YOUNGER THAN LAST_COMMON_ANCESTOR (%s)" % anc, txt)
	else:
		assert info[node]['Duplication'] == 0, ("DUPLICATION SHOULD BE 0", txt)
		assert info[node]['taxon_name'] in phylTree.listSpecies, ("UNKNOWN TAXON_NAME", txt)


print >> sys.stderr, "Verification des arbres ...",
allidsS = set()
allidsL = []
nb = 0
for (r,data,info) in utils.myProteinTree.loadTree(arguments["proteinTree"]):
	check(r)
	allidsS.update(info)
	allidsL.extend(info)
	nb += 1
assert len(allidsS) == len(allidsL), "NON UNIQUE IDS"
print >> sys.stderr, "%d arbres OK" % nb

