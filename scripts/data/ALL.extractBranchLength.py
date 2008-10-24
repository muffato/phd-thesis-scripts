#! /users/ldog/muffato/python

__doc__ = """
	Lit les arbres de proteines et cree les fichiers de genes ancestraux
"""

import sys
import utils.myTools
import utils.myPhylTree
import utils.myProteinTree

arguments = utils.myTools.checkArgs( [("phylTree.conf",file), ("proteinTree",file)], [], __doc__ )

phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])

lengths = utils.myTools.defaultdict(list)

# Parcours recursif de la famille de genes
def do(node):
	if node in data:
		for (g,d) in data[node]:
			# Une distance ne peut etre prise qu'entre deux noeuds de speciation
			# Les tests de duplication sont par rapport a 0 parce que l'arbre n'a pas encore ete modifie
			if (info[node]['Duplication'] == 0) and (info[g]['Duplication'] == 0):
				t1 = info[node]['taxon_name']
				t2 = info[g]['taxon_name']
				# Les deux noeuds doivent etre strictement consecutifs
				if (phylTree.parent[t2][0] == t1) and (d != 0):
					lengths[(t1,t2)].append(d)
					print utils.myFile.myTSV.printLine([t1,t2,d])
			do(g)

for (r,data,info) in utils.myProteinTree.loadTree(arguments["proteinTree"]):
	do(r)

# On trie les listes des longueurs
for l in lengths.itervalues():
	l.sort()

# Parcourt recursivement l'arbre et l'ecrit au format avec des parentheses, avec les longueurs de branche medianes
def convertToFlatFile(anc):

	a = phylTree.fileName[anc]
	if anc in phylTree.listSpecies:
		# On est arrive sur une feuille
		return a
	else:
		# On est sur un noeud, on construit la liste des distances
		l = []
		for (e,_) in phylTree.items[anc]:
			if len(lengths[(anc,e)]) == 0:
				# Par defaut, on revient a 0
				l.append(0)
			else:
				# La mediane est plus valable que la moyenne
				l.append( lengths[(anc,e)][len(lengths[(anc,e)])/2] )
		# Constuction de la chaine finale
		return "(" + ",".join([convertToFlatFile(e) + ":" + str(l) for ((e,_),l) in zip(phylTree.items[anc],l) ]) + ")%s|%d" % (a,phylTree.ages[anc])


print convertToFlatFile(phylTree.root), ";"


