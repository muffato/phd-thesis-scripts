#!/usr/bin/env python2

__doc__ = """
A partir de toutes les diagonales extraites entre les especes,
  reconstruit les chromosomes (ou scaffold) de chaque ancetre.
"""

import sys
import utils.myTools
import utils.myGenomes
import utils.myPhylTree

it = utils.myTools.myIterator.tupleOnStrictUpperList

###############################################################################
# La fonction de calcul du score de presence sur le meme chromosome ancestral #
###############################################################################
#@utils.myTools.memoize
def calcProba(comparedEsp, communEsp):
	
	# On renvoie une valeur non nulle uniquement si au moins deux branches sont couvertes par des especes OK
	if len([l for l in lstEspParNoeudsFils if len(l.intersection(communEsp)) > 0]) < 2:
		return 0

	# Table des valeurs
	values = {}
	for e in comparedEsp:
		values[e] = espIncertitude[e]
	for e in communEsp:
		values[e] = espCertitude[e]
	
	if arguments["scoringMethod"] == 0:
		# Methode initiale: poids pre-calcules de chaque espece et bi-produits
		prop = [sum([x*dicPoidsEspeces[f] for (f,x) in values.iteritems() if f in l]) for l in lstEspParNoeudsFils]
		return sum( [f1 * f2 for (f1,f2) in it(prop)] )

	if arguments["scoringMethod"] == 3:
		# Methode initiale: poids pre-calcules de chaque espece et bi-produits
		prop = [sum([x*dicPoidsEspeces[f] for (f,x) in values.iteritems() if f in l]) for l in lstEspParNoeudsFils]
		return sum( prop )

	if arguments["scoringMethod"] == 1:
		# Methode intermediaire: interpolation sur les sous arbres et bi-produits
		prop = [phylTree.calcDist(values, f) for f in lstNoeudsFils]
		return sum( [f1 * f2 for (f1,f2) in it([x for x in prop if x != None])] )

	if arguments["scoringMethod"] == 4:
		# Methode intermediaire: interpolation sur les sous arbres et bi-produits
		prop = [phylTree.calcDist(values, f) for f in lstNoeudsFils]
		return sum( [x for x in prop if x != None] )

	# Nouvelle methode: interpolation sur l'arbre entier
	return phylTree.calcWeightedValue(values, -1, None)[phylTree.indNames[anc]]



arguments = utils.myTools.checkArgs( \
	[("phylTree.conf",file), ("signatures",file),("ancestr",str)], \
	[("scoringMethod",int,[0,1,2,3,4]), ("weightNbChr+",bool,False), ("weightNbChr-",bool,False)], \
	__doc__ \
)


#  Chargement et initialisation
################################
phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])
anc = phylTree.officialName[arguments["ancestr"]]
lstNoeudsFils = [a for (a,_) in phylTree.items[anc]]
lstEspParNoeudsFils = [phylTree.allDescendants[f] for f in lstNoeudsFils]
if anc in phylTree.parent:
	lstNoeudsFils.append( phylTree.parent[anc][0] )
	lstEspParNoeudsFils.append( phylTree.allDescendants[phylTree.root].difference(phylTree.allDescendants[anc]) )

# Les poids et taux de precisions de chaque espece
###################################################

# L'apport en terme de couverture de l'arbre phylogenetique
def calcPoidsFils(node, calc):
	dicPoidsEspeces[node] = calc
	if len(phylTree.tmpItems[node]) > 0:
		poids = calc / float(len(phylTree.tmpItems[node]))
		for (f,_) in phylTree.tmpItems[node]:
			calcPoidsFils(f, poids)
dicPoidsEspeces = dict.fromkeys(phylTree.listSpecies, 0.)
phylTree.initCalcDist(anc, True)
calcPoidsFils(anc, float(len(phylTree.tmpItems[0])))

espCertitude = dict.fromkeys(phylTree.commonNames, 1.)
espIncertitude = dict.fromkeys(phylTree.commonNames, 0.)

print >> sys.stderr, "%-25s\t%s\t%s\t%s" % ("Ancetre", "Poids", "Cert", "Incert")
for e in dicPoidsEspeces:
	print >> sys.stderr, "%-25s\t%.3f\t%.3f\t%.3f" % (e, dicPoidsEspeces[e], espCertitude[e], espIncertitude[e])

maxS = [3,3,1,3,3]

for (s,_,_) in utils.myFile.myTSV.readTabular(arguments["signatures"], (str,int,int)):
	comparedEsp = set()
	communEsp = set()
	for (c,e) in zip(s, phylTree.listSpecies):
		if c == '+':
			communEsp.add(e)
		elif c == '-':
			comparedEsp.add(e)
	x = calcProba(comparedEsp, communEsp)
	print int(x>0.5), x, maxS[arguments["scoringMethod"]]-x

