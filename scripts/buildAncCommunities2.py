#! /users/ldog/muffato/python

__doc__ = """
	Reconstruit les chromosomes d'un ancetre a partir de toutes les diagonales extraites entre les especes
"""

import sys
import operator
import collections

import utils.myTools
import utils.myMaths
import utils.myPhylTree
import utils.walktrap

it = utils.myTools.myIterator.tupleOnStrictUpperList

###############################################################################
# La fonction de calcul du score de presence sur le meme chromosome ancestral #
###############################################################################
def calcProba(comparedEsp, communEsp):
	
	# On renvoie une valeur non nulle uniquement si au moins deux branches sont couvertes par des especes OK
	if len([l for l in lstEspParNoeudsFils if len(l.intersection(communEsp)) > 0]) < 2:
		return None

	# Table des valeurs
	values = {}
	for e in comparedEsp:
		values[e] = espIncertitude[e]
	for e in communEsp:
		values[e] = espCertitude[e]
	
	if arguments["scoringMethod"] == 0:
		# Methode initiale: poids pre-calcules de chaque espece et bi-produits
		prop = [sum(x*dicPoidsEspeces[f] for (f,x) in values.iteritems() if f in l) for l in lstEspParNoeudsFils]
		return sum( f1 * f2 for (f1,f2) in it(prop) )

	if arguments["scoringMethod"] == 3:
		# Methode initiale modifiee: poids pre-calcules de chaque espece et somme de tout le monde
		return sum(x*dicPoidsEspeces[f] for (f,x) in values.iteritems())

	if arguments["scoringMethod"] == 1:
		# Methode intermediaire: interpolation sur les sous arbres et bi-produits
		prop = [phylTree.calcDist(values, f) for f in lstNoeudsFils]
		return sum( f1 * f2 for (f1,f2) in it([x for x in prop if x != None]) )

	if arguments["scoringMethod"] == 4:
		# Methode intermediaire modifiee: interpolation sur les sous arbres et somme
		prop = [phylTree.calcDist(values, f) for f in lstNoeudsFils]
		return sum( x for x in prop if x != None )

	#if arguments["scoringMethod"] == 2:
	# Nouvelle methode: interpolation sur l'arbre entier
	return phylTree.calcWeightedValue(values, -1, None)[phylTree.indNames[anc]]


arguments = utils.myTools.checkArgs( \
	[("phylTree.conf",file), ("diagsList",file)], \
	[("minDiagLength",int,2), ("withMemoize",bool,False), \
	("weightNbChr+",bool,False), ("weightNbChr-",bool,False), ("weightFile",file,""), \
	("scoringMethod",int,[0,1,2,3,4]), ("scoreThresholdMin",float,-1.), ("scoreThresholdMax",float,10.), \
	("walktrapLength",int,5), ("cutChoice",str,["earliestCut","latestCut","highestRelevance","lowestRelevance"])], \
	__doc__ \
)

phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])

# Les diagonales et les scores possibles
#########################################
print >> sys.stderr, "Chargement des diagonales de %s ..." % arguments["diagsList"],
lstDiags = []
for (anc,_,diag,strands,e1,e2) in utils.myFile.myTSV.readTabular(arguments["diagsList"], [str,int,(str,4)]):
	
	# La diagonale
	d = [int(x) for x in diag.split(' ')]
	if len(d) < arguments["minDiagLength"]:
		continue
	
	# Les chromosomes en syntenie avec leurs scores
	tmp1 = {}
	for (e,c,n) in [y.split("/") for y in e1.split("|") if len(y) > 0]:
		tmp1[(e,c)] = n
	
	# Les localisations sur les chromosomes
	esp = {}
	for (e,c,n) in [y.split("/") for y in e2.split("|") if len(y) > 0]:
		if e in phylTree.listSpecies:
			n = 0
			# Le resultat sera pour chaque espece le chromosome d'appartenance majoritaire, departage par la syntenie en cas d'egalite
			esp[e] = max(esp.get(e, (0,None,None)), (n, tmp1.get((e,c),0), c) )
	espChr = frozenset((e,esp[e][2]) for e in esp)
	#print [c for (e,c) in espChr if e == "Homo sapiens"], [c for (e,c) in espChr if e == "Pan troglodytes"], [c for (e,c) in espChr if e == "Pongo pygmaeus"]
	lstDiags.append( (d,espChr,frozenset(esp),diag,strands) )

print >> sys.stderr, "OK (%d diagonales)" % len(lstDiags)

#sys.exit(0)

#  Chargement et initialisation
################################
lstNoeudsFils = [a for (a,_) in phylTree.items[anc]]
lstEspParNoeudsFils = [phylTree.allDescendants[f] for f in lstNoeudsFils]
if anc in phylTree.parent:
	lstNoeudsFils.append( phylTree.parent[anc][0] )
	lstEspParNoeudsFils.append( phylTree.allDescendants[phylTree.root].difference(phylTree.allDescendants[anc]) )

# Les taux de precisions de chaque espece
##########################################
espCertitude = dict.fromkeys(phylTree.listSpecies, 1.)
espIncertitude = dict.fromkeys(phylTree.listSpecies, 0.)
dicPoidsEspeces = dict.fromkeys(phylTree.listSpecies, 0.)

if arguments["weightNbChr+"] or arguments["weightNbChr-"]:
	# Chargement du fichier
	weights = dict(utils.myFile.myTSV.readTabular(arguments["weightFile"], [str,float]))
	# Poids de certitude
	if arguments["weightNbChr+"]:
		for e in phylTree.listSpecies:
			espCertitude[e] = 1. - weights[e]
	# Poids d'incertitude
	if arguments["weightNbChr-"]:
		alphaIncertitude = min(weights.values()) * max(weights.values())
		for e in phylTree.listSpecies:
			espIncertitude[e] = alphaIncertitude / weights[e]

# Initialisation de chaque methode de calcul de score
######################################################
if arguments["scoringMethod"] != 2:
	phylTree.initCalcDist(anc, True)

if arguments["scoringMethod"] in [0, 3]:

	# L'apport en terme de couverture de l'arbre phylogenetique
	def calcPoidsFils(node, calc):
		dicPoidsEspeces[node] = calc
		if len(phylTree.tmpItems[node]) > 0:
			poids = calc / float(len(phylTree.tmpItems[node]))
			for (f,_) in phylTree.tmpItems[node]:
				calcPoidsFils(f, poids)
	calcPoidsFils(anc, float(len(phylTree.tmpItems[0])))

# Les scores maxima
x = float(len(lstEspParNoeudsFils))
maxS = [x*(x-1)/2,x*(x-1)/2,1.,x,x][arguments["scoringMethod"]]

# Resume
print >> sys.stderr, "%-25s\t%s\t%s\t%s" % ("Espece", "Poids", "Cert", "Incert")
for e in phylTree.listSpecies:
	print >> sys.stderr, "%-25s\t%.3f\t%.3f\t%.3f" % (e, dicPoidsEspeces[e], espCertitude[e], espIncertitude[e])

if arguments["withMemoize"]:
	calcProba = utils.myTools.memoize(calcProba)

# On calcule les scores
print >> sys.stderr, "Calcul de la matrice ...",
edges = collections.defaultdict(dict)
for i1 in xrange(len(lstDiags)):
	(d1,ec1,e1,_,_) = lstDiags[i1]
	for i2 in xrange(i1):
		(d2,ec2,e2,_,_) = lstDiags[i2]
		x = calcProba( e1.intersection(e2), frozenset(e for (e,_) in ec1.intersection(ec2)) )
		if x is None:
			continue
		assert 0 < x <= maxS
		if x < arguments["scoreThresholdMin"]:
			continue
		if x >= arguments["scoreThresholdMax"]:
			x = maxS
		edges[i1][i2] = edges[i2][i1] = x
print >> sys.stderr, "OK"

res = utils.walktrap.doWalktrap(edges, showProgress=True, randomWalksLength=arguments["walktrapLength"])

clusters = []
relev = []
# Chaque composante connexe
for (nodes,cuts,_,dend) in res:
	print >> sys.stderr, cuts,

	if arguments["cutChoice"] == "earliestCut":
		cuts = sorted(cuts)
	elif arguments["cutChoice"] == "latestCut":
		cuts = sorted(cuts, reverse=True)
	elif arguments["cutChoice"] == "lowestRelevance":
		cuts = sorted(cuts, key=operator.itemgetter(1))
	elif arguments["cutChoice"] == "highestRelevance":
		cuts = sorted(cuts, key=operator.itemgetter(1), reverse=True)
	
	for (alpha,score) in sorted(cuts):
		# Un score de relevance > 0.1
		if score > 0.1:
			(clust,lonely) = dend.cut(alpha)
			# Les noeuds seuls doivent representer < de la moitie de l'ensemble des noeuds
			if len(lonely) < len(nodes)/2:
				print >> sys.stderr, "+%d/%d/%f/%f" % (len(clust), len(lonely), alpha, score)
				clusters.extend(clust)
				relev.append(score)
				break
	else:
		print >> sys.stderr, "-"
		clusters.append(nodes)
		relev.append(0)

print >> sys.stderr, "Impression des %d chromosomes ancestraux (score = %f)..." % (len(clusters),utils.myMaths.myStats.mean(relev)),

notused = set(xrange(len(lstDiags)))
for (i,chrom) in enumerate(clusters):
	for d in chrom:
		notused.remove(d)
		(_,_,_,d,s) = lstDiags[d]
		print "%d\t%s\t%s" % (i+1,d,s)

for d in notused:
	(_,_,_,d,s) = lstDiags[d]
	print "-\t%s\t%s" % (d,s)

print >> sys.stderr, "OK"

