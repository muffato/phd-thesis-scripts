#! /users/ldog/muffato/python

__doc__ = """
A partir de toutes les diagonales extraites entre les especes,
  reconstruit les chromosomes (ou scaffold) de chaque ancetre.
"""

import sys
import numpy
import itertools

import utils.myTools
import utils.myMaths
import utils.myGenomes
import utils.myPhylTree
import utils.walktrap

it = utils.myTools.myIterator.tupleOnStrictUpperList

###############################################################################
# La fonction de calcul du score de presence sur le meme chromosome ancestral #
###############################################################################
#@utils.myTools.memoize
def calcProba(comparedEsp, communEsp, espIncertitude, espCertitude, scoringMethod):
	
	# On renvoie une valeur non nulle uniquement si au moins deux branches sont couvertes par des especes OK
	if len([l for l in lstEspParNoeudsFils if len(l.intersection(communEsp)) > 0]) < 2:
		return 0

	# Table des valeurs
	values = {}
	for e in comparedEsp:
		values[e] = espIncertitude[e]
	for e in communEsp:
		values[e] = espCertitude[e]
	
	if scoringMethod == 0:
		# Methode initiale: poids pre-calcules de chaque espece et bi-produits
		prop = [sum([x*dicPoidsEspeces[f] for (f,x) in values.iteritems() if f in l]) for l in lstEspParNoeudsFils]
		return sum( [f1 * f2 for (f1,f2) in it(prop)] )

	if scoringMethod == 3:
		# Methode initiale: poids pre-calcules de chaque espece et bi-produits
		prop = [sum([x*dicPoidsEspeces[f] for (f,x) in values.iteritems() if f in l]) for l in lstEspParNoeudsFils]
		return sum( prop )

	if scoringMethod == 1:
		# Methode intermediaire: interpolation sur les sous arbres et bi-produits
		prop = [phylTree.calcDist(values, f) for f in lstNoeudsFils]
		return sum( [f1 * f2 for (f1,f2) in it([x for x in prop if x != None])] )

	if scoringMethod == 4:
		# Methode intermediaire: interpolation sur les sous arbres et bi-produits
		prop = [phylTree.calcDist(values, f) for f in lstNoeudsFils]
		return sum( [x for x in prop if x != None] )

	# Nouvelle methode: interpolation sur l'arbre entier
	return phylTree.calcWeightedValue(values, -1, None)[phylTree.indNames[anc]]


#############################################################################################################
# Rajoute les genes non presents dans des diagonales, en les considerant comme des diagonales de longueur 1 #
#############################################################################################################
def checkLonelyGenes():

	# Charge les genomes des ancetres outgroup
	genesAnc = {}
	a = anc
	while a in phylTree.parent:
		(a,_) = phylTree.parent[a]
		genesAnc[a] = utils.myGenomes.Genome(arguments["ancGenesFile"] % phylTree.fileName[a])
	
	print >> sys.stderr, "Ajout des genes solitaires ...",
	
	# Liste des genes solitaires
	genesSeuls = set(xrange(len(lstGenesAnc)))
	for (d,_,_,_,_) in lstDiags:
		genesSeuls.difference_update(d)
	
	nb = 0
	new = []
	# Les genes seuls vont devenir des diagonales de 1
	for i in genesSeuls:

		lst = []
		for e in phylTree.listSpecies:
			if e in phylTree.outgroupSpecies[anc]:
				a = phylTree.dicParents[anc][e]
				names = genesAnc[a].getOtherNames(lstGenesAnc[i].names)
			else:
				names = lstGenesAnc[i].names
			tmp = [dicGenes[x] for x in names if x in dicGenes]
			tmp = set( c for (x,c) in tmp if x == e )

			# Si les orthologues sont sur un unique chromosome
			if len(tmp) == 1:
				c = str(tmp.pop()).replace("_random","")
				if 'Un' not in c:
					lst.append( (e,c) )
		
		esp = frozenset([e for (e,_) in lst])
		# Filtre: il faut que le gene soit dans 2 groupes parmi (fils1 + ... + filsN + outgroup)
		if len([l for l in lstEspParNoeudsFils if len(l.intersection(esp)) > 0]) >= 2:
			new.append( ([i],frozenset(lst),esp,str(i),"0") )
			nb += 1
	
	print >> sys.stderr, "%d genes OK" % len(new)
	return new





arguments = utils.myTools.checkArgs( \
	[("phylTree.conf",file), ("diagsList",file), ("ancestr",str)], \
	[("nbiter",int,10), ("useLonelyGenes",bool,True), \
	("outputAncestralGenomes",str,""), \
	("genesFile",str,"~/work/data/genes/genes.%s.list.bz2"), \
	("ancGenesFile",str,"~/work/data/ancGenes/ancGenes.%s.list.bz2")], \
	__doc__ \
)


#  Chargement et initialisation
################################

# Donnees phylogenetique
phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])
anc = phylTree.officialName[arguments["ancestr"]]
genesAnc = utils.myGenomes.Genome(arguments["ancGenesFile"] % phylTree.fileName[anc])
lstGenesAnc = genesAnc.lstGenes[None]
lstNoeudsFils = [a for (a,_) in phylTree.items[anc]]
lstEspParNoeudsFils = [phylTree.allDescendants[f] for f in lstNoeudsFils]
if anc in phylTree.parent:
	lstNoeudsFils.append( phylTree.parent[anc][0] )
	lstEspParNoeudsFils.append( phylTree.allDescendants[phylTree.root].difference(phylTree.allDescendants[anc]) )

# Diagonales
lstDiags = loadDiagsFile(arguments["diagsList"], [arguments["ancestr"]], phylTree.officialName)[arguments["ancestr"]]

# Genomes pour les genes singletons
lengths = []
dicGenes = {}
for e in phylTree.listSpecies:
	genome = utils.myGenomes.Genome(arguments["genesFile"] % phylTree.fileName[e])
	lengths.append(len(genome.lstChr))
	for (g,(c,_)) in genome.dicGenes.iteritems():
		dicGenes[g] = (e,c)
if arguments["useLonelyGenes"]:
	lstDiags.extend( checkLonelyGenes() )
dicGenes.clear()

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
maxScores = range(5)
maxScores[0] = maxScores[3] = float(len(lstEspParNoeudsFils)*(len(lstEspParNoeudsFils)-1))/2.
maxScores[1] = maxScores[4] = float(len(lstEspParNoeudsFils))
maxScores[2] = 1.

# Pour stocker les numeros de chromosomes de chaque gene et de chaque solution
positions = numpy.empty( (len(lstDiags),160) )
lastPos = range(len(lstDiags))
lastGrp = []


def doIter(niter):

	print >> sys.stderr, "\n> Iteration ", niter
	nsol = 0
	unused = 0
	positions.fill(0)

	for (scoringMethod,weightM,weightP) in itertools.product([0,1,2,3,4],[False,True],[False,True]):

		print >> sys.stderr, ">> Parametres:", scoringMethod, weightM, weightP,

		espCertitude = dict.fromkeys(phylTree.commonNames, 1.)
		espIncertitude = dict.fromkeys(phylTree.commonNames, 0.)
		# Poids de certitude
		if weightM:
			for (i,e) in enumerate(phylTree.listSpecies):
				if lengths[i] > 0:
					espCertitude[e] = 1. - 1. / float(lengths[i])
		# Poids d'incertitude
		if weightP:
			alphaIncertitude = 1. / float(min([x for x in lengths if x > 0]) * max(lengths))
			for (i,e) in enumerate(phylTree.listSpecies):
				espIncertitude[e] = lengths[i] * alphaIncertitude

		edges = utils.myTools.defaultdict(dict)

		# On les stocke directement dans l'instance de walktrap
		nbedges = 0
		for i1 in xrange(len(lstDiags)):
			(d1,ec1,e1,_,_) = lstDiags[i1]
			for i2 in xrange(i1):
				if lastPos[i1] == lastPos[i2]:
					edges[i1][i2] = edges[i2][i1] = maxScores[scoringMethod]
					nbedges += 1
				else:
					(d2,ec2,e2,_,_) = lstDiags[i2]
					x = calcProba( e1.intersection(e2), frozenset([e for (e,_) in ec1.intersection(ec2)]), espIncertitude, espCertitude, scoringMethod )
					if x > 0:
						edges[i1][i2] = edges[i2][i1] = x
						nbedges += 1
		print >> sys.stderr, nbedges

		for randomWalksLength in [2,3,4,5,7,10,15,25]:

			print >> sys.stderr, ">>> WalktrapLength", randomWalksLength
			res = utils.walktrap.doWalktrap(edges.copy(), showProgress=True, randomWalksLength=randomWalksLength)

			clusters = []
			relev = []
			# Chaque composante connexe
			for (nodes,cuts,_,dend) in res:
				print >> sys.stderr, cuts,
				# Un score de relevance > 0.1
				interessant = [(alpha,score,dend.cut(alpha)) for (alpha,score) in cuts if score > 0.1]
				# Les noeuds seuls doivent representer < de la moitie de l'ensemble des noeuds
				interessant = [(alpha,score,clust,lonely) for (alpha,score,(clust,lonely)) in interessant if len(lonely) < len(nodes)/2]
				if len(interessant) == 0:
					print >> sys.stderr, "-"
					clusters.append(nodes)
					relev.append(0)
				else:
					(alpha,score,clust,lonely) = interessant[-1]
					# Au choix, on prend la version la moins fusionnee
					print >> sys.stderr, "+%d/%d/%f/%f" % (len(clust), len(lonely), alpha, score)
					clusters.extend(clust)
					#clusters.append(lonely)
					relev.append(score)
			# -> clusters contient la repartition des diagonales
			print >> sys.stderr

			print >> sys.stderr, ">>> %d chromosomes (score = %s)..." % (len(clusters),utils.myMaths.myStats.txtSummary(relev)),
			d = {False:"-", True:"+"}
			f = utils.myTools.myOpenFile(arguments["outputAncestralGenomes"] + "Iter%d.Length%d.Scoring%d.Genes%s.WeightM%s.WeightP%s"%(niter,randomWalksLength,scoringMethod,d[arguments["useLonelyGenes"]],d[weightM],d[weightP]), "w")
			for (i,chrom) in enumerate(clusters):
				for d in chrom:
					positions[d][nsol] = i+1
					(_,_,_,d,s) = lstDiags[d]
					print >> f, "%d\t%s\t%s" % (i+1,d,s)

			for d in xrange(len(lstDiags)):
				if positions[d][nsol] == 0:
					unused -= 1
					positions[d][nsol] = unused
					(_,_,_,d,s) = lstDiags[d]
					print >> f, "-\t%s\t%s" % (d,s)
			f.close()
			print >> sys.stderr, "OK"
			nsol += 1

	f = utils.myTools.myOpenFile(arguments["outputAncestralGenomes"] + "Common%d"%niter, "w")
	groups = utils.myTools.defaultdict(set)
	for d in xrange(len(lstDiags)):
		lastPos[d] = x = hash(tuple(positions[d]))
		groups[x].add(d)
		(_,_,_,d,s) = lstDiags[d]
		print >> f, "%d\t%s\t%s" % (x,d,s)
	f.close()

	print >> sys.stderr, "> Associations OK:", utils.myMaths.myStats.txtSummary([len(g) for g in groups.itervalues() if len(g) > 1])
	tmp = set([frozenset(x) for x in groups.itervalues()])
	if tmp == lastGrp:
		print >> sys.stderr, "> Reconstruction stable"
		return True
	lastGrp = tmp
	return False

# On lance les iterations
for niter in xrange(arguments["nbiter"]):
	if doIter(niter + 1):
		break

