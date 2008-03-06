#! /users/ldog/muffato/python -OO

__doc__ = """
A partir de toutes les diagonales extraites entre les especes,
  reconstruit les chromosomes (ou scaffold) de chaque ancetre.
"""

import sys
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
	
	if options["scoringMethod"] == 0:
		# Methode initiale: poids pre-calcules de chaque espece et bi-produits
		prop = [sum([x*dicPoidsEspeces[f] for (f,x) in values.iteritems() if f in l]) for l in lstEspParNoeudsFils]
		return sum( [f1 * f2 for (f1,f2) in it(prop)] )

	if options["scoringMethod"] == 1:
		# Methode intermediaire: interpolation sur les sous arbres et bi-produits
		prop = [phylTree.calcDist(values, f) for f in lstNoeudsFils]
		return sum( [f1 * f2 for (f1,f2) in it([x for x in prop if x != None])] )

	# Nouvelle methode: interpolation sur l'arbre entier
	return phylTree.calcWeightedValue(values, -1, None)[phylTree.indNames[anc]]


####################################
# Charge le fichier des diagonales #
####################################
def loadDiagsFile(nom, ancName):
	
	print >> sys.stderr, "Chargement du fichier de diagonales ...",
	f = utils.myTools.myOpenFile(nom, 'r')
	lst = []
	for l in f:
		# TODO non necessaire dans l'avenir
		# Selection de l'ancetre
		if not l.startswith(ancName):
			continue
		# On enleve les "_random" et on extrait chaque colonne
		ct = l.replace('\n', '').replace('_random', '').split('\t')
		
		# La diagonale
		d = [int(x) for x in ct[2].split(' ')]
		
		# On joint les especes qui ont vu la diagonale et celles qui n'apportent que le chromosome
		tmp = [y.split("/") for y in "|".join([x for x in ct[4:] if len(x) > 0]).split("|")]
		# Les chromosomes de ces especes
		espChr = frozenset( [(phylTree.officialName[e],c) for (e,c) in tmp if ('Un' not in c) and (e in phylTree.officialName)] )
		esp = frozenset( [e for (e,_) in espChr] )
		
		lst.append( (d,espChr,esp,ct[2],ct[3]) )

	f.close()
	print >> sys.stderr, "OK (%d diagonales)" % len(lst)
	return lst


#############################################################################################################
# Rajoute les genes non presents dans des diagonales, en les considerant comme des diagonales de longueur 1 #
#############################################################################################################
def checkLonelyGenes():

	# Charge les genomes des ancetres outgroup
	genesAnc = {}
	a = anc
	while a in phylTree.parent:
		(a,_) = phylTree.parent[a]
		genesAnc[a] = utils.myGenomes.Genome(options["ancGenesFile"] % phylTree.fileName[a])
	
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
		
		# Petit test: ne peuvent etre utilises que les genes avec au moins 2 especes
		# Pour etre exact, il faut avoir 2 groupes parmi (fils1 + ... + filsN + outgroup)
		if len(lst) >= 2:
			new.append( ([i],frozenset(lst),frozenset([e for (e,_) in lst]),str(i),"0") )
			nb += 1
	
	print >> sys.stderr, "%d genes OK" % len(new)
	return new


################################################################################
# Charge les ancetres deja construits et rajoute les infos dans les diagonales #
################################################################################
def checkAlreadyBuildAnc():
	genAlready = {}
	for f in lstNoeudsFils:
		s = options["alreadyBuiltAnc"] % phylTree.fileName[f]
		print >> sys.stderr, "Checking %s ..." % f,
		if utils.myTools.fileAccess(s):
			genAlready[f] = utils.myGenomes.Genome(s)
			s = len(genAlready[f].lstChr)
			if options["weightNbChr+"] and (s > 0):
				espCertitude[f] = 1. - 1. / float(s)
			if options["weightNbChr-"]:
				espIncertitude[f] = s * alphaIncertitude
		else:
			print >> sys.stderr, "Not found"

	print >> sys.stderr, "Mise a jour des chromosomes des diagonales ...",
	for (d,espC,esp,_,_) in lstDiags:
		g = utils.myMaths.flatten([lstGenesAnc[i].names for i in d])
		for f in genAlready:
			espC.update([(f,genAlready[f].dicGenes[s][0]) for s in g if s in genAlready[f].dicGenes])
			esp.update([f for s in g if s in genAlready[f].dicGenes])
	print >> sys.stderr, "OK"





(noms_fichiers, options) = utils.myTools.checkArgs( \
	["phylTree.conf", "diagsList"], \
	[("ancestr",str,""), ("scoringMethod",int,[0,1,2]), \
	("useOutgroups",bool,True), ("alreadyBuiltAnc",str,""), \
	("useLonelyGenes",bool,False), ("weightNbChr+",bool,False), ("weightNbChr-",bool,False), ("walktrapLength",int,5), \
	("genesFile",str,"~/work/data/genes/genes.%s.list.bz2"), \
	("ancGenesFile",str,"~/work/data/ancGenes/ancGenes.%s.list.bz2")], \
	__doc__ \
)


#  Chargement et initialisation
################################
phylTree = utils.myPhylTree.PhylogeneticTree(noms_fichiers["phylTree.conf"])
anc = phylTree.officialName[options["ancestr"]]
genesAnc = utils.myGenomes.Genome(options["ancGenesFile"] % phylTree.fileName[anc])
lstGenesAnc = genesAnc.lstGenes[None]
lstNoeudsFils = [a for (a,_) in phylTree.items[anc]]
lstEspParNoeudsFils = [phylTree.allDescendants[f] for f in lstNoeudsFils]
if (anc in phylTree.parent) and options["useOutgroups"]:
	lstNoeudsFils.append( phylTree.parent[anc][0] )
	lstEspParNoeudsFils.append( phylTree.allDescendants[phylTree.root].difference(phylTree.allDescendants[anc]) )

# Chargement des genomes si necessaire
#######################################
lengths = []
dicGenes = {}
if options["weightNbChr+"] or options["weightNbChr-"] or options["useLonelyGenes"]:
	for e in phylTree.listSpecies:
		genome = utils.myGenomes.Genome(options["genesFile"] % phylTree.fileName[e])
		lengths.append(len(genome.lstChr))
		if options["useLonelyGenes"]:
			for (g,(c,_)) in genome.dicGenes.iteritems():
				dicGenes[g] = (e,c)


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
phylTree.initCalcDist(anc, options["useOutgroups"])
calcPoidsFils(anc, float(len(phylTree.tmpItems[0])))

espCertitude = dict.fromkeys(phylTree.commonNames, 1.)
espIncertitude = dict.fromkeys(phylTree.commonNames, 0.)
# Poids de certitude
if options["weightNbChr+"]:
	for (i,e) in enumerate(phylTree.listSpecies):
		if lengths[i] > 0:
			espCertitude[e] = 1. - 1. / float(lengths[i])
# Poids d'incertitude
if options["weightNbChr-"]:
	alphaIncertitude = 1. / float(min([x for x in lengths if x > 0]) * max(lengths))
	for (i,e) in enumerate(phylTree.listSpecies):
		espIncertitude[e] = lengths[i] * alphaIncertitude

# Les diagonales et les scores possibles
#########################################
lstDiags = loadDiagsFile(noms_fichiers["diagsList"], anc)

# On doit rajouter les genes non presents dans des diagonales
if options["useLonelyGenes"]:
	lstDiags.extend(checkLonelyGenes())

# On doit noter les chromosomes des diagonales sur ces ancetres deja construits
if options["alreadyBuiltAnc"] != "":
	checkAlreadyBuildAnc()

dicGenes.clear()

print >> sys.stderr, "%-25s\t%s\t%s\t%s" % ("Ancetre", "Poids", "Cert", "Incert")
for e in dicPoidsEspeces:
	print >> sys.stderr, "%-25s\t%.3f\t%.3f\t%.3f" % (e, dicPoidsEspeces[e], espCertitude[e], espIncertitude[e])

# On calcule les scores
print >> sys.stderr, "Calcul de la matrice ...",

walktrapInstance = utils.walktrap.WalktrapLauncher(showProgress=True, randomWalksLength=options["walktrapLength"])

# On les stocke directement dans l'instance de walktrap
for i1 in xrange(len(lstDiags)):
	(d1,ec1,e1,_,_) = lstDiags[i1]
	for i2 in xrange(i1):
		(d2,ec2,e2,_,_) = lstDiags[i2]
		x = calcProba( e1.intersection(e2), frozenset([e for (e,_) in ec1.intersection(ec2)]) )
		if x > 0:
			walktrapInstance.edges[i1][i2] = walktrapInstance.edges[i2][i1] = x

print >> sys.stderr, "OK"
walktrapInstance.doWalktrap()

clusters = []
# Chaque composante connexe
for (nodes,cuts,_,dend) in walktrapInstance.res:
	print >> sys.stderr, cuts,
	# Un score de relevance > 0.1
	interessant = [(alpha,score,dend.cut(alpha)) for (alpha,score) in cuts if score > 0.1]
	# Les noeuds seuls doivent representer < de la moitie de l'ensemble des noeuds
	interessant = [(alpha,score,clust,lonely) for (alpha,score,(clust,lonely)) in interessant if len(lonely) < len(nodes)/2]
	if len(interessant) == 0:
		print >> sys.stderr, "-"
		clusters.append(nodes)
	else:
		(alpha,score,clust,lonely) = interessant[-1]
		# Au choix, on prend la version la moins fusionnee
		print >> sys.stderr, "+%d/%d/%f/%f" % (len(clust), len(lonely), alpha, score)
		clusters.extend(clust)
		#clusters.append(lonely)
# -> clusters contient la repartition des diagonales
print >> sys.stderr

print >> sys.stderr, "Impression des %d chromosomes ancestraux ..." % len(clusters),

used = [False] * len(lstDiags)
for (i,chrom) in enumerate(clusters):
	for d in chrom:
		used[d] = True
		(_,_,_,d,s) = lstDiags[d]
		print "%d\t%s\t%s" % (i+1,d,s)

for (d,x) in enumerate(used):
	if not x:
		(_,_,_,d,s) = lstDiags[d]
		print "-\t%s\t%s" % (d,s)
