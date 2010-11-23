#!/usr/bin/env python2
# Copyright (c) 2010 IBENS/Dyogen : Matthieu MUFFATO, Alexandra LOUIS, Hugues ROEST CROLLIUS
# Mail : hrc@ens.fr or alouis@biologie.ens.fr
# Licences GLP v3 and CeCILL v2

__doc__ = """
A partir de toutes les diagonales extraites entre les especes,
  reconstruit les chromosomes (ou scaffold) de chaque ancetre.
"""

import sys
import utils.myTools
import utils.myDiags
import utils.myMaths
import utils.myGenomes
import utils.myPhylTree
import utils.walktrap
import utils.svm

it = utils.myTools.myIterator.tupleOnStrictUpperList


##############################
# L'appel au classifieur SVM #
##############################
#@utils.myTools.memoize
def calcProba(comparedEsp, communEsp):

	def transform(e):
		if e in communEsp:
			return 1
		elif e in comparedEsp:
			return -1
		else:
			return 0
	
	data = [transform(e) for e in phylTree.listSpecies]

	#return svm_model.predict(data)
	x = float(svm_model.predict_probability(data)[1][1])
	return x
	if x > 0.5:
		return 1.
		#return 0.75 + (x-0.5)/2.
	else:
		return 0.
		#return x
	#return float(svm_model.predict_probability(data)[1][1])


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
				names = genesAnc[a].getOtherNames(lstGenesAnc[i])
			else:
				names = lstGenesAnc[i]
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
	print "#", " ".join([str(x[0][0]) for x in new])
	return new




arguments = utils.myTools.checkArgs( \
	[("phylTree.conf",file), ("diagsList",file), ("ancestr",str), ("svm.model",file)], \
	[("useOutgroups",bool,True), ("useLonelyGenes",bool,False), ("walktrapLength",int,5), \
	("genesFile",str,"~/work/data/genes/genes.%s.list.bz2"), \
	("ancGenesFile",str,"~/work/data/ancGenes/ancGenes.%s.list.bz2")], \
	__doc__ \
)


#  Chargement et initialisation
################################
phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])
anc = phylTree.officialName[arguments["ancestr"]]
lstGenesAnc = [g.names for g in utils.myGenomes.Genome(arguments["ancGenesFile"] % phylTree.fileName[anc]).lstGenes[None]]
lstNoeudsFils = [a for (a,_) in phylTree.items[anc]]
lstEspParNoeudsFils = [phylTree.allDescendants[f] for f in lstNoeudsFils]
if (anc in phylTree.parent) and arguments["useOutgroups"]:
	lstNoeudsFils.append( phylTree.parent[anc][0] )
	lstEspParNoeudsFils.append( phylTree.allDescendants[phylTree.root].difference(phylTree.allDescendants[anc]) )

# Chargement des genomes si necessaire
#######################################
lengths = []
dicGenes = {}
if arguments["useLonelyGenes"]:
	for e in phylTree.listSpecies:
		genome = utils.myGenomes.Genome(arguments["genesFile"] % phylTree.fileName[e])
		lengths.append(len(genome.lstChr))
		if arguments["useLonelyGenes"]:
			for (g,(c,_)) in genome.dicGenes.iteritems():
				dicGenes[g] = (e,c)


# Les diagonales et les scores possibles
#########################################
lstDiags = utils.myDiags.loadDiagsFile(arguments["diagsList"], [anc], phylTree.officialName)[anc]

# On doit rajouter les genes non presents dans des diagonales
if arguments["useLonelyGenes"]:
	lstDiags.extend(checkLonelyGenes())

dicGenes.clear()

svm_model = utils.svm.svm_model(arguments["svm.model"])

# On calcule les scores
print >> sys.stderr, "Calcul de la matrice ...",
edges = collections.defaultdict(dict)
for i1 in xrange(len(lstDiags)):
	(d1,ec1,e1,_,_) = lstDiags[i1]
	for i2 in xrange(i1):
		(d2,ec2,e2,_,_) = lstDiags[i2]
		x = calcProba( e1.intersection(e2), frozenset([e for (e,_) in ec1.intersection(ec2)]) )
		if x > 0: #arguments["scoreThreshold"]:
			edges[i1][i2] = edges[i2][i1] = x
print >> sys.stderr, "OK"

res = utils.walktrap.doWalktrap(edges, showProgress=True, randomWalksLength=arguments["walktrapLength"])

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
		#(alpha,score,clust,lonely) = interessant[-1]
		(alpha,score,clust,lonely) = interessant[0]
		# Au choix, on prend la version la moins fusionnee
		print >> sys.stderr, "+%d/%d/%f/%f" % (len(clust), len(lonely), alpha, score)
		clusters.extend(clust)
		#clusters.append(lonely)
		relev.append(score)
# -> clusters contient la repartition des diagonales
print >> sys.stderr

print >> sys.stderr, "Impression des %d chromosomes ancestraux (score = %f)..." % (len(clusters),utils.myMaths.myStats.mean(relev)),

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

print >> sys.stderr, "OK"

