#! /users/ldog/muffato/python

__doc__ = """
A partir de toutes les diagonales extraites entre les especes,
  reconstruit les chromosomes (ou scaffold) de chaque ancetre.
"""

import sys
import random
import itertools

import utils.svm
import utils.myTools
import utils.myDiags
import utils.myMaths
import utils.myGenomes
import utils.myPhylTree


it = utils.myTools.myIterator.tupleOnStrictUpperList

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
	[("phylTree.conf",file), ("diagsList",file), ("realGenome",file), ("ancestr",str)], \
	[("useOutgroups",bool,True), ("useLonelyGenes",bool,False), \
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
phylTree.initCalcDist(anc, arguments["useOutgroups"])
calcPoidsFils(anc, float(len(phylTree.tmpItems[0])))


# Les diagonales et les scores possibles
#########################################
lstDiags = utils.myDiags.loadDiagsFile(arguments["diagsList"], [anc], phylTree.officialName)[anc]

# On doit rajouter les genes non presents dans des diagonales
if arguments["useLonelyGenes"]:
	lstDiags.extend(checkLonelyGenes())


dicGenes.clear()

# On charge le genome pour savoir le chromosome de chaque diagonale
chromAssocGene = {}
gen = utils.myGenomes.Genome(arguments["realGenome"])
for (i,g) in enumerate(lstGenesAnc):
	x = gen.getPosition(g)
	if len(x) != 1:
		print >> sys.stderr, g, i, x
	else:
		(c,_) = x.pop()
		chromAssocGene[i] = c

# On charge les diagonales

numDiagsOK = set()
chromAssocDiag = {}
for (i,(d,_,_,_,_)) in enumerate(lstDiags):
	if len(set([chromAssocGene[g] for g in d])) == 1:
		numDiagsOK.add(i)
		chromAssocDiag[i] = chromAssocGene[d[0]]


print >> sys.stderr, "Calcul de la matrice ...", len(numDiagsOK)
edges = utils.myTools.defaultdict(dict)
nfacile = 0
goods = []
bads = []
nn = [0] * 4

comparedEsp = frozenset(phylTree.listSpecies)
for distr in itertools.product(*([[True,False]] * len(phylTree.listSpecies))):
	continue
	communEsp = frozenset([e for (b,e) in zip(distr,phylTree.listSpecies) if b])
	values = dict.fromkeys(phylTree.listSpecies, -1)
	for e in communEsp:
		values[e] = 1
	s = phylTree.calcWeightedValue(values, -1, None)[phylTree.indNames[anc]]
	print int(s>0), " ".join(["%d:%d" % (phylTree.indNames[e],v) for (e,v) in values.iteritems()])
	
#sys.exit(0)


def constr():
	return [0,0]
results = utils.myTools.defaultdict(constr)
for i1 in xrange(len(lstDiags)):
	#continue
	if i1 not in numDiagsOK:
		continue
	(d1,ec1,e1,_,_) = lstDiags[i1]
	for i2 in xrange(i1):
		if i2 not in numDiagsOK:
			continue
		(d2,ec2,e2,_,_) = lstDiags[i2]
		comparedEsp = e1.intersection(e2)
		communEsp = frozenset([e for (e,_) in ec1.intersection(ec2)])

		#x = [len(l.intersection(communEsp)) for l in lstEspParNoeudsFils]
		#x = [len(l.intersection(comparedEsp)) for l in lstEspParNoeudsFils]
		#x = [int(chromAssocDiag[i1] == chromAssocDiag[i2]), len(x)-x.count(0)] + x
		#print utils.myFile.myTSV.printLine(x, ' ')
		#continue
		#x = [int(e in comparedEsp) for e in phylTree.listSpecies] + [int(e in communEsp) for e in phylTree.listSpecies]

		results[(communEsp,comparedEsp)][int(chromAssocDiag[i1] == chromAssocDiag[i2])] += 1
		continue
		
		if chromAssocDiag[i1] == chromAssocDiag[i2]:
			results[(communEsp,comparedEsp)] += 1
		else:
			results[(communEsp,comparedEsp)] -= 1
		continue
		
		if len([l for l in lstEspParNoeudsFils if len(l.intersection(communEsp)) > 0]) < 2:
			if chromAssocDiag[i1] == chromAssocDiag[i2]:
				nn[0] += 1
				#results[(communEsp,comparedEsp)] += 1
			else:
				nn[1] += 1
				results[(communEsp,comparedEsp)] -= 1
		else:
			if chromAssocDiag[i1] == chromAssocDiag[i2]:
				#results[(communEsp,comparedEsp)] += 1
				nn[2] += 1
			else:
				results[(communEsp,comparedEsp)] -= 1
				nn[3] += 1

print >> sys.stderr, "OK", nn
#sys.exit(0)

def transform((communEsp,comparedEsp), e):
	if e in communEsp:
		return 1
	elif e in comparedEsp:
		return -1
	else:
		return 0

tr = "- +"
for (x,n) in results.iteritems():
	print utils.myFile.myTSV.printLine( ["".join([tr[transform(x,e)+1] for e in phylTree.listSpecies]), n[1]] )
	#print "".join([tr[transform(x,e)+1] for e in phylTree.listSpecies])

sys.exit(0)

count = utils.myTools.defaultdict(int)
for x in results.itervalues():
	count[x] += 1
for (x,y) in count.iteritems():
	print x, y

sys.exit(0)

for (x,c) in results.iteritems():
	if c == 0:
		continue
	print int(c>0),
	for e in phylTree.listSpecies:
		print "%d:%d" % (phylTree.indNames[e],transform(x,e)),
	print



