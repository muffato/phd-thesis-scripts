#! /users/ldog/muffato/python -OO

__doc__ = """
A partir de toutes les diagonales extraites entre les especes,
  reconstruit les chromosomes (ou scaffold) de chaque ancetre.
"""


##################
# INITIALISATION #
##################

# Librairies
import sys
import math
import operator
import utils.myGenomes
import utils.myTools
import utils.myMaths
import utils.myCommunities

#############
# FONCTIONS #
#############

# Charge le fichier de toutes les diagonales (supposees non chevauchantes)
def loadDiagsFile(nom, diagEntry):
	
	print >> sys.stderr, "Chargement du fichier de diagonales ...",
	f = utils.myTools.myOpenFile(nom, 'r')
	for l in f:

		ct = l.split('\t')
		anc = ct[0]
		#l = int(ct[1])
		d = [int(x) for x in ct[1].split(' ')]
		l = len(d)
		esp = set([tuple(x.split('/')) for x in ct[2].split()])
		if len(ct) == 5:
			esp.update( set([tuple(x.split('/')) for x in ct[4].split()]) )
		esp = set([(phylTree.officialName[e],c) for (e,c) in esp])
		diagEntry[anc].append( (l, d, esp) )

	f.close()
	print >> sys.stderr, "OK"
	

########
# MAIN #
########

# Arguments
(noms_fichiers, options) = utils.myTools.checkArgs( \
	["phylTree.conf", "diagsList"], \
	[("ancestr",str,""), ("ancGenesFile",str,"~/work/data/ancGenes/ancGenes.%s.list.bz2")], \
	__doc__ \
)

# L'arbre phylogenetique
phylTree = utils.myBioObjects.PhylogeneticTree(noms_fichiers["phylTree.conf"])

# Les genes ancestraux
genesAnc = utils.myGenomes.AncestralGenome(options["ancGenesFile"] % options["ancestr"])
lstGenesAnc = genesAnc.lstGenes[utils.myGenomes.Genome.defaultChr]

# On separe les especes en trois
(fils1, fils2) = phylTree.branchesSpecies[options["ancestr"]]
outgroup = phylTree.outgroupSpecies[options["ancestr"]]

print >> sys.stderr, fils1, fils2, outgroup

# Les genomes modernes
#geneBank = utils.myGenomes.GeneBank(noms_fichiers["genesList.conf"], phylTree.species[options["ancestr"]])

# Les diagonales
diagEntry = {}
for anc in phylTree.items:
	diagEntry[anc] = []
loadDiagsFile(noms_fichiers["diagsList"], diagEntry)
lstDiags = diagEntry[options["ancestr"]]
del diagEntry


print >> sys.stderr, len(lstDiags)
# Les genes seuls vont devenir des diagonales de 1
#genesSeuls = set(xrange(len(lstGenesAnc)))
#for (_,d,_) in lstDiags:
#	genesSeuls.difference_update(d)
#for i in genesSeuls:
#	esp = [geneBank.dicGenes[s] for s in lstGenesAnc[i].names]
#	if len(esp) >= 2:
#		lstDiags.append( (1,[i],set([(e,c) for (e,c,_) in esp])) )
#	#print "ajout de ", lstDiags[-1]
#
#print >> sys.stderr, len(lstDiags)
#sys.exit(0)


# Il faut calculer l'apport de chaque espece
dicPoidsEspeces = {}
def calcPoidsFils(node, calc):
	if node in phylTree.listSpecies:
		dicPoidsEspeces[node] = calc
	else:
		poids = calc / float(len(phylTree.items[node]))
		for f in phylTree.branches[node]:
			calcPoidsFils(f, poids)

def calcPoids(node):
	# Les fils a egalite avec un poids de 1
	for f in phylTree.branches[node]:
		calcPoidsFils(f, 1.)
	outgroup = []
	anc = node
	while anc in phylTree.parent:
		par = phylTree.parent[anc]
		outgroup.extend([(e,2*phylTree.ages[par]-phylTree.ages[node]) for (e,_) in phylTree.items[par] if e != anc])
		anc = par
	s = sum([1./math.log(a) for (_,a) in outgroup])
	for (e,a) in outgroup:
		calcPoidsFils(e, 1. / (math.log(a)*s))

calcPoids(options["ancestr"])

print >> sys.stderr, dicPoidsEspeces
#dicPoidsEspeces = dict([(phylTree.commonNames[esp][0],dicPoidsEspeces[esp]) for esp in dicPoidsEspeces])
#print >> sys.stderr, dicPoidsEspeces

def calcScore(i1, i2):

	global lstDiags, dicPoidsEspeces, fils1, fils2, outgroup
	
	(_,_,e1) = lstDiags[i1]
	(_,_,e2) = lstDiags[i2]
	communEsp = set([e for (e,_) in e1.intersection(e2)])

	propF1 = sum([dicPoidsEspeces[e] for e in communEsp.intersection(fils1)])
	propF2 = sum([dicPoidsEspeces[e] for e in communEsp.intersection(fils2)])
	propOut = sum([dicPoidsEspeces[e] for e in communEsp.intersection(outgroup)])
	
	#print propF1*propF2 + propF1*propOut + propF2*propOut
	return propF1*propF2 + propF1*propOut + propF2*propOut
	

(_,clusters) = utils.myCommunities.launchCommunitiesBuild(len(lstDiags), calcScore, keepLonelyNodes = False, minRelevance = 0.3, minCoverage = 0, bestRelevance = True)

#lstCommunities.sort(key = operator.itemgetter(1), reverse = True)
#if len(lstCommunities) == 0:
#	clusters = [range(len(lstDiags))]
#elif lstCommunities[0][1] > 0.3:
#	clusters = lstCommunities[0][2]
#else:
#	clusters = [range(len(lstDiags))]

print >> sys.stderr, "Impression des chromosomes ancestraux ...",
lstChr = []
for c in clusters:
	lst = set([])
	for i in c:
		lst.update(lstDiags[i][1])
	lstChr.append(lst)

for (i1,i2) in utils.myTools.myMatrixIterator(len(lstChr), len(lstChr), utils.myTools.myMatrixIterator.StrictUpperMatrix):
	inter = lstChr[i1].intersection(lstChr[i2])
	lstChr[i1].difference_update(inter)
	lstChr[i2].difference_update(inter)

chrIndex = 0
for c in lstChr:
	chrIndex += 1
	for i in c:
		print chrIndex, " ".join(lstGenesAnc[i].names)


print >> sys.stderr, "OK"
