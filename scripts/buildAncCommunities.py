#! /users/ldog/muffato/python -OO

__doc__ = """
A partir de toutes les diagonales extraites entre les especes,
  reconstruit les chromosomes (ou scaffold) de chaque ancetre.
"""


##################
# INITIALISATION #
##################

# Librairies
import math
import sys
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
		l = int(ct[1])
		d = [int(x) for x in ct[2].split(' ')]
		esp = set([tuple(x.split('/')) for x in ct[3].split()])
		if len(ct) == 5:
			esp.update( set([tuple(x.split('/')) for x in ct[4].split()]) )
		diagEntry[anc].append( (l, d, esp) )

	f.close()
	print >> sys.stderr, "OK"
	

########
# MAIN #
########

# Arguments
(noms_fichiers, options) = utils.myTools.checkArgs( \
	["phylTree.conf", "genesList.conf", "diagsList"], \
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

print >> sys.stderr, "debut", dicPoidsEspeces
print >> sys.stderr, fils1, fils2, outgroup

def calcScore(i1, i2):

	global lstDiags, dicPoidsEspeces, fils1, fils2, outgroup
	
	(_,_,e1) = lstDiags[i1]
	(_,_,e2) = lstDiags[i2]
	communEsp = set([e for (e,_) in e1.intersection(e2)])

	propF1 = sum([dicPoidsEspeces[e] for e in communEsp.intersection(fils1)])
	propF2 = sum([dicPoidsEspeces[e] for e in communEsp.intersection(fils2)])
	propOut = sum([dicPoidsEspeces[e] for e in communEsp.intersection(outgroup)])
	
	return propF1*propF2 + propF1*propOut + propF2*propOut
	

clusters = utils.myCommunities.launchCommunitiesBuild(len(lstDiags), calcScore)

for c in clusters:
	print c
