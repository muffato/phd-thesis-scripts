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
		esp = set([tuple(x.split('/')) for x in ct[3].split()])
		diagEntry[anc].append( (l, esp) )

	f.close()
	print >> sys.stderr, "OK"
	

########
# MAIN #
########

# Arguments
(noms_fichiers, options) = utils.myTools.checkArgs( \
	["phylTree.conf", "diagsList"], \
	[("ancestr",str,"")], \
	__doc__ \
)

# L'arbre phylogenetique
phylTree = utils.myBioObjects.PhylogeneticTree(noms_fichiers["phylTree.conf"])
listEspeces = phylTree.species[phylTree.root] + ['Opossum']

# Les genes ancestraux
diagEntry = {}
for anc in phylTree.items:
	diagEntry[anc] = []
loadDiagsFile(noms_fichiers["diagsList"], diagEntry)


# Impression du resultat
nn = max([len(anc) for anc in phylTree.items])
for anc in phylTree.items:
	continue
	chr50 = []
	chr250 = []
	for (l,d,esp) in diagEntry[anc]:
		#print "%s\t%d\t%s" % (anc, len(d), " ".join([str(x) for x in d]))
		if l >= 50:
			chr50.append(l)
		if l >= 250:
			chr250.append(l)
	s = utils.myMaths.fullStats(chr50)
	print ("%" + ("%d" % nn) + "s\t>50\t%3d chromosomes\t%5d genes\tmin=%4d\tmax=%5d\tmean=%.2f\tmedian=%s\tstddev=%.2f") % \
	(anc, s[5], sum(chr50), s[3], s[4], s[0], s[2], s[1])
	s = utils.myMaths.fullStats(chr250)
	print ("%" + ("%d" % nn) + "s\t>250\t%3d chromosomes\t%5d genes\tmin=%4d\tmax=%5d\tmean=%.2f\tmedian=%s\tstddev=%.2f") % \
	(anc, s[5], sum(chr250), s[3], s[4], s[0], s[2], s[1])


anc = options["ancestr"]
(fils1, fils2) = phylTree.branchesSpecies[anc]
outgroup = set(phylTree.species[phylTree.root]).difference(phylTree.species[anc])
lst = diagEntry[anc]
nb = len(lst)
combin = utils.myTools.myCombinator([])
dicAretes = dict([(i,{}) for i in xrange(nb)])
for (i1,i2) in utils.myTools.myMatrixIterator(nb, nb, utils.myTools.myMatrixIterator.StrictUpperMatrix):
	(_,e1) = lst[i1]
	(_,e2) = lst[i2]
	commun = e1.intersection(e2)
	communEsp = set([e for (e,_) in commun])
	
	propF1 = float(len(communEsp.intersection(fils1))) / float(len(fils1))
	propF2 = float(len(communEsp.intersection(fils2))) / float(len(fils2))
	propOut = float(len(communEsp.intersection(outgroup))) / float(len(outgroup))
	propOut = 0
	
	score = propF1*propF2 + propF1*propOut + propF2*propOut
	
	if score > 0:
		combin.addLink([i1,i2])
		dicAretes[i1][i2] = score
		dicAretes[i2][i1] = score
		#tout.discard(i1)
		#tout.discard(i2)
		#print i1, i2, round(score, 2)
		#print i1, i2, round(math.exp(score), 2)


# On traite chaque composante connexe
indComp = 0
for g in combin:
	indComp += 1
	nb = len(g)
	f = open('/users/ldog/muffato/work/temp/walktrap/%s/nodes.%d' % (anc, indComp), 'w')
	for i in xrange(nb):
		print >> f, i, g[i]
	f.close()
	f = open('/users/ldog/muffato/work/temp/walktrap/%s/graph.%d' % (anc, indComp), 'w')
	for (i1,i2) in utils.myTools.myMatrixIterator(nb, nb, utils.myTools.myMatrixIterator.StrictUpperMatrix):
		if g[i2] in dicAretes[g[i1]]:
			print >> f, i1, i2, dicAretes[g[i1]][g[i2]]
	f.close()
	
