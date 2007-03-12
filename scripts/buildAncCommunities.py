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
def loadDiagsFile(nom, ancName):
	
	print >> sys.stderr, "Chargement du fichier de diagonales ...",
	f = utils.myTools.myOpenFile(nom, 'r')
	lst = []
	for l in f:

		ct = l[:-1].split('\t')
		anc = ct[0]
		if anc != ancName:
			continue
		l = int(ct[1])
		d = [int(x) for x in ct[2].split(' ')]
		esp = set([])
		if len(ct[3]) > 0:
			esp.update( set([tuple(x.split('/')) for x in ct[3].split('|')]) )
		if len(ct) == 5 and len(ct[4]) > 0:
			#print >> sys.stderr, ct[4], ct[4].split('|')
			esp.update( set([tuple(x.split('/')) for x in ct[4].split('|')]) )
		#print >> sys.stderr, esp
		esp = set([(phylTree.officialName[e],c) for (e,c) in esp])
		lst.append( (l, d, esp) )

	f.close()
	print >> sys.stderr, "OK (%d diagonales)" % len(lst)
	return lst


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
genesAnc = utils.myGenomes.AncestralGenome(options["ancGenesFile"] % phylTree.fileName[options.ancestr])
lstGenesAnc = genesAnc.lstGenes[utils.myGenomes.Genome.defaultChr]

# On separe les especes en trois
(fils1, fils2) = phylTree.branchesSpecies[options["ancestr"]]
outgroup = phylTree.outgroupSpecies[options["ancestr"]]

print >> sys.stderr, fils1, fils2, outgroup

# Les genomes modernes
#geneBank = utils.myGenomes.GeneBank(noms_fichiers["genesList.conf"], phylTree.species[options["ancestr"]])

# Les diagonales
lstDiags = loadDiagsFile(noms_fichiers["diagsList"], options["ancestr"])


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
	#s = sum([1./math.log(a) for (_,a) in outgroup])
	s = sum([1./float(a) for (_,a) in outgroup])
	for (e,a) in outgroup:
		#calcPoidsFils(e, 1. / (math.log(a)*s))
		calcPoidsFils(e, 1. / (float(a)*s))

calcPoids(options["ancestr"])

print >> sys.stderr, dicPoidsEspeces

def calcScore(i1, i2):

	global lstDiags, dicPoidsEspeces, fils1, fils2, outgroup
	
	(_,_,e1) = lstDiags[i1]
	(_,_,e2) = lstDiags[i2]
	communEsp = set([e for (e,_) in e1.intersection(e2)])

	propF1 = sum([dicPoidsEspeces[e] for e in communEsp.intersection(fils1)])
	propF2 = sum([dicPoidsEspeces[e] for e in communEsp.intersection(fils2)])
	propOut = sum([dicPoidsEspeces[e] for e in communEsp.intersection(outgroup)])
	
	return propF1*propF2 + propF1*propOut + propF2*propOut
	

lstLstComm = utils.myCommunities.launchCommunitiesBuild(items = range(len(lstDiags)), scoreFunc = calcScore)
clusters = []

# Chaque composante connexe
for lst in lstLstComm:
	# relevance >= 0.3 && noeudsOublies = 0
	#interessant = [comm for comm in lst if (len(comm[3]) == 0) and (comm[1] >= 0.3)]
	interessant = lst
	interessant.sort(key = operator.itemgetter(1), reverse = True)
	if len(interessant) == 0:
		clusters.append(utils.myMaths.flatten(lst[0][2])+lst[0][-1])
	else:
		clusters.extend(interessant[0][2])


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
