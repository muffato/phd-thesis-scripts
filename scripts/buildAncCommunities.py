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
		d = [int(x) for x in ct[2].split(' ')]
		esp = set()
		if len(ct[3]) > 0:
			esp.update( set([tuple(x.split('/')) for x in ct[3].split('|')]) )
		if len(ct) == 5 and len(ct[4]) > 0:
			esp.update( set([tuple(x.split('/')) for x in ct[4].split('|')]) )
		esp = set([(phylTree.officialName[e],c) for (e,c) in esp if e in phylTree.officialName])
		lst.append( (d, esp) )

	f.close()
	print >> sys.stderr, "OK (%d diagonales)" % len(lst)
	return lst


#
# Rajoute les genes non presents dans des diagonales, en les considerant comme des diagonales de 1
#
def checkLonelyGenes():

	# Les genomes modernes
	phylTree.loadAllSpeciesSince(None, options["genesFile"])
	genesAnc = {}
	for a in phylTree.items:
		if phylTree.getFirstParent(options["ancestr"], a) == a:
			genesAnc[a] = utils.myGenomes.AncestralGenome(options["ancGenesFile"] % phylTree.fileName[a])
	
	# Les genes seuls vont devenir des diagonales de 1
	genesSeuls = set(xrange(len(lstGenesAnc)))
	for (d,_) in lstDiags:
		genesSeuls.difference_update(d)
	nb = 0
	new = {}
	for i in genesSeuls:
		esp = [phylTree.dicGenes[s] for s in lstGenesAnc[i].names]
		lst = [(e,c) for (e,c,_) in esp]
		# Les outgroup
		
		for e in phylTree.outgroupSpecies[options["ancestr"]]:
			# L'ancetre commun
			a = phylTree.getFirstParent(options["ancestr"], e)
			# On va intersecter les chromosomes des orthologues de chaque gene
			poss = set()
			# Le gene chez l'ancetre commun
			g = genesAnc[a].getPosition(lstGenesAnc[i])
			# Cas d'un gene specifique de la lignee
			if len(g) == 0:
				continue
			# Le gene dans l'autre espece
			tmp = [c for (c,_) in phylTree.dicGenomes[e].getPosition(genesAnc[a].lstGenes[utils.myGenomes.Genome.defaultChr][g[0][1]])]
			if len(tmp) == 1:
				lst.append( (e,tmp[0]) )
		
		if len(lst) >= 2:
			lst = tuple(sorted(set(lst)))
			new[lst] = new.get(lst,[]) + [i]
			nb += 1
	
	del phylTree.dicGenomes

	lst = [(d,set(e)) for (e,d) in new.iteritems()]
	print >> sys.stderr, "Ajout des %d (%d) genes solitaires" % (nb,len(new))
	return lst






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
	if not options["useOutgroups"]:
		return
	outgroup = []
	anc = node
	while anc in phylTree.parent:
		par = phylTree.parent[anc]
		outgroup.extend([(e,phylTree.ages[par]-phylTree.ages[node]) for (e,_) in phylTree.items[par] if e != anc])
		#outgroup.extend([(e,2*phylTree.ages[par]-phylTree.ages[node]) for (e,_) in phylTree.items[par] if e != anc])
		anc = par
	#par = phylTree.parent[anc]
	#outgroup.extend([(e,phylTree.ages[par]-phylTree.ages[node]) for (e,_) in phylTree.items[par] if e != anc])
	#s = sum([1./math.log(a) for (_,a) in outgroup])
	s = sum([1./float(a) for (_,a) in outgroup])
	for (e,a) in outgroup:
		#calcPoidsFils(e, 1. / (math.log(a)*s))
		#calcPoidsFils(e, 0)
		calcPoidsFils(e, 1. / (float(a)*s))

########
# MAIN #
########

# Arguments
(noms_fichiers, options) = utils.myTools.checkArgs( \
	["phylTree.conf", "diagsList"], \
	[("useOutgroups",bool,True), ("useLonelyGenes",bool,False), ("ancestr",str,""), ("alreadyBuildAnc",str,""), \
	("genesFile",str,"~/work/data/genes/genes.%s.list.bz2"), \
	("ancGenesFile",str,"~/work/data/ancGenes/ancGenes.%s.list.bz2")], \
	__doc__ \
)

#  Chargement et initialisation
phylTree = utils.myBioObjects.PhylogeneticTree(noms_fichiers["phylTree.conf"])
genesAnc = utils.myGenomes.AncestralGenome(options["ancGenesFile"] % phylTree.fileName[options.ancestr])
lstDiags = loadDiagsFile(noms_fichiers["diagsList"], options["ancestr"])
lstGenesAnc = genesAnc.lstGenes[utils.myGenomes.Genome.defaultChr]
filsAnc = phylTree.branches[options["ancestr"]]
filsEsp = phylTree.branchesSpecies[options["ancestr"]]
outgroup = phylTree.outgroupSpecies[options["ancestr"]]


# On doit rajouter les genes non presents dans des diagonales
if options["useLonelyGenes"]:
	checkLonelyGenes()


# On doit noter les chromosomes des diagonales sur ces ancetres deja construits
if options["alreadyBuildAnc"] != "":
	genAlready = {}
	for f in filsAnc:
		if f in phylTree.items:
			genAlready[f] = utils.myGenomes.AncestralGenome(options["alreadyBuildAnc"] % phylTree.fileName[f])
	print >> sys.stderr, "Mise a jour des chromosomes des diagonales ...",
	for (d,e) in lstDiags:
		g = utils.myMaths.flatten([lstGenesAnc[i].names for i in d])
		for f in genAlready:
			e.update([(f,genAlready[f].dicGenes[s][0]) for s in g if s in genAlready[f].dicGenes])
	print >> sys.stderr, "OK"


# Il faut calculer l'apport de chaque espece
dicPoidsEspeces = dict.fromkeys(phylTree.listSpecies, 0.)
calcPoids(options["ancestr"])

print >> sys.stderr, dicPoidsEspeces

def calcScore(i1, i2):

	(_,e1) = lstDiags[i1]
	(_,e2) = lstDiags[i2]
	communEsp = set([e for (e,_) in e1.intersection(e2)])

	propF = [max(sum([dicPoidsEspeces[e] for e in communEsp.intersection(filsEsp[i])]), filsAnc[i] in communEsp) for i in xrange(len(filsEsp))]
	propOut = sum([dicPoidsEspeces[e] for e in communEsp.intersection(outgroup)])
	
	s = 0
	for i in xrange(len(filsEsp)):
		s += propF[i] * propOut
		for j in xrange(i+1, len(filsEsp)):
			s += propF[i] * propF[j]
	return s
	

lstLstComm = utils.myCommunities.launchCommunitiesBuild(items = range(len(lstDiags)), scoreFunc = calcScore)
clusters = []

# Chaque composante connexe
for lst in lstLstComm:
	# relevance >= 0.3 && noeudsOublies = 0
	#interessant = [comm for comm in lst if (len(comm[3]) == 0) and (comm[1] >= 0.3)]
	#interessant = lst
	#interessant = [comm for comm in lst if (len(comm[3]) == 0) and (comm[1] >= 0.2)]
	interessant = [comm for comm in lst if (len(comm[3]) == 0)]
	interessant.sort(key = operator.itemgetter(1), reverse = True)
	if len(interessant) == 0:
		print >> sys.stderr, "souci"
		clusters.append(utils.myMaths.flatten(lst[0][2])+lst[0][-1])
	else:
		clusters.extend(interessant[0][2])


print >> sys.stderr, "Impression des chromosomes ancestraux ...",
lstChr = []
for c in clusters:
	lst = set()
	for i in c:
		lst.update(lstDiags[i][0])
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
