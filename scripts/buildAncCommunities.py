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
		
		if not l.startswith(ancName):
			continue
		ct = l[:-1].replace("_random", "").split('\t')
		d = [int(x) for x in ct[2].split(' ')]
		esp = set()
		if len(ct[3]) > 0:
			esp.update( set([tuple(x.split('/')) for x in ct[3].split('|')]) )
		if len(ct) == 5 and len(ct[4]) > 0:
			esp.update( set([tuple(x.split('/')) for x in ct[4].split('|')]) )
		esp = set([(phylTree.officialName[e],c) for (e,c) in esp if (e in phylTree.officialName) and ('Un' not in c)])
		lst.append( (d, esp) )

	f.close()
	print >> sys.stderr, "OK (%d diagonales)" % len(lst)
	return lst


#
# Rajoute les genes non presents dans des diagonales, en les considerant comme des diagonales de 1
#
def checkLonelyGenes():


	if len(phylTree.dicGenomes) == 0:
		phylTree.loadAllSpeciesSince(None, options["genesFile"])
	genesAnc = {}
	for a in phylTree.items:
		if phylTree.getFirstParent(options["ancestr"], a) == a:
			genesAnc[a] = utils.myGenomes.AncestralGenome(options["ancGenesFile"] % phylTree.fileName[a])
	
	print >> sys.stderr, "Ajout des genes solitaires ...",
	# Les genes seuls vont devenir des diagonales de 1
	genesSeuls = set(xrange(len(lstGenesAnc)))
	for (d,_) in lstDiags:
		genesSeuls.difference_update(d)
	nb = 0
	new = {}
	for i in genesSeuls:

		lst = []
		# Les especes filles
		esp = [phylTree.dicGenes[s] for s in lstGenesAnc[i].names if s in phylTree.dicGenes]

		for e in phylTree.listSpecies:
			if e in outgroup:
				# Le gene chez l'ancetre commun
				a = phylTree.getFirstParent(options["ancestr"], e)
				g = genesAnc[a].getPosition(lstGenesAnc[i])
				# Cas d'un gene specifique de la lignee
				if len(g) == 0:
					continue
				# Le gene dans l'autre espece
				tmp = [c for (c,_) in phylTree.dicGenomes[e].getPosition(genesAnc[a].lstGenes[utils.myGenomes.Genome.defaultChr][g[0][1]])]
			else:
				# Le gene dans l'autre espece
				tmp = [c for (x,c,_) in esp if x == e]
			
			# Au final
			if len(tmp) == 1:
				c = str(tmp[0]).replace("_random","")
				if 'Un' not in c:
					lst.append( (e,c) )
		
		if len(lst) >= 2:
			lst = tuple(sorted(set(lst)))
			new[lst] = new.get(lst,[]) + [i]
			nb += 1
	
	lstDiags.extend([(d,set(e)) for (e,d) in new.iteritems()])
	print >> sys.stderr, "%d (%d) OK" % (nb,len(new))




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
		outgroup.extend([(e,1./float(2*phylTree.ages[par]-phylTree.ages[node])) for (e,_) in phylTree.items[par] if e != anc])
		anc = par
	s = sum([a for (_,a) in outgroup])
	for (e,a) in outgroup:
		calcPoidsFils(e, a/s)


########
# MAIN #
########

# Arguments
(noms_fichiers, options) = utils.myTools.checkArgs( \
	["phylTree.conf", "diagsList"], \
	[("ancestr",str,""), ("alreadyBuiltAnc",str,""), \
	("useOutgroups",bool,True), ("useLonelyGenes",bool,False), ("weightNbChromosomes",bool,True), \
	("genesFile",str,"~/work/data/genes/genes.%s.list.bz2"), \
	("ancGenesFile",str,"~/work/data/ancGenes/ancGenes.%s.list.bz2")], \
	__doc__ \
)

#  Chargement et initialisation
phylTree = utils.myBioObjects.PhylogeneticTree(noms_fichiers["phylTree.conf"])
genesAnc = utils.myGenomes.AncestralGenome(options["ancGenesFile"] % phylTree.fileName[options.ancestr])
lstDiags = loadDiagsFile(noms_fichiers["diagsList"], options["ancestr"])
lstGenesAnc = genesAnc.lstGenes[utils.myGenomes.Genome.defaultChr]

lstToutesEspeces = set(phylTree.listSpecies)
lstEspOutgroup = set(phylTree.outgroupSpecies[options["ancestr"]])
lstEspFilles = set(phylTree.species[options["ancestr"]])
lstNoeudsFils = phylTree.branches[options["ancestr"]]
lstEspParNoeudsFils = [set(s) for s in phylTree.branchesSpecies[options["ancestr"]]]

filsAnc = phylTree.branches[options["ancestr"]]
filsEsp = [set(s) for s in phylTree.branchesSpecies[options["ancestr"]]]
outgroup = set(phylTree.outgroupSpecies[options["ancestr"]])


# L'apport en terme de couverture de l'arbre phylogenetique
dicPoidsEspeces = dict.fromkeys(phylTree.listSpecies, 0.)
def calcPoidsFils(node, calc):
	if node in phylTree.listSpecies:
		dicPoidsEspeces[node] = calc
	else:
		poids = calc / float(len(phylTree.items[node]))
		for f in phylTree.branches[node]:
			calcPoidsFils(f, poids)
phylTree.travelFunc(options["ancestr"], calcPoidsFils, options["useOutgroups"])

# Le taux de precision de chaque espece
dicNbChr = dict.fromkeys(phylTree.listSpecies, 0.)
dicNbChr2 = dict.fromkeys(phylTree.listSpecies, 0.)
if options["weightNbChromosomes"]:
	phylTree.loadAllSpeciesSince(None, options["genesFile"])
	for e in phylTree.listSpecies:
		s = len(phylTree.dicGenomes[e].lstChr)
		if s > 0:
			dicNbChr[e] = 1. / float(s)
	alpha = min(dicNbChr.values()) * max(dicNbChr.values())
	for e in phylTree.listSpecies:
		dicNbChr2[e] = len(phylTree.dicGenomes[e].lstChr) * alpha

# On doit rajouter les genes non presents dans des diagonales
if options["useLonelyGenes"]:
	checkLonelyGenes()

phylTree.dicGenomes.clear()

# On doit noter les chromosomes des diagonales sur ces ancetres deja construits
if options["alreadyBuiltAnc"] != "":
	genAlready = {}
	for f in filsAnc:
		if f in phylTree.items:
			try:
				genAlready[f] = utils.myGenomes.AncestralGenome(options["alreadyBuiltAnc"] % phylTree.fileName[f], chromPresents=True)
				dicPoidsEspeces[f] = 1. / (1. - float(phylTree.ages[f])/float(phylTree.ages[options["ancestr"]]))
			except IOError:
				genAlready.pop(f, None)

	print >> sys.stderr, "Mise a jour des chromosomes des diagonales ...",
	for (d,e) in lstDiags:
		g = utils.myMaths.flatten([lstGenesAnc[i].names for i in d])
		for f in genAlready:
			e.update([(f,genAlready[f].dicGenes[s][0]) for s in g if s in genAlready[f].dicGenes])
	genAlready = {}
	print >> sys.stderr, "OK"


print >> sys.stderr, dicPoidsEspeces
print >> sys.stderr, dicNbChr


def calcScore(i1, i2):

	(_,e1) = lstDiags[i1]
	(_,e2) = lstDiags[i2]
	comparedEsp = set([e for (e,_) in e1]).intersection([e for (e,_) in e2])
	communEsp = set([e for (e,_) in e1.intersection(e2)])

	propF = []
	for i in xrange(len(filsEsp)):
		# Chez l'ancetre du dessous - meme chr
		if filsAnc[i] in communEsp:
			propF.append( dicPoidsEspeces[filsAnc[i]] )
		# Chez l'ancetre du dessous - diff chr
		elif filsAnc[i] in comparedEsp:
			propF.append(0.)
		# Sinon, on revient aux genomes modernes
		else:
			s = sum([dicPoidsEspeces[e]*(1.-dicNbChr[e]) for e in filsEsp[i].intersection(communEsp)])
			if s > 0:
				s += sum([dicPoidsEspeces[e]*dicNbChr2[e] for e in filsEsp[i].difference(communEsp)])
			propF.append(s)

		
	propOut = sum([dicPoidsEspeces[e]*(1.-dicNbChr[e]) for e in communEsp.intersection(outgroup)])
	if propOut > 0:
		propOut += sum([dicPoidsEspeces[e]*dicNbChr2[e] for e in communEsp.difference(outgroup)])
	
	s = sum(propF) * propOut
	for (i1,i2) in utils.myTools.myMatrixIterator(len(filsEsp), len(filsEsp), utils.myTools.myMatrixIterator.StrictUpperMatrix):
		s += propF[i1] * propF[i2]
	
	return s
	
lstLstComm = utils.myCommunities.launchCommunitiesBuild(items = range(len(lstDiags)), scoreFunc = calcScore)
clusters = []

# Chaque composante connexe
for lst in lstLstComm:
	#interessant = [comm for comm in lst if (len(comm[3]) == 0) and (comm[1] >= 0.3)]
	interessant = lst
	#interessant = [comm for comm in lst if (len(comm[3]) == 0)]
	interessant.sort(key = operator.itemgetter(1), reverse = True)
	#interessant = []
	if len(interessant) == 0:
		sys.stderr.write('-')
		clusters.append(utils.myMaths.flatten(lst[0][2])+lst[0][-1])
	else:
		sys.stderr.write('+')
		clusters.extend(interessant[0][2])
print >> sys.stderr

print >> sys.stderr, "Impression des chromosomes ancestraux ...",
lstChr = []
ind = 0
for c in clusters:
	ind += 1
	lst = set()
	for i in c:
		lst.update(lstDiags[i][0])
	lstChr.append(lst)

inter = set()
for (i1,i2) in utils.myTools.myMatrixIterator(len(lstChr), len(lstChr), utils.myTools.myMatrixIterator.StrictUpperMatrix):
	inter.update(lstChr[i1].intersection(lstChr[i2]))

chrIndex = 0
for c in lstChr:
	chrIndex += 1
	for i in c:
		if i not in inter:
			print chrIndex, " ".join(lstGenesAnc[i].names)


print >> sys.stderr, "OK"
