#! /users/ldog/muffato/python -OO

__doc__ = """
A partir de toutes les diagonales extraites entre les especes,
  reconstruit les chromosomes (ou scaffold) de chaque ancetre.
"""


##################
# INITIALISATION #
##################

# Librairies
import os
import sys
import operator
import utils.myGenomes
import utils.myTools
import utils.myMaths
import utils.myPhylTree
import utils.walktrap

#############
# FONCTIONS #
#############

# Charge le fichier de toutes les diagonales (supposees non chevauchantes)
def loadDiagsFile(nom, ancName):
	
	print >> sys.stderr, "Chargement du fichier de diagonales ...",
	f = utils.myTools.myOpenFile(nom, 'r')
	lst = []
	for l in f:
		# Selection de l'ancetre
		if not l.startswith(ancName):
			continue
		# On enleve les "_random" et on extrait chaque colonne
		ct = l[:-1].replace("_random", "").split('\t')
		# La diagonale
		d = [int(x) for x in ct[2].split(' ')]
		# Les especes qui l'ont vue (on enleve les "Un")
		esp = set()
		if len(ct[3]) > 0:
			esp.update( set([tuple(x.split('/')) for x in ct[3].split('|')]) )
		if len(ct) == 5 and len(ct[4]) > 0:
			esp.update( set([tuple(x.split('/')) for x in ct[4].split('|')]) )
		esp = set([(phylTree.officialName[e],c) for (e,c) in esp if (e in phylTree.officialName) and ('Un' not in c)])
		lst.append( (d,esp,set([e for (e,_) in esp])) )

	f.close()
	print >> sys.stderr, "OK (%d diagonales)" % len(lst)
	return lst


#
# Rajoute les genes non presents dans des diagonales, en les considerant comme des diagonales de 1
#
def checkLonelyGenes():

	# Charge les genomes des ancetres outgroup
	genesAnc = {}
	for a in phylTree.listAncestr:
		if (phylTree.dicParents[options["ancestr"]][a] == a) and (a != options["ancestr"]):
			genesAnc[a] = utils.myGenomes.Genome(options["ancGenesFile"] % phylTree.fileName[a])
	
	print >> sys.stderr, "Ajout des genes solitaires ...",
	# Les genes seuls vont devenir des diagonales de 1
	genesSeuls = set(xrange(len(lstGenesAnc)))
	for (d,_,_) in lstDiags:
		genesSeuls.difference_update(d)
	
	nb = 0
	new = utils.myTools.defaultdict(list)
	for i in genesSeuls:

		lst = []
		for e in phylTree.listSpecies:
			if e in outgroup:
				a = phylTree.dicParents[options["ancestr"]][e]
				names = genesAnc[a].getOtherNames(lstGenesAnc[i].names)
			else:
				names = lstGenesAnc[i].names
			tmp = [dicGenes[x] for x in names if x in dicGenes]
			tmp = set([c for (x,c) in tmp if x == e])

			# Au final
			if len(tmp) == 1:
				c = str(tmp.pop()).replace("_random","")
				if 'Un' not in c:
					lst.append( (e,c) )
		
		if len(lst) >= 2:
			lst = tuple(sorted(set(lst)))
			new[lst].append(i)
			nb += 1
	
	if nbDiagsIni + len(new) < 50000:
		lstDiags.extend([(d,set(esp),set([e for (e,_) in esp])) for (esp,d) in new.iteritems()])
		print >> sys.stderr, "%d (%d) OK" % (nb,len(new))
	else:
		print >> sys.stderr, "too many genes !"


# Charge les ancetres deja construits et rajoute les infos dans les diagonales
def checkAlreadyBuildAnc():
	genAlready = {}
	for f in filsAnc:
		s = options["alreadyBuiltAnc"] % phylTree.fileName[f]
		print >> sys.stderr, "Checking %s ..." % f,
		if os.access(s, os.R_OK):
			genAlready[f] = utils.myGenomes.Genome(s, withChr=True)
			s = len(genAlready[f].lstChr)
			if options["weightNbChr+"] and (s > 0):
				espCertitude[f] = 1. - 1. / float(s)
			if options["weightNbChr-"]:
				espIncertitude[f] = s * alphaIncertitude
		else:
			print >> sys.stderr, "Not found"

	print >> sys.stderr, "Mise a jour des chromosomes des diagonales ...",
	for (d,esp,esp2) in lstDiags:
		g = utils.myMaths.flatten([lstGenesAnc[i].names for i in d])
		for f in genAlready:
			esp.update([(f,genAlready[f].dicGenes[s][0]) for s in g if s in genAlready[f].dicGenes])
			esp2.update([f for s in g if s in genAlready[f].dicGenes])
	print >> sys.stderr, "OK"




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
		outgroup.extend([(e,1./float(2*phylTree.ages[par]-phylTree.ages[node])) for e in phylTree.branches[par] if e != anc])
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
	[("ancestr",str,""), ("alreadyBuiltAnc",str,""), ("printDiags",bool,False),\
	("useOutgroups",bool,True), ("useLonelyGenes",bool,False), ("weightNbChr+",bool,False), ("weightNbChr-",bool,False), ("newScoring",bool,False), ("walktrapLength",int,5), \
	("genesFile",str,"~/work/data/genes/genes.%s.list.bz2"), \
	("ancGenesFile",str,"~/work/data/ancGenes/ancGenes.%s.list.bz2")], \
	__doc__ \
)

#  Chargement et initialisation
phylTree = utils.myPhylTree.PhylogeneticTree(noms_fichiers["phylTree.conf"])
genesAnc = utils.myGenomes.Genome(options["ancGenesFile"] % phylTree.fileName[options.ancestr])
lstGenesAnc = genesAnc.lstGenes[None]
lstDiags = loadDiagsFile(noms_fichiers["diagsList"], options["ancestr"])
nbDiagsIni = len(lstDiags)

lstToutesEspeces = set(phylTree.listSpecies)
lstEspOutgroup = set(phylTree.outgroupSpecies[options["ancestr"]])
lstEspFilles = set(phylTree.species[options["ancestr"]])
lstNoeudsFils = phylTree.branches[options["ancestr"]]
lstEspParNoeudsFils = [set(s) for s in phylTree.branchesSpecies[options["ancestr"]]]

filsAnc = phylTree.branches[options["ancestr"]]
filsEsp = [set(s) for s in phylTree.branchesSpecies[options["ancestr"]]]
outgroup = set(phylTree.outgroupSpecies[options["ancestr"]])


# L'apport en terme de couverture de l'arbre phylogenetique
def calcPoidsFils(node, calc):
	dicPoidsEspeces[node] = calc
	if node in phylTree.tmpItems:
		poids = calc / float(len(phylTree.tmpItems[node]))
		for (f,_) in phylTree.tmpItems[node]:
			calcPoidsFils(f, poids)
dicPoidsEspeces = dict.fromkeys(phylTree.listSpecies, 0.)
phylTree.initCalcDist(options["ancestr"], options["useOutgroups"])
calcPoidsFils(options["ancestr"], float(len(phylTree.tmpItems[0])))

# Chargement des genomes si necessaire
lengths = []
dicGenes = {}
if options["weightNbChr+"] or options["weightNbChr-"] or options["useLonelyGenes"]:
	for e in phylTree.listSpecies:
		genome = utils.myGenomes.Genome(options["genesFile"] % phylTree.fileName[e])
		lengths.append(len(genome.lstChr))
		if options["useLonelyGenes"]:
			for (g,(c,_)) in genome.dicGenes.iteritems():
				dicGenes[g] = (e,c)

# Le taux de precision de chaque espece
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

# On doit rajouter les genes non presents dans des diagonales
if options["useLonelyGenes"]:
	checkLonelyGenes()
nbDiags = len(lstDiags)

# On doit noter les chromosomes des diagonales sur ces ancetres deja construits
if options["alreadyBuiltAnc"] != "":
	checkAlreadyBuildAnc()

print >> sys.stderr, "%-25s\t%s\t%s\t%s" % ("Ancetre", "Poids", "Cert", "Incert")
for e in dicPoidsEspeces:
	print >> sys.stderr, "%-25s\t%.3f\t%.3f\t%.3f" % (e, dicPoidsEspeces[e], espCertitude[e], espIncertitude[e])


def calcScore(i1, i2):
	
	(_,ec1,e1) = lstDiags[i1]
	(_,ec2,e2) = lstDiags[i2]
	comparedEsp = e1.intersection(e2)
	communEsp = set([e for (e,_) in ec1.intersection(ec2)])

	propF = range(len(filsAnc))
		
	for (i,f) in enumerate(filsAnc):
		# Chez l'ancetre du dessous - meme chr
		if f in communEsp:
			propF[i] = dicPoidsEspeces[f]*espCertitude[f]
		# Chez l'ancetre du dessous - diff chr
		elif f in comparedEsp:
			propF[i] = dicPoidsEspeces[f]*espIncertitude[f]
		# Sinon, on revient aux genomes modernes
		else:
			f = filsEsp[i]
			s = sum([dicPoidsEspeces[e]*espCertitude[e] for e in f.intersection(communEsp)])
			# On peut rajouter l'incertitude des autres especes
			#   - si au moins une espece a valide la fusion
			#   - si la branche est constituee d'une unique espece
			if (s > 0) or (len(f) == 1):
				s += sum([dicPoidsEspeces[e]*espIncertitude[e] for e in f.difference(communEsp)])
			propF[i] = s
	
	s = sum(propF)
	if s == 0:
		return 0

	propOut = sum([dicPoidsEspeces[e]*espCertitude[e] for e in communEsp.intersection(outgroup)])
	# On rajoute l'incertitude
	if (propOut > 0):
		propOut += sum([dicPoidsEspeces[e]*espIncertitude[e] for e in communEsp.difference(outgroup)])
	
	s *= propOut
	for (f1,f2) in utils.myTools.myIterator.tupleOnStrictUpperList(propF):
		s += f1 * f2

	return s
	
def calcScore2(i1, i2):

	(_,ec1,e1) = lstDiags[i1]
	(_,ec2,e2) = lstDiags[i2]

	values = {}
	for e in e1.intersection(e2):
		values[e] = espIncertitude[e]
	for (e,_) in ec1.intersection(ec2):
		values[e] = espCertitude[e]

	prop = [phylTree.calcDist(values, f) for f in filsAnc]
	
	s = sum( [f1 * f2 for (f1,f2) in it([x for x in prop if x != None])] )

	return s

it = utils.myTools.myIterator.tupleOnStrictUpperList

if options["newScoring"]:
	if len(outgroup) != 0:
		filsAnc.append(phylTree.parent[options["ancestr"]])
	f = calcScore2
else:
	f = calcScore

walktrapInstance = utils.walktrap.WalktrapLauncher(showProgress=True, randomWalksLength=options["walktrapLength"])
print >> sys.stderr, "Calcul de la matrice ...",
walktrapInstance.updateFromFunc(range(nbDiags), f)
walktrapInstance.doWalktrap()

clusters = []
# Chaque composante connexe
for (nodes,cuts,_,dend) in walktrapInstance.res:
	print >> sys.stderr, cuts
	# Un score de relevance > 0.2
	interessant = [(alpha,score,dend.cut(alpha)) for (alpha,score) in cuts if score > 0.1]
	# Les noeuds seuls representent < de la moitie de l'ensemble des noeuds
	interessant = [(alpha,score,clust) for (alpha,score,(clust,lonely)) in interessant if len(lonely) < len(nodes)/2]
	if len(interessant) == 0:
		print >> sys.stderr, "-",
		clusters.append(nodes)
	else:
		# Au choix, on prend la version la moins fusionnee
		print >> sys.stderr, "+%d/%f/%f" % (len(interessant[-1][-1]), interessant[-1][0], interessant[-1][1]),
		clusters.extend(interessant[-1][-1])
# -> clusters contient la repartition des diagonales
print >> sys.stderr

print >> sys.stderr, "Impression des %d chromosomes ancestraux ..." % len(clusters),
lstChr = []
for clust in clusters:
	lst = set()
	for i in clust:
		lst.update(lstDiags[i][0])
	lstChr.append(lst)
# -> lstChr contient la repartition des genes

inter = set()
for (l1,l2) in utils.myTools.myIterator.tupleOnStrictUpperList(lstChr):
	inter.update(l1.intersection(l2))
# -> inter contient les genes presents dans deux chromosomes a la fois

for indChr in xrange(len(lstChr)):
	if options["printDiags"]:
		for d in clusters[indChr]:
			if d < nbDiagsIni:
				print indChr+1, " ".join([str(g) for g in lstDiags[d][0] if g not in inter])
			else:
				for g in lstDiags[d][0]:
					if g not in inter:
						print indChr+1, g
	else:
		for g in lstChr[indChr]:
			if g not in inter:
				print indChr+1, " ".join(lstGenesAnc[g].names)

print >> sys.stderr, "OK"

