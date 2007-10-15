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
		# On joint les especes qui ont vu la diagonale et celles qui n'apportent que le chromosome
		tmp = [y.split("/") for y in "|".join([x for x in ct[4:] if len(x) > 0]).split("|")]
		# Les chromosomes de ces especes
		espChr = frozenset( (phylTree.officialName[e],c) for (e,c) in tmp if ('Un' not in c) )
		# On la garde en memoire
		lst.append( (d,espChr,ct[2],ct[3]) )

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
	for (d,_,_,_) in lstDiagsIni:
		genesSeuls.difference_update(d)
	
	nb = 0
	new = []
	for i in genesSeuls:

		lst = []
		for e in phylTree.listSpecies:
			if e in lstEspOutgroup:
				a = phylTree.dicParents[options["ancestr"]][e]
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
		# Pour etre exact, il faudrait avec 2 groupes parmi (fils1 + ... + filsN + outgroup)
		if len(lst) >= 2:
			new.append( ([i],frozenset(lst),str(i),"1") )
			nb += 1
	
	print >> sys.stderr, "%d genes OK" % len(new)
	return new


# Charge les ancetres deja construits et rajoute les infos dans les diagonales
def checkAlreadyBuildAnc():
	genAlready = {}
	for f in lstNoeudsFils:
		s = options["alreadyBuiltAnc"] % phylTree.fileName[f]
		print >> sys.stderr, "Checking %s ..." % f,
		if utils.myTools.fileAccess(s):
			genAlready[f] = utils.myGenomes.Genome(s, withChr=True)
			s = len(genAlready[f].lstChr)
			if options["weightNbChr+"] and (s > 0):
				espCertitude[f] = 1. - 1. / float(s)
			if options["weightNbChr-"]:
				espIncertitude[f] = s * alphaIncertitude
		else:
			print >> sys.stderr, "Not found"

	print >> sys.stderr, "Mise a jour des chromosomes des diagonales ...",
	for (d,esp,_,_) in lstDiagsIni:
		g = utils.myMaths.flatten([lstGenesAnc[i].names for i in d])
		for f in genAlready:
			esp.update([(f,genAlready[f].dicGenes[s][0]) for s in g if s in genAlready[f].dicGenes])
	print >> sys.stderr, "OK"



########
# MAIN #
########

# Arguments
############
(noms_fichiers, options) = utils.myTools.checkArgs( \
	["phylTree.conf", "diagsList"], \
	[("ancestr",str,""), ("alreadyBuiltAnc",str,""), ("useOutgroups",bool,True), \
	("mergeDiags",bool,False), ("useLonelyGenes",bool,False), ("weightNbChr+",bool,False), ("weightNbChr-",bool,False), ("newScoring",bool,False), ("walktrapLength",int,5), \
	("genesFile",str,"~/work/data/genes/genes.%s.list.bz2"), \
	("ancGenesFile",str,"~/work/data/ancGenes/ancGenes.%s.list.bz2")], \
	__doc__ \
)


#  Chargement et initialisation
################################
phylTree = utils.myPhylTree.PhylogeneticTree(noms_fichiers["phylTree.conf"])
genesAnc = utils.myGenomes.Genome(options["ancGenesFile"] % phylTree.fileName[options.ancestr])
lstGenesAnc = genesAnc.lstGenes[None]

lstNoeudsFils = phylTree.branches[options["ancestr"]]
lstEspParNoeudsFils = [set(s) for s in phylTree.branchesSpecies[options["ancestr"]]]
lstEspOutgroup = set(phylTree.outgroupSpecies[options["ancestr"]])


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
	if node in phylTree.tmpItems:
		poids = calc / float(len(phylTree.tmpItems[node]))
		for (f,_) in phylTree.tmpItems[node]:
			calcPoidsFils(f, poids)
dicPoidsEspeces = dict.fromkeys(phylTree.listSpecies, 0.)
phylTree.initCalcDist(options["ancestr"], options["useOutgroups"])
calcPoidsFils(options["ancestr"], float(len(phylTree.tmpItems[0])))

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
lstDiagsIni = loadDiagsFile(noms_fichiers["diagsList"], options["ancestr"])

# On doit rajouter les genes non presents dans des diagonales
if options["useLonelyGenes"]:
	lstDiagsIni.extend(checkLonelyGenes())

# On doit noter les chromosomes des diagonales sur ces ancetres deja construits
if options["alreadyBuiltAnc"] != "":
	checkAlreadyBuildAnc()

# On reduit les diagonales
print >> sys.stderr, "Reduction de %d elements a ..." % len(lstDiagsIni),
red = utils.myTools.defaultdict(list)
for (i,(_,e,_,_)) in enumerate(lstDiagsIni):
	red[e].append(i)
lstDiags = [ (d,e,frozenset(x for (x,_) in e)) for (e,d) in red.iteritems() ]
nbDiags = len(lstDiags)
dicGenes.clear()
print >> sys.stderr, nbDiags

print >> sys.stderr, "%-25s\t%s\t%s\t%s" % ("Ancetre", "Poids", "Cert", "Incert")
for e in dicPoidsEspeces:
	print >> sys.stderr, "%-25s\t%.3f\t%.3f\t%.3f" % (e, dicPoidsEspeces[e], espCertitude[e], espIncertitude[e])
	if not options["newScoring"]:
		espCertitude[e] *= dicPoidsEspeces[e]
		espIncertitude[e] *= dicPoidsEspeces[e]

# On calcule les scores
print >> sys.stderr, "Calcul de la matrice ...",

#scores = numpy.zeros( (len(lstDiags),len(lstDiags)) )
if options["newScoring"] and (len(lstEspOutgroup) != 0):
	lstNoeudsFils.append(phylTree.parent[options["ancestr"]])
it = utils.myTools.myIterator.tupleOnStrictUpperList

# On les stocke directement dans l'instance de walktrap
walktrapInstance = utils.walktrap.WalktrapLauncher(showProgress=True, randomWalksLength=options["walktrapLength"])
edges = walktrapInstance.edges

for i1 in xrange(nbDiags):
	(d1,ec1,e1) = lstDiags[i1]
	for i2 in xrange(i1):
		(d2,ec2,e2) = lstDiags[i2]
		
		if options["newScoring"]:
			values = {}
			for e in e1.intersection(e2):
				values[e] = espIncertitude[e]
			for (e,_) in ec1.intersection(ec2):
				values[e] = espCertitude[e]
			prop = [phylTree.calcDist(values, f) for f in lstNoeudsFils]
			s = sum( [f1 * f2 for (f1,f2) in it([x for x in prop if x != None])] )
		else:
			comparedEsp = e1.intersection(e2)
			communEsp = set([e for (e,_) in ec1.intersection(ec2)])

			propF = range(len(lstNoeudsFils))
				
			for (i,f) in enumerate(lstNoeudsFils):
				# Chez l'ancetre du dessous - meme chr
				if f in communEsp:
					propF[i] = espCertitude[f]
				# Chez l'ancetre du dessous - diff chr
				elif f in comparedEsp:
					propF[i] = espIncertitude[f]
				# Sinon, on revient aux genomes modernes
				else:
					f = lstEspParNoeudsFils[i]
					s = sum([espCertitude[e] for e in f.intersection(communEsp)])
					# On peut rajouter l'incertitude des autres especes
					#   - si au moins une espece a valide la fusion
					#   - si la branche est constituee d'une unique espece
					if (s > 0) or (len(f) == 1):
						s += sum([espIncertitude[e] for e in f.difference(communEsp)])
					propF[i] = s
			
			s = sum(propF)
			if s != 0:
				propOut = sum([espCertitude[e] for e in communEsp.intersection(lstEspOutgroup)])
				# On rajoute l'incertitude
				if (propOut > 0):
					propOut += sum([espIncertitude[e] for e in communEsp.difference(lstEspOutgroup)])
				
				s *= propOut
				for (f1,f2) in it(propF):
					s += f1 * f2
		
		if s > 0:
			if options["mergeDiags"]:
				edges[i1][i2] = s
				edges[i2][i1] = s
			else:
				for x1 in d1:
					for x2 in d2:
						edges[x1][x2] = s
						edges[x2][x1] = s

print >> sys.stderr, "OK"

walktrapInstance.doWalktrap()

clusters = []
# Chaque composante connexe
for (nodes,cuts,_,dend) in walktrapInstance.res:
	print >> sys.stderr, cuts
	# Un score de relevance > 0.1
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
# Chaque cluster
for clust in clusters:
	lstD = []
	lstG = set()
	for t in clust:
		# Chaque element du cluster est le numero d'un ensemble de diagonales initiales
		if options["mergeDiags"]:
			t = lstDiags[t][0]
		else:
			t = [t]
		# lstD contient la liste des diagonales du chromosome
		lstD.extend(t)
		for i in t:
			# lstG contient la liste des genes de ces diagonales
			lstG.update(lstDiagsIni[i][0])
	lstChr.append( (lstD,lstG) )
# -> lstChr contient la repartition des genes

for indChr in xrange(len(lstChr)):
	for d in lstChr[indChr][0]:
		(_,_,d,s) = lstDiagsIni[d]
		print "%d\t%s\t%s" % (indChr+1,d,s)

print >> sys.stderr, "OK"

