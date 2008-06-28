#! /users/ldog/muffato/python -OO

__doc__ = """
Ce script scanne chaque genome d'espece non dupliquee en le comparant a chaque genome duplique et cree les DCS.
Un gene est inclus dans le DCS courant si son orthologue est proche d'un gene paralogue d'un gene proche des genes inclus dans le DCS.
On regroupe les genes en faisant des suites qui ont une origine commune.
Chaque gene ancestral recoit son chromosome ancestral par un vote a la majorite en fonction l'annotation de chaque tetrapode
"""

# INITIALISATION #

# Librairies
import sys
import operator
import itertools
import utils.myPhylTree
import utils.myGenomes
import utils.myTools
import utils.myMaths
defaultdict = utils.myTools.defaultdict


# FONCTIONS #


#
# Charge le fichier qui contient les alternances que l'on doit observer pour les especes utilisees
# Format de chaque ligne: 
# A	Tetraodon nigroviridis|2/3 5	Gasterosteus aculeatus|16/1	Oryzias latipes|21/2	Takifugu rubripes|/	Danio rerio|9/6 1
# Renvoie un dictionnaire des chromosomes ancestraux, qui contient les dictionnaires (par especes) des alternances
#
def loadChrAncIni(nom, especesDup):

	print >> sys.stderr, "Chargement des alternances predefinies ...",
	
	chrAnc = {}
	f = utils.myTools.myOpenFile(nom, 'r')
	for ligne in f:

		c = ligne.replace('\n', '').split('\t')
		dic = {}
		for x in c[1:]:
			(e,x) = x.split('|')
			(c1,c2) = x.split('/')
			dic[phylTree.officialName[e]] = (set([utils.myGenomes.commonChrName(x) for x in c1.split()]), set([utils.myGenomes.commonChrName(x) for x in c2.split()]))
		chrAnc[c[0]] = [dic.get(phylTree.officialName[e], (set(),set())) for e in especesDup]
	f.close()
	print >> sys.stderr, "OK"

	return chrAnc


#
# Construit les tables d'associations pour chaque especes dupliquee
#  des paralogues et des orthologues
#
def buildParaOrtho():

	print >> sys.stderr, "Formattage des listes d'orthologues et de paralogues ...",

	para = dict([(e,{}) for e in especesDup])
	ortho = dict([(e,{}) for e in especesDup])
	
	# On parcourt les genes ancestraux pour les especes en dessous de l'ancetre
	for g in lstGenesAnc:
		for e in especesDup:
			genomeDup = phylTree.dicGenomes[e]
			# Les genes de l'espece dupliquee
			gT = [x for x in g.names if x in genomeDup.dicGenes and genomeDup.dicGenes[x][0] not in genomeDup.lstRand]
			if len(gT) == 0:
				continue
			# On construit les paralogues
			for x in gT:
				para[e][x] = [y for y in gT if y != x]
			# On construit les orthologues
			gT = [genomeDup.dicGenes[y] for y in gT]
			for x in g.names:
				if x not in genomeDup.dicGenes:
					ortho[e][x] = gT
	
	print >> sys.stderr, "OK"

	# On rajoute les outgroup
	
	# Les ancetres correspondant a chaque outgroup
	toLoad = defaultdict(list)
	for e in especesNonDup:
		par = phylTree.dicParents[e][arguments["target"]]
		toLoad[par].append(e)
	toLoad.pop(arguments["target"], None)

	for (anc,outgroups) in toLoad.iteritems():
		print >> sys.stderr, "Rajout de %s ..." % "/".join(outgroups),
		outgroups =  frozenset([phylTree.officialName[e] for e in outgroups])
		for g in utils.myGenomes.Genome(arguments["ancGenesFile"] % phylTree.fileName[anc]):
			# On trie les genes ancestraux
			dicGenes = utils.myTools.defaultdict(list)
			newPos = None
			for x in g.names:
				# La position dans le genome que l'on reconstruit
				newPos = genesAnc.dicGenes.get(x, newPos)
				if x not in phylTree.dicGenes:
					continue
				(esp,_,_) = phylTree.dicGenes[x]
				if esp in dicCommonEspNames:
					dicGenes[dicCommonEspNames[esp]].append(x)
			if newPos == None:
				continue
			gNT = []
			for e in especesNonDup:
				if e in dicGenes:
					gNT.extend(dicGenes[e])

			# Mise a jour du dictionnaire pour qu'on puisse les utiliser
			for x in gNT:
				genesAnc.dicGenes[x] = newPos
			
			# Enregistrement des orthologues
			for e in especesDup:
				if e in dicGenes:
					gT = [phylTree.dicGenomes[e].dicGenes[x] for x in dicGenes[e]]
					for x in gNT:
						ortho[e][x] = gT
	return (para, ortho)


#
# Cree les DCS en parcourant un genome non duplique (eND) face a un genome duplique (eD)
#
def colorAncestr(eND, eD, phylTree, para, orthos):

	print >> sys.stderr, "Decoupage de", eND, "avec", eD, "",

	nbOrthos = 0
	genome = phylTree.dicGenomes[eND]
	genomeDup = phylTree.dicGenomes[eD]
	orthosDup = orthos[eD]
	parasDup = para[eD]
	lstBlocs = []

	# On parcourt les chromosomes de l'espece
	for c in genome.lstChr + genome.lstScaff:
		
		# Eviter d'essayer de creer des DCS sur des scaffolds trop petits
		if len(genome.lstGenes[c]) < arguments["minChrLen"]:
			continue
		
		bloc = None
		lastCT = []
		lastGT = set()
	
		# On parcourt les genes du chromosomes
		for tg in genome.lstGenes[c]:
			g = tg.names[0]

			# Il faut un orthologue avec l'espece dupliquee
			if g not in orthosDup:
				continue
			nbOrthos += 1

			# On parcourt les orthologues
			for (cT,i) in orthosDup[g]:
				# La region environnante (chez l'espece dupliquee)
				gTn = [gT.names[0] for gT in genomeDup.getGenesNear(cT, i, arguments["windowSize"])]
				
				# Si on reste sur le meme chromosome, on continue le DCS
				if cT in lastCT:
					break
				# Autre solution, on revient dans le voisinage d'une region deja visitee
				elif len(lastGT.intersection(gTn)) > 0:
					break
				# Sinon, il faut qu'il y ait un paralogue qui justifie le saut de chromosome
				elif len(lastGT.intersection(utils.myMaths.flatten([parasDup.get(s,[]) for s in gTn]))) > 0:
					break
				
			# Si on n'a pas fait de break, c'est qu'aucun orthologue ne convient, il faut arreter le DCS
			else:
				# On l'enregistre et on repart de zero
				if bloc != None:
					lstBlocs.append(bloc)
				bloc = []
				lastGT = set()
				
			# On rajoute les infos du gene qu'on vient de lire
			bloc.append( g )
			lastCT = [cT for (cT,i) in orthosDup[g]]
			lastGT.update(gTn)
		
		# Ne pas oublier le bloc courant
		if bloc != None:
			lstBlocs.append(bloc)
		sys.stderr.write(".")
	
	print >> sys.stderr, "", len(lstBlocs), "blocs, pour", nbOrthos, "genes orthologues"

	return lstBlocs



#
# Pour une espece non dupliquee donnee, fait la synthese de tous les DCS chevauchants
#
def doSynthese(combin, eND, orthos, col, dicGenesAnc, chrAnc):
	
	print >> sys.stderr, "Synthese des decoupages de", eND, "...",
	
	# combin a ete mis a jour avec tous les DCS
	# Il faut le lire pour avoir chaque DCS final (qu'il faut reformatter)
	lstBlocs = []
	for gr in combin:
		# On fait la liste des positions et des chromosomes des orthologues
		l = [ (phylTree.dicGenomes[eND].dicGenes[g], g, [[c for (c,_) in orthos[eD].get(g,[])] for eD in especesDup]) for g in gr]
		l.sort()
		lstBlocs.append(l)
	lstBlocs.sort()
	
	print >> sys.stderr, len(lstBlocs), "blocs pour", sum([len(x) for x in lstBlocs]), "orthologues",

	if arguments["showDCS"]:
		print utils.myTools.printLine( [eND, "", "", ""] + especesDup )
	
	nbDCS = 0
	DCSlen = 0
	eND = phylTree.officialName[eND]

	# On va assigner un chromosome ancestral a chaque DCS en fonction des alternances predefinies
	for gr in lstBlocs:
		cc = addDCS(gr, col, dicGenesAnc, chrAnc, eND)
		# Alternance detectee
		if cc != None:
			nbDCS += 1
			DCSlen += len(gr)
		if arguments["showDCS"]:
			for ((c,i),g,a) in gr:
				#fishContent = ["/".join([str(x) for x in set(y)]) for y in a]
				fishContent = ["/".join(["%s|%s" % (phylTree.dicGenomes[eD].lstGenes[cT][iT].names[0],cT) for (cT,iT) in orthos[eD].get(g,[])]) for eD in especesDup]
				print utils.myTools.printLine( [c, i, g, ""] + fishContent + [cc] )
			print "---"

	print >> sys.stderr, "/", nbDCS, "DCS pour", DCSlen, "orthologues"


#
# Construit les tables d'association "nom de gene" -> couleur
# Prend la couleur la plus probable
#
def addDCS(bloc, col, dicGenesAnc, chrAnc, eNonDup):

	score = {}
	# On calcule le score de chaque chromosome ancestral
	for c in chrAnc:
		(nb1,nb2) = (0,0)
		for (_,_,a) in bloc:
			(flag1,flag2) = (False,False)
			# On parcourt les especes
			for ((expected1,expected2),observed) in itertools.izip(chrAnc[c],a):
				flag1 |= len(expected1.intersection(observed)) > 0
				flag2 |= len(expected2.intersection(observed)) > 0
			# On met a jour nb1 et nb2: le nombre de genes vus a gauche / a droite
			nb1 += flag1
			nb2 += flag2
		
		# Est-ce qu'on observe bien une alternance sur ce chromosome
		#   <-> Est-ce qu'on a vu au moins un gene a gauche et un gene a droite
		if (nb1 == 0) or (nb2 == 0):
			# Non -> score = 0
			score[c] = 0
		else:
			# Oui -> score = nb de genes qui appartiennent au chromosome
			score[c] = nb1 + nb2

	(c,s) = max(score.items(), key = operator.itemgetter(1))
	
	if s == 0:
		return None
	
	for g in bloc:
		col[dicGenesAnc[g[1]][1]].append( (len(bloc),c,eNonDup) )
	
	return c

#
# Range chaque gene ancestral dans son chromosome
#
def buildChrAnc(genesAncCol, chrAncGenes):

	# Renvoie un score (~pourcentage d'especes) qui soutiennent l'attribution d'un gene a un chromosome
	def calcChrAncScore(col, ch):
		
		# Une tetrapode rapporte 1 si un de ses representants a vote pour le chromosome ancestral
		values = utils.myTools.hashabledict()
		for (_,c,eND) in col:
			values[eND] = max(values.get(eND,0), float(c == ch))

		# Soit on fait un calcul phylogenetique
		if arguments["usePhylTreeScoring"]:
			return probaMemory(values, 0, rootNonDup)[2]

		# On fait les groupes
		rTot = [[values[e] for e in gr if e in values] for gr in especesNonDupGrp]
		# La moyenne des moyennes de chaque groupe non vide
		return utils.myMaths.mean([utils.myMaths.mean(r) for r in rTot if len(r) > 0])

	probaMemory = utils.myTools.memoize(phylTree.calcWeightedValue)
	chrNames = sorted(chrAncGenes)
	for (i,col) in enumerate(genesAncCol):
	
		if len(col) == 0:
			# Certains genes n'ont pas de chance !
			continue
		nb = [(calcChrAncScore(col,x), x) for x in chrNames]
		
		# On verifie les egalites
		s = sorted(nb)
		if (s[-1][0] == s[-2][0]) and not arguments["keepUncertainGenes"]:
			continue

		genesAncCol[i] = nb
		chrAncGenes[s[-1][1]].append(i)


#
# Affichage des chromosomes ancestraux
#
def printColorAncestr(genesAnc, chrAncGenes):
	
	print >> sys.stderr, "Impression des associations genes / chromosomes ancestraux ... ",

	chrNames = sorted(chrAncGenes)
	
	if arguments["showQuality"]:
		print utils.myTools.printLine( [""] + chrNames)

	for (j,c) in enumerate(chrNames):
		nb = 0
		for i in chrAncGenes[c]:
			nb += 1
			if arguments["showQuality"]:
				print utils.myTools.printLine( [c, nb, utils.myTools.printLine(["%.2f" % (100*x) for (x,_) in col[i]]), 100*col[i][j][0]] )
			if arguments["showAncestralGenome"]:
				print c, " ".join(genesAnc[i].names)
		
	print >> sys.stderr, sum([len(chrAncGenes[c]) for c in chrNames]), "genes dans le genome ancestral"


########
# MAIN #
########


# Arguments
arguments = utils.myTools.checkArgs( \
	[("draftPreDupGenome.conf",file), ("phylTree.conf",file)],
	[("minChrLen",int,20), ("windowSize",int,25), ("usePhylTreeScoring",bool,True), ("keepUncertainGenes",bool,False), \
	("especesNonDup",str,""), ("especesDup",str,""), ("target",str,""), \
	("genesFile",str,"~/work/data/genes/genes.%s.list.bz2"), \
	("ancGenesFile",str,"~/work/data/ancGenes/ancGenes.%s.list.bz2"), \
	("showDCS",bool,False), ("showQuality",bool,False), ("showAncestralGenome",bool,True)], \
	__doc__ \
)

# Chargement des fichiers
phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])

def proceedList(l):
	l1 = [s[1:] for s in l if s[0] == "."]
	l2 = [s for s in l if s[0] != "."]
	return utils.myMaths.flatten([phylTree.species[s] for s in l2]) + l1
especesDup = proceedList(arguments["especesDup"].split(','))
especesNonDupGrp = [proceedList(x.split('+')) for x in arguments["especesNonDup"].split(',')]
especesNonDup = utils.myMaths.flatten(especesNonDupGrp)
rootNonDup = especesNonDup[0]
for e in especesNonDup[1:]:
	rootNonDup = phylTree.dicParents[rootNonDup][e]
dicCommonEspNames = dict([(phylTree.officialName[esp],esp) for esp in especesDup+especesNonDup])
phylTree.loadSpeciesFromList(especesNonDup+especesDup, arguments["genesFile"])
genesAnc = utils.myGenomes.Genome(arguments["ancGenesFile"] % phylTree.fileName[arguments["target"]])
lstGenesAnc = genesAnc.lstGenes[None]
(para,orthos) = buildParaOrtho()
chrAnc = loadChrAncIni(arguments["draftPreDupGenome.conf"], especesDup)

col = [[] for i in xrange(len(lstGenesAnc))]

# Decoupage de chaque tetrapode
for eND in especesNonDup:
	combin = utils.myTools.myCombinator([])
	# Avec chaque poisson
	for eD in especesDup:
		lstBlocs = colorAncestr(eND, eD, phylTree, para, orthos)
		for bloc in lstBlocs:
			combin.addLink(bloc)
	
	doSynthese(combin, eND, orthos, col, genesAnc.dicGenes, chrAnc)


# On construit les chromosomes ancestraux
print >> sys.stderr, "Synthese des genes ancestraux ...",
chrAncGenes = dict([ (x,[]) for x in chrAnc])
buildChrAnc(col, chrAncGenes)
print >> sys.stderr, "OK"
	
# On affiche le resultat
printColorAncestr(lstGenesAnc, chrAncGenes)

