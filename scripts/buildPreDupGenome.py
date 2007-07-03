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



# FONCTIONS #


#
# Charge le fichier qui contient les alternances que l'on doit observer pour les especes utilisees
# Format de chaque ligne: 
# A	Tetraodon nigroviridis|2/3 5	Gasterosteus aculeatus|16/1	Oryzias latipes|21/2	Takifugu rubripes|/	Danio rerio|9/6 1
# Renvoie un dictionnaire des chromosomes ancestraux, qui contient les dictionnaires (par especes) des alternances
#
def loadChrAncIni(nom, especesDup):

	# On convertit en entier le nom du chromosome si possible
	def getEasierName(x):
		try:
			return int(x)
		except Exception:
			return x

	print >> sys.stderr, "Chargement des alternances predefinies ...",
	
	chrAnc = {}
	f = utils.myTools.myOpenFile(nom, 'r')
	for ligne in f:

		c = ligne[:-1].split('\t')
		dic = {}
		for x in c[1:]:
			(e,x) = x.split('|')
			(c1,c2) = x.split('/')
			dic[phylTree.officialName[e]] = (set([getEasierName(x) for x in c1.split()]), set([getEasierName(x) for x in c2.split()]))
		chrAnc[c[0]] = [dic.get(phylTree.officialName[e], (set(),set())) for e in especesDup]
	f.close()
	print >> sys.stderr, "OK"

	return chrAnc


#
# Construit les tables d'associations pour chaque especes dupliquee
#  des paralogues et des orthologues
#
def buildParaOrtho(lstGenesAnc):

	print >> sys.stderr, "Formattage des listes d'orthologues et de paralogues ...",
	para = dict([(e,{}) for e in especesDup])
	ortho = dict([(e,{}) for e in especesDup])
	
	# On parcourt les genes ancestraux
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
		if len(genome.lstGenes[c]) < options["minChrLen"]:
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
				gTn = [gT.names[0] for gT in genomeDup.getGenesNear(cT, i, options["windowSize"])]
				
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

	if options["showDCS"]:
		print "%s\t\t\t\t%s\t" % (eND, "\t".join(especesDup))
	
	nbDCS = 0
	DCSlen = 0

	# On va assigner un chromosome ancestral a chaque DCS en fonction des alternances predefinies
	for gr in lstBlocs:
		cc = addDCS(gr, col, dicGenesAnc, chrAnc, eND)
		# Alternance detectee
		if cc != None:
			nbDCS += 1
			DCSlen += len(gr)
		if options["showDCS"]:
			for ((c,i),g,a) in gr:
				print "%s\t%d\t%s\t\t%s\t%s" % \
				(c, i, g, "\t".join(["/".join([str(x) for x in set(y)]) for y in a]), cc)
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

	#
	# Renvoie un score (~pourcentage d'especes) qui soutiennent l'attribution d'un gene a un chromosome
	#
	def calcChrAncScore(col, ch):
		
		# phylTree.calcDist requiert les noms latins
		if options["usePhylTreeScoring"]:
			values = phylTree.newCommonNamesMapperInstance()
		else:
			values = {}
		
		# Une tetrapode rapporte 1 si un de ses representants a vote pour le chromosome ancestral
		for (_,c,eND) in col:
			values[eND] = max(values.get(eND,0), float(c == ch))

		# Soit on fait un calcul phylogenetique
		if options["usePhylTreeScoring"]:
			return phylTree.calcDist(values)

		# On fait les groupes
		rTot = [[values[e] for e in gr if e in values] for gr in especesNonDupGrp]
		# La moyenne des moyennes de chaque groupe non vide
		return utils.myMaths.mean([utils.myMaths.mean(r) for r in rTot if len(r) > 0])

	
	if options["usePhylTreeScoring"]:
		phylTree.initCalcDist(rootNonDup, False)

	chrNames = sorted(chrAncGenes)
	for (i,col) in enumerate(genesAncCol):
	
		if len(col) == 0:
			# Certains genes n'ont pas de chance !
			continue
	
		nb = [(calcChrAncScore(col,x), x) for x in chrNames]
		
		# On verifie les egalites
		s = sorted(nb, reverse=True)
		if (s[0][0] == s[1][0]) and not options["keepUncertainGenes"]:
			continue

		genesAncCol[i] = nb
		chrAncGenes[s[0][1]].append(i)


	#for i in xrange(len(genesAncCol)):
	#
	#	if len(genesAncCol[i]) == 0:
	#		# Certains genes n'ont pas de chance !
	#		continue
	#
	#	nb = [(calcChrAncScore(genesAncCol[i],x), x) for x in chrNames]
	#	
	#	# On verifie les egalites
	#	s = sorted(nb, reverse=True)
	#	if (s[0][0] == s[1][0]) and not options["keepUncertainGenes"]:
	#		continue
	#
	#	genesAncCol[i] = nb
	#	chrAncGenes[s[0][1]].append(i)


#
# Affichage des chromosomes ancestraux
#
def printColorAncestr(genesAnc, chrAncGenes):
	
	print >> sys.stderr, "Impression des associations genes / chromosomes ancestraux ... ",

	chrNames = sorted(chrAncGenes)
	
	if options["showQuality"]:
		print "\t\t%s" % "\t".join([str(c) for c in chrNames])

	#for j in xrange(len(chrNames)):
	#	c = chrNames[j]
	for (j,c) in enumerate(chrNames):
		nb = 0
		for i in chrAncGenes[c]:
			nb += 1
			if options["showQuality"]:
				print "%s\t%d\t%s\t%.2f" % (c, nb, "\t".join(["%.2f" % (100*x) for (x,_) in col[i]]), 100*col[i][j][0])
			if options["showAncestralGenome"]:
				print c, " ".join(genesAnc[i].names)
		
	print >> sys.stderr, sum([len(chrAncGenes[c]) for c in chrNames]), "genes dans le genome ancestral"


########
# MAIN #
########


# Arguments
(noms_fichiers, options) = utils.myTools.checkArgs( \
	["genesAncestraux.list", "draftPreDupGenome.conf", "phylTree.conf"],
	[("minChrLen",int,20), ("windowSize",int,25), ("usePhylTreeScoring",bool,False), ("keepUncertainGenes",bool,False), \
	("especesNonDup",str,""), ("especesDup",str,""), \
	("genesFile",str,"~/work/data/genes/genes.%s.list.bz2"), \
	("showDCS",bool,False), ("showQuality",bool,False), ("showAncestralGenome",bool,True)], \
	__doc__ \
)

# Chargement des fichiers
phylTree = utils.myPhylTree.PhylogeneticTree(noms_fichiers["phylTree.conf"])
especesDup = options["especesDup"].split(',')
especesNonDupGrp = [x.split('+') for x in options["especesNonDup"].split(',')]
especesNonDup = utils.myMaths.flatten(especesNonDupGrp)
for x in phylTree.branches[phylTree.root]:
	if phylTree.dicParents[x][especesDup[0]] != x:
		rootNonDup = x
phylTree.loadSpeciesFromList(especesNonDup+especesDup, options["genesFile"])
genesAnc = utils.myGenomes.AncestralGenome(noms_fichiers["genesAncestraux.list"])
lstGenesAnc = genesAnc.lstGenes[utils.myGenomes.AncestralGenome.defaultChr]
(para,orthos) = buildParaOrtho(lstGenesAnc)
chrAnc = loadChrAncIni(noms_fichiers["draftPreDupGenome.conf"], especesDup)

# On colorie les matrices actuelles
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

