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
import utils.myPhylTree
import utils.myGenomes
import utils.myTools
import utils.myMaths
import utils.walktrap

defaultdict = utils.myTools.defaultdict


# FONCTIONS #

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
def doSynthese(combin, eND, orthos):
	
	print >> sys.stderr, "Synthese des decoupages de", eND, "...",
	
	# combin a ete mis a jour avec tous les DCS
	# Il faut le lire pour avoir chaque DCS final (qu'il faut reformatter)
	lstBlocs = []
	for gr in combin:
		# On fait la liste des positions et des chromosomes des orthologues
		l = []
		for g in gr:
			p = phylTree.dicGenomes[eND].dicGenes[g]
			a = {}
			for eD in especesDup:
				a[eD] = [c for (c,_) in orthos[eD].get(g,[])]
			l.append( (p,g,a) )
		l.sort()
		lstBlocs.append(l)
	lstBlocs.sort()
	
	print >> sys.stderr, len(lstBlocs), "blocs pour", sum([len(x) for x in lstBlocs]), "orthologues",

	if arguments["showDCS"]:
		print utils.myTools.printLine( [eND, "", "", ""] + especesDup )
	
	res = []
	DCSlen = 0

	# On extrait le profil d'alternance de chaque DCS et on ne garde que les DCS qui alternent
	for gr in lstBlocs:
		# Alternance detectee
		(ok,alt) = addDCS(gr)
		if ok:
			res.append( (gr,alt) )
			DCSlen += len(gr)
		if arguments["showDCS"]:
			for ((c,i),g,a) in gr:
				print utils.myTools.printLine( [c, i, g, ""] + [utils.myTools.printLine(frozenset(a[eD]), "/") for eD in especesDup] )
			print "---"

	print >> sys.stderr, "/", len(res), "DCS pour", DCSlen, "orthologues"
	return res


#
# Construit les tables d'association "nom de gene" -> couleur
# Prend la couleur la plus probable
#
def addDCS(dcs):

	#
	# Compte le score d'alternance de chaque paire de chromosomes de l'espece eD
	#
	def countAltern(lstDCS, eD):
		
		# La liste des chromosomes de l'alternance
		lst = [x[eD] for (_,_,x) in lstDCS.__reversed__()]

		# Compte le nombre d'occurrences de c dans la liste courante
		def countChr(c):
			nb = 0
			for x in lst:
				if c not in x:
					break
				nb += 1
			return nb
		
		# Le compte final
		count = defaultdict(int)
		# La derniere ligne lue
		last = defaultdict(int)
		# On parcourt la liste
		while len(lst) > 0:
			curr = lst.pop()
			for x in curr:
				# Les alternances sont mesurees entre deux positions consecutives
				for y in last:
					if y == x:
						continue
					count[(x,y)] += (countChr(x)+1) * last[y]
					count[(y,x)] = count[(x,y)]
				# Et aussi entre les paralogues
				for y in curr:
					if y >= x:
						continue
					count[(x,y)] += 1
					count[(y,x)] = count[(x,y)]
			
			# On met a jour last
			for y in last:
				if y not in curr:
					last[y] = 0
			for x in curr:
				last[x] += 1

		return count

	# On parcourt les DCS en ne gardant que ceux qui alternent
	res = {}
	altern = False
	for eD in especesDup:
		res[eD] = countAltern(dcs, eD)
		if len(res[eD]) > 0 and max(res[eD].values()) > 0:
			altern = True
	
	return (altern,res)


#
# Range chaque gene ancestral dans son chromosome
#
def buildChrAnc(genesAncCol, chrAncGenes):

	#
	# Renvoie un score (~pourcentage d'especes) qui soutiennent l'attribution d'un gene a un chromosome
	#
	def calcChrAncScore(col, ch):
		
		# Une tetrapode rapporte 1 si un de ses representants a vote pour le chromosome ancestral
		values = utils.myTools.hashabledict()
		for (_,c,eND) in col:
			values[eND] = max(values.get(eND,0), float(c == ch))

		# Soit on fait un calcul phylogenetique
		if arguments["usePhylTreeScoringNonDup"]:
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

	for c in chrNames:
		nb = 0
		for i in chrAncGenes[c]:
			nb += 1
			if arguments["showQuality"]:
				print utils.myTools.printLine( [c, nb, utils.myTools.printLine(["%.2f" % (100*x) for (x,_) in col[i]]), 100*col[i][c-1][0]] )
			if arguments["showAncestralGenome"]:
				print c, " ".join(genesAnc[i].names)
		
	print >> sys.stderr, sum([len(chrAncGenes[c]) for c in chrNames]), "genes dans le genome ancestral"


########
# MAIN #
########


# Arguments
arguments = utils.myTools.checkArgs( \
	[("phylTree.conf",file)],
	[("minChrLen",int,20), ("windowSize",int,25), ("usePhylTreeScoringDup",bool,False), ("usePhylTreeScoringNonDup",bool,False), ("keepUncertainGenes",bool,False), \
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
rootDup = especesDup[0]
for e in especesDup[1:]:
	rootDup = phylTree.dicParents[rootDup][e]
dicCommonEspNames = dict([(phylTree.officialName[esp],esp) for esp in especesDup+especesNonDup])
phylTree.loadSpeciesFromList(especesNonDup+especesDup, arguments["genesFile"])
genesAnc = utils.myGenomes.Genome(arguments["ancGenesFile"] % phylTree.fileName[arguments["target"]])
lstGenesAnc = genesAnc.lstGenes[None]
(para,orthos) = buildParaOrtho()

# De quoi stocker les DCS
col = [[] for i in xrange(len(lstGenesAnc))]
allDCS = []

# Decoupage de chaque tetrapode
for eND in especesNonDup:
	combin = utils.myTools.myCombinator([])
	# Avec chaque poisson
	for eD in especesDup:
		lstBlocs = colorAncestr(eND, eD, phylTree, para, orthos)
		for bloc in lstBlocs:
			combin.addLink(bloc)
	
	allDCS.extend( doSynthese(combin, eND, orthos))
		
allDCSe2 = [dict([(e,frozenset(alt[e])) for e in alt]) for (_,alt) in allDCS]

print >> sys.stderr, "Comparaison des %d DCS ..." % len(allDCS),

if arguments["usePhylTreeScoringDup"]:
	dicCoeff = {}
	for e in especesDup:
		values = dict.fromkeys(especesDup, 0)
		values[e] = 1
		dicCoeff[e] = phylTree.calcWeightedValue(values, 0, rootDup)[2]
else:
	dicCoeff = dict.fromkeys(especesDup, 1./len(especesDup))


edges = utils.myTools.defaultdict(dict)
for i1 in xrange(len(allDCS)):
	(_,alt1) = allDCS[i1]
	esp1 = allDCSe2[i1]
	
	for i2 in xrange(i1):
		(_,alt2) = allDCS[i2]

		score = 0
		for e in especesDup:
			val = 0
			for x in (esp1[e] & allDCSe2[i2][e]):
				val += min(alt1[e][x], alt2[e][x])
			score += val*dicCoeff[e]
		
		if score > 0:
			edges[i1][i2] = edges[i2][i1] = score

print >> sys.stderr, "Lancement de walktrap ...",
res = utils.walktrap.doWalktrap(edges)

# A partir d'ici, on a une association DCS <-> chromosomes, on retombe sur la regle de base, le vote a la majorite

chrInd = 0
for (nodes,cuts,_,dend) in res:
	print >> sys.stderr, "Communaute de %d noeuds:" % len(nodes)
	# Un point de coupure virtuel si il n'y en a pas
	cuts.append( (1,0) )
	(alpha,relevance,(clusters,lonely)) = utils.walktrap.askPartitionChoice(dend, cuts)
	# On enregistre les resultats
	for cl in clusters:
		chrInd += 1
		for g in cl:
			for (_,s,_) in allDCS[g][0]:
				(_,anc) = genesAnc.dicGenes[s]
				(e,_,_) = phylTree.dicGenes[s]
				col[anc].append( (None, chrInd, e) )

print >> sys.stderr, "%d chromosomes ancestraux" % chrInd

# On construit les chromosomes ancestraux
print >> sys.stderr, "Synthese des genes ancestraux ...",
chrAncGenes = dict([ (x,[]) for x in xrange(1, chrInd) ])
buildChrAnc(col, chrAncGenes)
print >> sys.stderr, "OK"
	
# On affiche le resultat
printColorAncestr(lstGenesAnc, chrAncGenes)

