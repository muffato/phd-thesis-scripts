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
import utils.myPhylTree
import utils.myGenomes
import utils.myTools
import utils.myMaths
import utils.walktrap
from collections import defaultdict

# FONCTIONS #

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
def doSynthese(combin, eND, orthos):
	
	print >> sys.stderr, "Synthese des decoupages de", eND, "...",
	
	# combin a ete mis a jour avec tous les DCS
	# Il faut le lire pour avoir chaque DCS final (qu'il faut reformatter)
	lstBlocs = []
	for gr in combin:
		l = []
		# On fait la liste des positions et des chromosomes des orthologues
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

	if options["showDCS"]:
		print "%s\t\t\t\t%s\t" % (eND, "\t".join(especesDup))
	
	res = []
	DCSlen = 0

	# On extrait le profil d'alternance de chaque DCS et on ne garde que les DCS qui alternent
	for gr in lstBlocs:
		# Alternance detectee
		if addDCS(gr):
			res.append(gr)
			DCSlen += len(gr)
		if options["showDCS"]:
			for ((c,i),g,a) in gr:
				print "%s\t%d\t%s\t\t%s" % \
				(c, i, g, "\t".join(["/".join([str(x) for x in set(a[eD])]) for eD in especesDup]))
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
		lst = [x[eD] for (_,_,x) in lstDCS]

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
			curr = lst.pop(0)
			for x in curr:
				# Les alternances sont mesurees entre deux positions consecutives
				for y in last:
					if y == x:
						continue
					count[(x,y)] += (countChr(x)+1) * last[y]
					count[(y,x)] = count[(x,y)]
					#s = (countChr(x)+1) * last[y] + count.get((x,y), 0)
					#count[(x,y)] = count[(y,x)] = s
				# Et aussi entre les paralogues
				for y in curr:
					if y >= x:
						continue
					count[(x,y)] += 1
					count[(y,x)] = count[(x,y)]
					#s = 1 + count.get((x,y), 0)
					#count[(x,y)] = count[(y,x)] = s
			
			# On met a jour last
			for y in last:
				if y not in curr:
					last[y] = 0
			for x in curr:
				last[x] += 1

		return count

	global allDCSe1, allDCSe2
	# On parcourt les DCS en ne gardant que ceux qui alternent
	res = {}
	altern = False
	for eD in especesDup:
		res[eD] = countAltern(dcs, eD)
		if len(res[eD]) > 0 and max(res[eD].values()) > 0:
			altern = True
	if altern:
		for eD in res:
			allDCSe1[eD].append(res[eD])
			allDCSe2[eD].append(set(res[eD]))

	return altern


#
# Range chaque gene ancestral dans son chromosome
#
def buildChrAnc(genesAncCol, chrAncGenes):

	#
	# Renvoie un score (~pourcentage d'especes) qui soutiennent l'attribution d'un gene a un chromosome
	#
	def calcChrAncScore(col, ch):
		espOK = [eND for (_,c,eND) in col if c == ch]
		espNO = [eND for (_,c,eND) in col if c != ch]
		
		def recCalc(node):
			if node in espALL:
				return float(espOK.count(node))/float(espOK.count(node)+espNO.count(node))
			r = []
			for fils in phylTree.branches[node]:
				if len(espALL.intersection(phylTree.species[fils])) != 0:
					r.append(recCalc(fils))
			return utils.myMaths.mean(r)

		if options["usePhylTreeScoring"]:
			espALL = set(espOK + espNO)
			return recCalc(rootNonDup)

		rTot = []
		for gr in espNames:
			r = []
			for e in gr:
				if e in espOK:
					r.append(1)
				elif e in espNO:
					r.append(0)
			if len(r) > 0:
				rTot.append( utils.myMaths.mean(r) )
		return utils.myMaths.mean(rTot)
		
	
	espNames = [[phylTree.officialName[x] for x in g] for g in especesNonDupGrp]
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


#
# Affichage des chromosomes ancestraux
#
def printColorAncestr(genesAnc, chrAncGenes):
	
	print >> sys.stderr, "Impression des associations genes / chromosomes ancestraux ... ",

	chrNames = sorted(chrAncGenes)
	
	if options["showQuality"]:
		print "\t\t%s" % "\t".join([str(c) for c in chrNames])

	for c in chrNames:
		nb = 0
		for i in chrAncGenes[c]:
			nb += 1
			if options["showQuality"]:
				print "%s\t%d\t%s\t%.2f" % (c, nb, "\t".join(["%.2f" % (100*x) for (x,_) in col[i]]), 100*col[i][c-1][0])
			if options["showAncestralGenome"]:
				print c, " ".join(genesAnc[i].names)
		
	print >> sys.stderr, sum([len(chrAncGenes[c]) for c in chrNames]), "genes dans le genome ancestral"


########
# MAIN #
########


# Arguments
(noms_fichiers, options) = utils.myTools.checkArgs( \
	["genesAncestraux.list", "phylTree.conf"],
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
	else:
		rootDup = x
phylTree.loadSpeciesFromList(especesNonDup+especesDup, options["genesFile"])
genesAnc = utils.myGenomes.AncestralGenome(noms_fichiers["genesAncestraux.list"])
lstGenesAnc = genesAnc.lstGenes[utils.myGenomes.AncestralGenome.defaultChr]
(para,orthos) = buildParaOrtho(lstGenesAnc)

# De quoi stocker les DCS
col = [[] for i in xrange(len(lstGenesAnc))]
allDCSe1 = dict( [(eD,[]) for eD in especesDup] )
allDCSe2 = dict( [(eD,[]) for eD in especesDup] )
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
		
def scorePaireDCS(i1, i2):

	val = 0
	nb = 0
	for e in especesDup:
		
		if (len(allDCSe2[e][i1]) != 0) and (len(allDCSe2[e][i2]) != 0):
			nb += 1
		for x in (allDCSe2[e][i1] & allDCSe2[e][i2]):
			val += min(allDCSe1[e][i1][x], allDCSe1[e][i2][x])
	
	if nb == 0:
		return 0
	return val / nb

	# On calcule par une moyenne les autres distances
	scores = phylTree.newCommonNamesMapperInstance()
	for e in especesDup:
		if len(allDCSe2[e][i1]) == 0 or len(allDCSe2[e][i2]) == 0:
			continue
		val = 0
		for x in (allDCSe2[e][i1] & allDCSe2[e][i2]):
			val += min(allDCSe1[e][i1][x], allDCSe1[e][i2][x])
		scores[e] = val

	# On calcule par une moyenne les autres distances
	return phylTree.calcDist(scores)


print >> sys.stderr, "Comparaison des %d DCS ..." % len(allDCS),
phylTree.initCalcDist(rootDup, False)
walktrapInstance = utils.walktrap.WalktrapLauncher()
walktrapInstance.updateFromFunc(range(len(allDCS)), scorePaireDCS)
print >> sys.stderr, "Lancement de walktrap ...",
walktrapInstance.doWalktrap()

# A partir d'ici, on a une association DCS <-> chromosomes, on retombe sur la regle de base, le vote a la majorite

chrInd = 0
for (nodes,cuts,_,dend) in walktrapInstance.res:
	print >> sys.stderr, "Communaute de %d noeuds:" % len(nodes)
	# Un point de coupure virtuel si il n'y a pas
	cuts.append( (1,0) )
	# Les clusterings
	res = [(alpha,relevance,dend.cut(alpha)) for (alpha,relevance) in cuts]
	# Le choix par defaut
	x = 0
	if utils.myTools.stdinInput and (len(res) > 1):
		# Si on peut, on propose a l'utilisateur de choisir
		for (alpha,relevance,(clusters,lonely)) in res:
			print >> sys.stderr, "> alpha=%f relevance=%f clusters=%d size=%d lonely=%d sizes={%s}" % \
				(alpha,relevance,len(clusters),sum([len(c) for c in clusters]),len(lonely),utils.myMaths.myStats([len(c) for c in clusters]))
		while True:
			try:
				print >> sys.stderr, "Choix ? ",
				x = int(raw_input())
				break
			except Exception:
				pass
	(alpha,relevance,(clusters,lonely)) = res[x]
	print >> sys.stderr, "Choix de alpha=%f relevance=%f clusters=%d size=%d lonely=%d sizes={%s}" % \
		(alpha,relevance,len(clusters),sum([len(c) for c in clusters]),len(lonely),utils.myMaths.myStats([len(c) for c in clusters]))
	# On enregistre les resultats
	for cl in clusters:
		chrInd += 1
		for g in cl:
			for (_,s,_) in allDCS[g]:
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
