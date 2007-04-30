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
import utils.myBioObjects
import utils.myGenomes
import utils.myTools
import utils.myMaths
import utils.myCommunities


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
		
		bloc = []
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
				gTn = [gT.names[0] for gT in genomeDup.getGenesNearN(cT, i, options["precisionChrAnc"]) if gT.names[0] in parasDup]
				
				# Si on reste sur le meme chromosome, on continue le DCS
				if cT in lastCT:
					break
				# Sinon, il faut qu'il y ait un paralogue qui nous fasse revenir vers une region environnante deja visitee
				else:
					gTp = utils.myMaths.flatten([parasDup[s] for s in gTn])
					if len(lastGT.intersection(gTp)) > 0:
						break
				
			# Si on n'a pas fait de break, c'est qu'aucun orthologue ne convient, il faut arreter le DCS
			else:
				# On l'enregistre et on repart de zero
				if len(bloc) != 0:
					lstBlocs.append(bloc)
				bloc = []
				lastGT = set()
				
			# On rajoute les infos du gene qu'on vient de lire
			bloc.append( g )
			lastCT = [cT for (cT,i) in orthosDup[g]]
			lastGT.update(gTn)
		
		# Ne pas oublier le bloc courant
		if len(bloc) != 0:
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
		count = {}
		# La derniere ligne lue
		last = {}
		# On parcourt la liste
		while len(lst) > 0:
			curr = lst.pop(0)
			for x in curr:
				# Les alternances sont mesurees entre deux positions consecutives
				for y in last:
					if y == x:
						continue
					s = (countChr(x)+1) * last[y] + count.get((x,y), 0)
					count[(x,y)] = count[(y,x)] = s
				# Et aussi entre les paralogues
				for y in curr:
					if y >= x:
						continue
					s = 1 + count.get((x,y), 0)
					count[(x,y)] = count[(y,x)] = s
			
			# On met a jour last
			for y in last:
				if y not in curr:
					last[y] = 0
			for x in curr:
				last[x] = last.get(x,0) + 1

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
		for gr in especesNonDupGrp:
			r = []
			for e in gr:
				if e in espOK:
					r.append(1)
				elif e in espNO:
					r.append(0)
			if len(r) > 0:
				rTot.append( utils.myMaths.mean(r) )
		return utils.myMaths.mean(rTot)
		
	
	for i in xrange(len(genesAncCol)):
	
		if len(genesAncCol[i]) == 0:
			# Certains genes n'ont pas de chance !
			continue
	
		nb = dict([(x,calcChrAncScore(genesAncCol[i], x)) for x in chrAncGenes])

		tmp = utils.myMaths.sortDict(nb)
		c = tmp[0]
		genesAncCol[i] = nb
		
		chrAncGenes[c].append(i)


#
# Affichage des chromosomes ancestraux
#
def printColorAncestr(genesAnc, chrAncGenes):
	
	print >> sys.stderr, "Impression des associations genes / chromosomes ancestraux ... ",

	lstChr = sorted(chrAncGenes)
	
	if options["showQuality"]:
		print "\t\t%s" % "\t".join([str(c) for c in lstChr])
	
	for c in lstChr:
		nb = 0
		for i in chrAncGenes[c]:
			nb += 1
			if options["showQuality"]:
				print "%s\t%d\t%s\t%.2f" % (c, nb, "\t".join(["%.2f" % (100*col[i][x]) for x in lstChr]), 100*col[i][c])
			if options["showAncestralGenome"]:
				print c, " ".join(genesAnc[i].names)
		
	print >> sys.stderr, sum([len(chrAncGenes[c]) for c in lstChr]), "genes dans le genome ancestral"



# MAIN #


# Arguments
(noms_fichiers, options) = utils.myTools.checkArgs( \
	["genesAncestraux.list", "phylTree.conf"],
	[("minChrLen",int,20), ("precisionChrAnc",int,25), ("usePhylTreeScoring",bool,False), \
	("especesNonDup",str,""), ("especesDup",str,""), \
	("genesFile",str,"~/work/data/genes/genes.%s.list.bz2"), \
	("showDCS",bool,False), ("showQuality",bool,False), ("showAncestralGenome",bool,True)], \
	__doc__ \
)

# Chargement des fichiers
phylTree = utils.myBioObjects.PhylogeneticTree(noms_fichiers["phylTree.conf"])
especesDup = options["especesDup"].split(',')
especesNonDupGrp = [x.split('+') for x in options["especesNonDup"].split(',')]
especesNonDup = utils.myMaths.flatten(especesNonDupGrp)
for x in phylTree.branches[phylTree.root]:
	if phylTree.getFirstParent(x, especesDup[0]) == x:
		rootDup = x
	else:
		rootNonDup = x
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


print >> sys.stderr, "Lancement des communautes sur les %d DCS :" % len(allDCS),
phylTree.initCalcDist(rootDup, False)
lstLstComm = utils.myCommunities.launchCommunitiesBuild(items = range(len(allDCS)), scoreFunc = scorePaireDCS)

# A partir d'ici, on a une association DCS <-> chromosomes, on retombe sur la regle de base, le vote a la majorite

chrInd = 0
for lstComm in lstLstComm:
	if len(lstComm) == 0:
		continue
	(alpha,relevance,clusters,lonely) = lstComm[0]
	print >> sys.stderr, "Resultat alpha=%f relevance=%f clusters=%d size=%d lonely=%d" % \
		(alpha,relevance,len(clusters),sum([len(c) for c in clusters]),len(lonely))
	for i in range(len(clusters)):
		chrInd += 1
		for g in clusters[i]:
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

