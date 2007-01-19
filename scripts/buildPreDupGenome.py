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



# FONCTIONS #


#
# Cette fonction renvoie les chromosomes ancestraux tels qu'on peut les definir
# grace aux paralogues
#
def loadChrAncIni(nom):

	chrAnc = {}
	f = utils.myTools.myOpenFile(nom, 'r')
	for ligne in f:
		c = ligne.split()
		dic = {}
		for x in c[1:]:
			if x[0] == '*':
				e = x[1:]
				dic[e] = set([])
			else:
				try:
					x = int(x)
				except Exception:
					pass
				dic[e].add(x)
		chrAnc[c[0]] = dic
	f.close()
	return chrAnc


#
# Construit les tables d'associations pour chaque especes dupliquee
#  des paralogues et des orthologues
#
def buildParaOrtho(lstGenesAnc, geneBank):
	para = dict([(e,{}) for e in geneBank.lstEspecesDup])
	ortho = dict([(e,{}) for e in geneBank.lstEspecesDup])
	
	for g in lstGenesAnc:
		for e in geneBank.lstEspecesDup:
			genomeDup = geneBank.dicEspeces[e]
			gT = [x for x in g.names if x in genomeDup.dicGenes]
			if len(gT) == 0:
				continue
			for x in gT:
				para[e][x] = [y for y in gT if y != x]
			gNT = [x for x in g.names if x not in genomeDup.dicGenes]
			for x in gNT:
				ortho[e][x] = [genomeDup.dicGenes[y] for y in gT]
	return (para, ortho)


#
# Cree les regions de chromosomes ancestraux dans la classe d'orthologues
#
def colorAncestr(esp, e, geneBank, para, orthos):

	print >> sys.stderr, "Decoupage de", esp, "avec", e, "",

	nbOrthos = 0
	nbBlocs = 0
	genome = geneBank.dicEspeces[esp]
	genomeDup = geneBank.dicEspeces[e]
	orthosDup = orthos[e]
	parasDup = para[e]
	lstBlocs = []

	# On parcourt les chromosomes de l'espece
	for c in genome.lstChr:
		
		# Eviter d'essayer de creer des DCS sur des scaffolds trop petits
		if len(genome.lstGenes[c]) < options["minChrLen"]:
			continue
		
		bloc = []
		lastCT = []
		lastGT = set([])
		
		for tg in genome.lstGenes[c]:
			g = tg.names[0]

			# Il faut un orthologue avec Tetraodon
			if g not in orthos[e]:
				continue
			nbOrthos += 1
	
			for (cT,i) in orthosDup[g]:
				gTn = [gT.names[0] for gT in genomeDup.getGenesNearN(cT, i, options["precisionChrAnc"]) if gT.names[0] in parasDup]
				
				if cT in lastCT:
					ok = True
				else:
					gTp = utils.myMaths.flatten([parasDup[s] for s in gTn])
					ok = (len(lastGT.intersection(gTp)) > 0)
				
				if ok:
					break
		
			if not ok:
				if len(bloc) != 0:
					nbBlocs += 1
					lstBlocs.append(bloc)
				bloc = []
				lastGT = set([])
				
			bloc.append( g )
			lastCT = [cT for (cT,i) in orthosDup[g]]
			lastGT.update(gTn)

		if len(bloc) != 0:
			nbBlocs += 1
			lstBlocs.append(bloc)
		sys.stderr.write(".")
	
	print >> sys.stderr, "", nbBlocs, "blocs, pour", nbOrthos, "genes orthologues"

	return lstBlocs



# Synthese des DCS

def doSynthese(combin, eND, geneBank, orthos, col, dicGenesAnc, chrAnc):
	
	print >> sys.stderr, "Synthese des decoupages de", eND, "...",
	
	# On retrouve les groupes de genes pour chaque tetrapode
	lstBlocs = []
	for gr in combin:
		l = []
		for g in gr:
			p = geneBank.dicEspeces[eND].dicGenes[g]
			a = {}
			for eD in geneBank.lstEspecesDup:
				a[eD] = [c for (c,_) in orthos[eD].get(g,[])]
			l.append( (p,g,a) )
		l.sort()
		lstBlocs.append(l)
	lstBlocs.sort()
	
	print >> sys.stderr, len(lstBlocs), "blocs pour", sum([len(x) for x in lstBlocs]), "orthologues",

	if options["showDCS"]:
		print "%s\t\t\t\t%s\t" % (eND, "\t".join(geneBank.lstEspecesDup))
	
	nbDCS = 0
	DCSlen = 0

	for gr in lstBlocs:
		cc = addDCS(gr, col, dicGenesAnc, chrAnc, eND)
		if cc != "":
			nbDCS += 1
			DCSlen += len(gr)
		if options["showDCS"]:
			for ((c,i),g,a) in gr:
				print "%s\t%d\t%s\t\t%s\t%s" % \
				(c, i, g, "\t".join(["/".join([str(x) for x in set(a[eD])]) for eD in geneBank.lstEspecesDup]), cc)
			print "---"

	print >> sys.stderr, "/", nbDCS, "DCS pour", DCSlen, "orthologues"


#
# Construit les tables d'association "nom de gene" -> couleur
# Prend la couleur la plus probable
#
def addDCS(bloc, col, dicGenesAnc, chrAnc, eNonDup):

	score = dict([(c,0) for c in chrAnc])
	for eDup in geneBank.lstEspecesDup:
		l = [a[eDup] for (_,_,a) in bloc if len(a[eDup]) > 0]
		# On veut une veritable alternance
		if len(set(utils.myMaths.flatten(l))) < 2:
			continue
		for c in chrAnc:
			score[c] += float(len([x for x in l if len(chrAnc[c][eDup].intersection(x))>0 ])) / float(len(l))
			
	s = max(score.values())

	if s == 0:
		return ""
	
	for c in chrAnc:
		if score[c] == s:
			cc = c
			for g in bloc:
				col[dicGenesAnc[g[1]][1]].append( (len(bloc),c,eNonDup) )
	return cc


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
			return utils.myMaths.moyenne(r)

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
				rTot.append( utils.myMaths.moyenne(r) )
		return utils.myMaths.moyenne(rTot)
		
	
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
		print "\t\t%s" % "\t".join(lstChr)
	
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
	["genesList.conf", "genesAncestraux.list", "draftPreDupGenome.conf", "phylTree.conf"],
	[("minChrLen",int,20), ("precisionChrAnc",int,10), ("usePhylTreeScoring",bool,False), \
	("especesNonDup",str,""), ("especesDup",str,""), \
	("showDCS",bool,False), ("showQuality",bool,False), ("showAncestralGenome",bool,True)], \
	__doc__ \
)

# Chargement des fichiers
especesDup = options["especesDup"].split(',')
especesNonDupGrp = [x.split('+') for x in options["especesNonDup"].split(',')]
especesNonDup = utils.myMaths.flatten(especesNonDupGrp)
phylTree = utils.myBioObjects.PhylogeneticTree(noms_fichiers["phylTree.conf"])
rootNonDup = [x for x in phylTree.branches[phylTree.root] if len(set(phylTree.species[x]).intersection(especesDup)) == 0][0]
geneBank = utils.myGenomes.GeneBank(noms_fichiers["genesList.conf"], especesNonDup+especesDup)
genesAnc = utils.myGenomes.AncestralGenome(noms_fichiers["genesAncestraux.list"])
lstGenesAnc = genesAnc.lstGenes[utils.myGenomes.AncestralGenome.defaultChr]
(para,orthos) = buildParaOrtho(lstGenesAnc, geneBank)
chrAnc = loadChrAncIni(noms_fichiers["draftPreDupGenome.conf"])

# On colorie les matrices actuelles
col = [[] for i in xrange(len(lstGenesAnc))]

# Decoupage de chaque tetrapode
for eND in geneBank.lstEspecesNonDup:
	combin = utils.myTools.myCombinator([])
	# Avec chaque poisson
	for eD in geneBank.lstEspecesDup:
		lstBlocs = colorAncestr(eND, eD, geneBank, para, orthos)
		for bloc in lstBlocs:
			combin.addLink(bloc)
	
	doSynthese(combin, eND, geneBank, orthos, col, genesAnc.dicGenes, chrAnc)


# On construit les chromosomes ancestraux
print >> sys.stderr, "Synthese des genes ancestraux ...",
chrAncGenes = dict([ (x,[]) for x in chrAnc])
buildChrAnc(col, chrAncGenes)
print >> sys.stderr, "OK"
	
# On affiche le resultat
printColorAncestr(lstGenesAnc, chrAncGenes)

