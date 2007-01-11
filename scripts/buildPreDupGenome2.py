#! /users/ldog/muffato/python -OO

__doc__ = """
Ce script scanne chaque genome d'espece non dupliquee en le comparant a chaque genome duplique.
Pour chaque gene, on dresse une liste de ses origines possibles en utilisant son orthologue chez l'espece dupliquee puis en etudiant les paralogues autour.
On regroupe les genes en faisant des suites qui ont une origine commune.
Chaque gene ancestral recoit son chromosome ancestral par un vote a la majorite en fonction de chaque annotation.
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
				ortho[e][x] = genomeDup.dicGenes[gT[0]]
	return (para, ortho)

#
# Cree les regions de chromosomes ancestraux dans la classe d'orthologues
#
def colorAncestr(esp, e, geneBank, para, orthos, col, dicGenesAnc, chrAnc):

	print >> sys.stderr, "Decoupage de", esp, "avec", e, "",

	nbOrthos = 0
	nbBlocs = 0
	genome = geneBank.dicEspeces[esp]
	genomeDup = geneBank.dicEspeces[e]
	orthosDup = orthos[e]
	parasDup = para[e]

	# On parcourt les chromosomes de l'espece
	for c in genome.lstChr:
		
		# Eviter d'essayer de creer des DCS sur des scaffolds trop petits
		if len(genome.lstGenes[c]) < options["minChrLen"]:
			continue
		
		lastOrig = set([])
		bloc = []
		(lastCT,lastI) = ("",0)
		
		for tg in genome.lstGenes[c]:
			g = tg.names[0]

			# Il faut un orthologue avec Tetraodon
			if g not in orthos[e]:
				continue
			nbOrthos += 1
			
			(cT,i) = orthosDup[g]
			ok = False
			if cT == lastCT:
				ok = True
			else:
				gTn1 = [gT.names[0] for gT in genomeDup.getGenesNearN(lastCT, lastI, options["precisionChrAnc"])]
				gTn2 = [gT.names[0] for gT in genomeDup.getGenesNearN(cT, i, options["precisionChrAnc"])]
				gTp = set(utils.myMaths.flatten([parasDup[s] for s in gTn2 if s in parasDup]))
				print ":", len(gTn2), len(gTp)
				if len(gTp.intersection(gTn1)) > 0:
					ok = True
				
			if not ok:
				if len(bloc) != 0:
					nbBlocs += 1
					print "----"
					addDCS(bloc, col, dicGenesAnc, chrAnc, esp, e)
				bloc = []
			print g, c, genomeDup.lstGenes[cT][i].names[0], cT
			bloc.append( (g,cT) )
			(lastCT,lastI) = (cT,i)

		if len(bloc) != 0:
			nbBlocs += 1
			print "----"
			addDCS(bloc, col, dicGenesAnc, chrAnc, esp, e)
		sys.stderr.write(".")
	
	print >> sys.stderr, "", nbBlocs, "blocs, pour", nbOrthos, "genes orthologues"
	sys.exit(0)

#
# Construit les tables d'association "nom de gene" -> couleur
# Prend la couleur la plus probable
#
def addDCS(bloc, col, dicGenesAnc, chrAnc, eNonDup, eDup):

	# Renvoie le chromosome ancestral qui a le meilleur score pour la region consideree
	#print len(bloc)
	score = {}
	for c in chrAnc:
		score[c] = len([cT for (_,cT) in bloc if cT in chrAnc[c][eDup]])
	s = max(score.values())
	
	if s == 0:
		sys.stderr.write("!")
		return
	
	r = float(s-1) / float(len(bloc))
	for c in chrAnc:
		if score[c] == s:
			for g in bloc:
				col[dicGenesAnc[g[0]][1]].append( (r,len(bloc),c,(eNonDup,eDup)) )
		

#
# Affichage des chromosomes ancestraux
#
def printColorAncestr(genesAnc, chrAncGenes):
	
	print >> sys.stderr, "Impression des associations genes / chromosomes ancestraux ... ",
	nb = 0	

	lstChr = sorted(chrAncGenes)
	
	if options["showStats"]:
		print "\t\t%s" % "\t".join(lstChr)
	
	for c in lstChr:
		for i in chrAncGenes[c]:
			nb += 1
			if options["showStats"]:
				print "%s\t%d\t%s" % (c, nb, "\t".join(["%.2f" % (100*col[i][x]) for x in lstChr]))
			else:
				print c, " ".join(genesAnc[i].names)
		
	print >> sys.stderr, nb, "genes dans le genome ancestral"


#
# Range chaque gene ancestral dans son chromosome
#
def buildChrAnc(genesAncCol, chrAncGenes):

	#
	# Renvoie un score (~pourcentage d'especes) qui soutiennent l'attribution d'un gene a un chromosome
	#
	def calcChrAncScore(col, ch):
		espOK = [eND for (_,_,c,(eND,eD)) in col if c == ch]
		espNO = [eND for (_,_,c,(eND,eD)) in col if c != ch]
		espALL = set(espOK + espNO)
		
		def recCalc(node):
			if node in espALL:
				return float(espOK.count(node))/float(espOK.count(node)+espNO.count(node))
			r = []
			for fils in phylTree.branches[node]:
				if len(espALL.intersection(phylTree.species[fils])) != 0:
					r.append(recCalc(fils))
			return utils.myMaths.moyenne(r)

		return recCalc(rootNonDup)

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
				dic[e] = []
			else:
				try:
					x = int(x)
				except Exception:
					pass
				dic[e].append(x)
		chrAnc[c[0]] = dic
	f.close()
	return chrAnc

# MAIN #


# Arguments
(noms_fichiers, options) = utils.myTools.checkArgs( \
	["genesList.conf", "genesAncestraux.list", "draftPreDupGenome.conf", "phylTree.conf"],
	[("minChrLen",int,20), ("precisionChrAnc",int,10), ("especesNonDup",str,""), ("especesDup",str,""), ("showStats",bool,False)], \
	__doc__ \
)

# Chargement des fichiers
especesDup = options["especesDup"].split(',')
especesNonDup = options["especesNonDup"].split(',')
phylTree = utils.myBioObjects.PhylogeneticTree(noms_fichiers["phylTree.conf"])
rootNonDup = [x for x in phylTree.branches[phylTree.root] if len(set(phylTree.species[x]).intersection(especesDup)) == 0][0]
geneBank = utils.myGenomes.GeneBank(noms_fichiers["genesList.conf"], especesNonDup+especesDup)
genesAnc = utils.myGenomes.AncestralGenome(noms_fichiers["genesAncestraux.list"])
lstGenesAnc = genesAnc.lstGenes[utils.myGenomes.AncestralGenome.defaultChr]
(para,orthos) = buildParaOrtho(lstGenesAnc, geneBank)
chrAnc = loadChrAncIni(noms_fichiers["draftPreDupGenome.conf"])

# On colorie les matrices actuelles
col = [[] for i in xrange(len(lstGenesAnc))]
for eND in geneBank.lstEspecesNonDup:
	for eD in geneBank.lstEspecesDup:
		colorAncestr(eND, eD, geneBank, para, orthos, col, genesAnc.dicGenes, chrAnc)

# On construit les chromosomes ancestraux
print >> sys.stderr, "Synthese des genes ancestraux ...",
chrAncGenes = dict([ (x,[]) for x in chrAnc])
buildChrAnc(col, chrAncGenes)
print >> sys.stderr, "OK"
	
# On affiche le resultat
printColorAncestr(lstGenesAnc, chrAncGenes)

