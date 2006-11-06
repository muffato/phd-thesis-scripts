#! /usr/bin/python2.4

__doc__ = """
Ce script scanne chaque genome d'espece non dupliquee en le comparant a chaque genome duplique.
Pour chaque gene, on dresse une liste de ses origines possibles en utilisant son orthologue chez l'espece dupliquee puis en etudiant les paralogues autour.
On regroupe les genes en faisant des suites qui ont une origine commune.
Chaque gene ancestral recoit son chromosome ancestral par un vote a la majorite en fonction de chaque annotation.
"""

# INITIALISATION #

# Librairies
import sys
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
				for y in gT:
					if x != y:
						para[e][x] = y
			gNT = [x for x in g.names if x not in genomeDup.dicGenes]
			for x in gNT:
				ortho[e][x] = genomeDup.dicGenes[gT[0]]
	return (para, ortho)

#
# Cree les regions de chromosomes ancestraux dans la classe d'orthologues
#
def colorAncestr(esp, geneBank, para, orthos):

	lstBlocs = {}
	genome = geneBank.dicEspeces[esp]

	for e in geneBank.lstEspecesDup:

		print >> sys.stderr, "Decoupage de", esp, "avec", e, "",

		nbOrthos = 0
		lstBlocs[e] = []
		genomeDup = geneBank.dicEspeces[e]

		# On parcourt les chromosomes de l'espece
		for c in genome.lstGenes:
			lastOrig = set([])
			bloc = []
			for tg in genome.lstGenes[c]:
				g = tg.names[0]

				# Il faut un orthologue avec Tetraodon
				if g not in orthos[e]:
					continue
				nbOrthos += 1
				(cT,i) = orthos[e][g]
				x1 = genomeDup.lstGenes[cT][i].beginning
				x2 = genomeDup.lstGenes[cT][i].end
				orig = set([cT])
				for gT in genomeDup.getGenesAt(cT, x1-options["precisionChrAnc"], x2+options["precisionChrAnc"]):
					s = gT.names[0]
					if s in para[e]:
						(cT2,_) = genomeDup.dicGenes[para[e][s]]
						orig.add(cT2)
				
				nouvOrig = lastOrig & orig
				if len(nouvOrig) == 0:
					if len(bloc) != 0:
						lstBlocs[e].append(bloc)
					bloc = []
					nouvOrig = orig
				bloc.append( (g,cT) )
				lastOrig = nouvOrig

			lstBlocs[e].append(bloc)
			sys.stderr.write(".")
		print >> sys.stderr, "", len(lstBlocs[e]), "blocs, pour", nbOrthos, "genes orthologues"

	return lstBlocs

#
# Construit les tables d'association "nom de gene" -> couleur
# Prend la couleur la plus probable
#
def buildColorTable(lstBlocs, col, dicGenesAnc, chrAnc):

	#
	# Renvoie le chromosome ancestral qui a le meilleur score pour la region consideree
	#
	def getMaxScore(bloc, especeDup):
		score = {}
		for c in chrAnc:
			score[c] = 0
			for (_,cT) in bloc:
				if cT in chrAnc[c][especeDup]:
					score[c] += 1
		m = max(score.values())
		return (m, [c for c in chrAnc if score[c] == m])

	for e in lstBlocs:
		for b in lstBlocs[e]:
			if len(b) == 0:
				continue
			(s, l) = getMaxScore(b, e)
			r = float(s-1) / float(len(b))
			for c in l:
				for g in b:
					col[dicGenesAnc[g[0]][1]].append( (r,len(b),c,e) )
		

#
# Affichage des chromosomes ancestraux
#
def printColorAncestr(genesAnc, chrAncGenes):
	
	print >> sys.stderr, "Impression des associations genes / chromosomes ancestraux ... ",
	nb = 0	

	lstChr = chrAncGenes.keys()
	lstChr.sort()
	if options["showStats"]:
		print "\t\t%s" % "\t".join(lstChr)
	
	for c in lstChr:
		for i in chrAncGenes[c]:
			nb += 1
			if options["showStats"]:
				print "%s\t%d\t%s\t%d" % (c, nb, "\t".join([str(col[i][x][0]) for x in lstChr]), (100*col[i][c][0])/sum([col[i][x][0] for x in lstChr]))
			else:
				print c, " ".join(genesAnc[i].names)
		
	print >> sys.stderr, nb, "genes dans le genome ancestral"


#
# Range chaque gene ancestral dans son chromosome
#
def buildChrAnc(genesAncCol, chrAncGenes):

	lstChr = chrAncGenes.keys()
	lstChr.sort()
	
	for i in xrange(len(genesAncCol)):
	
		if len(genesAncCol[i]) == 0:
			# Certains genes n'ont pas de chance !
			continue
		
		# Le chromosome le plus frequent
		nb = dict([(x,(0,0,0,[])) for x in chrAncGenes])
		for (s,l,c,e) in genesAncCol[i]:
			nb[c] = (nb[c][0]+1, max(nb[c][1], l), max(nb[c][2], s), nb[c][3]+[e])
		tmp = utils.myMaths.sortDict(nb)
		c = tmp[0]
		#genesAncCol[i] = nb
		
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
	["genesList.conf", "genesAncestraux.list", "draftPreDupGenome.conf"],
	[("precisionChrAnc", int, 1000000), ("especesNonDup",str,""), ("especesDup",str,""), ("showStats",bool,False)], \
	__doc__ \
)

# Chargement des fichiers
especesNonDup = options["especesNonDup"].split(',')
especesDup = options["especesDup"].split(',')
geneBank = utils.myGenomes.GeneBank(noms_fichiers["genesList.conf"], especesNonDup+especesDup)
genesAnc = utils.myGenomes.loadGenome(noms_fichiers["genesAncestraux.list"])
lstGenesAnc = genesAnc.lstGenes[utils.myGenomes.AncestralGenome.defaultChr]
(para,orthos) = buildParaOrtho(lstGenesAnc, geneBank)
chrAnc = loadChrAncIni(noms_fichiers["draftPreDupGenome.conf"])

# On colorie les matrices actuelles
blocs = {}
for e in geneBank.lstEspecesNonDup:
	blocs[e] = colorAncestr(e, geneBank, para, orthos)

# On colorie les genes ancestraux
print >> sys.stderr, "Synthese des genes ancestraux ",
col = [[] for i in xrange(len(lstGenesAnc))]
for e in blocs:
	buildColorTable(blocs[e], col, genesAnc.dicGenes, chrAnc)
	sys.stderr.write(".")

# On construit les chromosomes ancestraux
print >> sys.stderr, " ...",
chrAncGenes = dict([ (x,[]) for x in chrAnc])
buildChrAnc(col, chrAncGenes)
print >> sys.stderr, "OK"
	
# On affiche le resultat
printColorAncestr(lstGenesAnc, chrAncGenes)

