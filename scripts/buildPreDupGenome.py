#! /usr/bin/python2.4

# Ce script prend en entree les paralogues tetraodon, les orthologues d'une
#    espece avec tetraodon et assigne pour chaque orthologue ses origines
#    possibles chez tetraodon en regardant a quels paralogues il correspond

# On scanne le genome de l'espece
# Pour chaque orthologue, on dresse une liste des origines possibles des
#    genes en fonction des paralogues
# On regroupe les orthologues en faisant des suites des genes qui ont une 
#    origine commune
# Une autre solution est de rassembler les suites de genes qui s'autorisent
#    mutuellement


# INITIALISATION #

# Librairies
import string
import sys
import math
import random
import os

sys.path.append(os.environ['HOME'] + "/work/scripts/utils")
import myOrthos
import myTools
import myMaths

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
def colorAncestr(esp, geneBank, para, orthos, dupEspLst):

	lstBlocs = {}
	genome = geneBank.dicEspeces[esp]

	for e in geneBank.lstEspecesDup:

		if e not in dupEspLst:
			continue
	
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
				lstT = genomeDup.getGenesAt(cT, x1-options["precisionChrAnc"], x2+options["precisionChrAnc"])
				for gT in lstT:
					if gT[0] in para[e]:
						(cT2,_) = genomeDup.dicGenes[para[e][gT[0]]]
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
					col[dicGenesAnc[g[0]][1]].append( (r,len(b),c) )
		

#
# Affichage des chromosomes ancestraux
#
def printColorAncestr(genesAnc, chrAncGenes):
	
	print >> sys.stderr, "Impression des associations genes / chromosomes ancestraux ... ",
	nb = 0	

	for c in chrAncGenes:
		for i in chrAncGenes[c]:
			print c, " ".join(genesAnc[i].names)
		nb += len(chrAncGenes[c])
		
	print >> sys.stderr, nb, "genes dans le genome ancestral"


#
# Range chaque gene ancestral dans son chromosome
#
def buildChrAnc(genesAncCol, chrAncGenes):

	for i in range(len(genesAncCol)):
	
		if len(genesAncCol[i]) == 0:
			# Ce sont des groupes de genes uniquement Tetraodon
			#  ou sans genes Tetraodon
			continue
		
		# Le chromosome le plus frequent
		nb = dict([(x,(0,0,0)) for x in chrAncGenes])
		for (s,l,c) in genesAncCol[i]:
			nb[c] = (nb[c][0]+1, max(nb[c][1], l), max(nb[c][2], s))
		c = myMaths.sortDict(nb)[0]
		
		chrAncGenes[c].append(i)

#
# Cette fonction renvoie les chromosomes ancestraux tels qu'on peut les definir
# grace aux paralogues
#
def loadChrAncIni(nom):

	chrAnc = {}
	espUtil = set([])
	f = open(nom, 'r')
	for ligne in f:
		c = ligne.split()
		dic = {}
		for x in c[1:]:
			if x[0] == '*':
				e = x[1:]
				dic[e] = []
				espUtil.add(e)
			else:
				try:
					x = int(x)
				except Exception:
					pass
				dic[e].append(x)
		chrAnc[c[0]] = dic
	f.close()
	return (chrAnc, espUtil)

# MAIN #

# Arguments
(noms_fichiers, options) = myTools.checkArgs( \
	["genesList.conf", "genesAncestraux.list", "draftPreDupGenome.conf"],
	[("precisionChrAnc", int, 1000000)], \
	"Retrouve le genome pre-duplication grace aux alternances predites" \
)

# Chargement des fichiers
geneBank = myOrthos.GeneBank(noms_fichiers[0])
genesAnc = myOrthos.AncestralGenome(noms_fichiers[1], False)
lstGenesAnc = genesAnc.lstGenes[myOrthos.AncestralGenome.defaultChr]
(para,orthos) = buildParaOrtho(lstGenesAnc, geneBank)
(chrAnc, especesDupliqueesUtilisees) = loadChrAncIni(noms_fichiers[2])

# On colorie les matrices actuelles
blocs = {}
for e in geneBank.lstEspecesNonDup:
	blocs[e] = colorAncestr(e, geneBank, para, orthos, especesDupliqueesUtilisees)

# On colorie les genes ancestraux
print >> sys.stderr, "Synthese des genes ancestraux ",
col = [[] for i in range(len(lstGenesAnc))]
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

