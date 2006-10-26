#! /usr/bin/python2.4

__doc__ = """
Ce script scanne les paralogues d'une espece et essaie de reconstruire les chromosomes pre-duplication de maniere appoximative.
C'est a dire savoir sur quels chromosomes actuels ils se repartissent, et donc les alternances que l'on devra observer.

TODO: Reparer en reprenant la version du M2 !!!
"""


##################
# INITIALISATION #
##################

# Librairies
import string
import sys
import math
import random
import os
import utils.myGenomes
import utils.myTools
import utils.myMaths



#############
# FONCTIONS #
#############

#
# Construit les tables d'associations:
#  - des paralogues de Tetraodon
#  - des orthologues avec Tetraodon
#
def buildParaOrtho(lstGenesAnc, genomeDup):
	para = {}
	ortho = {}
	for g in lstGenesAnc:
		gT = [x for x in g.names if x in genomeDup.dicGenes]
		if len(gT) == 0:
			continue
		for x in gT:
			for y in gT:
				if x != y:
					para[x] = y
		gNT = [x for x in g.names if x not in genomeDup.dicGenes]
		for x in gNT:
			ortho[x] = genomeDup.dicGenes[gT[0]]
	return (para, ortho)



#
# Compte le nombre de passages d'un chromosome a l'autre
#
def countAltern3(genome, orthos, count):

	for c in genome.lstGenes:
		grp3Tet = []
		for tg in genome.lstGenes[c]:
			
			g = tg.names[0]
			
			# Il faut un orthologue avec Tetraodon
			if g not in orthos:
				continue
		
			(KT,_) = orthos[g]
			
			if KT in grp3Tet:
				# On deplace le numero de chromosome en fin de groupe
				grp3Tet.remove(KT)
				grp3Tet.append(KT)
				
			elif len(grp3Tet) < 3:
				# On construit l'alternance
				grp3Tet.append(KT)
			else:
				# On enleve le plus vieux et on rajoute le nouveau
				grp3Tet.pop(0)
				grp3Tet.append(KT)
			
			if len(grp3Tet) == 3:
				# On ajoute 1 au score
				t = grp3Tet[:]
				t.sort()
				t = tuple(t)
				if t in count:
					count[t] += 1
				else:
					count[t] = 1

#
# Renvoie le nombre de paralogues qui lie chaque couple de chromosomes Tetraodon
#
def countPara():
	
	count = {}
	for g1 in para:
		(c1,_) = genomeDup.dicGenes[g1]
		(c2,_) = genomeDup.dicGenes[para[g1]]
		if c1 < c2:
			count[ (c1,c2) ] = count.get( (c1,c2), 0) + 1
	return count

#
# Affiche sous forme de table les quantites de paralogues liant les chromosomes
#
def printNbParaTable():
	dicNb = {}
	for ((c1,c2), v) in countP.items():
		if v in dicNb:
			dicNb[v].append( (c1,c2) )
		else:
			dicNb[v] = [(c1,c2)]
	
	print "NbPara\tNbCouples\tCouples..."
	lst = dicNb.keys()
	lst.sort(reverse = True)
	for v in lst:
		print "%d\t%d\t" % (v,len(dicNb[v])), dicNb[v]


########
# MAIN #
########

# Arguments
(noms_fichiers, options) = utils.myTools.checkArgs( \
	["genesList.conf", "genesAncestraux.list"],
	[("especeDupliquee", str, ""), ("nbMinParaloguesCoupleChr", int, -1), ("coefMinAltern3", float, -1.)], \
	__doc__ \
)

# Chargement des fichiers
geneBank = utils.myGenomes.GeneBank(noms_fichiers["genesList.conf"], [options["especeDupliquee"]])
if options["especeDupliquee"] not in geneBank.lstEspecesDup:
	print >> sys.stderr, "-ERREUR- Pas de -%(especeDupliquee)s- dans la liste des especes dupliquees ..." % options
	sys.exit(1)
genomeDup = geneBank.dicEspeces[options["especeDupliquee"]]
genesAnc = utils.myGenomes.AncestralGenome(noms_fichiers["genesAncestraux.list"], False)
lstGenesAnc = genesAnc.lstGenes[utils.myGenomes.AncestralGenome.defaultChr]
(para,orthos) = buildParaOrtho(lstGenesAnc, genomeDup)


# On initialise les calculs

print >> sys.stderr, "Calculs des 3-alternances ",
count3 = {}
for e in geneBank.lstEspecesNonDup:
	countAltern3(geneBank.dicEspeces[e], orthos, count3)
	sys.stderr.write(".")

for c in count3:
	count3[c] /= len(geneBank.lstEspecesNonDup)
print >> sys.stderr, " OK"


countP = countPara()

if options["nbMinParaloguesCoupleChr"] < 0:
	printNbParaTable()
	sys.exit(0)

# On filtre les couples de paralogues
r = []
for ((c1,c2), v) in countP.items():
	if v >= options["nbMinParaloguesCoupleChr"]:
		r.append( (v, c1, c2) )
r.sort()
r.reverse()
countP = r

chrom = []
nbPara = []
rep = dict([(k,[]) for k in genomeDup.lstGenes])

# 1ere etape on extrait tous les couples qui ne se chevauchent pas
#  (priorite aux mieux representes)
newCountP = []
for (v, c1, c2) in countP:
	if len(rep[c1]) == 0 and len(rep[c2]) == 0:
	
		# Les deux chromosomes sont inutilises: nouveau chromosome ancestral
		print >> sys.stderr, "Nouveau chromosome:", c1, c2
		rep[c1].append(len(chrom))
		rep[c2].append(len(chrom))
		chrom.append( set([c1,c2]) )
		nbPara.append( v )
	else:
		newCountP.append( (v,c1,c2) )
print >> sys.stderr

# 2eme etape, on rajoute le reste
for (v, c1, c2) in newCountP:

	if len(set(rep[c1]).intersection(rep[c2])) > 0:
		continue
	
	best = 0
	bestC = -1
	
	for i in rep[c1] + rep[c2]:
		for c in chrom[i]:
			if c == c1 or c == c2:
				continue
			t = [c,c1,c2]
			t.sort()
			x = float(count3.get(tuple(t), 0)) / float(nbPara[i]+v)
			print >> sys.stderr, "Test de", c1, c2, "sur", chrom[i], "(score %f)" % x
			if x > best:
				best = x
				bestC = i

	if options["coefMinAltern3"] < 0:
		continue

	if best >= options["coefMinAltern3"]:
		print >> sys.stderr, "Ajout de", c1, c2, "sur", chrom[bestC], "(meilleur choix, score %f)" % best
		chrom[bestC].add( c1 )
		chrom[bestC].add( c2 )
		rep[c1].append(bestC)
		rep[c2].append(bestC)
		nbPara[bestC] += v
	else:
		print >> sys.stderr, "Nouveau chromosome (aucune 3-alternance satisfaisante):", c1, c2
		rep[c1].append(len(chrom))
		rep[c2].append(len(chrom))
		chrom.append( set([c1,c2]) )
		nbPara.append( v )

	print >> sys.stderr, "\n"

if options["coefMinAltern3"] < 0:
	sys.exit(0)

# Affichage des resultats
print >> sys.stderr, "%d chromosomes trouves :" % len(chrom)
for i in range(len(chrom)):
	print >> sys.stderr, "Chromosome %c (%d):" % (chr(65+i), nbPara[i]),
	for c in chrom[i]:
		print >> sys.stderr, " ", c,
	print >> sys.stderr

