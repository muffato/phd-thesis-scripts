#! /users/ldog/muffato/python

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

sys.path.append(os.environ['HOME'] + "/M2/scripts/utils")
import myOrthos
import myTools
import myMaths

# FONCTIONS #

#
# Construit les tables d'associations:
#  - des paralogues de Tetraodon
#  - des orthologues avec Tetraodon
#
def buildParaOrtho(lstGenesAnc, genomeDup):
	para = {}
	ortho = {}
	for g in lstGenesAnc:
		gT = [x for x in g if x in genomeDup.dicGenes]
		if len(gT) == 0:
			continue
		for x in gT:
			for y in gT:
				if x != y:
					para[x] = y
		gNT = [x for x in g if x not in genomeDup.dicGenes]
		for x in gNT:
			ortho[x] = genomeDup.dicGenes[gT[0]]
	return (para, ortho)

#
# Cree les regions de chromosomes ancestraux dans la classe d'orthologues
#
def colorAncestr(genome, genomeDup, para, orthos):

	print >> sys.stderr, "Decoupage de", genome.nom, "",

	lstBlocs = []
	nbOrthos = 0

	# On parcourt les chromosomes de l'espece
	for c in genome.lstGenes:
		lastOrig = set([])
		bloc = []
		for (_,_,_,g) in genome.lstGenes[c]:
			# Il faut un orthologue avec Tetraodon
			if g not in orthos:
				continue
			nbOrthos += 1
			(cT,i) = orthos[g]
			(x1,x2,_,_) = genomeDup.lstGenes[cT][i]
			orig = set([cT])
			for gT in genomeDup.getGenesAt(cT, x1-options["precisionChrAnc"], x2+options["precisionChrAnc"]):
				s = gT.names[0]
				if s in para:
					(cT2,_) = genomeDup.dicGenes[para[s]]
					orig.add(cT2)
			
			nouvOrig = lastOrig & orig
			if len(nouvOrig) == 0:
				if len(bloc) != 0:
					lstBlocs.append(bloc)
				bloc = []
				nouvOrig = orig
			bloc.append( (g,cT) )
			lastOrig = nouvOrig

		lstBlocs.append(bloc)
		sys.stderr.write(".")
	print >> sys.stderr, "", len(lstBlocs), "blocs, pour", nbOrthos, "genes orthologues"
	return lstBlocs

#
# Construit les tables d'association "nom de gene" -> couleur
# Prend la couleur la plus probable
#
def buildColorTable(lstBlocs, col, dicGenesAnc, chrAnc):

	#
	# Renvoie le chromosome ancestral qui a le meilleur score pour la region consideree
	#
	def getMaxScore(bloc):
		score = []
		for i in range(len(chrAnc)):
			s = 0
			for (_,cT) in bloc:
				if cT in chrAnc[i]:
					s += 1
			score.append( s )
		m = max(score)
		return (m, [i for i in range(len(chrAnc)) if score[i] == m])


	print >> sys.stderr, "Assignation des chromosomes ancestraux sur", "... ",
	
	for b in lstBlocs:
		(s, l) = getMaxScore(b)
		r = float(s-1) / float(len(b))
		for c in l:
			for g in b:
				col[dicGenesAnc[g[0]][1]].append( (r,len(b),c) )
	
	print >> sys.stderr, "OK"

#
# Affichage des chromosomes ancestraux
#
def printColorAncestr(genesAnc, chrAncGenes):
	
	print >> sys.stderr, "Impression des associations genes / chromosomes ancestraux ... ",
	nb = 0	

	for c in range(len(chrAncGenes)):
		for i in chrAncGenes[c]:
			r = chr(65+c)
			for s in genesAnc[i]:
				r += " " + s
			print r
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
		
		if options["bestScore"]:
			# Le meilleur score de purete
			c = max(genesAncCol[i])[2]
		else:
			# Le chromosome le plus frequent
			nb = [[0,0,j] for j in range(len(chrAncGenes))]
			for (s,l,c) in genesAncCol[i]:
				nb[c][1] = max(nb[c][1], s)
				nb[c][0] += 1
			c = max(nb)[2]
		
		chrAncGenes[c].append(i)


#
# Scanne les chromosomes Tetraodon et recherche des translocations des
# chromosomes ancestraux qu'on n'aurait pas decele lors de leur definition
#
def findNewChrAncTetra():

	print >> sys.stderr, "Etude de la nouvelle distribution des chromosomes ancestraux ... ",
	# Tout d'abord, on construit le profil de chaque chromosome Tetraodon
	profils = dict([(k,[]) for k in para.lstChr1])
	
	for i in range(len(chrAncGenes)):
		for x in chrAncGenes[i]:
			for (g, e) in genesAnc[x]:
				if e == options["especeDupliquee"]:
					profils[g[0]].append( (g,i) )

	# On prend chaque chromosome Tetraodon
	for j in para.lstChr1:
		
		profils[j].sort()
		
		# On extrait les groupes d'au moins 'nbGenesTransloc' genes
		p = []
		lastC = -1
		groupe = []
		for (g,c) in profils[j]:
			if c != lastC:
				if len(groupe) >= options["nbGenesTransloc"]:
					p.append( (lastC,groupe) )
				groupe = []
			lastC = c
			groupe.append(g)

		# On rajoute chacun de ces groupes s'il n'est pas present
		for (c,t) in p:
			(x1,x2) = (min([g[1] for g in t])-options["precisionChrAnc"], max([g[2] for g in t])+options["precisionChrAnc"])
			x1 -= options["precisionChrAnc"]
			x2 += options["precisionChrAnc"]
			flag = False
			for x in chrAncTetra[c]:
				if x[0] != j:
					continue
				if len(x) == 1:
					flag = True
				elif x[1] < x1 or x[2] > x2:
					flag = True
			if not flag:
				#sys.stderr.write("Ajout de %c sur %d (%d-%d)\n" % (chr(65+c), j, x1, x2))
				#chrAncTetra[c].add( (j,) )
				chrAncTetra[c].add( (j,x1,x2) )
			
	print >> sys.stderr, "OK"
				

#
# Cette fonction renvoie les chromosomes ancestraux tels qu'on peut les definir
# grace aux paralogues
#
def getChrAncIni(genomeDup, geneBank, para, orthos):

	#
	# Compte le nombre de passages d'un chromosome a l'autre
	#
	def countAltern3(genome, orthos, count):
	
		for c in genome.lstGenes:
			grp3Tet = []
			for (_,_,_,g) in genome.lstGenes[c]:
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
	def countPara(nbMin):
		
		count = {}
		for g1 in para:
			(c1,_) = genomeDup.dicGenes[g1]
			(c2,_) = genomeDup.dicGenes[para[g1]]
			if c1 < c2:
				count[ (c1,c2) ] = count.get( (c1,c2), 0) + 1
		r = []
		for (c,v) in count.items():
			if v >= nbMin:
				r.append( (v, c[0], c[1]) )
	
		r.sort()
		r.reverse()
	
		return r

	print >> sys.stderr, "Recherche des chromosomes ancestraux ... ",
	
	count3 = {}
	for e in geneBank.dicEspeces:
		countAltern3(geneBank.dicEspeces[e], orthos, count3)

	for c in count3:
		count3[c] /= len(geneBank.dicEspeces)
	
	# On va parcourir les couples de chromosomes dupliques des plus au
	# moins presents en construisant les chromosomes au fur et a mesure
	chrom = []
	nbPara = []
	rep = dict([(k,[]) for k in genomeDup.lstGenes])

	for (v, c1, c2) in countPara(options["nbMinParaCouple"]):
		if len(rep[c1]) == 0 and len(rep[c2]) == 0:
		
			# Les deux chromosomes sont inutilises: nouveau chromosome ancestral
			#print >> sys.stderr, "CHR", c1, c2
			rep[c1].append(len(chrom))
			rep[c2].append(len(chrom))
			chrom.append( set([c1,c2]) )
			nbPara.append( v )
			
		elif len(rep[c1]) == 0 or len(rep[c2]) == 0:
			
			# Juste un des deux est utilise, on le parcourt pour
			# voir si on peut creer un trio
			if len(rep[c2]) == 0:
				(x, y) = (c1, c2)
			else:
				(x, y) = (c2, c1)
				
			# On scanne les representations de x
			for i in rep[x]:
				
				# On prend le troisieme chromosome de l'alternance
				for c in chrom[i]:
					if c != c1 and c != c2:
						break
				t = [c,c1,c2]
				t.sort()
				#print >> sys.stderr, "test de %s: %d / %d" % (t, count3.get(tuple(t), 0), nbPara[i])
				
				# Si l'alternance est vue assez de fois, on rajoute le chromosome
				if count3.get(tuple(t), 0) > options["nbMinAltern3"]*nbPara[i]:
					t = chrom[i].union(set([y]))
					#print >> sys.stderr, "CHRN", t
					for j in t:
						rep[j].append(len(chrom))
					chrom.append( t )
					nbPara.append( nbPara[i] + v)
					break
			else:
				# Si aucune triple alternance trouvee, on cree
				# un nouveau chromosome ancestral
				#print >> sys.stderr, "CHRP", c1, c2
				rep[c1].append(len(chrom))
				rep[c2].append(len(chrom))
				chrom.append( set([c1,c2]) )
				nbPara.append( v )
		else:
			# Les deux sont deja pris, on va choisir celui dans
			# lequel on a la meilleure alternance de 3
			best = 1
			bestC = -1
			for i in rep[c1] + rep[c2]:
				for c in chrom[i]:
					if c == c1 or c == c2:
						continue
					t = [c,c1,c2]
					t.sort()
					x = count3.get(tuple(t), 0)
					if x < best:
						continue
					if x > best or len(chrom[i]) < len(chrom[bestC]):
						best = count3[tuple(t)]
						bestC = i
			if bestC != -1:
				#print >> sys.stderr, "FUS", c1, c2, chrom[bestC]
				chrom[bestC].add( c1 )
				chrom[bestC].add( c2 )
				rep[c1].append(bestC)
				rep[c2].append(bestC)
				nbPara[bestC] += v
	
	# Affichage des resultats
	print >> sys.stderr, "%d chromosomes trouves :" % len(chrom)
	for i in range(len(chrom)):
		print >> sys.stderr, "Chromosome %c (%d):" % (chr(65+i), nbPara[i]),
		for c in chrom[i]:
			print >> sys.stderr, " %d" % c,
		print >> sys.stderr
	
	return chrom



# MAIN #

# Arguments
(noms_fichiers, options) = myTools.checkArgs(["genes_list.conf", "GENES_ANC"], [("especeDupliquee", str, 'T'), ("bestScore", bool, False), ("precisionChrAnc", int, 1000000), ("nbGenesTransloc", int, 4), ("nbMinParaCouple", int, 4), ("nbMinAltern3", float, 1.), ("loopChrAnc", bool, False)] )

# Chargement des fichiers
geneBank = myOrthos.MyGeneBank(noms_fichiers[0])
genomeDup = geneBank.dicEspeces[options["especeDupliquee"]]
del geneBank.dicEspeces[options["especeDupliquee"]]
genesAnc = myOrthos.AncestralGenome(noms_fichiers[1], False)
lstGenesAnc = genesAnc.lstGenes[myOrthos.AncestralGenome.defaultChr]
(para,orthos) = buildParaOrtho(lstGenesAnc, genomeDup)

# On colorie les matrices actuelles
blocs = {}
for e in geneBank.dicEspeces:
	blocs[e] = colorAncestr(geneBank.dicEspeces[e], genomeDup, para, orthos)

# Les chromosomes ancestraux
chrAncTetra = getChrAncIni(genomeDup, geneBank, para, orthos)
fin = False
lastCol = []
nbR = 0

# On rentre dans la boucle
while not fin:
	nbR += 1
	print >> sys.stderr, "Passage", nbR
	
	# On colorie les genes ancestraux
	col = [[] for i in range(len(lstGenesAnc))]
	for e in blocs:
		buildColorTable(blocs[e], col, genesAnc.dicGenes, chrAncTetra)

	# On construit les chromosomes ancestraux
	chrAncGenes = [[] for x in chrAncTetra]
	buildChrAnc(col, chrAncGenes)
	
	if not options["loopChrAnc"]:
		break

	# A t'on le meme resultat qu'avant
	if col == lastCol:
		fin = True
	else:
		# Non, on continue
		findNewChrAncTetra()
		lastCol = col

# On affiche le resultat
printColorAncestr(lstGenesAnc, chrAncGenes)

