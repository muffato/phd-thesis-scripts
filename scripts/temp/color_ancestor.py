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

sys.path.append(os.environ['HOME'] + "/M2/scripts/utils")
import myOrthos
import myTools
import myMaths

# FONCTIONS #


#
# Chargement du fichier de configuration dees donnees
#  "EF = data/orthoEF" pour le fichier des orthologues entre E et F
# On a besoin de connaitre l'initiale de l'espece qui a subi la duplication
# Des lors, on separe:
#   - les paralogues de cette espece
#   - les orthologues avec cette espece (qui seront colories)
#   - les autres orthologues (pour avoir plus d'information)
#
def loadConfigDataFile(nom, t):

	tabOrthos = []
	tabOrthosT = []
	f = open(nom, 'r')
	
	for ligne in f:
		# Les noms des deux especes
		e1 = ligne[0]
		e2 = ligne[1]
		s = ligne[ligne.index('=')+1:ligne.index('\n')]
		ortho = myOrthos.Ensembl2SpeciesOrtho(s.strip())
		if e1 == t and e2 == t:
			para = ortho
		elif e1 == t:
			ortho.swapGene1Gene2()
			tabOrthosT.append( (ortho, e2, t) )
		elif e2 == t:
			tabOrthosT.append( (ortho, e1, t) )
		else:
			tabOrthos.append( (ortho, e1, e2) )

	f.close()
	return (para, tabOrthos, tabOrthosT)


#
# Creation des genes ancestraux, comme des listes de genes modernes
#
def buildAncGenes(para, tabOrthos, tabOrthosT):

	#
	# Construit les genes ancestraux, c'est a dire des listes de genes modernes
	# On rajoute chaque couple d'orthologues, en regroupant avec d'autres couples deja lus
	#
	def buildAncestrGenes(ortho, e1, e2):

		# Parcours des orthologues
		for (g1, g2, _, _) in ortho:
			
			if g1[3] in dicGenesAnc:
				i = dicGenesAnc[g1[3]]
				if g2[3] in dicGenesAnc:
				
					# Les deux genes sont deja dans la liste, on les relie
					j = dicGenesAnc[g2[3]]
					
					if i == j:
						continue

					genesAnc[i].extend(genesAnc[j])
					for (g,_) in genesAnc[j]:
						dicGenesAnc[g[3]] = i
					del genesAnc[j][:]
					
				else:
					# Un seul des deux, on rajoute a cette liste
					genesAnc[i].append( (g2,e2) )
					dicGenesAnc[g2[3]] = i
			
			elif g2[3] in dicGenesAnc:
				# Un seul des deux, on rajoute a cette liste
				j = dicGenesAnc[g2[3]]
				genesAnc[j].append( (g1,e1) )
				dicGenesAnc[g1[3]] = j
			
			else:
				# Aucun des deux, on cree un nouveau groupe
				dicGenesAnc[g1[3]] = len(genesAnc)
				dicGenesAnc[g2[3]] = len(genesAnc)
				genesAnc.append([(g1,e1), (g2,e2)])



	# Tout d'abord, on cree les genes ancestraux, comme des listes de genes modernes
	print >> sys.stderr, "Construction des genes ancestraux ",
	genesAnc = []
	dicGenesAnc = {}
	
	# Insertion de chaque fichier d'orthologues
	for (ortho, e1, e2) in tabOrthosT:
		buildAncestrGenes(ortho, e1, e2)
		sys.stderr.write(".")
	for (ortho, e1, e2) in tabOrthos:
		buildAncestrGenes(ortho, e1, e2)
		sys.stderr.write(".")
	buildAncestrGenes(para, options["especeDupliquee"], options["especeDupliquee"])
	sys.stderr.write(".")
	
	esp = [set([e for (g,e) in gene]) for gene in genesAnc]
	genesAnc = [genesAnc[i] for i in range(len(genesAnc)) if options["especeDupliquee"] in esp[i] and len(esp[i]) > 1]
	
	# On associe noms et genes
	for i in range(len(genesAnc)):
		genesAnc[i].sort()
		for (g,e) in genesAnc[i]:
			dicGenesAnc[g[3]] = i
	
	print >> sys.stderr, "", len(genesAnc), "dans le genome ancestral"

	return (genesAnc, dicGenesAnc)

#
# Recherche des nouveaux orthologues qu'on pourrait utiliser
#
def searchNewOrtho(ortho, e1, e2, genesAnc):

	print >> sys.stderr, "Recherche de nouveaux orthologues dans", ortho.nom, "",
	nb = 0

	for t in genesAnc:

		t1 = [x for (x,e) in t if e == e1 and x[3] not in ortho.dictGenes1]
		t2 = [x for (x,e) in t if e == e2 and x[3] not in ortho.dictGenes2]
	
		while len(t1) >= 1 and len(t2) >= 1:
			ortho.addOrtho(t1.pop(), t2.pop())
			sys.stderr.write(".")
			nb += 1
	
	print >> sys.stderr, "", nb, "ajouts"

#
# Recherche des nouveaux paralogues qu'on pourrait utiliser
#
def searchNewPara(ortho, e, genesAnc):

	print >> sys.stderr, "Recherche de nouveaux paralogues dans", ortho.nom,
	nb = 0

	for t in genesAnc:

		tab = [x for (x,f) in t if e == f and x[3] not in ortho.dictGenes1]

		while len(tab) >= 2:
			g1 = tab.pop()
			for i in range(len(tab)):
				g2 = tab[i]
				if g1[0] == g2[0]:
					continue
				ortho.addOrtho(g1, g2)
				ortho.addOrtho(g2, g1)
				sys.stderr.write(".")
				nb += 1
				tab.pop(i)
				break
	
	print >> sys.stderr, "", nb, "ajouts"


#
# Cree les regions de chromosomes ancestraux dans la classe d'orthologues
#
def colorAncestr(ortho, para):

	print >> sys.stderr, "Decoupage de", ortho.nom, "...",

	# La premiere etape est de faire une liste de chromosome Tetraodon pour
	# chaque orthologue
	ortho.orig = [ set([]) for i in range(ortho.nbGenes) ]

	# On parcourt les orthologues
	for (_, gT, _, i) in ortho:
		
		# On construit une liste avec le numero de chromosome Tetraodon
		ortho.orig[i].add(gT[0])

		for g in para.getGenes2(para.getGene1Filtered(gT[0], gT[1]-options["precisionChrAnc"], gT[2]+options["precisionChrAnc"])):
			# On rajoute les numeros de chromosome trouves
			ortho.orig[i].add(g[0])
	

	# On rassemble en regions lorsqu'il y a une origine commune
	bloc = []
	ortho.blocs = []
	lastOrig = ortho.orig[ortho.tabGenes1[0][4]]
	lastChrE = 1
	
	# Parcours des genes poulets
	for i in range(ortho.nbGenes):

		ind = ortho.tabGenes1[i][4]
		gT = ortho.tabGenes2[ortho.indGenes2[ind]]
		
		# On intersecte les origines possibles des genes precedents et du courant
		ori = ortho.orig[ind]
		nouvOrig = lastOrig & ori
		
		# On arrete le bloc si:
		#   - aucune origine en commun
		#   - on change de chromosome sur l'espece E
		if ortho.tabGenes1[i][0] != lastChrE or len(nouvOrig) == 0:
			ortho.blocs.append( (lastChrE, bloc) )
			bloc = []
			nouvOrig = ori
		
		bloc.append( ind )
		lastOrig = nouvOrig
		lastChrE = ortho.tabGenes1[i][0]

	# Le dernier bloc
	ortho.blocs.append( (lastChrE, bloc) )

	print >> sys.stderr, len(ortho.blocs), "blocs"

#
# Renvoie le chromosome ancestral qui a le meilleur score pour la region consideree
#
def getMaxScore(ortho, bloc, chrAnc):
	score = [0 for i in range(len(chrAnc))]
	for gi in bloc:
		g = ortho.tabGenes2[ortho.indGenes2[gi]]
		for i in range(len(chrAnc)):
			for c in chrAnc[i]:
				if c[0] == g[0]:
					if len(c) == 1:
						score[i] += 1
						break
					elif g[1] <= c[2] and g[2] >= c[1]:
						score[i] += 1
						break
	m = max(score)
	return (m, score.index(m))


#
# Construit les tables d'association "nom de gene" -> couleur
# Prend la couleur la plus probable
#
def buildColorTable(ortho, col, dicGenesAnc, chrAnc):


	print >> sys.stderr, "Assignation des chromosomes ancestraux sur", ortho.nom, "... ",
	
	ortho.tabAncestr = range(ortho.nbGenes)
	for i in range(len(ortho.blocs)):
		(_, b) = ortho.blocs[i]
		(s, c) = getMaxScore(ortho, b, chrAnc)
		for j in b:
			ortho.tabAncestr[j] = ((float(s-1)/float(len(b)), len(b)), c)
	
	for i in range(ortho.nbGenes):

		g1 = ortho.tabGenes1[ortho.indGenes1[i]]
		g2 = ortho.tabGenes2[ortho.indGenes2[i]]

		j = dicGenesAnc[g1[3]]
		col[j].append(ortho.tabAncestr[i])
	
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
				r += " " + s[0][3]
			print r
		nb += len(chrAncGenes[c])
		
	print >> sys.stderr, nb, "genes dans le genome ancestral"


#
# Range chaque gene ancestral dans son chromosome, et trie les chromosomes
#
def buildChrAnc(genesAncCol, chrAncGenes):

	for i in range(len(genesAncCol)):
	
		if len(genesAncCol[i]) == 0:
			# Ce sont des groupes de genes uniquement Tetraodon
			#  ou sans genes Tetraodon
			continue
		
		if options["bestScore"]:
			# Le meilleur score de purete
			c = max(genesAncCol[i])[1]
		else:
			# Le chromosome le plus frequent
			nb = [[0,0,j] for j in range(len(chrAncGenes))]
			for ((s,l),c) in genesAncCol[i]:
				nb[c][1] += s
				nb[c][0] += 1
			for j in range(len(chrAncGenes)):
				if nb[j][0] != 0:
					nb[j][1] /= nb[j][0]
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
def getChrAncIni():

	#
	# Compte le nombre de passages d'un chromosome a l'autre
	#
	def countAltern3(ortho, count):
		
		# Contient le liste des trois derniers chromosomes Tetraodon vus
		grp3Tet = []
		lastK = 0
		for g in ortho.tabGenes1:
			
			# On verifie qu'on n'a pas change de chromosome dans l'espece
			if g[0] != lastK:
				lastK = g[0]
				grp3Tet = []
			
			KT = ortho.tabGenes2[ortho.indGenes2[g[4]]][0]
			
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
		for i in range(para.nbGenes):
			g1 = para.tabGenes1[i]
			g2 = para.tabGenes2[para.indGenes2[g1[4]]]
			if g1[0] < g2[0]:
				count[ (g1[0],g2[0]) ] = count.get( (g1[0],g2[0]), 0) + 1
		r = []
		for (c,v) in count.items():
			if v >= nbMin:
				r.append( (v, c[0], c[1]) )
	
		r.sort()
		r.reverse()
	
		return r

	print >> sys.stderr, "Recherche des chromosomes ancestraux ... ",
	
	count3 = {}
	for (ortho, e1, e2) in tabOrthosT:
		countAltern3(ortho, count3)

	for c in count3:
		count3[c] /= len(tabOrthosT)
	
	# On va parcourir les couples de chromosomes dupliques des plus au
	# moins presents en construisant les chromosomes au fur et a mesure
	chrom = []
	nbPara = []
	rep = dict([(k,[]) for k in para.lstChr1])

	for (v, c1, c2) in countPara(options["nbMinParaCouple"]):
		if len(rep[c1]) == 0 and len(rep[c2]) == 0:
		
			# Les deux chromosomes sont inutilises: nouveau chromosome ancestral
			#print >> sys.stderr, "CHR", c1, c2
			rep[c1].append(len(chrom))
			rep[c2].append(len(chrom))
			chrom.append( set([ (c1,), (c2,)]) )
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
				for (c,) in chrom[i]:
					if c != c1 and c != c2:
						break
				t = [c,c1,c2]
				t.sort()
				# print >> sys.stderr, "test de %s: %d" % (t, count3.get(tuple(t), 0))
				
				# Si l'alternance est vue assez de fois, on rajoute le chromosome
				if count3.get(tuple(t), 0) > options["nbMinAltern3"]*nbPara[i]:
					t = chrom[i].union(set([(y,)]))
					#print >> sys.stderr, "CHRN", t
					for (j,) in t:
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
				chrom.append( set([ (c1,), (c2,)]) )
				nbPara.append( v )
		else:
			# Les deux sont deja pris, on va choisir celui dans
			# lequel on a la meilleure alternance de 3
			best = 1
			bestC = -1
			for i in rep[c1] + rep[c2]:
				for (c,) in chrom[i]:
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
				chrom[bestC].add( (c1,) )
				chrom[bestC].add( (c2,) )
				rep[c1].append(bestC)
				rep[c2].append(bestC)
				nbPara[bestC] += v
	
	# Affichage des resultats
	print >> sys.stderr, "%d chromosomes trouves :" % len(chrom)
	for i in range(len(chrom)):
		print >> sys.stderr, "Chromosome %c (%d):" % (chr(65+i), nbPara[i]),
		for (c,) in chrom[i]:
			print >> sys.stderr, " %d" % c,
		print >> sys.stderr
	
	return chrom



# MAIN #

# Arguments
(noms_fichiers, options) = myTools.checkArgs(["color_ancestor.data.conf"], [("especeDupliquee", str, 'T'), ("bestScore", bool, False), ("precisionChrAnc", int, 1000000), ("nbGenesTransloc", int, 4), ("nbMinParaCouple", int, 5), ("nbMinAltern3", float, 1.), ("loopChrAnc", bool, False)] )

# Chargement du fichier de configuration
# On charge ici tous les fichiers d'orthologues
(para, tabOrthos, tabOrthosT) = loadConfigDataFile(noms_fichiers[0], options["especeDupliquee"])

# Tout d'abord, on cree les genes ancestraux, comme des listes de genes modernes
(genesAnc, dicGenesAnc) = buildAncGenes(para, tabOrthos, tabOrthosT)

# Pour etablir de nouveaux orthologues, il suffit de parcourir cette liste
for (ortho, e1, e2) in tabOrthosT:
	searchNewOrtho(ortho, e1, e2, genesAnc)
searchNewPara(para, options["especeDupliquee"], genesAnc)

# On colorie les matrices actuelles
for (ortho, e1, e2) in tabOrthosT:
	colorAncestr(ortho, para)

# Les chromosomes ancestraux
chrAncTetra = getChrAncIni()
fin = False
lastCol = []
nbR = 0

# On rentre dans la boucle
while not fin:
	nbR += 1
	print >> sys.stderr, "Passage", nbR
	
	# On colorie les genes ancestraux
	col = [[] for i in range(len(genesAnc))]
	for (ortho, e1, e2) in tabOrthosT:
		buildColorTable(ortho, col, dicGenesAnc, chrAncTetra)

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
printColorAncestr(genesAnc, chrAncGenes)

