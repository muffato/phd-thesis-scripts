#! /usr/bin/python2.4

#
# Fonctions communes de traitement des diagonales
#

import sys
import myGenomes
import myMaths


#
# Fait la liste de tous les genes de tab1 en diagonale avec ceux de tab2
#  (eventuellement des singletons)
# Pour optimiser, on demande 
#   tab1 qui est la liste des numeros des genes ancestraux du chromosome 1 (cf translateGenome)
#   dic2 qui associe a un numero de gene ancestral ses positions sur le genome 2
#

def extractDiags(tab1, dic2, largeurTrou):
	
	diag = []
	lastI2 = []
	lastC2 = []
	listI1 = []
	listI2 = []
	
	# Parcours du genome 1
	for i1 in xrange(len(tab1)):
		if tab1[i1] < 0:
			presI2 = []
		else:
			presI2 = dic2[tab1[i1]]

		# On regarde chaque orthologue du gene
		for (c2,i2) in presI2:
			# Est-ce qu'on est dans le prolongement d'une diagonale
			if c2 in lastC2 and (((i2+1) in lastI2) or ((i2-1) in lastI2)):
				listI2.append(i2)
				lastI2 = [i2]
				lastC2 = [c2]
				break
		else:
			# On n'a pas trouve de i2 satisfaisant, c'est la fin de la diagonale
			if len(listI1) > 0 and len(listI2) > 0:
				diag.append( (listI1,listI2,lastC2[0],(deb1,fin1),myMaths.getMinMax(listI2)) )
			deb1 = i1
			listI1 = []
			lastI2 = [i for (_,i) in presI2]
			listI2 = [i for (_,i) in presI2]
			lastC2 = [c for (c,_) in presI2]
		listI1.append(i1)
		fin1 = i1
	
	if len(listI1) > 0 and len(listI2) > 0:
		diag.append( (listI1,listI2,lastC2[0],(deb1,fin1),myMaths.getMinMax(listI2)) )
	elif len(diag) == 0:
		return []
	
	#for d in diag:
	#	print d
	#sys.exit(0)

	# On rassemble des diagonales separees par une espace pas trop large
	combiDiag = []
	while len(diag) > 0:
		#print "on est sur", diag[0]
		(d1,d2,c2,(deb1,fin1),(deb2,fin2)) = diag.pop(0)
		i = 0
		while i < len(diag):
			(dd1,dd2,cc2,(debb1,finn1),(debb2,finn2)) = diag[i]
			#print "test de", diag[i]
			# Aucune chance de poursuivre la diagonale
			if debb1 > (fin1+largeurTrou+1):
				#print "meme pas la peine"
				break
			a = min(abs(deb2-finn2), abs(debb2-fin2))
			if a <= (largeurTrou+1) and c2==cc2:
				#print "OK"
				d1.extend(dd1)
				d2.extend(dd2)
				fin1 = finn1
				deb2 = min(deb2,finn2)
				fin2 = max(deb2,finn2)
				del diag[i]
			else:
				#print "non"
				i += 1
		combiDiag.append( (d1,d2,c2,(deb1,fin1),(deb2,fin2)) )
	
	return combiDiag


#
# Convertit un genome en liste de genes ancestraux 
#
def translateGenome(genome, genesAnc):
	newGenome = {}
	for c in genome.lstChr:
		genes = genome.lstGenes[c]
		newGenome[c] = [genesAnc.dicGenes.get(g.names[0], (0,-1))[1] for g in genes]

	return newGenome

#
# Construit la liste des positions des genes ancestraux dans les especes de la banque
#
def buildAncGenesLocations(geneBank, genesAnc):
	lstGenesAnc = genesAnc.lstGenes[myGenomes.AncestralGenome.defaultChr]
	locations = dict( [(e,[[] for x in lstGenesAnc]) for e in geneBank.dicEspeces] )
	for ianc in xrange(len(lstGenesAnc)):
		for g in lstGenesAnc[ianc].names:
			if g not in geneBank.dicGenes:
				continue
			(e,c,i) = geneBank.dicGenes[g]
			locations[e][ianc].append( (c,i) )
	
	return locations


def iterateDiags(genome1, dic2, threshold, callBackFunc):

	for c1 in genome1:
	
		diags = extractDiags(genome1[c1], dic2, threshold)
		for (d1,d2,c2,_,_) in diags:
			callBackFunc(c1,c2,d1,d2)


