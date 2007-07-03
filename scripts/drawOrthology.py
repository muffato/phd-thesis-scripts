#! /users/ldog/muffato/python -OO

__doc__ = """
	Compare deux genomes:
		- Dessine la matrice des genes orthologues entre deux genomes.
		- Dessine le karyotype d'un genome face a l'autre
		- Renvoie les couples de chromosomes orthologues avec des stats sur l'evolution des genomes (nb genes, rearrangements)
		- Renvoie la liste des couples de genes orthologues avec les details sur leurs positions
		- Reordonne le genome 1 pour qu'il soit plus ressemblant au genome 2
		- Renvoie les diagonales entre les deux genomes
"""

##################
# INITIALISATION #
##################

# Librairies
import sys
import random
import utils.myGenomes
import utils.myTools
import utils.myDiags
import utils.myPsOutput


#############
# FONCTIONS #
#############

# Fabrique la liste des orthologues entre les deux genomes
def buildOrthosTable(genome1, chr1, genome2, chr2):

	# Tous les orthologues entre pour les chromosomes OK du genome 1
	res = {}
	for c1 in chr1:
		res[c1] = {}
		for g1 in genome1.lstGenes[c1]:
			tmp = set()
			for s in g1.names:
				if s in genome2.dicGenes:
					tmp.add(genome2.dicGenes[s])
				if options["orthologuesList"] == "":
					continue
				if s in genesAnc.dicGenes:
					(c,i) = genesAnc.dicGenes[s]
					for ss in genesAnc.lstGenes[c][i].names:
						if ss in genome2.dicGenes:
							tmp.add(genome2.dicGenes[ss])
			# On ne garde que les chromsomes OK du genome 2
			tmp = [(c,i) for (c,i) in tmp if c in chr2]
			# +/- includeGaps
			if (not options["includeGaps"]) and (len(tmp) == 0):
				continue
			res[c1][genome1.dicGenes[s][1]] = [(c,i) for (c,i) in tmp if c in chr2]
			

	return res


# Renvoie la liste des chromosomes orthologues
def getOrthosChr(table, chr):
	
	res = {}
	for c1 in chr:

		# On ne garde que les chromosomes orthologues majoritaires
		lst = [x for (x,_) in utils.myMaths.flatten(table[c1].itervalues())]
		count = [(lst.count(x),x) for x in set(lst)]
		count.sort()
		nb = (len(lst)*options["minHomology"])/100
		tmp = []
		for (n,c2) in count.__reversed__():
			tmp.append( (c2,n) )
			nb -= n
			if nb <= 0:
				break
		res[c1] = tmp
	return res




########
# MAIN #
########

# Arguments
modeMatrix = "Matrix"
modeKaryo = "Karyotype"
modeGenomeEvol = "GenomeEvolution"
modeOrthos = "OrthosGenes"
modeReindexedChr = "ReindexedChr"
modeDiags = "Diags"
modes = [modeMatrix, modeKaryo, modeOrthos, modeGenomeEvol, modeReindexedChr, modeDiags]
(noms_fichiers, options) = utils.myTools.checkArgs( \
	["studiedGenome", "referenceGenome"], \
	[("orthologuesList",str,""), ("includeGaps",bool,False), ("includeScaffolds",bool,False), ("includeRandoms",bool,False), \
	("output",str,modes), ("reverse",bool,False), ("scaleY",bool,False), ("minHomology",int,90), \
	("pointSize",float,-1), ("colorFile",str,""), ("defaultColor",str,"black"), ("penColor",str,"black"), ("backgroundColor",str,"")], \
	__doc__
)


# Chargement des fichiers
genome1 = utils.myGenomes.loadGenome(noms_fichiers["studiedGenome"])
genome2 = utils.myGenomes.loadGenome(noms_fichiers["referenceGenome"])
if options["reverse"]:
	x = genome1
	genome1 = genome2
	genome2 = x
if options["orthologuesList"] != "":
	genesAnc = utils.myGenomes.loadGenome(options["orthologuesList"])
else:
	genesAnc = genome2
try:
	colors = utils.myGenomes.loadGenome(options["colorFile"])
except Exception:
	colors = options["defaultColor"]

# Les chromosomes a etudier
chr1 = genome1.lstChr
chr2 = genome2.lstChr
if options["includeScaffolds"]:
	chr1.extend(genome1.lstScaff)
	chr2.extend(genome2.lstScaff)
if options["includeRandoms"]:
	chr1.extend(genome1.lstRand)
	chr2.extend(genome2.lstRand)


if (options["output"] != modeDiags):
	table12 = buildOrthosTable(genome1, chr1, genome2, chr2)

# Matrice
if (options["output"] == modeMatrix):

	print >> sys.stderr, "Affichage ",
	
	utils.myPsOutput.printPsHeader()
	if len(options["backgroundColor"]) > 0:
		utils.myPsOutput.drawBox(0,0, 21,29.7, options["backgroundColor"], options["backgroundColor"])
	sys.stderr.write('.')

	# Initialisations
	table21 = buildOrthosTable(genome2, chr2, genome1, chr1)
	nb = sum([len(table12[c]) for c in table12])
	scaleX = 19. / float(nb)
	if options["scaleY"]:
		scaleY = 19. / float(sum([len(table21[c]) for c in table21]))
	else:
		scaleY = scaleX
	if options["pointSize"] < 0:
		dp = scaleX
	else:
		dp = options["pointSize"]
	sys.stderr.write('.')

	def prepareGenome(dicOrthos, func):
		i = 0
		y = 0
		lstNum = {}
		func(y)
		for c in sorted(dicOrthos):
			y += len(dicOrthos[c])
			func(y)
			for gene in dicOrthos[c]:
				lstNum[(c,gene)] = i
				i += 1
		return lstNum

	lstNum1 = prepareGenome(table12, lambda x: utils.myPsOutput.drawLine(1 + x*scaleX, 1, 0, float(sum([len(table21[c]) for c in table21]))*scaleY, options["penColor"]))
	sys.stderr.write('.')
	lstNum2 = prepareGenome(table21, lambda y: utils.myPsOutput.drawLine(1, 1 + y*scaleY, 19, 0, options["penColor"]))
	sys.stderr.write('.')

	print "0 setlinewidth"

	for c1 in table12:
		for i1 in table12[c1]:
			xx = 1 + float(lstNum1[(c1,i1)]) * scaleX
			for (c2,i2) in table12[c1][i1]:

				if type(colors) == str:
					coul = colors
				else:
					tmp = colors.getPosition(genome1.lstGenes[c1][i1]) + colors.getPosition(genome2.lstGenes[c2][i2])
					if len(tmp) == 0:
						continue
					coul = tmp[0][0]
				
				yy = 1 + float(lstNum2[(c2,i2)]) * scaleY
				#coul = utils.myPsOutput.getColor(coul, options["defaultColor"])
				utils.myPsOutput.drawBox( xx, yy, dp, dp, coul, coul)

	utils.myPsOutput.printPsFooter()
	print >> sys.stderr, " OK"


# Caryotype
elif (options["output"] == modeKaryo):

	print >> sys.stderr, "Affichage ...",

	utils.myPsOutput.printPsHeader()

	# On dessine
	dx = (19.*3.)/(5.*len(chr1)-2.)
	dy = float(max([len(x) for x in table12.values()])) / 26.
	xx = 1
	y0 = 1.

	for c in chr1:
		
		if options["scaleY"]:
			dy = float(len(table12[c])) / 26.

		utils.myPsOutput.drawText(xx, y0, str(c), options["penColor"])
		y = y0 + 1
		
		last = ""
		nb = 0
		for i in table12[c]:
			tmp = table12[c][i]
			if len(tmp) == 0:
				col = options["defaultColor"]
			else:
				col = tmp[0][0]
			if col == last:
				nb += 1
			else:
				if nb != 0:
					utils.myPsOutput.drawBox(xx, y, dx, nb/dy, last, last)
					y += nb/dy
				last = col
				nb = 1
		if nb != 0:
			utils.myPsOutput.drawBox(xx, y, dx, nb/dy, last, last)
		xx += (5.*dx)/3.

	utils.myPsOutput.printPsFooter()
	print >> sys.stderr, "OK"


# Fichier avec les noms des paires de genes orthologues et leurs coordonnees
elif (options["output"] == modeOrthos):

	for c1 in chr1:
		for i1 in table12[c1]:
			g1 = genome1.lstGenes[c1][i1]
			for (c2,i2) in table12[c1][i1]:
				g2 = genome2.lstGenes[c2][i2]
				r = (c1, g1.beginning, g1.end, g1.strand, "/".join(g1.names), \
					c2, g2.beginning, g2.end, g2.strand, "/".join(g2.names))
				print "\t".join([str(x) for x in r])


# Liste de chromosomes orthologues
elif (options["output"] == modeGenomeEvol):

	# On a besoin des equivalences d'un genome a l'autre
	res1 = getOrthosChr(table12, chr1)
	table21 = buildOrthosTable(genome2, chr2, genome1, chr1)
	res2 = getOrthosChr(table21, chr2)
	
	# Et d'en tirer les best-hits reciproques chromosomiques
	besthits1 = {}
	besthits2 = {}
	for c1 in chr1:
		for (c2,_) in res1[c1]:
			if len([(c,x) for (c,x) in res2[c2] if c == c1]) > 0:
				besthits1[c1] = besthits1.get(c1,[]) + [c2]
				besthits2[c2] = besthits2.get(c2,[]) + [c1]
	
	# Initialisation
	cassures = 0
	nbChr1 = 0
	pertes = 0
	nb1 = 0
	dupliques = {}
	tmp2 = {}
	
	# Parcours du genome 1
	for c1 in chr1:
		
		# Les calculs sur les nombres de genes
		nb1 += len(genome1.lstGenes[c1])
		pertes += len(genome1.lstGenes[c1]) - len(table12[c1])
		for (g1,g2) in table12[c1].iteritems():
			# Les duplications du genome 2 vers le genome 1
			for x in g2:
				tmp2[x] = tmp2.get(x,0) + 1
			# Des pertes du genome 1 vers le genome 2
			if len(g2) == 0:
				pertes += 1
			# Les duplications du genome 1 vers le genome 2
			elif len(g2) > 1:
				dupliques[len(g2)] = dupliques.get(len(g2),0) + 1
		
		# Les calculs sur les rearrangements
		if len(besthits1.get(c1,[])) == 0:
			continue
		nbChr1 += 1
		cassures += len(besthits1[c1]) - 1
			
		
	# Parcours du genome 2
	nb2 = 0
	nbChr2 = 0
	fusions = 0
	#tt = 0
	for c2 in chr2:
		# On met a jour le nombre de genes du genome 1
		nb2 += len(genome2.lstGenes[c2])
		
		# Les calculs sur les rearrangements
		if len(besthits2.get(c2,[])) == 0:
			continue
			
		nbChr2 += 1
		fusions += len(besthits2[c2]) - 1
	
	nouveaux = nb2 - len(tmp2)
	
	# Les duplications dans l'autre sens
	dupliques2 = {}
	for x in tmp2.itervalues():
		dupliques2[x] = dupliques2.get(x,0) + 1
	dupliques2.pop(1, None)
	
	print >> sys.stderr, "Evolution de A (%d chr / %d genes) vers B (%d chr / %d genes)" % (nbChr1,nb1, nbChr2,nb2)
	print >> sys.stderr, "Nb nouveaux genes", nouveaux
	print >> sys.stderr, "Nb duplications", sum(dupliques.values()), dupliques
	print >> sys.stderr, "Nb duplications (inverse)", sum(dupliques2.values()), dupliques2
	print >> sys.stderr, "Nb genes perdus", pertes
	print >> sys.stderr, "Nb points de cassures", cassures
	print >> sys.stderr, "Nb points de fusions", fusions

	for c1 in chr1:
		print "%s\t%s" % (c1, "\t".join(["%s (%d)" % (c2,n) for (c2,n) in res1[c1]]))

# On echange l'ordre des chromosomes pour que les deux genomes paraissent plus colineaires
elif (options["output"] == modeReindexedChr):

	# D'abord on fait la liste des paires de chromosomes homologues
	res1 = getOrthosChr(table12, chr1)
	besthits = []
	for c1 in chr1:
		if len(res1[c1]) == 0:
			besthits.append( (100000,0,c1) )
		else:
			(c2,nb) = res[c1][0]
			besthits.append( (c2,-nb,c1) )

	# On renvoie les chromosomes du genome 1 dans l'ordre des best hits avec le genome 2
	besthits.sort()
	#for i in xrange(len(besthits)):
	#	(c2,nb,c1) = besthits[i]
	for (i,(c2,nb,c1)) in enumerate(besthits):
		# Faut-il retourner le chromosome ?
		memeSens = 0
		# On restreint c1 a ses orthologues avec c2
		tmp = []
		for i1 in xrange(len(genome1.lstGenes[c1])):
			tg2 = [i2 for (cc2,i2) in table12[c1].get(i1,[]) if cc2 == c2]
			if len(tg2) > 0:
				tmp.append(tg2)
		# Le test de colinearite
		for j in xrange(len(tmp)-1):
			if max(tmp[j]) < min(tmp[j+1]):
				memeSens += 1
			elif min(tmp[j]) > max(tmp[j+1]):
				memeSens -= 1
		# Le resultat
		if memeSens > 0:
			res = genome1.lstGenes[c1].__iter__()
		else:
			res = genome1.lstGenes[c1].__reversed__()
		for g in res:
			print i+1, " ".join(g.names)


# Les diagonales entre les deux genomes
elif (options["output"] == modeDiags):
	def printDiag( ((e1,c1,d1), (e2,c2,d2), _) ):
		print '\t'.join([str(len(d1)), str(c1)," ".join([genome1.lstGenes[c1][i1].names[0] for i1 in d1]), str(c2)," ".join([genome2.lstGenes[c2][i2].names[0] for i2 in d2])])
	genome1.lstChr = chr1
	genome1.lstScaff = []
	genome2.lstChr = chr2
	genome2.lstScaff = []

	utils.myDiags.calcDiags(genome1.nom, genome2.nom, genome1, genome2, genesAnc, printDiag, \
		options["minHomology"], -1, options["scaleY"], not options["includeGaps"])


