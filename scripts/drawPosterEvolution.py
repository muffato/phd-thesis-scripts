#! /users/ldog/muffato/python -OO

__doc__ = """
	Compare deux genomes:
		- Dessine la matrice des genes orthologues entre deux genomes.
		- Dessine le karyotype d'un genome face a l'autre
		- Renvoie les couples de chromosomes orthologues
"""

##################
# INITIALISATION #
##################

# Librairies
import sys
import utils.myGenomes
import utils.myTools
import utils.myPsOutput
import utils.myPhylTree


#############
# FONCTIONS #
#############

# Fabrique la liste des orthologues utilises pour le dessin
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
				continue
				if s in genesAnc.dicGenes:
					(c,i) = genesAnc.dicGenes[s]
					for ss in genesAnc.lstGenes[c][i].names:
						if ss in genome2.dicGenes:
							tmp.add(genome2.dicGenes[ss])
			# On ne garde que les chromsomes OK du genome 2
			tmp = [(c,i) for (c,i) in tmp if c in chr2]
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


# On echange l'ordre des chromosomes pour que les deux genomes paraissent plus colineaires
def reorderGenome(genome1, genome2):

	# D'abord on fait la liste des paires de chromosomes homologues
	table12 = buildOrthosTable(genome1, genome1.lstChr, genome2, genome2.lstChr)
	res1 = getOrthosChr(table12, genome1.lstChr)
	besthits = []
	for c1 in genome1.lstChr:
		if len(res1[c1]) == 0:
			besthits.append( (100000,0,c1) )
		else:
			(c2,nb) = res1[c1][0]
			besthits.append( (c2,-nb,c1) )

	newGenome = utils.myGenomes.Genome(genome2)
	
	# On renvoie les chromosomes du genome 1 dans l'ordre des best hits avec le genome 2
	besthits.sort()
	for (i,(c2,nb,c1)) in enumerate(besthits):
	#for i in xrange(len(besthits)):
	#	(c2,nb,c1) = besthits[i]
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
		if memeSens < 0:
			genome1.lstGenes[c1].reverse()
		for j in xrange(len(genome1.lstGenes[c1])):
			g = genome1.lstGenes[c1][j]
			newGenome.addGene( utils.myPhylTree.Gene(g.names, i+1, j, j, 0) )
	newGenome.lstChr = range(1,len(besthits)+1)
	newGenome.sortGenome()
	return newGenome

########
# MAIN #
########


# Arguments
(noms_fichiers, options) = utils.myTools.checkArgs( ["phylTree.conf"], [("penColor",str,"black"), ("minHomology",int,90), ("ancGenesFile",str,"")], __doc__)


phylTree = utils.myPhylTree.PhylogeneticTree(noms_fichiers["phylTree.conf"])

order = ["Euteleostomi", "Tetrapoda", "Amniota", "Mammalia", "Theria", "Boreoeutheria", "Euarchontoglires", "Catarrhini", "Homo/Pan/Gorilla group"]
#order = order[:2]

genomesAnc = []
for anc in order:
	gen = utils.myGenomes.AncestralGenome( options["ancGenesFile"] % phylTree.fileName[anc], chromPresents=True)
	if len(genomesAnc) > 0:
		#gen = reorderGenome(genomesAnc[-1], gen)
		gen = reorderGenome(gen, genomesAnc[-1])
	genomesAnc.append(gen)

utils.myPsOutput.printPsHeader(0.001)
print "1 cm 1 cm translate"
print "0.1 0.1 scale"

for (ia1,ia2) in utils.myTools.myIterator.tupleOnWholeList(range(len(order))):

	print >> sys.stderr, "%s _/\_ %s " % (order[ia1], order[ia2]),
	print "%f cm %f cm translate" % (20*ia1, 20*ia2)

	genome1 = genomesAnc[ia1]
	genome2 = genomesAnc[ia2]

	table12 = buildOrthosTable(genome1, genome1.lstChr, genome2, genome2.lstChr)
	sys.stderr.write('.')

	# Initialisations
	table21 = buildOrthosTable(genome2, genome2.lstChr, genome1, genome1.lstChr)
	nb = sum([len(table12[c]) for c in table12])
	scaleX = 19. / float(nb)
	scaleY = 19. / float(sum([len(table21[c]) for c in table21]))
	dp = scaleX
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
		#func(y)
		return lstNum

	lstNum1 = prepareGenome(table12, lambda x: utils.myPsOutput.drawLine(1 + x*scaleX, 1, 0, float(sum([len(table21[c]) for c in table21]))*scaleY, options["penColor"]))
	sys.stderr.write('.')
	lstNum2 = prepareGenome(table21, lambda y: utils.myPsOutput.drawLine(1, 1 + y*scaleY, 19, 0, options["penColor"]))
	sys.stderr.write('.')

	for c1 in table12:
		for i1 in table12[c1]:
			
			xx = 1 + float(lstNum1[(c1,i1)]) * scaleX
			
			for (c2,i2) in table12[c1][i1]:

				yy = 1 + float(lstNum2[(c2,i2)]) * scaleY
				utils.myPsOutput.drawLine( xx, yy, dp, dp, None)

	sys.stderr.write('+\n')

	print "%f cm %f cm translate" % (-20*ia1, -20*ia2)

utils.myPsOutput.printPsFooter()


