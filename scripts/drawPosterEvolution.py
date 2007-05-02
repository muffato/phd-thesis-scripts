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


########
# MAIN #
########


# Arguments
(noms_fichiers, options) = utils.myTools.checkArgs( \
	["phylTree.conf"], \
	[("pointSize",float,-1), ("colorFile",str,""), ("defaultColor",str,"black"), ("penColor",str,"black"), ("minHomology",int,90), \
	("ancGenesFile",str,"")], \
	__doc__
)


phylTree = utils.myBioObjects.PhylogeneticTree(noms_fichiers["phylTree.conf"])

#order = ["Euteleostomi", "Tetrapoda", "Amniota", "Mammalia", "Theria", "Boreoeutheria", "Euarchontoglires", "Catarrhini", "Homo/Pan/Gorilla group"]
#order = ["Euteleostomi", "Tetrapoda", "Amniota"]
order = ["Euteleostomi", "Homo/Pan/Gorilla group"]

genomesAnc = []
for anc in order:
	genomesAnc.append( utils.myGenomes.AncestralGenome( options["ancGenesFile"] % phylTree.fileName[anc], chromPresents=True) )

utils.myPsOutput.printPsHeader(0.001)
print "1 cm 1 cm translate"
print "0.1 0.1 scale"

for ia1 in xrange(len(order)):
	for ia2 in xrange(len(order)):

		print >> sys.stderr, "%s _/\_ %s " % (order[ia1], order[ia2]),
		print "%f cm %f cm translate" % (20*ia1, 20*ia2)

		genome1 = genomesAnc[ia1]
		genome2 = genomesAnc[ia2]

		# Les chromosomes a etudier
		chr1 = genome1.lstChr
		chr2 = genome2.lstChr
		chr1.extend(genome1.lstScaff)
		chr2.extend(genome2.lstScaff)


		table12 = buildOrthosTable(genome1, chr1, genome2, chr2)

		sys.stderr.write('.')

		# Initialisations
		table21 = buildOrthosTable(genome2, chr2, genome1, chr1)
		nb = sum([len(table12[c]) for c in table12])
		scaleX = 19. / float(nb)
		scaleY = 19. / float(sum([len(table21[c]) for c in table21]))
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


