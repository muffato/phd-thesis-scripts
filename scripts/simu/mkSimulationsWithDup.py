#! /users/ldog/muffato/python -OO

__doc__ = """
Prend un arbre phylogenetique et simule l'evolution d'un genome ancestral.
Genere des fichiers similaires a ceux d'Ensembl
-> Variante pour ressembler a la vraie evolution:
    - Contraintes de rearrangements sur poulet/opossum/rongeurs
    - Especes en cours d'assemblage: orny/xenope
    - Especes a 2X de couverture
"""


##################
# INITIALISATION #
##################

# Librairies
import sys
import math
import random
import utils.myMaths
import utils.myTools
import utils.myPhylTree

#############
# FONCTIONS #
#############


# +1 ou -1 au hasard
def randomStrand():
	return random.choice([-1,1])


# Un taux specifique compris entre rate^-1 et rate^1
def randomAccel():
	return math.pow(arguments["eventAccel"], random.vonmisesvariate(0, arguments["vonMisesKappa"]) / math.pi)


# Un chromosome au hasard (en tenant compte des tailles)
def randomChromosome(genome):
	x = utils.myMaths.randomValue()
	x.initBisect([len(_) for _ in genome])
	return x.getRandomBisectPos()

# Un endroit du genome au hasard (mettre includeEnd=1 si on autorise la fin du chromosome comme position)
# Un chromosome plus long a plus de chance d'etre choisi
def randomPlace(genome, includeEnd = 0):
	c = randomChromosome(genome)
	return (c,random.randint(0,len(genome[c])-1+includeEnd))


# Une region du genome au hasard
def randomSlice(genome):

	# Le chromosome
	c = randomChromosome(genome)
	# On passe de 0-1 a 0-len
	l = int(len(genome[c]) * utils.myMaths.randomValue().myVonMises())
	x1 = random.randint(0, len(genome[c])-1-l)
	return (c,x1,x1+l)


# Casse un chromosome 
def doChrBreak(genome):
	(c,x) = randomPlace(genome)
	genome.append(genome[c][x:])
	genome[c] = genome[c][:x]


# Retourne la region genomique si necessaire
def applyStrand(chr, strand):
	if strand > 0:
		return chr
	else:
		return [(gene,-strand) for (gene,strand) in chr.__reversed__()]



#
# La construction recursive des genomes
#
def buildGenomes(node):

	global nbTotalGenes

	# Pour chaque fils
	for (fils,dist) in phylTree.items.get(node,[]):
		
		# On construit un nouveau genome selon la distance qui les separe avec les taux de rearrangements fixes
		print >> sys.stderr, "Computing the %s -> %s branch ..." % (node,fils),
		newGenome = [list(x) for x in genomes[node]]
		
		# CONSTRAINT-SET
		if (fils == "Murinae") or (node == "Murinae"):
			geneAccel = 3.
		elif fils == "Gallus gallus":
			geneAccel = 0.3
		else:
			geneAccel = randomAccel()
		print >> sys.stderr, "genes {accel=%.3f" % geneAccel,

		rates = [arguments[x]*randomAccel() for x in ["geneLossWeight", "geneGainWeight", "geneDuplicationWeight"]]
		nbGeneEvents = (arguments["geneEventRate"] * dist * geneAccel) / sum(rates)
		nbPertes = int(nbGeneEvents * rates[0])
		nbGains = int(nbGeneEvents * rates[1])
		nbDup = int(nbGeneEvents * rates[2])

		# CONSTRAINT-SET
		if phylTree.dicParents[fils].get("Tetraodon nigroviridis","Euteleostomi") != "Euteleostomi":
			nbDup *= 2
			nbGains /= 2
	
		# Duplications
		dup = utils.myTools.defaultdict(list)
		for i in xrange(nbDup):
			(c,x) = randomPlace(newGenome, 0)
			dup[newGenome[c][x][0]].append(nbTotalGenes+i)
		dupGenes[fils] = dup
		print >> sys.stderr, "*%d" % nbDup,
		
		# On enleve les genes perdus
		for _ in xrange(nbPertes):
			(c,x) = randomPlace(newGenome, 0)
			del newGenome[c][x]
		print >> sys.stderr, "-%d" % nbPertes,
		
		# On place les nouveaux genes
		for i in xrange(nbTotalGenes, nbTotalGenes+nbDup+nbGains):
			(c,x) = randomPlace(newGenome, 1)
			newGenome[c].insert(x, (i,randomStrand()))
		nbTotalGenes += nbGains+nbDup
		print >> sys.stderr, "+%d}" % nbGains,

		# CONSTRAINT-SET
		if fils == "Monodelphis domestica":
			rates = (85,5,8,2)
		elif fils == "Gallus gallus":
			eventAccel = 0.3
			rates = (85,9,1,5)
		elif (fils == "Murinae") or (node == "Murinae"):
			eventAccel = 3.
		else:
			rates = [arguments[x]*randomAccel() for x in ["chrInvertWeight","chrTranslocWeight","chrFusionWeight","chrBreakWeight"]]
		
		# Rearrangements
		eventAccel = randomAccel()
		nbEvents = int(arguments["chrEventRate"] * dist * eventAccel)
		nb = [0,0,0,0]
		ratesmgr = utils.myMaths.randomValue()
		ratesmgr.initBisect(rates)
		for _ in xrange(nbEvents):
			evt = ratesmgr.getRandomBisectPos()
			
			# C'est une inversion
			if evt == 0:
				# La region qui s'inverse
				(c,x1,x2) = randomSlice(newGenome)
				# Le nouveau chromosome avec la region qui s'inverse au milieu
				newGenome[c]= newGenome[c][:x1] + applyStrand(newGenome[c][x1:x2], -1) + newGenome[c][x2:]
			
			# C'est une translocation
			elif evt == 1:
				# La region qui se deplace
				(c,x1,x2) = randomSlice(newGenome)
				# On l'enleve
				r = newGenome[c][x1:x2]
				del newGenome[c][x1:x2]
				# On choisit une destination parmi le reste du genome
				(newC,newX) = randomPlace(newGenome)
				# On l'insere (eventuellement en le retournant)
				newGenome[newC] = newGenome[newC][:newX] + applyStrand(r, randomStrand()) + newGenome[newC][newX:]
			
			# C'est une fusion
			elif evt == 2:
				if len(newGenome) < 2:
					continue
				# Les deux chromosomes a fusionnerr
				(c1,c2) = random.sample(xrange(len(newGenome)), 2)
				# On met l'un au bout de l'autre (eventuellement en le retournant)
				newGenome[c1].extend( applyStrand(newGenome[c2], randomStrand()) )
				# Et on supprime l'original
				del newGenome[c2]
			
			# Finalement, c'est une cassure
			else:
				doChrBreak(newGenome)
			
			nb[evt] += 1
		
		print >> sys.stderr, "chr {accel=%.3f <>%d ^%d +%d /%d} %d chromosomes" % (eventAccel,nb[0],nb[1],nb[2],nb[3],len(newGenome))

		# On boucle recursivement dessus
		genomes[fils] = newGenome
		buildGenomes(fils)

	# On rajoute les couvertures des sequencages
	if node in genomes2x:
		print >> sys.stderr, "Applying 2x coverage on %s ..." % node,
		tmp = genomes[node]
		for _ in xrange(random.randint(12000,18000)):
			doChrBreak(tmp)
		# On n'a que 2/3 du genome
		random.shuffle(tmp)
		genomes[node] = tmp[:(len(tmp)*2)/3]
		print >> sys.stderr, "OK (%d chromosomes)" % len(genomes[node])
	
	if node in genomesScaffolds:
		print >> sys.stderr, "Applying 6x coverage on %s ..." % node,
		tmp = genomes[node]
		for _ in xrange(random.randint(2000,4000)):
			doChrBreak(tmp)
		print >> sys.stderr, "OK (%d chromosomes)" % len(genomes[node])




# Imprime les donnees en tenant compte des duplications successives
def printData(node):

	# On ecrit le veritable genome ancestral et on prepare les familles
	# Initialement les familles de genes contiennent les genes eux-memes
	familles = {}
	print >> sys.stderr, "Writing %s genome (nbChr=%d) ..." % (node,len(genomes[node])),
	s = phylTree.fileName[node]
	f = utils.myTools.myOpenFile(arguments["genomeFile"] % s, 'w')
	for (c,lst) in enumerate(genomes[node]):
		for (i,(gene,strand)) in enumerate(lst):
			nom = "%s.%d" % (s,gene)
			familles[gene] = nom + " "
			print >> f, utils.myTools.printLine( (c+1,i,i,strand,nom) )
	f.close()
	print >> sys.stderr, "OK"

	# On ecrit les familles de genes ancestraux correspondantes
	if node in phylTree.items:
	
		for (fils,_) in phylTree.items[node]:
			subFam = printData(fils)

			# Quelques fausses assignations
			falseOrthologs = random.sample(subFam.keys(), int(len(subFam)*(100-arguments["orthologyQuality"])/100.))
			falseAssociations = [(x,subFam[x]) for x in falseOrthologs]
			random.shuffle(falseAssociations)
			for (old,new) in falseAssociations:
				subFam[old] = new

			# On rajoute les noms des genes dans les especes filles
			for g in familles:
				familles[g] += subFam.get(g,"")
			# ... en tenant compte des duplications
			for (ini,new) in dupGenes[fils].iteritems():
				for g in new:
					familles[ini] += subFam.get(g,"")

		print >> sys.stderr, "Writing %s ancestral genes ..." % node,
		f = utils.myTools.myOpenFile(arguments["ancGenesFile"] % s, 'w')
		for fam in familles.itervalues():
			print >> f, fam
		f.close()
		print >> sys.stderr, "OK"
	
	return familles


########
# MAIN #
########

# Arguments
arguments = utils.myTools.checkArgs( \
	[("phylTree.conf",file), ("root",str)], \
	# Genome initial
	[("nbOrigGenes",int,20000), ("nbOrigChr",int,20), \
	# Parametres d'evolution
	#("rearrRateAccel",float,1.732), ("vonMisesMean",float,.2), ("vonMisesKappa",float,2), \
	("eventAccel",float,3.), ("vonMisesMean",float,.2), ("vonMisesKappa",float,2.), \
	# Repertoire de genes
	#("geneLossRate",float,6), ("geneGainRate",float,6), ("geneDuplicationRate",float,3), \
	("geneEventRate",float,50.), ("geneLossWeight",float,40.), ("geneGainWeight",float,30.), ("geneDuplicationWeight",float,30.), \
	# Rearrangements de chromosomes
	("chrEventRate",float,2.), ("chrInvertWeight",float,91.), ("chrTranslocWeight",float,4.), ("chrFusionWeight",float,2.5), ("chrBreakWeight",float,2.5), \
	# Couvertures de sequencage
	("2Xspecies",str,'Loxodonta^africana/Echinops^telfairi/Dasypus^novemcinctus/Felis^catus/Erinaceus^europaeus/Myotis^lucifugus/Tupaia^belangeri/Otolemur^garnettii/Oryctolagus^cuniculus/Cavia^porcellus/Spermophilus^tridecemlineatus/Sorex^araneus'), \
	("6Xspecies",str,'Xenopus^tropicalis/Ornithorhynchus^anatinus/Takifugu^rubripes'), \
	# Qualite d'Ensembl
	("orthologyQuality",float,98.), \
	# Fichiers
	("genomeFile",str,"simu/genes/genes.%s.list.bz2"), \
	("ancGenesFile",str,"simu/ancGenes/ancGenes.%s.list.bz2")], \
	__doc__ \
)


# L'arbre phylogenetique
phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])
if arguments["root"] not in phylTree.listAncestr:
	print >> sys.stderr, "Unknown root '%s'" % arguments["root"]
	sys.exit(1)

utils.myTools.mkDir(arguments["genomeFile"])
utils.myTools.mkDir(arguments["ancGenesFile"])

genomes2x = arguments["2Xspecies"].replace('^', ' ').split('/')
genomesScaffolds = arguments["6Xspecies"].replace('^', ' ').split('/')

dupGenes = {}
genomes = {}

utils.myMaths.randomValue.vonmisesmean = arguments["vonMisesMean"]
utils.myMaths.randomValue.vonmiseskappa = arguments["vonMisesKappa"]

# Le genome original
tmp = [random.random() for _ in xrange(arguments["nbOrigChr"])]
facteur = arguments["nbOrigGenes"]/sum(tmp)
genome = []
nbTotalGenes = 0
for chrlen in tmp:
	nb = int(chrlen * facteur)
	genome.append( [(x+nbTotalGenes,randomStrand()) for x in xrange(nb)] )
	nbTotalGenes += nb
genomes[arguments["root"]] = genome

# On construit le scenario evolutif
buildGenomes(arguments["root"])

# On ecrit les familles de genes
printData(arguments["root"])

