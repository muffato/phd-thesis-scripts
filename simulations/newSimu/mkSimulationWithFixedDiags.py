#! /users/ldog/muffato/python

__doc__ = """
Prend un arbre phylogenetique et simule l'evolution d'un genome ancestral.
Genere des fichiers similaires a ceux d'Ensembl
-> Variante pour ressembler a la vraie evolution:
    - Contraintes de rearrangements sur poulet/opossum/rongeurs
    - Especes en cours d'assemblage: orny/xenope
    - Especes a 2X de couverture
"""

import sys
import math
import random
import utils.myMaths
import utils.myTools
import utils.myPhylTree


# Un taux specifique compris entre rate^-1 et rate^1
def randomAccel():
	return math.pow(arguments["eventAccel"], random.vonmisesvariate(0, arguments["vonMisesKappa"]) / math.pi)

# Un booleen au hasard
def randomBool():
	return random.choice([True,False])

# Un chromosome au hasard (en tenant compte des tailles)
def randomChromosome(genome):
	x = utils.myMaths.randomValue()
	x.initBisect([len(_) for _ in genome])
	while True:
		c = x.getRandomBisectPos()
		if len(genome[c]) >= 2:
			return c

# Un endroit du genome au hasard (mettre includeEnd=1 si on autorise l'extremite du chromosome comme position)
# Un chromosome plus long a plus de chance d'etre choisi
def randomPlace(genome, includeEnd):
	c = randomChromosome(genome)
	return (c,random.randint(1-int(includeEnd), len(genome[c])-1+int(includeEnd)))

# Une region du genome au hasard
def randomSlice(genome):
	# Le chromosome
	c = randomChromosome(genome)
	while True:
		# La taille de la region inversee suit la distribution
		l = int(len(genome[c]) * utils.myMaths.randomValue().myVonMises())
		x1 = random.randint(0, len(genome[c])-l)
		if (x1 != 0) or (x1+l != len(genome[c])):
			return (c,x1,x1+l)


# Casse un chromosome 
def doChrBreak(genome):
	(c,x) = randomPlace(genome, False)
	genome.append(genome[c][x:])
	genome[c] = genome[c][:x]


#
# La construction recursive des genomes
#
def buildGenomes(node, genome):

	# On rajoute les couvertures des sequencages
	if node in genomesScaffolds:
		print >> sys.stderr, "Applying 6x coverage on %s ..." % node,
		for _ in xrange(random.randint(2000,4000)):
			doChrBreak(genome)
		print >> sys.stderr, "OK (%d chromosomes)" % len(genome)

	# On ecrit le veritable genome ancestral en organisation de diagonales
	# On cree le dictionnaire diagonale -> chromosome du genome original
	print >> sys.stderr, "Writing %s genome ..." % node,
	s = phylTree.fileName[node]
	f = utils.myTools.myOpenFile(arguments["genomeFile"] % s, 'w')
	dicDiags = {}
	for (c,lst) in enumerate(genome):
		for i in lst:
			dicDiags[i] = c
			d = lstDiags[node][i]
			print >> f, utils.myTools.printLine( (c+1,utils.myTools.printLine(d," "),utils.myTools.printLine([0]*len(d)," ")) )
	f.close()
	print >> sys.stderr, "OK"

	if node in phylTree.items:

		# Le dictionnaire qui ucontiendra le futur de chaque diagonale
		# Pour chaque fils
		for (fils,dist) in phylTree.items[node]:
			
			# On construit un nouveau genome selon la distance qui les separe avec les taux de rearrangements fixes
			print >> sys.stderr, "Computing the %s -> %s branch ..." % (node,fils),
	
			# 1. On traduit les diagonales de node en diagonales de fils

			# Les diagonales (ou genes) qui apparaissent
			newDiags = set()
			# Dictionnaire ancien_id -> [nouveau_id] (y compris genes et duplications)
			translate = utils.myTools.defaultdict(list)
			future = eval(diagsEvolution.readline())
			for (i,corresp) in enumerate(future):
				if len(corresp) == 0:
					newDiags.add(i)
				else:
					# On choisit le chromosome majoritaire
					count = utils.myTools.defaultdict(int)
					for (d,c) in corresp:
						count[dicDiags[d]] += c
					newC = max(count, key=count.__getitem__)
					# On en deduit la diagonale ancestrale qui a engendre la nouvelle
					translate[ max([(c,d) for (d,c) in corresp if dicDiags[d] == newC])[1]].append( i )
			# On remplace 
			newGenome = [utils.myMaths.flatten([translate[d] for d in l if d in translate]) for l in genome]
			# On ajoute les nouvelles diagonales/genes
			for d in newDiags:
				(c,x) = randomPlace(newGenome, True)
				newGenome[c].insert(x, d)
			
			nbEvents = [int(round(arguments[x]*dist*randomAccel())) for x in ["chrInvertRate","chrTranslocRate","chrFusionRate","chrBreakRate"]]

			# CONSTRAINT-SET
			if fils == "Monodelphis domestica":
				# L'opossum a fusionne ses chromosomes
				nbEvents[2] *= 2
				nbEvents[3] /= 2
			elif fils == "Gallus gallus":
				# Le poulet a un genome proche de la version ancestrale, avec une preference pour les micro-chromosomes
				nbEvents[0] /= 2
				nbEvents[1] /= 2
				nbEvents[2] /= 2
				nbEvents[3] *= 2
			elif (fils == "Murinae") or (node == "Murinae"):
				# Les rongeurs ont diverge tres vite
				nbEvents[0] *= 2
				nbEvents[1] *= 2
				nbEvents[2] *= 2
				nbEvents[3] *= 2

			# La liste des evenements a appliquer
			events = []
			for i in xrange(4):
				events.extend( [i] * nbEvents[i] )
			random.shuffle(events)
			for evt in events:
				
				# C'est une inversion
				if evt == 0:
					# La region qui s'inverse
					(c,x1,x2) = randomSlice(newGenome)
					# Le nouveau chromosome avec la region qui s'inverse au milieu
					newGenome[c]= newGenome[c][:x1] + list(newGenome[c][x1:x2].__reversed__()) + newGenome[c][x2:]
				
				# C'est une translocation
				elif evt == 1:
				
					while True:
						# Les chromosomes a echanger
						c1 = randomChromosome(newGenome)
						c2 = randomChromosome(newGenome)
						
						# Differents si possible
						if c2 == c1:
							continue

						# Relativement l'un a l'autre, ils peuvent etre inverses
						if randomBool():
							newGenome[c2].reverse()

						# Des tailles aleatoires
						x1 = random.randint(1, len(newGenome[c1])-1)
						x2 = random.randint(1, len(newGenome[c2])-1)

						# On procede aux echanges
						chr1 = newGenome[c1][:x1] + newGenome[c2][x2:]
						chr2 = newGenome[c2][:x2] + newGenome[c1][x1:]
					
						# On verifie que les chromosomes crees ne sont pas trop petits
						if sum([diagsLength[fils][i] for i in chr1]) < arguments["minChrSize"]:
							continue
						if sum([diagsLength[fils][i] for i in chr2]) < arguments["minChrSize"]:
							continue
						break

					newGenome[c1] = chr1
					newGenome[c2] = chr2
					

				# C'est une fusion
				elif evt == 2:

					if len(newGenome) < 2:
						continue
					# Les deux chromosomes a fusionnerr
					(c1,c2) = random.sample(xrange(len(newGenome)), 2)
					# Quelles extremites se fusionnent
					if randomBool():
						newGenome[c1].reverse()
					if randomBool():
						newGenome[c2].reverse()
					# On met l'un au bout de l'autre (eventuellement en le retournant)
					newGenome[c1].extend(newGenome[c2])
					# Et on supprime l'original
					del newGenome[c2]
				
				# C'est une cassure
				elif evt == 3:

					while True:
						c = randomChromosome(newGenome)
						l = random.randint(1, len(newGenome[c])-1)
						# On verifie que les chromosomes crees ne sont pas trop petits
						if sum([diagsLength[fils][i] for i in newGenome[c][:l]]) < arguments["minChrSize"]:
							continue
						if sum([diagsLength[fils][i] for i in newGenome[c][l:]]) < arguments["minChrSize"]:
							continue
						break

					newGenome.append(newGenome[c][l:])
					newGenome[c] = newGenome[c][:l]
				
				# On ne devrait pas arriver ici ...
				else:
					print >> sys.stderr, "?%s?" % evt,
			
			print >> sys.stderr, "chr {<>%d ^%d +%d /%d}" % tuple(nbEvents), "-> %d chr" % len(newGenome)

			# On continue le parcours de l'arbre
			buildGenomes(fils, newGenome)


# Arguments
arguments = utils.myTools.checkArgs( \
	[("phylTree.conf",file), ("diagStatsFile",str), ("chrSizes",str)], \
	# Parametres d'evolution
	[("eventAccel",float,3.), ("vonMisesMean",float,.2), ("vonMisesKappa",float,2.), ("minChrSize",int,150), \
	# Rearrangements de chromosomes
	("chrInvertRate",float,1.7), ("chrTranslocRate",float,0.1), ("chrFusionRate",float,0.05), ("chrBreakRate",float,0.05), \
	# Couvertures de sequencage
	("6Xspecies",str,'Xenopus^tropicalis/Ornithorhynchus^anatinus/Takifugu^rubripes'), \
	# Fichiers
	("genomeFile",str,"simu/genes/genome.%s.bz2")], \
	__doc__ \
)

# L'arbre phylogenetique
phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])

genomesScaffolds = arguments["6Xspecies"].replace('^', ' ').split('/')

utils.myMaths.randomValue.vonmisesmean = arguments["vonMisesMean"]
utils.myMaths.randomValue.vonmiseskappa = arguments["vonMisesKappa"]

# Chargement des donnees
print >> sys.stderr, "Chargement des diagonales ",
diagsEvolution = utils.myTools.myOpenFile(arguments["diagStatsFile"], "r")
lstDiags = {}
diagsLength = {}
for _ in xrange(len(phylTree.allNames)):
	# Pour chaque espece, la liste des diagonales (y compris les singletons) en fonction des numeros des genes ancestraux
	(esp,lst) = eval(diagsEvolution.readline())
	lstDiags[esp] = lst
	diagsLength[esp] = [len(x) for x in lst]
	sys.stderr.write(".")
print >> sys.stderr, "OK"

tmp = [float(x) for x in arguments["chrSizes"].split()]
facteur = len(lstDiags[phylTree.root]) / sum(tmp)
chrlen = [int(x*facteur) for x in tmp]

# Le genome original
while True:
	lst = range(len(lstDiags[phylTree.root]))
	random.shuffle(lst)
	# Des longueurs de chromosomes aleatoires
	# On place les genes dans les chromosomes
	genome = []
	for nb in chrlen:
		if sum([diagsLength[phylTree.root][i] for i in lst[:nb]]) < arguments["minChrSize"]:
			break
		genome.append(lst[:nb])
		lst = lst[nb:]
	else:
		# Pour les erreurs d'arrondis ...
		genome[-1].extend(lst)
		break

# On construit le scenario evolutif
buildGenomes(phylTree.root, genome)

