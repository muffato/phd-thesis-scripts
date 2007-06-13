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
import utils.myPhylTree
import utils.myGenomes
import utils.myTools
import utils.myMaths


#############
# FONCTIONS #
#############



# +1 ou -1 au hasard
def randomStrand():
	return random.randint(0,1)*2-1

# Un taux specifique compris entre 1/rate et rate
def randomRate():
	return math.pow(options["rearrRateAccel"], random.vonmisesvariate(0, options["vonMisesKappa"]) / math.pi)

# Retourne la region genomique si necessaire
def applyStrand(chr, strand):
	if strand == 1:
		return list(chr)
	else:
		return [(gene,-strand) for (gene,strand) in chr.__reversed__()]


# Un endroit du genome au hasard (mettre includeEnd=1 si on autorise la fin du chromosome comme position)
# Un chromosome plus long a plus de chance d'etre choisi
def randomPlace(genome, includeEnd=0):
	s = sum([len(x) for x in genome])
	r = random.uniform(0, s)
	for c in xrange(len(genome)):
		if r < len(genome[c]):
			return (c, int(r))
		r -= len(genome[c])


# Une region du genome au hasard
def randomSlice(genome):

	while True:
		# Le chromosome au hasard
		(c,_) = randomPlace(genome)
		# La longueur de la region (plutot petite)
		l = int(abs(random.vonmisesvariate(0, 1)) * len(genome[c]) / math.pi)
		# On ne peut pas renvoyer une region d'un chromosome trop petit
		if len(genome[c]) > l:
			x1 = random.randint(0, len(genome[c])-l)
			x2 = x1 + l
			#(x1,x2) = sorted(random.sample(xrange(len(genome[c])), 2))
			return (c,x1,x2)


# Imprime un genome, en le transformant en scaffolds si necessaire
def printGenome(name, genome):

	if options["realLifeConstraints"]:
		if name in ["Loxodonta africana", "Echinops telfairi", "Dasypus novemcinctus", "Felis catus", "Erinaceus europaeus", "Myotis lucifugus", "Tupaia belangeri", "Otolemur garnettii", "Oryctolagus cuniculus", "Cavia porcellus", "Spermophilus tridecemlineatus"]:
			# especes en 2X
			for i in xrange(15000):
				(c,x) = randomPlace(genome)
				genome.append(genome[c][x:])
				genome[c] = genome[c][:x]
			random.shuffle(genome)
			genome = genome[:(len(genome)*2)/3]

		
		elif name in ["Xenopus tropicalis", "Ornithorhynchus anatinus", "Takifugu rubripes"]:
			# Larges scaffolds
			for i in xrange(3000):
				(c,x) = randomPlace(genome)
				genome.append(genome[c][x:])
				genome[c] = genome[c][:x]
	
	print >> sys.stderr, "Writing %s genome (nbChr=%d) ..." % (name,len(genome)),
	s = phylTree.fileName[name]
	f = utils.myTools.myOpenFile(options["genomeFile"] % s, 'w')
	for c in xrange(len(genome)):
		for i in xrange(len(genome[c])):
			(gene,strand) = genome[c][i]
			print >> f, "%d\t%d\t%d\t%d\t%s.%d" % (c+1,i,i,strand,s,gene)
	f.close()
	print >> sys.stderr, "OK"


# La construction recursive des genomes
def launchRecSimu(node, genomeIni):

	global nbTotalGenes

	# On ecrit le genome en question
	printGenome(node, genomeIni)

	# Si il n'y a pas de fils, c'est qu'on est arrive sur une espece moderne, on peut s'arreter
	if node not in phylTree.items:
		return

	# On ecrit les familles de genes ancestraux correspondantes
	print >> sys.stderr, "Writing %s ancestral genes ..." % node,
	f = utils.myTools.myOpenFile(options["ancGenesFile"] % phylTree.fileName[node], 'w')
	tmp = set()
	for (e1,e2) in utils.myTools.myMatrixIterator(phylTree.species[node], None, utils.myTools.myMatrixIterator.StrictUpperMatrix):
		tmp.update(phylTree.dicLinks[e1][e2])
	tmp = [phylTree.fileName[x] for x in tmp]
	for c in xrange(len(genomeIni)):
		for (gene,_) in genomeIni[c]:
			print >> f, '\t'.join(["%s.%d" % (x,gene) for x in tmp])
	f.close()
	print >> sys.stderr, "OK"


	# Pour chaque fils
	for (fils,dist) in phylTree.items[node]:
		
		# On construit un nouveau genome selon la distance qui les separe avec les taux de rearrangements fixes
		
		print >> sys.stderr, "Computing the %s -> %s branch ..." % (node,fils),
		
		newGenome = [list(x) for x in genomeIni]
		evolRate = randomRate()
		# CONSTRAINT-SET1
		if options["realLifeConstraints"]:
			if (fils == "Murinae") or (node == "Murinae"):
				evolRate = 4.
			elif fils == "Gallus gallus":
				evolRate = 0.2
		print >> sys.stderr, "rate=%.3f" % evolRate,

		# Perte de genes
		nbPertes = int(options["geneLossRate"] * dist * evolRate * randomRate())
		for i in xrange(nbPertes):
			(c,x) = randomPlace(newGenome)
			del newGenome[c][x]
		print >> sys.stderr, "lostGenes=%d" % nbPertes,
		
		# Gain de genes
		nbGains = int(options["geneGainRate"] * dist * evolRate * randomRate())
		for i in xrange(nbGains):
			(c,x) = randomPlace(newGenome, includeEnd=1)
			newGenome[c].insert(x, (nbTotalGenes,randomStrand()))
			nbTotalGenes += 1
		print >> sys.stderr, "newGenes=%d" % nbGains,

		# Rearrangements
		nbEvents = int(options["chrEventRate"] * dist * evolRate * randomRate())
		(InvertRate,TranslocRate,FusionRate,BreakRate) = [options[x]*randomRate() for x in ["chrInvertWeight","chrTranslocWeight","chrFusionWeight","chrBreakWeight"]]
		# CONSTRAINT-SET2
		if options["realLifeConstraints"]:
			if fils == "Monodelphis domestica":
				(InvertRate,TranslocRate,FusionRate,BreakRate)=(75,10,15,0)
			elif fils == "Gallus gallus":
				(InvertRate,TranslocRate,FusionRate,BreakRate)=(90,9,1,5)
		
		s = InvertRate + TranslocRate + FusionRate + BreakRate
		nb = [0,0,0,0]
		for i in xrange(nbEvents):
			r = random.uniform(0, s)
			
			# C'est une inversion
			if r < InvertRate:
				# La region qui s'inverse
				(c,x1,x2) = randomSlice(newGenome)
				# Le nouveau chromosome avec la region qui s'inverse au milieu
				newGenome[c]= newGenome[c][:x1] + applyStrand(newGenome[c][x1:x2], -1) + newGenome[c][x2:]
				nb[0] += 1
				continue
			
			r -= InvertRate
			# C'est une translocation
			if r < TranslocRate:
				# La region qui se deplace
				(c,x1,x2) = randomSlice(newGenome)
				# On l'enleve
				r = newGenome[c][x1:x2]
				del newGenome[c][x1:x2]
				# On choisit une destination parmi le reste du genome
				(newC,newX) = randomPlace(newGenome)
				# On l'insere (eventuellement en le retournant)
				newGenome[newC] = newGenome[newC][:newX] + applyStrand(r, randomStrand()) + newGenome[newC][newX:]
				nb[1] += 1
				continue
			
			r -= TranslocRate
			# C'est une fusion
			if r < FusionRate:
				# Les deux chromosomes a fusionnerr
				(c1,c2) = random.sample(xrange(len(newGenome)), 2)
				# On met l'un au bout de l'autre (eventuellement en le retournant)
				newGenome[c1].extend( applyStrand(newGenome[c2], randomStrand()) )
				# Et on supprime l'original
				del newGenome[c2]
				nb[2] += 1
				continue
			
			
			# Finalement, c'est une cassure
			# Le point de cassure
			(c,x) = randomPlace(newGenome)
			# On cree le nouveau chromosome
			newGenome.append(newGenome[c][x:])
			# On enleve l'ancienne partie
			newGenome[c] = newGenome[c][:x]
			nb[3] += 1
		
		print >> sys.stderr, "events=%s OK" % nb

		# On boucle recursivement dessus
		launchRecSimu(fils, newGenome)

########
# MAIN #
########

# Arguments
(noms_fichiers, options) = utils.myTools.checkArgs( \
	["phylTree.conf"], \
	[("root",str,""), ("realLifeConstraints",bool,False), \
	("nbOrigGenes",int,20000), ("nbMaxChr",int,20), \
	("geneLossRate",float,10), ("geneGainRate",float,10), ("chrEventRate",float,1), ("rearrRateAccel",float,1.7), ("vonMisesKappa",float,2), \
	("chrInvertWeight",float,90), ("chrTranslocWeight",float,5), ("chrFusionWeight",float,2.5), ("chrBreakWeight",float,2.5), \
	("genomeFile",str,"~/work/simu/genes/genes.%s.list.bz2"), \
	("ancGenesFile",str,"~/work/simu/ancGenes/ancGenes.%s.list.bz2")], \
	__doc__ \
)


# L'arbre phylogenetique
phylTree = utils.myPhylTree.PhylogeneticTree(noms_fichiers["phylTree.conf"])
if options["root"] not in phylTree.listAncestr:
	print >> sys.stderr, "Unknown root '%s'" % options["root"]
	sys.exit(1)

# Le genome original
nbTotalGenes = 0
genome = []
while nbTotalGenes < options["nbOrigGenes"]:
	if len(genome) == options["nbMaxChr"]-1:
		genome.append(range(nbTotalGenes, options["nbOrigGenes"]))
		break
	tmp = (options["nbOrigGenes"]-nbTotalGenes)/(options["nbMaxChr"]-len(genome))
	new = int(random.gauss(tmp, tmp/2))
	if new > 0:
		genome.append(range(nbTotalGenes, nbTotalGenes+new))
		nbTotalGenes += new



nbTotalGenes = options["nbOrigGenes"]
launchRecSimu(options["root"], [[(x,randomStrand()) for x in c] for c in genome])



