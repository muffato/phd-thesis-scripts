#! /users/ldog/muffato/python -OO

__doc__ = """
Prend un arbre phylogenetique et simule l'evolution d'un genome ancestral.
Genere des fichiers similaires a ceux d'Ensembl
"""


##################
# INITIALISATION #
##################

# Librairies
import sys
import math
import random
import utils.myBioObjects
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
	return math.pow(options["rearrRateAccel"], random.vonmisesvariate(0, 1) / math.pi)

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


# La construction recursive des genomes
def launchRecSimu(node, genomeIni):

	global nbTotalGenes

	# On ecrit le genome en question
	print >> sys.stderr, "Writing %s genome (nbChr=%d) ..." % (node,len(genomeIni)),
	s = phylTree.fileName[node]
	f = utils.myTools.myOpenFile(options["genomeFile"] % s, 'w')
	for c in xrange(len(genomeIni)):
		for i in xrange(len(genomeIni[c])):
			(gene,strand) = genomeIni[c][i]
			print >> f, "%d\t%d\t%d\t%d\t%s.%d" % (c+1,i,i,strand,s,gene)
	f.close()
	print >> sys.stderr, "OK"

	# Si il n'y a pas de fils, c'est qu'on est arrive sur une espece moderne, on peut s'arreter
	if node not in phylTree.items:
		return

	# On ecrit les familles de genes ancestraux correspondantes
	print >> sys.stderr, "Writing %s ancestral genes ..." % node,
	f = utils.myTools.myOpenFile(options["ancGenesFile"] % s, 'w')
	tmp = [phylTree.fileName[x] for x in phylTree.species[node]] + [s]
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
		s = (options["chrInvertWeight"] * randomRate()) + (options["chrTranslocWeight"] * randomRate()) \
		+ (options["chrFusionWeight"] * randomRate()) + (options["chrBreakWeight"] * randomRate())
		nb = [0,0,0,0]
		for i in xrange(nbEvents):
			r = random.uniform(0, s)
			
			# C'est une inversion
			if r < options["chrInvertWeight"]:
				# La region qui s'inverse
				(c,x1,x2) = randomSlice(newGenome)
				# Le nouveau chromosome avec la region qui s'inverse au milieu
				newGenome[c]= newGenome[c][:x1] + applyStrand(newGenome[c][x1:x2], -1) + newGenome[c][x2:]
				nb[0] += 1
				continue
			
			r -= options["chrInvertWeight"]
			# C'est une translocation
			if r < options["chrTranslocWeight"]:
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
			
			r -= options["chrTranslocWeight"]
			# C'est une fusion
			if r < options["chrFusionWeight"]:
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
	[("root",str,""), \
	("nbOrigGenes",int,20000), ("nbMaxChr",int,20), \
	("geneLossRate",float,10), ("geneGainRate",float,10), ("chrEventRate",float,1), ("rearrRateAccel",float,2), \
	("chrInvertWeight",float,90), ("chrTranslocWeight",float,5), ("chrFusionWeight",float,2.5), ("chrBreakWeight",float,2.5), \
	("genomeFile",str,"~/work/simu/genes/genes.%s.list.bz2"), \
	("ancGenesFile",str,"~/work/simu/ancGenes/ancGenes.%s.list.bz2")], \
	__doc__ \
)


# L'arbre phylogenetique
phylTree = utils.myBioObjects.PhylogeneticTree(noms_fichiers["phylTree.conf"])
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



