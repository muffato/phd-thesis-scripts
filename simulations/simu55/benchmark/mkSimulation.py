#!/usr/bin/env python2

__doc__ = """
Prend un arbre phylogenetique et simule l'evolution d'un genome ancestral.
Genere des fichiers similaires a ceux d'Ensembl
-> Variante pour ressembler a la vraie evolution:
    - Contraintes de rearrangements sur poulet/opossum/rongeurs
    - Especes en cours d'assemblage: orny/xenope
    - Especes a 2X de couverture
"""

import os
import sys
import math
import random
import itertools
import collections

import utils.myMaths
import utils.myTools
import utils.myPhylTree


# +1 ou -1 au hasard
def randomStrand():
	return random.choice([-1,1])


# Un taux specifique compris entre rate^-1 et rate^1
def randomAccel():
	return math.pow(arguments["rate:eventMaxAccel"], random.vonmisesvariate(0, arguments["rate:vonMisesKappa"]) / math.pi)


# Un chromosome au hasard (en tenant compte des tailles)
def randomChromosome(genome, minLength):
	r = utils.myMaths.randomValue.bisectChooser([len(x) for x in genome])
	while True:
		c = r()
		if len(genome[c]) >= minLength:
			return c


# Une region du genome au hasard
def randomSlice(genome):
	# Le chromosome
	c = randomChromosome(genome, 2)
	# On passe de [0,1] a [[1,len-1]]
	l = int((len(genome[c])-1) * utils.myMaths.randomValue.myVonMises(arguments["chr:vonMisesMean"], arguments["chr:vonMisesKappa"]))
	l = min(l+1, len(genome[c])-1)
	x1 = random.randint(0, len(genome[c])-l)
	return (c,x1,x1+l)


# Casse un chromosome 
def doChrBreak(genome):
	c = randomChromosome(genome, 2)
	x = random.randint(1, len(genome[c])-1)
	genome.append(genome[c][x:])
	del genome[c][x:]


# Retourne la region genomique si necessaire
def applyStrand(l, strand):
	return l if strand > 0 else [(gene,-strand) for (gene,strand) in l.__reversed__()]


# Liste de n tailles dont la somme vaut s, et valant chacune au minimum m
def randomSizes(n, s, m):
	while True:
		tmp = [random.random() for _ in xrange(n)]
		facteur = s / sum(tmp)
		lengths = [max(int(x*facteur),m) for x in tmp]
		# Pour les erreurs d'arrondis ...
		lengths[-1] += s-sum(lengths)
		if lengths[-1] >= m:
			return lengths


# On ecrit le veritable genome ancestral et on prepare les familles
# Initialement les familles de genes contiennent les genes eux-memes
def writeGenome(genome, node):
	print >> sys.stderr, "Writing %s genome ..." % node,
	f = utils.myFile.openFile(arguments["out:genomeFile"] % phylTree.fileName[node], 'w')
	if node in phylTree.listSpecies:
		t = node.split()
		s = t[0][0] + t[1][:3]
	else:
		s = node[:7]
	print >> f, ">%s" % s
	for (c,lst) in enumerate(genome):
		print >> f, "#%d" % c
		assert len(lst) >= 1, [len(x) for x in genome]
		print >> f, ' '.join("%s%d" % ("+" if strand > 0 else "-", gene) for (gene,strand) in lst), "$"
	f.close()
	print >> sys.stderr, "OK"


# Ecrit la liste des genes ancestraux
def writeAncGenes(node):
	print >> sys.stderr, "Writing %s ancestral genes ..." % node,
	os.link(arguments["out:ancGenesFile"], arguments["out:ancGenesFile"] % phylTree.fileName[node])
	print >> sys.stderr, "OK"


# Decompose un segment en sous-segments en suivant la distribution de Pareto
def segmentsLengths(n, alpha):
	res = []
	while n > 0:
		x = min(n, int(random.paretovariate(alpha)))
		res.append(x)
		n -= x
	return res


# Definit le nombre de de rearrangements a effectuer (au hasard)
def choiceChromEventsRandomly(genome, fils):
	
	dist = phylTree.parent[fils][1]
	nbEvents = [int(round(arguments[x]*arguments["chr:rateMultiplier"]*dist*randomAccel())) for x in ["chr:invertRate","chr:translocRate","chr:fusionRate","chr:breakRate"]]

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
	elif phylTree.dicParents[fils]["Rodentia"] == "Rodentia":
		# Les rongeurs ont diverge tres vite
		for i in xrange(4):
			nbEvents[i] *= 2

	return nbEvents


# Definit le nombre de de rearrangements a effectuer (selon un taux lu dans un arbre phylogenetique)
def choiceChromEventsFromFile(genome, fils):

	nbGenes = sum(len(x) for x in genome)
	dist = phylTreeRate.parent[fils][1]
	dist *= nbGenes * arguments["chr:rateMultiplier"] / (2*arguments["chr:invertRate"]+3*arguments["chr:translocRate"]+arguments["chr:fusionRate"]+arguments["chr:breakRate"])
	nbEvents = [int(round(arguments[x]*dist*randomAccel())) for x in ["chr:invertRate","chr:translocRate","chr:fusionRate","chr:breakRate"]]

	# CONSTRAINT-SET
	if fils == "Monodelphis domestica":
		# L'opossum a fusionne ses chromosomes
		nbEvents[2] *= 2
		nbEvents[3] /= 2
	elif fils == "Gallus gallus":
		# Le poulet a une preference pour les micro-chromosomes
		nbEvents[2] /= 2
		nbEvents[3] *= 2

	return nbEvents


# Evenements sur les chromosomes
def performChromEvents(newGenome, nbEvents):
	
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
			newGenome[c][x1:x2] = applyStrand(newGenome[c][x1:x2], -1)
		
		# C'est une translocation
		elif evt == 1:
			# La region qui se deplace
			(c,x1,x2) = randomSlice(newGenome)
			# On l'enleve
			r = newGenome[c][x1:x2]
			del newGenome[c][x1:x2]
			# On choisit une destination parmi le reste du genome
			newC = randomChromosome(newGenome, 1)
			newX = random.randint(0, len(newGenome[newC]))
			# On l'insere (eventuellement en le retournant)
			newGenome[newC][newX:newX] = applyStrand(r, randomStrand())
		
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
		
		# C'est une cassure
		elif evt == 3:
			doChrBreak(newGenome)

	print >> sys.stderr, "chr: <>%d ^%d +%d /%d" % tuple(nbEvents), "##%d" % len(newGenome),


# La construction recursive des genomes
def buildGenomes(genome, node):
	
	writeGenome(genome, node)

	if node in phylTree.items:
	
		for (fils,dist) in phylTree.items[node]:
			
			# On construit un nouveau genome selon la distance qui les separe avec les taux de rearrangements fixes
			print >> sys.stderr, "Computing the %s -> %s branch ..." % (node,fils),

			newGenome = [list(x) for x in genome]

			chromEvents = chooseChromEvents(newGenome, fils)
			performChromEvents(newGenome, chromEvents)
			
			print >> sys.stderr, "OK"

			# On continue le parcours de l'arbre
			buildGenomes(newGenome, fils)

		writeAncGenes(node)


# Arguments
arguments = utils.myTools.checkArgs(
	[("phylTree.conf",file), ("root",str), ("iniChrNumber",int)], [
	
	# Taux de rearrangements
	("rate:eventMaxAccel",float,1.5), ("rate:vonMisesKappa",float,2.),

	# Parametres des genes
	("gene:initialNumber",int,22000),
	
	# Rearrangements de chromosomes
	("chr:ratesPhylTree",str,""), ("chr:rateMultiplier",float,1.),
	("chr:invertRate",float,1.), ("chr:translocRate",float,0.3), ("chr:fusionRate",float,0.05), ("chr:breakRate",float,0.05),
	("chr:vonMisesMean",float,.2), ("chr:vonMisesKappa",float,2.),
	
	# Fichiers
	("out:genomeFile",str,"simu/genes/genes.%s.list.bz2"),
	("out:ancGenesFile",str,"simu/ancGenes/ancGenes.%s.list.bz2")],
	__doc__
)

# L'arbre phylogenetique
phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])

# Des genes
genes = [(x,randomStrand()) for x in xrange(arguments["gene:initialNumber"])]
random.shuffle(genes)

f = utils.myFile.openFile(arguments["out:ancGenesFile"], 'w')
for i in xrange(arguments["gene:initialNumber"]):
	print >> f, i
f.close()

# Des longueurs de chromosomes aleatoires
chrlen = randomSizes(arguments["iniChrNumber"], len(genes), 1)

# On place les genes dans les chromosomes
genome = []
for nb in chrlen:
	genome.append(genes[:nb])
	del genes[:nb]

if utils.myFile.hasAccess(arguments["chr:ratesPhylTree"]):
	phylTreeRate = utils.myPhylTree.PhylogeneticTree(arguments["chr:ratesPhylTree"])
	chooseChromEvents = choiceChromEventsFromFile
else:
	chooseChromEvents = choiceChromEventsRandomly

# On lance l'evolution
buildGenomes(genome, arguments["root"])

os.unlink(arguments["out:ancGenesFile"])

