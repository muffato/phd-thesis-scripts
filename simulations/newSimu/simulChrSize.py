#! /users/ldog/muffato/python

__doc__ = """
Teste differents modeles de cassures/fusions de chromosomes
"""

import sys
import math
import random
import utils.myMaths
import utils.myTools

# Arguments
arguments = utils.myTools.checkArgs( \
	[("iniChrNumber",int), ("iniGeneNumber",int), ("nbSamples",int), ("nbIter",int)], \
	[("chrTranslocRate",int,3), ("chrFusionRate",int,1), ("chrBreakRate",int,1), ("minChrSize",int,200)], \
	__doc__, \
)

utils.myMaths.randomValue.vonmisesmean = 0.2
utils.myMaths.randomValue.vonmiseskappa = 2.

count = [[0] * arguments["iniChrNumber"] for _ in xrange(arguments["nbIter"])]

def randomChromosome(genome, mode):
	while True:
		if mode == "r":
			c = random.randint(0, len(genome)-1)
		elif mode == "s":
			x = utils.myMaths.randomValue()
			x.initBisect(genome)
			c = x.getRandomBisectPos()
		elif mode == "i":
			x = utils.myMaths.randomValue()
			x.initBisect([1./x for x in genome])
			c = x.getRandomBisectPos()
		if genome[c] >= 2:
			return c

def randomLength(l, mode):
	if mode == "r":
		return random.randint(1, l-1)
	elif mode == "v":
		return 1 + int((l-1) * utils.myMaths.randomValue().myVonMises())



for j in xrange(arguments["nbSamples"]):

	nbGenes = arguments["iniGeneNumber"]

	# Des longueurs de chromosomes aleatoires
	while True:
		tmp = [random.random() for _ in xrange(arguments["iniChrNumber"])]
		facteur = nbGenes / sum(tmp)
		genome = [int(x*facteur) for x in tmp]
		genome[-1] += nbGenes-sum(genome)
		if (min(genome) >= arguments["minChrSize"]) and (sum(genome) == nbGenes):
			break

	events = ([0] * arguments["chrTranslocRate"]) + ([1] * arguments["chrFusionRate"]) + ([2] * arguments["chrBreakRate"])
	for i in xrange(arguments["nbIter"]):
		random.shuffle(events)
		for evt in events:

			if evt == 0:

				while True:
					# TRANSLOCATION
					c1 = randomChromosome(genome, "s")
					x1 = randomLength(genome[c1], "v")
					
					c2 = randomChromosome(genome, "s")
					x2 = randomLength(genome[c2], "v")
					
					if ((genome[c1]+x2-x1) < arguments["minChrSize"]) or ((genome[c2]+x1-x2) < arguments["minChrSize"]):
						continue

					genome[c1] += x2-x1
					genome[c2] += x1-x2
					break

			elif evt == 1:

				# FUSION
				(c1,c2) = random.sample(range(len(genome)), 2)
				genome[c1] += genome[c2]
				del genome[c2]

			elif evt == 2:

				while True:
					# CASSURE
					c = randomChromosome(genome, "s")
					x = randomLength(genome[c], "r")

					if ((genome[c]-x) < arguments["minChrSize"]) or (x < arguments["minChrSize"]):
						continue

					genome.append(x)
					genome[c] -= x
					break

			else:
				print >> sys.stderr, "Evenement inconnu", evt

		for (c,x) in enumerate(sorted(genome)):
			count[i][c] += x

for i in xrange(arguments["nbIter"]):
	print utils.myFile.myTSV.printLine(count[i], " ")

