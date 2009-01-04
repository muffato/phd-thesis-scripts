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
import collections

import utils.myMaths
import utils.myTools
import utils.myPhylTree


# +1 ou -1 au hasard
def randomStrand():
	return random.choice([-1,1])


# Un taux specifique compris entre rate^-1 et rate^1
def randomAccel():
	return math.pow(arguments["eventAccel"], random.vonmisesvariate(0, arguments["vonMisesKappa"]) / math.pi)


statsMgr = utils.myMaths.randomValue()
# Un chromosome au hasard (en tenant compte des tailles)
def randomChromosome(genome, minLength):
	statsMgr.initBisect([len(x) for x in genome])
	while True:
		c = statsMgr.getRandomBisectPos()
		if len(genome[c]) >= minLength:
			return c

# Un endroit du genome au hasard (regler includeEnd selon les extremites qu'on autorise)
# Un chromosome plus long a plus de chance d'etre choisi
def randomPlace(genome, includeEnd):
	n = 0
	while True:
		n += 1
		c = randomChromosome(genome, 1-includeEnd)
		return (c,random.randint(0, len(genome[c])-1+includeEnd))

# Une region du genome au hasard
def randomSlice(genome):
	# Le chromosome
	c = randomChromosome(genome, 0)
	# On passe de 0-1 a 0-len
	l = int(len(genome[c]) * utils.myMaths.randomValue().myVonMises())
	x1 = random.randint(0, len(genome[c])-1-l)
	return (c,x1,x1+l)


# Casse un chromosome 
def doChrBreak(genome):
	(c,x) = randomPlace(genome, -1)
	genome.append(genome[c][x+1:])
	genome[c] = genome[c][:x+1]


# Retourne la region genomique si necessaire
def applyStrand(chr, strand):
	if strand > 0:
		return chr
	else:
		return [(gene,-strand) for (gene,strand) in chr.__reversed__()]


#
# La construction recursive des genomes
#
def buildGenomes(genome, node, dupPairs):
	
	if node in phylTree.listSpecies:
		# On rajoute les couvertures des sequencages
		if node in phylTree.lstEsp2X:
			print >> sys.stderr, "Applying 2x coverage on %s ..." % node,
			cov = random.uniform(1.1, 2.9)
			for _ in xrange(int(sum(len(c) for c in genome)/cov)):
				doChrBreak(genome)
			# On n'a que 2/3 du genome
			random.shuffle(genome)
			genomes = genome[:(len(genome)*2)/3]
			print >> sys.stderr, "OK (%d chromosomes)" % len(genomes)
		
		elif node in phylTree.lstEsp6X:
			print >> sys.stderr, "Applying 6x coverage on %s ..." % node,
			for _ in xrange(random.randint(2000,4000)):
				doChrBreak(genome)
			print >> sys.stderr, "OK (%d chromosomes)" % len(genome)
		
		else:
			assert node in phylTree.lstEspFull
	else:
		assert node in phylTree.listAncestr

	# On ecrit le veritable genome ancestral et on prepare les familles
	# Initialement les familles de genes contiennent les genes eux-memes
	familles = {}
	print >> sys.stderr, "Writing %s genome ..." % node,
	s = phylTree.fileName[node]
	f = utils.myFile.openFile(arguments["genomeFile"] % s, 'w')
	for (c,lst) in enumerate(genome):
		for (i,(gene,strand)) in enumerate(lst):
			nom = "%s.%s" % (s,gene)
			familles[gene] = nom + " "
			print >> f, utils.myFile.myTSV.printLine( (c+1,i,i,strand,nom) )
	f.close()
	print >> sys.stderr, "OK"

	if node in phylTree.items:
	
		# Pour chaque fils
		for (fils,dist) in phylTree.items[node]:
			
			# On construit un nouveau genome selon la distance qui les separe avec les taux de rearrangements fixes
			print >> sys.stderr, "Computing the %s -> %s branch ..." % (node,fils),

			# L'evolution du repertoire de genes est calque sur les donnees reelles
			(correspondances,newGenes) = eval(geneEvolution.readline())
			todelete = set()
			nbDup = 0
			# Analyse du devenir de chaque gene
			for (i,x) in enumerate(correspondances):
				if len(x) == 0:
					# Le gene disparait
					todelete.add(i)
				else:
					# Un peu de duplications
					nbDup += len(x)-1
					if not arguments["tandemDuplications"]:
						newGenes.update(x[1:])
			# On ecrit donc le nouveau genome
			newGenome = []
			for chrom in genome:
				tmp = []
				for (g,strand) in chrom:
					# En enlevant les genes a supprimer
					if g not in todelete:
						# En les renommant
						tmp.append( (correspondances[g][0],strand) )
						if arguments["tandemDuplications"]:
							# Et en les dupliquant
							for x in correspondances[g][1:]:
								tmp.append( (x,randomStrand()) )
				newGenome.append(tmp)
			# Et on ajoute les nouveaux
			for g in newGenes:
				(c,x) = randomPlace(newGenome, 1)
				newGenome[c].insert(x, (g,randomStrand()))
			# On construit egalement les nouveaux groupes de genes dupliques
			newDupPairs = [[i] for i in newGenes]
			for l in dupPairs:
				new = []
				for x in l:
					new.extend(correspondances[x])
				newDupPairs.append(new)

			nbEvents = [int(round(arguments[x]*dist*randomAccel())) for x in ["chrInvertRate","chrTranslocRate","chrFusionRate","chrBreakRate"]]
			nbEvents.extend( [len(newGenes),nbDup,len(todelete)] )

			# CONSTRAINT-SET
			if fils == "Monodelphis domestica":
				# L'opossum a fusionne ses chromosomes
				nbEvents[2] *= 2
				nbEvents[3] /= 2
			elif fils == "Gallus gallus":
				# Le poulet a un genome proche de la version ancestrale, avec une preference pour les micro-chromosomes
				for i in xrange(3):
					nbEvents[i] /= 2
				nbEvents[3] *= 2
			elif phylTree.dicParents[fils]["Rodentia"] == "Rodentia":
				# Les rongeurs ont diverge tres vite
				for i in xrange(4):
					nbEvents[i] *= 2

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
					newGenome[c]= newGenome[c][:x1] + applyStrand(newGenome[c][x1:x2], -1) + newGenome[c][x2:]
				
				# C'est une translocation
				elif evt == 1:
					# La region qui se deplace
					(c,x1,x2) = randomSlice(newGenome)
					# On l'enleve
					r = newGenome[c][x1:x2]
					del newGenome[c][x1:x2]
					# On choisit une destination parmi le reste du genome
					(newC,newX) = randomPlace(newGenome, 1)
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
				
				# C'est une cassure
				elif evt == 3:
					doChrBreak(newGenome)
				
			
			print >> sys.stderr, "chr {<>%d ^%d +%d /%d} genes {+%d *%d -%d}" % tuple(nbEvents), "-> %d chr / %d genes" % (len(newGenome),sum([len(x) for x in newGenome]))

			# On continue le parcours de l'arbre
			subFam = buildGenomes(newGenome, fils, newDupPairs)
	
			# On rajoute les noms des genes dans les especes filles
			for (g,lnew) in enumerate(correspondances):
				for x in lnew:
					familles[g] += subFam[x]

		print >> sys.stderr, "Writing %s ancestral genes ..." % node,
		f = utils.myFile.openFile(arguments["ancGenesFile"] % s, 'w')
		for fam in familles.itervalues():
			print >> f, fam[:-1]
		f.close()
		print >> sys.stderr, "OK"

	else:
		# On echange quelques genes dupliques
		print >> sys.stderr, "Swapping some duplicates genes ...",
		poss = []
		for grp in dupPairs:
			poss.extend( [(i1,i2) for (i1,i2) in utils.myTools.myIterator.tupleOnStrictUpperList(grp)] )
		toSwap = random.sample(poss, int(len(poss) * arguments["duplicationGeneSwap"]))
		for (i1,i2) in toSwap:
			(familles[i1],familles[i2]) = (familles[i2],familles[i1])
		print >> sys.stderr, len(toSwap), "OK"
	
	return familles


# Arguments
arguments = utils.myTools.checkArgs( \
	[("phylTree.conf",file), ("genesStatsFile",str), ("iniChrNumber",int)], \
	# Parametres d'evolution
	[("eventAccel",float,3.), ("vonMisesMean",float,.2), ("vonMisesKappa",float,2.), \
	# Rearrangements de chromosomes
	("chrInvertRate",float,1.), ("chrTranslocRate",float,0.3), ("chrFusionRate",float,0.1), ("chrBreakRate",float,0.1), \
	# Divers
	("duplicationGeneSwap",float,.02), ("tandemDuplications",bool,False), \
	# Fichiers
	("genomeFile",str,"simu/genes/genes.%s.list.bz2"), \
	("ancGenesFile",str,"simu/ancGenes/ancGenes.%s.list.bz2")], \
	__doc__ \
)

# L'arbre phylogenetique
phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])

utils.myMaths.randomValue.vonmisesmean = arguments["vonMisesMean"]
utils.myMaths.randomValue.vonmiseskappa = arguments["vonMisesKappa"]

geneEvolution = utils.myFile.openFile(arguments["genesStatsFile"], "r")

# Le genome original
(root,nbGenes) = eval(geneEvolution.readline())
# Des genes
genes = [(x,randomStrand()) for x in xrange(nbGenes)]
random.shuffle(genes)
# Des longueurs de chromosomes aleatoires
tmp = [random.random() for _ in xrange(arguments["iniChrNumber"])]
facteur = nbGenes / sum(tmp)
chrlen = [int(x*facteur) for x in tmp]
# On place les genes dans les chromosomes
genome = []
for nb in chrlen:
	genome.append(genes[:nb])
	genes = genes[nb:]
# Pour les erreurs d'arrondis ...
genome[-1].extend(genes)

# On construit le scenario evolutif
buildGenomes(genome, root, [[i] for i in xrange(nbGenes)])

