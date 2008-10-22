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
def buildGenomes(genome, node):

	global nbTotalGenes
			
	# On rajoute les couvertures des sequencages
	if node in genomesScaffolds:
		print >> sys.stderr, "Applying 6x coverage on %s ..." % node,
		for _ in xrange(random.randint(2000,4000)):
			doChrBreak(genome)
		print >> sys.stderr, "OK (%d chromosomes)" % len(genome)

	# On ecrit le veritable genome ancestral et on prepare les familles
	# Initialement les familles de genes contiennent les genes eux-memes
	familles = {}
	print >> sys.stderr, "Writing %s genome ..." % node,
	s = phylTree.fileName[node]
	f = utils.myTools.myOpenFile(arguments["genomeFile"] % s, 'w')
	for (c,lst) in enumerate(genome):
		for (i,(gene,strand)) in enumerate(lst):
			nom = "%s.%d" % (s,gene)
			familles[gene] = nom + " "
			print >> f, utils.myTools.printLine( (c+1,i,i,strand,nom) )
	f.close()
	print >> sys.stderr, "OK"

	if node in phylTree.items:
	
		# Pour chaque fils
		for (fils,dist) in phylTree.items[node]:
			
			# On construit un nouveau genome selon la distance qui les separe avec les taux de rearrangements fixes
			print >> sys.stderr, "Computing the %s -> %s branch ..." % (node,fils),
			
			newGenome = [list(x) for x in genome]
			dup = utils.myTools.defaultdict(list)
			maxGeneID = nbTotalGenes
			nbEvents = [int(round(arguments[x]*dist*randomAccel())) for x in ["chrInvertRate","chrTranslocRate","chrFusionRate","chrBreakRate","geneGainRate","geneDuplicationRate","geneLossRate"]]
			
			# CONSTRAINT-SET
			if phylTree.dicParents[fils].get("Tetraodon nigroviridis","Euteleostomi") != "Euteleostomi":
				nbEvents[4] /= 2
				nbEvents[5] *= 2
			elif fils == "Monodelphis domestica":
				nbEvents[2] *= 2
				nbEvents[3] /= 2
			elif fils == "Gallus gallus":
				for i in xrange(7):
					nbEvents[i] /= 2
			elif len(set([fils,node]).intersection(set(["Murinae", "Boreoeutheria"]))) > 0:
				for i in xrange(7):
					nbEvents[i] *= 2
			
			# On cree la liste des evenements et on les execute
			events = []
			for i in xrange(7):
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
				
				# C'est une cassure
				elif evt == 3:
					doChrBreak(newGenome)
				
				# Gain de gene
				elif evt == 4:
					(c,x) = randomPlace(newGenome, 1)
					newGenome[c].insert(x, (nbTotalGenes,randomStrand()))
					nbTotalGenes += 1

				# Duplication de gene
				elif evt == 5:
					while True:
						(c1,x1) = randomPlace(newGenome, 0)
						i = newGenome[c1][x1][0]
						# le gene duplique doit etre un des genes d'origine
						if (i < maxGeneID) and (i not in dup):
							break
					# Insertion de la nouvelle copie
					(c2,x2) = randomPlace(newGenome, 1)
					newGenome[c2].insert(x2, (nbTotalGenes,randomStrand()))
					# Mise a jour des listes de genes dupliques
					dup[i].append(nbTotalGenes)
					allDup.addLink([i,nbTotalGenes])
					nbTotalGenes += 1

				# Perte de gene
				elif evt == 6:
					while True:
						(c,x) = randomPlace(newGenome, 0)
						i = newGenome[c][x][0]
						# Le gene perdu doit etre un des genes d'origine et ne peut etre un gene qui s'est duplique
						if (i < maxGeneID) and (i not in dup):
							break
					del newGenome[c][x]
			
			print >> sys.stderr, "chr {<>%d ^%d +%d /%d} genes {+%d *%d -%d}" % tuple(nbEvents), "-> %d chr / %d genes" % (len(newGenome),sum([len(x) for x in newGenome]))

			# On boucle recursivement dessus
			subFam = buildGenomes(newGenome, fils)
	
			# On rajoute les noms des genes dans les especes filles
			for g in familles:
				familles[g] += subFam.get(g,"")
			# ... en tenant compte des duplications
			for (ini,new) in dup.iteritems():
				for g in new:
					familles[ini] += subFam.get(g,"")

		print >> sys.stderr, "Writing %s ancestral genes ..." % node,
		f = utils.myTools.myOpenFile(arguments["ancGenesFile"] % s, 'w')
		for fam in familles.itervalues():
			print >> f, fam[:-1]
		f.close()
		print >> sys.stderr, "OK"

	else:
		# On echange quelques genes dupliques
		print >> sys.stderr, "Swapping some duplicates genes ...",
		poss = []
		for grp in allDup:
			poss.extend( [(i1,i2) for (i1,i2) in utils.myTools.myIterator.tupleOnStrictUpperList(grp) if (i1 in familles) and (i2 in familles)] )
		toSwap = random.sample(poss, int(len(poss) * arguments["duplicationGeneSwap"]))
		for (i1,i2) in toSwap:
			(familles[i1],familles[i2]) = (familles[i2],familles[i1])
		print >> sys.stderr, len(toSwap), "OK"
	
	return familles


# Arguments
arguments = utils.myTools.checkArgs( \
	[("phylTree.conf",file), ("root",str), ("iniGeneNumber",int), ("iniChrNumber",int)], \
	# Parametres d'evolution
	[("eventAccel",float,3.), ("vonMisesMean",float,.2), ("vonMisesKappa",float,2.), \
	# Repertoire de genes
	("geneLossRate",float,6.), ("geneGainRate",float,6.), ("geneDuplicationRate",float,3.), \
	# Rearrangements de chromosomes
	("chrInvertRate",float,1.), ("chrTranslocRate",float,0.3), ("chrFusionRate",float,0.1), ("chrBreakRate",float,0.1), \
	# Couvertures de sequencage
	("6Xspecies",str,'Xenopus^tropicalis/Ornithorhynchus^anatinus/Takifugu^rubripes'), \
	# Qualite d'Ensembl
	("duplicationGeneSwap",float,.02), \
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

genomesScaffolds = arguments["6Xspecies"].replace('^', ' ').split('/')

utils.myMaths.randomValue.vonmisesmean = arguments["vonMisesMean"]
utils.myMaths.randomValue.vonmiseskappa = arguments["vonMisesKappa"]

allDup = utils.myTools.myCombinator()

# Le genome original
tmp = [random.random() for _ in xrange(arguments["iniChrNumber"])]
facteur = arguments["iniGeneNumber"]/sum(tmp)
chrlen = [int(x*facteur) for x in tmp]
chrlen[-1] += arguments["iniGeneNumber"] - sum(chrlen)
nbTotalGenes = 0
genome = []
for nb in chrlen:
	genome.append( [(x+nbTotalGenes,randomStrand()) for x in xrange(nb)] )
	nbTotalGenes += nb

# On construit le scenario evolutif
buildGenomes(genome, arguments["root"])

