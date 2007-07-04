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
import itertools
import utils.myTools
import utils.myPhylTree

#############
# FONCTIONS #
#############

# +1 ou -1 au hasard
def randomStrand():
	return random.choice([-1,1])

# Un taux specifique compris entre 1/rate et rate
def randomRate():
	return math.pow(options["rearrRateAccel"], random.vonmisesvariate(0, options["vonMisesKappa"]) / math.pi)

# Retourne la region genomique si necessaire
def applyStrand(chr, strand):
	if strand > 0:
		return chr
	else:
		return [(gene,-strand) for (gene,strand) in chr.__reversed__()]


# Un endroit du genome au hasard (mettre includeEnd=1 si on autorise la fin du chromosome comme position)
# Un chromosome plus long a plus de chance d'etre choisi
def randomPlace(genome, includeEnd = 0):
	tmp = [len(x) for x in genome]
	tmp[-1] += includeEnd
	r = random.randint(0, sum(tmp)-1)
	for (c,l) in enumerate(tmp):
		if r < l:
			return (c, r)
		r -= l


# Une region du genome au hasard
def randomSlice(genome):

	(c,_) = randomPlace(genome)
	
	l = int(abs(random.vonmisesvariate(0, options["vonMisesKappa"])) * len(genome[c]) / math.pi)
	x1 = random.randint(0, len(genome[c])-1-l)
	return (c,x1,x1+l)


# Casse un chromosome 
def doChrBreak(genome):
	(c,x) = randomPlace(genome)
	genome.append(genome[c][x:])
	genome[c] = genome[c][:x]



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
			evolRate = 3.
		elif fils == "Gallus gallus":
			evolRate = 0.3
		else:
			evolRate = randomRate()
		print >> sys.stderr, "rate=%.3f" % evolRate,

		nbPertes = int(options["geneLossRate"] * dist * evolRate * randomRate())
		nbGains = int(options["geneGainRate"] * dist * evolRate * randomRate())
		nbDup = int(options["geneDuplicationRate"] * dist * evolRate * randomRate())

		# CONSTRAINT-SET
		if phylTree.dicParents[fils]["Tetraodon nigroviridis"] != "Euteleostomi":
			nbDup *= 2
			nbGains /= 2
	
		# Duplications
		dup = utils.myTools.defaultdict(list)
		for i in xrange(nbDup):
			(c,x) = randomPlace(newGenome, 0)
			dup[newGenome[c][x][0]].append(nbTotalGenes+i)
		dupGenes[fils] = dup
		print >> sys.stderr, "dupGenes=%d" % nbDup,
		
		# On enleve les genes perdus
		for i in xrange(nbPertes):
			(c,x) = randomPlace(newGenome, 0)
			del newGenome[c][x]
		print >> sys.stderr, "lostGenes=%d" % nbPertes,
		
		# On place les nouveaux genes
		for i in xrange(nbTotalGenes, nbTotalGenes+nbDup+nbGains):
			(c,x) = randomPlace(newGenome, 1)
			newGenome[c].insert(x, (i,randomStrand()))
		nbTotalGenes += nbGains+nbDup
		print >> sys.stderr, "newGenes=%d" % nbGains,

		# Rearrangements
		nbEvents = int(options["chrEventRate"] * dist * randomRate() * evolRate)
		# CONSTRAINT-SET
		if fils == "Monodelphis domestica":
			(InvertRate,TranslocRate,FusionRate,BreakRate) = (85,5,8,2)
		elif fils == "Gallus gallus":
			(InvertRate,TranslocRate,FusionRate,BreakRate) = (85,9,1,5)
		else:
			(InvertRate,TranslocRate,FusionRate,BreakRate) = [options[x]*randomRate() for x in ["chrInvertWeight","chrTranslocWeight","chrFusionWeight","chrBreakWeight"]]
		
		s = InvertRate + TranslocRate + FusionRate + BreakRate
		nb = [0,0,0,0]
		for i in xrange(nbEvents):
			r = random.uniform(0, s)
			
			# C'est une inversion
			r -= InvertRate
			if r < 0:
				# La region qui s'inverse
				(c,x1,x2) = randomSlice(newGenome)
				# Le nouveau chromosome avec la region qui s'inverse au milieu
				newGenome[c]= newGenome[c][:x1] + applyStrand(newGenome[c][x1:x2], -1) + newGenome[c][x2:]
				nb[0] += 1
				continue
			
			# C'est une translocation
			r -= TranslocRate
			if r < 0:
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
			
			# C'est une fusion
			r -= FusionRate
			if r < 0 and len(newGenome) >= 2:
				# Les deux chromosomes a fusionnerr
				(c1,c2) = random.sample(xrange(len(newGenome)), 2)
				# On met l'un au bout de l'autre (eventuellement en le retournant)
				newGenome[c1].extend( applyStrand(newGenome[c2], randomStrand()) )
				# Et on supprime l'original
				del newGenome[c2]
				nb[2] += 1
				continue
			
			# Finalement, c'est une cassure
			doChrBreak(newGenome)
			nb[3] += 1
		
		print >> sys.stderr, "events=%s / %d chromosomes" % (nb,len(newGenome))

		# On boucle recursivement dessus
		genomes[fils] = newGenome
		buildGenomes(fils)

	# On rajoute les couvertures des sequencages
	if node in genomes2x:
		print >> sys.stderr, "Applying 2x coverage on %s ..." % node,
		tmp = genomes[node]
		for i in xrange(random.randint(12000,18000)):
			doChrBreak(tmp)
		# On n'a que 2/3 du genome
		random.shuffle(tmp)
		genomes[node] = tmp[:(len(tmp)*2)/3]
		print >> sys.stderr, "OK (%d chromosomes)" % len(genomes[node])
	
	if node in genomesScaffolds:
		print >> sys.stderr, "Applying 6x coverage on %s ..." % node,
		tmp = genomes[node]
		for i in xrange(random.randint(2000,4000)):
			doChrBreak(tmp)
		print >> sys.stderr, "OK (%d chromosomes)" % len(genomes[node])




# Imprime les donnees en tenant compte des duplications successives
def printData(node):

	# On ecrit le veritable genome ancestral
	print >> sys.stderr, "Writing %s genome (nbChr=%d) ..." % (node,len(genomes[node])),
	s = phylTree.fileName[node]
	f = utils.myTools.myOpenFile(options["genomeFile"] % s, 'w')
	for (c,lst) in enumerate(genomes[node]):
		for (i,(gene,strand)) in enumerate(lst):
			print >> f, "%d\t%d\t%d\t%d\t%s.%d" % (c+1,i,i,strand,s,gene)
	f.close()
	print >> sys.stderr, "OK"

	# Initialement les familles de genes contiennent les genes eux-memes
	familles = {}
	for lst in genomes[node]:
		for (gene,_) in lst:
			familles[gene] = "%s.%d " % (s,gene)
	
	for fils in phylTree.branches[node]:
		subFam = printData(fils)
		falseOrthologs = random.sample(subFam.keys(), int(len(subFam)*(100-options["orthologyQuality"])/100.))
		falseAssociations = [subFam[x] for x in falseOrthologs]
		random.shuffle(falseAssociations)
		for (old,new) in itertools.izip(falseOrthologs,falseAssociations):
			subFam[old] = new

		# On rajoute les noms des genes dans les especes filles
		for g in familles:
			familles[g] += subFam.get(g,"")
		# ... en tenant compte des duplications
		for (ini,new) in dupGenes[fils].iteritems():
			for g in new:
				familles[ini] += subFam.get(g,"")

	# On ecrit les familles de genes ancestraux correspondantes
	if node in phylTree.items:
		print >> sys.stderr, "Writing %s ancestral genes ..." % node,
		f = utils.myTools.myOpenFile(options["ancGenesFile"] % s, 'w')
		for fam in familles.itervalues():
			print >> f, fam
		f.close()
		print >> sys.stderr, "OK"
	
	return familles


########
# MAIN #
########

# Arguments
(noms_fichiers, options) = utils.myTools.checkArgs( \
	["phylTree.conf"], \
	[("root",str,""), ("orthologyQuality",float,99), \
	("nbOrigGenes",int,20000), ("nbOrigChr",int,20), \
	("chrEventRate",float,1.6), ("rearrRateAccel",float,1.5), ("vonMisesKappa",float,2), \
	("geneLossRate",float,6), ("geneGainRate",float,6), ("geneDuplicationRate",float,3), \
	("chrInvertWeight",float,91), ("chrTranslocWeight",float,4), ("chrFusionWeight",float,2.5), ("chrBreakWeight",float,2.5), \
	("genomeFile",str,"~/work/simu/genes/genes.%s.list.bz2"), \
	("ancGenesFile",str,"~/work/simu/ancGenes/ancGenes.%s.list.bz2")], \
	__doc__ \
)


# L'arbre phylogenetique
phylTree = utils.myPhylTree.PhylogeneticTree(noms_fichiers["phylTree.conf"])
if options["root"] not in phylTree.listAncestr:
	print >> sys.stderr, "Unknown root '%s'" % options["root"]
	sys.exit(1)


genomes2x = ["Loxodonta africana", "Echinops telfairi", "Dasypus novemcinctus", "Felis catus", "Erinaceus europaeus", "Myotis lucifugus", "Tupaia belangeri", "Otolemur garnettii", "Oryctolagus cuniculus", "Cavia porcellus", "Spermophilus tridecemlineatus", "Sorex araneus"]
genomesScaffolds = ["Xenopus tropicalis", "Ornithorhynchus anatinus", "Takifugu rubripes"]

dupGenes = {}
genomes = {}

# Le genome original
tmp = [random.random() for i in xrange(options["nbOrigChr"])]
facteur = options["nbOrigGenes"]/sum(tmp)
genome = []
nbTotalGenes = 0
for chrlen in tmp:
	nb = int(chrlen * facteur)
	genome.append( [(x+nbTotalGenes,randomStrand()) for x in xrange(nb)] )
	nbTotalGenes += nb
genomes[options["root"]] = genome

# On construit le scenario evolutif
buildGenomes(options["root"])

# On ecrit les familles de genes
printData(options["root"])



"""
# ESSAIS POUR randomSlice


# 5000 - 7500 - 2500 / 27"
def v1():
	x1 = random.randint(0, length-1)
	x2 = random.randint(x1, length-1)
	return (x1,x2)


# 1= 3400 - 6600 - 3200 / 36"
# 2= 3900 - 6100 - 2200 / 36"
def v2():
	l = int(abs(random.vonmisesvariate(0, 2)) * length / math.pi)
	l = 0
	x1 = random.randint(0, length-1-l)
	return (x1,x1+l)


# 2500 - 7500 - 5000 / 27"
def v3():
	l = random.randint(0, length-1)
	x1 = random.randint(0, length-1-l)
	return (x1,x1+l)


# 3333 - 6666 - 3333 / 30"
def v4():
	x = random.sample(xrange(length), 2)
	return (min(x),max(x))


# 3333 - 6666 - 3333 / 30"
def v5():
	x1 = random.randint(0, length-1)
	x2 = random.randint(0, length-1)
	return (min(x1,x2),max(x1,x2))

res1 = []
res2 = []
resL = []
for i in xrange(1000000):
	(x1,x2) = v2()
	res1.append(x1)
	res2.append(x2)
	resL.append(x2-x1+1)
print utils.myMaths.myStats(res1)
print utils.myMaths.myStats(res2)
print utils.myMaths.myStats(resL)
"""
