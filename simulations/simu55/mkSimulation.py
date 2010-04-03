#! /users/ldog/muffato/python

__doc__ = """
Prend un arbre phylogenetique et simule l'evolution d'un genome ancestral.
Genere des fichiers similaires a ceux d'Ensembl
-> Variante pour ressembler a la vraie evolution:
    - Contraintes de rearrangements sur poulet/opossum/rongeurs
    - Especes en cours d'assemblage: orny/xenope
    - Especes a 2X de couverture

*/v4 Pas de positionnement aleatoire des genes, et pas de longueurs de blocs de genes crees/perdus/dupliques car deduit du taux de clustering
"""

import sys
import math
import random
import operator
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
	return l if strand > 0 else [(gene,-strand) for (gene,strand) in reversed(l)]


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


# Score de ressemblance des profils des deux genes
def scoreSameEvents(v1, v2, lstindices):
	s = 0.
	n = 0
	for i in lstindices:
		s += v1[i] & v2[i]
		n += v1[i] | v2[i]
	return s/n if n != 0 else -1


# Renvoie le score relatif aux voisins et le sens dans lequel rentrer le gene
# Version ou on fait varier la liste des blocs
def scoreWithNeighboursB(lblocks, chrom, x, lstindices):
	nl = byname[chrom[max(x-1,0)][0]]
	nr = byname[chrom[min(x,len(chrom)-1)][0]]
	for block in lblocks:
		bl = byname[block[0]]
		br = byname[block[-1]]
		sf = scoreSameEvents(nl, bl, lstindices) + scoreSameEvents(nr, br, lstindices)
		sr = scoreSameEvents(nl, br, lstindices) + scoreSameEvents(nr, bl, lstindices)
		yield ((sf,1) if sf > sr else (sr,-1), block)


# Renvoie le score relatif aux voisins et le sens dans lequel rentrer le gene
# Version ou on fait varier la liste des positions
def scoreWithNeighboursP(block, genome, lpos, lstindices):
	bl = byname[block[0]]
	br = byname[block[-1]]
	for (c,x) in lpos:
		chrom = genome[c]
		nl = byname[chrom[max(x-1,0)][0]]
		nr = byname[chrom[min(x,len(chrom)-1)][0]]
		sf = scoreSameEvents(nl, bl, lstindices) + scoreSameEvents(nr, br, lstindices)
		sr = scoreSameEvents(nl, br, lstindices) + scoreSameEvents(nr, bl, lstindices)
		yield ((sf,1) if sf > sr else (sr,-1), c, x)


# On ecrit le veritable genome ancestral et on prepare les familles
# Initialement les familles de genes contiennent les genes eux-memes
def writeGenome(genome, node):
	familles = {}
	print >> sys.stderr, "Writing %s genome ..." % node,
	f = utils.myFile.openFile(arguments["out:genomeFile"] % phylTree.fileName[node], 'w')
	s = phylTree.fileName[node]
	for (c,lst) in enumerate(genome):
		assert len(lst) >= 1, [len(x) for x in genome]
		for (i,(gene,strand)) in enumerate(lst):
			nom = gene if isinstance(gene, str) else ("%s.%d" % (s,gene))
			familles[gene] = set([nom])
			print >> f, utils.myFile.myTSV.printLine( (c+1,i,i,strand,nom) )
	f.close()
	print >> sys.stderr, "OK"
	return familles


# On rajoute les couvertures des sequencages pour les especes modernes quand necessaire
def applyCoverage(genome, node):

	if node in phylTree.lstEsp2X:
		print >> sys.stderr, "Applying 2x coverage on %s ..." % node,
		# Nombre de scaffolds en fonction des longueurs moyennes observees
		nbscaff = int( sum(len(c) for c in genome) / random.uniform(1.1, 2.9) )
		for _ in xrange(nbscaff - len(genome)):
			doChrBreak(genome)
		# On n'a qu'une partie du genome
		random.shuffle(genome)
		del genome[int(len(genome) * random.uniform(0.55, 0.85)):]
		print >> sys.stderr, "OK (%d chromosomes)" % len(genome)
	
	elif node in phylTree.lstEsp6X:
		print >> sys.stderr, "Applying 6x coverage on %s ..." % node,
		for _ in xrange(random.randint(2000,4000)):
			doChrBreak(genome)
		print >> sys.stderr, "OK (%d chromosomes)" % len(genome)
	

# Ecrit la liste des genes ancestraux
def writeAncGenes(node, familles):
	print >> sys.stderr, "Writing %s ancestral genes ..." % node,
	f = utils.myFile.openFile(arguments["out:ancGenesFile"] % phylTree.fileName[node], 'w')
	for x in familles:
		print >> f, " ".join(familles[x])
	f.close()
	print >> sys.stderr, "OK"


# Definit la liste des genes a modifier (en lisant depuis le fichier reference d'evolution des genes)
def chooseGeneEventsFromFile(genome, fils, dist):

	# L'evolution du repertoire de genes est calque sur les donnees reelles
	(correspondances,newGenes) = eval(geneStatsFile.readline())

	# Nouveaux noms
	rename = dict( (g,l) for (g,l) in correspondances.iteritems() if len(l) > 0)

	# Les pertes de genes des especes a 2X ne sont pas bien mesurees ici, on prefere supprimer 1/3 du genome plus tard
	if fils in phylTree.lstEsp2X:
		rename.update( (g,[g]) for (g,l) in correspondances.iteritems() if len(l) == 0)

	return (rename, list(newGenes))


# Segmente une liste de genes et les clusterise
def mkClusterByScore(lstGenes, indices):
	
	print "NEWGENESET", len(lstGenes)

	if len(indices) == 0:
		print "NOINDICES"
		return [[g] for g in lstGenes]

	refG = random.choice(lstGenes)
	allblocks = [[refG]]
	lstGenes.remove(refG)
	
	epsilon = 1e-100
	# Une liste aleatoire de genes en associant les cooccurrences d'evenements
	while len(lstGenes) > 0:
		refS = byname[refG]

		# Un gene selon un hasard proportionnel aux scores des comparaisons
		if arguments["gene:clusteredGenesRatio"] < 0:
			l = [max(scoreSameEvents(refS, byname[g], indices),epsilon) for g in lstGenes]
			if max(l) <= epsilon:
				i = random.choice(range(len(lstGenes)))
				newBlock = True
				print "LOWRANDOM", epsilon
			else:
				i = utils.myMaths.randomValue.bisectChooser(l)()
				newBlock = False
				print "PROPRANDOM", l[i]

		# Un gene strictement au hasard
		elif random.random() > arguments["gene:clusteredGenesRatio"]:
			i = random.choice(range(len(lstGenes)))
			newBlock = True
			print "RANDOM", scoreSameEvents(refS, byname[lstGenes[i]], indices)

		# Un des genes avec le meilleur score
		else:
			l = [(scoreSameEvents(refS, byname[g], indices),i) for (i,g) in enumerate(lstGenes)]
			m = max(l)[0]
			if m <= epsilon:
				i = random.choice(range(len(lstGenes)))
				newBlock = True
				print "LOWBEST", l[i][0], m
			else:
				newBlock = False
				i = random.choice([i for (s,i) in l if s == m])
				print "BEST", m

		# On construit une structure en blocs de genes choisis selon le score
		if newBlock:
			allblocks.append([])
		refG = lstGenes.pop(i)
		allblocks[-1].append(refG)
	print "BLOCKSLENGTHS", [len(x) for x in allblocks]

	return allblocks


# Evenements sur les genes uniques
def performGeneEvents(genome, rename, newGenes, indices):

	nini = sum(len(chrom) for chrom in genome)
	nnew = len(newGenes)

	# Blocs de nouveaux genes
	newBlocks = [tuple((g,randomStrand()) for g in b) for b in mkClusterByScore(newGenes, indices)]

	# Renommage, suppression et insertion des genes dupliques en tandem
	nbDup = 0
	for (c,chrom) in enumerate(genome):
		chrom = [(g,strand) for (g,strand) in chrom if g in rename]
		tmp = []
		for (i,(g,strand)) in enumerate(chrom):
			if len(rename[g]) == 1:
				tmp.append( (rename[g][0],strand) )
			elif len(rename[g]) >= 2:
				# Identification des genes a dupliquer et repartition en blocs
				lb = list(scoreWithNeighboursB(mkClusterByScore(rename[g], indices), chrom, i, indices))
				inpos = max(lb)[1]
				tmp.extend( (x,strand) for x in inpos )
				newBlocks.extend(tuple((g,randomStrand()) for g in b) for (_,b) in lb if b != inpos)
				nbDup += 1
		genome[c] = tmp

	# Insertion des blocs 
	allPos = set((c,x) for (c,chrom) in enumerate(genome) for x in xrange(len(chrom)))
	todo = []
	
	# Une partie au meilleur endroit
	clusteredBlocks = random.sample(newBlocks, int(len(newBlocks)*arguments["gene:clusteredGenesRatio"]))
	for b in clusteredBlocks:
		scores = list(scoreWithNeighboursP(b, genome, allPos, indices))
		m = max(scores)[0][0]
		((_,s),c,x) = random.choice([x for x in scores if x[0][0] == m])
		todo.append( (applyStrand(b,s),(c,x)) )
		allPos.remove( (c,x) )

	# Une partie au hasard
	randomPosBlocks = set(newBlocks).difference(clusteredBlocks)
	todo.extend(itertools.izip(randomPosBlocks, random.sample(allPos, len(randomPosBlocks))))
	
	# Insertion des nouveaux blocs
	todo.sort(reverse=True, key=operator.itemgetter(1))
	for (l,(c,x)) in todo:
		genome[c][x:x] = l

	print >> sys.stderr, "genes: +%d *%d -%d ##%d" % (nnew,nbDup,nini-len(rename),sum(len(chrom) for chrom in genome)),


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
	elif ("Rodentia" in phylTree.indNames) and (phylTree.dicParents[fils]["Rodentia"] == "Rodentia"):
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
			(c1,c2) = random.sample(range(len(newGenome)), 2)
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
	
	applyCoverage(genome, node)
	familles = writeGenome(genome, node)

	if node in phylTree.items:
	
		for (fils,dist) in phylTree.items[node]:
			
			# On construit un nouveau genome selon la distance qui les separe avec les taux de rearrangements fixes
			print >> sys.stderr, "Computing the %s -> %s branch ..." % (node,fils),

			newGenome = [list(x) for x in genome]

			(rename,newGenes) = chooseGeneEvents(newGenome, fils, dist)
			indices = [phylTree.indNames[e] for e in phylTree.allNames if phylTree.isChildOf(e, fils) and (fils != e)]
			performGeneEvents(newGenome, rename, newGenes, indices)
			newGenome = [x for x in newGenome if len(x) > 0]
			chromEvents = chooseChromEvents(newGenome, fils)
			performChromEvents(newGenome, chromEvents)
			
			print >> sys.stderr, "OK"

			# On continue le parcours de l'arbre
			subFam = buildGenomes(newGenome, fils)
	
			# On rajoute les noms des genes dans les especes filles
			for (ini,new) in rename.iteritems():
				for g in new:
					if g in subFam:
						familles[ini].update(subFam[g])

		writeAncGenes(node, familles)
	
	return familles


# Cree un genome aleatoire, eventuellement clusterise en fonction des evenements de chaque gene
def randomGenomesWithFixedIDs():
	
	geneNames = eval(geneStatsFile.readline())
	assert set(geneNames).issubset(byname)

	allindices = range(1, len(phylTree.indNames))
	allblocks = mkClusterByScore(geneNames, allindices)
	random.shuffle(allblocks)
	genes = []
	for l in allblocks:
		genes.extend((x,randomStrand()) for x in l)
	for ((g1,_),(g2,_)) in utils.myTools.myIterator.slidingTuple(genes):
		print "FINAL", scoreSameEvents(byname[g1], byname[g2], allindices)
	return genes

# Arguments
arguments = utils.myTools.checkArgs(
	[("phylTree.conf",file), ("root",str), ("iniChrNumber",int)], [
	
	# Graine
	("seed",str,""),

	# Taux de rearrangements
	("rate:eventMaxAccel",float,1.5), ("rate:vonMisesKappa",float,2.),

	# Parametres des genes
	("gene:statsFile",str,""), ("gene:treesSignatures",str,""), ("gene:clusteredGenesRatio",float,.8),
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

if len(arguments["seed"]) > 0:
	random.setstate(eval(arguments["seed"]))
	print >> sys.stderr, "Random number generator internal state initialized from arguments"
else:
	print >> sys.stderr, "Random number generator internal state:", random.getstate()

print >> sys.stderr, "Building %s gene list" % arguments["root"],
assert utils.myFile.hasAccess(arguments["gene:statsFile"])
chooseGeneEvents = chooseGeneEventsFromFile
geneStatsFile = utils.myFile.openFile(arguments["gene:statsFile"], "r")

# Charge les signatures
class defaultdict2(dict):
	def __init__(self, val):
		self.defaultvalue = val

	def __getitem__(self, key):
		return dict.__getitem__(self, key) if key in self else self.defaultvalue
allsign = {}
byname = defaultdict2([0] * len(phylTree.indNames))
f = utils.myFile.openFile(arguments["gene:treesSignatures"], "r")
for l in f:
	(name,_,val) = l.partition(' ')
	s = tuple(eval(val.replace('None', '0')))
	# Rassemble les signatures identiques -> sauve de la memoire
	byname[name] = allsign.setdefault(s, s)
f.close()
genes = randomGenomesWithFixedIDs()
print >> sys.stderr, "OK"

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

